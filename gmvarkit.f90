! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Core algorithms translated from the R gmvarkit package.
module gmvarkit_mod
   !! Gaussian and Student-t mixture vector autoregressive algorithms.
   use kind_mod, only: dp
   use linalg_mod, only: identity_matrix, invert_matrix, inverse_logdet, &
      kronecker_product, cholesky_lower, matrix_rank, symmetric_eigen, &
      general_eigenvalues
   use random_mod, only: random_uniform, random_standard_normal, random_gamma
   use stats_mod, only: sorted, median, quantile, normal_quantile
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd, &
      finite_difference_hessian
   use special_functions_mod, only: regularized_gamma_q, regularized_beta
   implicit none
   private

   type, public :: gmvarkit_model_t
      !! Reduced-form G-StMVAR model parameters arranged by regime.
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: ar(:, :, :, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: weight(:)
      real(dp), allocatable :: degrees_of_freedom(:)
      integer :: gaussian_regimes = 0
      integer :: info = 0
   end type gmvarkit_model_t

   type, public :: gmvarkit_regime_moments_t
      !! Stationary regime means and covariance matrices of lag vectors.
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: lag_covariance(:, :, :)
      integer :: info = 0
   end type gmvarkit_regime_moments_t

   type, public :: gmvarkit_location_parameters_t
      !! Equivalent intercept and stationary-mean parameterizations.
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: mean(:, :)
      integer :: info = 0
   end type gmvarkit_location_parameters_t

   type, public :: gmvarkit_evaluation_t
      !! Mixture likelihood terms, weights, and conditional moments.
      real(dp), allocatable :: mixing_weight(:, :)
      real(dp), allocatable :: regime_mean(:, :, :)
      real(dp), allocatable :: arch_scalar(:, :)
      real(dp), allocatable :: conditional_mean(:, :)
      real(dp), allocatable :: conditional_covariance(:, :, :)
      real(dp), allocatable :: log_likelihood_term(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: info = 0
   end type gmvarkit_evaluation_t

   type, public :: gmvarkit_simulation_t
      !! Recursive observations, regimes, and mixing weights by path.
      real(dp), allocatable :: series(:, :, :)
      integer, allocatable :: regime(:, :)
      real(dp), allocatable :: mixing_weight(:, :, :)
      integer :: info = 0
   end type gmvarkit_simulation_t

   type, public :: gmvarkit_forecast_t
      !! Monte Carlo point forecasts, quantiles, and regime-weight forecasts.
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: median(:, :)
      real(dp), allocatable :: quantile(:, :, :)
      real(dp), allocatable :: probability(:)
      real(dp), allocatable :: mixing_weight_mean(:, :)
      integer :: simulations = 0
      integer :: info = 0
   end type gmvarkit_forecast_t

   type, public :: gmvarkit_girf_t
      !! Generalized impulse responses for variables and regime probabilities.
      real(dp), allocatable :: response(:, :, :)
      real(dp), allocatable :: mixing_weight_response(:, :, :)
      integer :: simulations = 0
      integer :: info = 0
   end type gmvarkit_girf_t

   type, public :: gmvarkit_gfevd_t
      !! Generalized forecast-error variance decompositions by horizon.
      real(dp), allocatable :: decomposition(:, :, :)
      integer :: info = 0
   end type gmvarkit_gfevd_t

   type, public :: gmvarkit_girf_inference_t
      !! GIRF distributions across fixed, observed, or random initial histories.
      real(dp), allocatable :: point_response(:, :, :)
      real(dp), allocatable :: lower_response(:, :, :)
      real(dp), allocatable :: upper_response(:, :, :)
      real(dp), allocatable :: response_draw(:, :, :, :)
      real(dp), allocatable :: point_mixing_weight(:, :, :)
      real(dp), allocatable :: lower_mixing_weight(:, :, :)
      real(dp), allocatable :: upper_mixing_weight(:, :, :)
      real(dp), allocatable :: mixing_weight_draw(:, :, :, :)
      real(dp) :: confidence_level = 0.0_dp
      integer :: inner_simulations = 0
      integer :: outer_replications = 0
      integer :: info = 0
   end type gmvarkit_girf_inference_t

   type, public :: gmvarkit_gfevd_inference_t
      !! Mean and history-specific GFEVDs for variables and mixing weights.
      real(dp), allocatable :: decomposition(:, :, :)
      real(dp), allocatable :: individual(:, :, :, :)
      real(dp), allocatable :: mixing_weight_decomposition(:, :, :)
      real(dp), allocatable :: mixing_weight_individual(:, :, :, :)
      integer :: outer_replications = 0
      integer :: info = 0
   end type gmvarkit_gfevd_inference_t

   type, public :: gmvarkit_linear_irf_t
      !! Regime-specific linear impulse responses under recursive identification.
      real(dp), allocatable :: response(:, :, :)
      real(dp), allocatable :: impact(:, :)
      integer, allocatable :: cumulative_variable(:)
      integer :: regime = 0
      integer :: info = 0
   end type gmvarkit_linear_irf_t

   type, public :: gmvarkit_unconditional_moments_t
      !! Mixture mean and lag-zero-through-p second moments.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: autocovariance(:, :, :)
      real(dp), allocatable :: autocorrelation(:, :, :)
      real(dp), allocatable :: regime_autocovariance(:, :, :, :)
      integer :: info = 0
   end type gmvarkit_unconditional_moments_t

   type, public :: gmvarkit_pearson_residuals_t
      !! Raw or conditional-covariance-standardized Pearson residuals.
      real(dp), allocatable :: residual(:, :)
      logical :: standardized = .false.
      integer :: info = 0
   end type gmvarkit_pearson_residuals_t

   type, public :: gmvarkit_information_criteria_t
      !! Akaike, Hannan-Quinn, and Bayesian information criteria.
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: hqic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: info = 0
   end type gmvarkit_information_criteria_t

   type, public :: gmvarkit_structural_t
      !! Common impact matrix and regime-relative structural variances.
      real(dp), allocatable :: impact(:, :)
      real(dp), allocatable :: relative_variance(:, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp) :: reconstruction_error = huge(1.0_dp)
      integer :: reference_regime = 1
      integer :: info = 0
   end type gmvarkit_structural_t

   type, public :: gmvarkit_linear_irf_bootstrap_t
      !! Wild-bootstrap linear IRF paths and equal-tail confidence bounds.
      real(dp), allocatable :: point(:, :, :)
      real(dp), allocatable :: lower(:, :, :)
      real(dp), allocatable :: upper(:, :, :)
      real(dp), allocatable :: draw(:, :, :, :)
      real(dp) :: confidence_level = 0.0_dp
      integer :: requested_replications = 0
      integer :: successful_replications = 0
      integer :: info = 0
   end type gmvarkit_linear_irf_bootstrap_t

   type, public :: gmvarkit_structural_constraints_t
      !! Fixed, zero, sign, and positive-linear structural restrictions.
      logical, allocatable :: impact_fixed(:, :)
      real(dp), allocatable :: impact_value(:, :)
      integer, allocatable :: impact_sign(:, :)
      real(dp), allocatable :: lambda_mapping(:, :)
      real(dp), allocatable :: fixed_relative_variance(:, :)
      integer :: info = 0
   end type gmvarkit_structural_constraints_t

   type, public :: gmvarkit_fit_t
      !! Locally optimized reduced-form mixture VAR model and likelihood.
      type(gmvarkit_model_t) :: model
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: parameter_count = 0
      logical :: converged = .false.
      integer :: info = 0
   end type gmvarkit_fit_t

   type, public :: gmvarkit_genetic_fit_t
      !! Genetic-search diagnostics and locally refined best model.
      type(gmvarkit_fit_t) :: fit
      real(dp), allocatable :: best_objective(:)
      integer :: population_size = 0
      integer :: generations = 0
      integer :: evaluations = 0
      integer :: info = 0
   end type gmvarkit_genetic_fit_t

   type, public :: gmvarkit_multistart_fit_t
      !! Ranked local fits with canonical duplicate classification.
      type(gmvarkit_fit_t), allocatable :: fit(:)
      integer, allocatable :: order(:)
      logical, allocatable :: successful(:)
      logical, allocatable :: distinct(:)
      integer, allocatable :: duplicate_of(:)
      integer :: best_index = 0
      integer :: successful_count = 0
      integer :: distinct_count = 0
      integer :: info = 0
   end type gmvarkit_multistart_fit_t

   type, public :: gmvarkit_constraints_t
      !! Linear AR, shared-mean, and fixed-weight restrictions.
      real(dp), allocatable :: ar_mapping(:, :)
      integer, allocatable :: mean_group(:)
      real(dp), allocatable :: fixed_weight(:)
      integer :: info = 0
   end type gmvarkit_constraints_t

   type, public :: gmvarkit_inference_t
      !! Numerical score, Hessian covariance, and transformed-scale standard errors.
      real(dp), allocatable :: parameter(:)
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: standard_error(:)
      integer :: info = 0
   end type gmvarkit_inference_t

   type, public :: gmvarkit_structural_fit_t
      !! Direct structural likelihood fit and transformed-coordinate inference.
      type(gmvarkit_fit_t) :: fit
      type(gmvarkit_structural_t) :: structural
      type(gmvarkit_inference_t) :: inference
      integer :: info = 0
   end type gmvarkit_structural_fit_t

   type, public :: gmvarkit_structural_multistart_fit_t
      !! Ranked direct structural fits with observational duplicate links.
      type(gmvarkit_structural_fit_t), allocatable :: fit(:)
      integer, allocatable :: order(:)
      logical, allocatable :: successful(:)
      logical, allocatable :: distinct(:)
      integer, allocatable :: duplicate_of(:)
      integer :: best_index = 0
      integer :: successful_count = 0
      integer :: distinct_count = 0
      integer :: info = 0
   end type gmvarkit_structural_multistart_fit_t

   type, public :: gmvarkit_likelihood_profile_t
      !! Fixed-nuisance likelihood slices through transformed parameters.
      integer, allocatable :: parameter(:)
      real(dp), allocatable :: center(:)
      real(dp), allocatable :: value(:, :)
      real(dp), allocatable :: log_likelihood(:, :)
      logical, allocatable :: valid(:, :)
      logical :: structural = .false.
      integer :: info = 0
   end type gmvarkit_likelihood_profile_t

   type, public :: gmvarkit_regime_conversion_t
      !! Student-to-Gaussian regime conversion and optional refit.
      type(gmvarkit_fit_t) :: fit
      type(gmvarkit_structural_t) :: structural
      integer, allocatable :: old_to_new(:)
      integer, allocatable :: new_to_old(:)
      logical, allocatable :: converted(:)
      logical :: has_structural = .false.
      logical :: reestimated = .false.
      integer :: info = 0
   end type gmvarkit_regime_conversion_t

   type, public :: gmvarkit_companion_eigen_t
      !! Regimewise VAR companion roots and stability flags.
      complex(dp), allocatable :: eigenvalue(:, :)
      real(dp), allocatable :: modulus(:, :)
      real(dp), allocatable :: spectral_radius(:)
      logical, allocatable :: stationary(:)
      logical, allocatable :: near_unit_root(:)
      real(dp) :: tolerance = 0.0_dp
      integer :: info = 0
   end type gmvarkit_companion_eigen_t

   type, public :: gmvarkit_covariance_eigen_t
      !! Covariance and pairwise generalized eigenvalue diagnostics.
      real(dp), allocatable :: eigenvalue(:, :)
      logical, allocatable :: near_singular(:)
      integer, allocatable :: pair(:, :)
      real(dp), allocatable :: ratio_eigenvalue(:, :)
      real(dp), allocatable :: minimum_separation(:)
      logical, allocatable :: weakly_identified(:)
      real(dp) :: positive_definite_tolerance = 0.0_dp
      real(dp) :: identification_tolerance = 0.0_dp
      integer :: info = 0
   end type gmvarkit_covariance_eigen_t

   type, public :: gmvarkit_hypothesis_test_t
      !! Chi-square Wald or likelihood-ratio test result.
      real(dp) :: statistic = 0.0_dp
      integer :: degrees_of_freedom = 0
      real(dp) :: p_value = 1.0_dp
      integer :: info = 0
   end type gmvarkit_hypothesis_test_t

   type, public :: gmvarkit_score_t
      !! Observationwise transformed-parameter scores and OPG information.
      real(dp), allocatable :: parameter(:)
      real(dp), allocatable :: observation(:, :)
      real(dp), allocatable :: total(:)
      real(dp), allocatable :: opg(:, :)
      integer :: info = 0
   end type gmvarkit_score_t

   type, public :: gmvarkit_quantile_residuals_t
      !! Sequential multivariate quantile residuals and conditional CDF values.
      real(dp), allocatable :: residual(:, :)
      real(dp), allocatable :: cdf(:, :)
      integer :: info = 0
   end type gmvarkit_quantile_residuals_t

   type, public :: gmvarkit_residual_tests_t
      !! Quantile-residual normality, correlation, and variance tests.
      type(gmvarkit_hypothesis_test_t) :: normality
      type(gmvarkit_hypothesis_test_t), allocatable :: autocorrelation(:)
      type(gmvarkit_hypothesis_test_t), allocatable :: heteroskedasticity(:)
      integer, allocatable :: autocorrelation_lag(:)
      integer, allocatable :: heteroskedasticity_lag(:)
      logical :: parameter_corrected = .false.
      integer :: info = 0
   end type gmvarkit_residual_tests_t

   interface gmvarkit_profile_likelihood
      module procedure gmvarkit_profile_likelihood_reduced
      module procedure gmvarkit_profile_likelihood_structural
   end interface gmvarkit_profile_likelihood

   public :: gmvarkit_regime_moments
   public :: gmvarkit_location_parameters, gmvarkit_model_from_regime_means
   public :: gmvarkit_evaluate
   public :: gmvarkit_simulate, gmvarkit_forecast
   public :: gmvarkit_girf, gmvarkit_gfevd
   public :: gmvarkit_girf_inference, gmvarkit_gfevd_inference
   public :: gmvarkit_linear_irf
   public :: gmvarkit_unconditional_moments
   public :: gmvarkit_pearson_residuals, gmvarkit_information_criteria
   public :: gmvarkit_structural_covariances, gmvarkit_identify_structural
   public :: gmvarkit_reorder_structural, gmvarkit_swap_structural_signs
   public :: gmvarkit_rebase_structural
   public :: gmvarkit_linear_irf_bootstrap
   public :: gmvarkit_estimate_structural
   public :: gmvarkit_structural_multistart_estimate
   public :: gmvarkit_profile_likelihood
   public :: gmvarkit_convert_student_regimes
   public :: gmvarkit_companion_eigenvalues
   public :: gmvarkit_covariance_eigenvalues
   public :: gmvarkit_estimate
   public :: gmvarkit_multistart_estimate
   public :: gmvarkit_genetic_estimate
   public :: gmvarkit_estimate_constrained
   public :: gmvarkit_inference, gmvarkit_wald_test
   public :: gmvarkit_likelihood_ratio
   public :: gmvarkit_score_matrix, gmvarkit_rao_test
   public :: gmvarkit_quantile_residuals, gmvarkit_quantile_residual_tests

contains

   pure function gmvarkit_structural_covariances(impact, relative_variance, &
      reference_regime) result(out)
      !! Construct regime covariances from structural impact parameters.
      real(dp), intent(in) :: impact(:, :) !! Common structural impact matrix W.
      real(dp), intent(in) :: relative_variance(:, :) !! Positive lambda values by shock and regime.
      integer, intent(in), optional :: reference_regime !! Regime normalized to unit relative variances.
      type(gmvarkit_structural_t) :: out
      real(dp), allocatable :: scaled(:, :)
      integer :: variables, regimes, regime, selected_reference

      selected_reference = 1
      if (present(reference_regime)) selected_reference = reference_regime
      variables = size(impact, 1)
      regimes = size(relative_variance, 2)
      if (variables < 1 .or. size(impact, 2) /= variables .or. &
         size(relative_variance, 1) /= variables .or. regimes < 1 .or. &
         any(relative_variance <= 0.0_dp) .or. selected_reference < 1 .or. &
         selected_reference > regimes) then
         out%info = 1
         return
      end if
      out%impact = impact
      out%relative_variance = relative_variance
      allocate(out%covariance(variables, variables, regimes))
      do regime = 1, regimes
         scaled = impact*spread(sqrt(relative_variance(:, regime)), 1, variables)
         out%covariance(:, :, regime) = matmul(scaled, transpose(scaled))
      end do
      out%reconstruction_error = 0.0_dp
      out%reference_regime = selected_reference
   end function gmvarkit_structural_covariances

   pure function gmvarkit_identify_structural(model, reference_regime, &
      comparison_regime, tolerance) result(out)
      !! Identify common structural shocks through covariance heteroskedasticity.
      type(gmvarkit_model_t), intent(in) :: model !! Reduced-form G-StMVAR model.
      integer, intent(in), optional :: reference_regime !! Covariance normalized to WW transpose.
      integer, intent(in), optional :: comparison_regime !! Covariance fixing the shock ordering.
      real(dp), intent(in), optional :: tolerance !! Relative simultaneous-diagonalization tolerance.
      type(gmvarkit_structural_t) :: out
      real(dp), allocatable :: values(:), vectors(:, :), square_root(:, :)
      real(dp), allocatable :: inverse_square_root(:, :), inverse_impact(:, :)
      real(dp), allocatable :: normalized(:, :), transformed(:, :), lambda(:, :)
      real(dp), allocatable :: reconstruction(:, :, :)
      real(dp) :: selected_tolerance, scale, error
      integer :: variables, regimes, reference, comparison, regime, shock, info

      selected_tolerance = 1.0e-8_dp
      if (present(tolerance)) selected_tolerance = tolerance
      reference = 1
      if (present(reference_regime)) reference = reference_regime
      if (.not. valid_model(model) .or. selected_tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      comparison = merge(2, 1, regimes > 1)
      if (reference == comparison .and. regimes > 1) comparison = 1
      if (present(comparison_regime)) comparison = comparison_regime
      if (reference < 1 .or. reference > regimes .or. comparison < 1 .or. &
         comparison > regimes .or. (regimes > 1 .and. comparison == reference)) then
         out%info = 2
         return
      end if
      call symmetric_eigen(model%covariance(:, :, reference), values, vectors, info)
      if (info /= 0 .or. any(values <= 0.0_dp)) then
         out%info = 100 + max(1, info)
         return
      end if
      square_root = matmul(vectors, matmul( &
         spread(sqrt(values), 1, variables)*identity_matrix(variables), &
         transpose(vectors)))
      inverse_square_root = matmul(vectors, matmul( &
         spread(1.0_dp/sqrt(values), 1, variables)*identity_matrix(variables), &
         transpose(vectors)))
      if (regimes > 1) then
         normalized = matmul(inverse_square_root, matmul( &
            model%covariance(:, :, comparison), inverse_square_root))
         call symmetric_eigen(0.5_dp*(normalized + transpose(normalized)), &
            values, vectors, info)
         if (info /= 0 .or. any(values <= 0.0_dp)) then
            out%info = 200 + max(1, info)
            return
         end if
         out%impact = matmul(square_root, vectors)
      else
         out%impact = square_root
      end if
      do shock = 1, variables
         if (out%impact(shock, shock) < 0.0_dp) then
            out%impact(:, shock) = -out%impact(:, shock)
         end if
      end do
      call invert_matrix(out%impact, inverse_impact, info)
      if (info /= 0) then
         out%info = 300 + info
         return
      end if
      allocate(lambda(variables, regimes), reconstruction(variables, variables, regimes))
      error = 0.0_dp
      do regime = 1, regimes
         transformed = matmul(inverse_impact, matmul( &
            model%covariance(:, :, regime), transpose(inverse_impact)))
         lambda(:, regime) = [(transformed(shock, shock), shock=1, variables)]
         if (any(lambda(:, regime) <= 0.0_dp)) then
            out%info = 400 + regime
            return
         end if
         reconstruction(:, :, regime) = matmul( &
            out%impact*spread(sqrt(lambda(:, regime)), 1, variables), &
            transpose(out%impact*spread(sqrt(lambda(:, regime)), 1, variables)))
         scale = max(1.0_dp, maxval(abs(model%covariance(:, :, regime))))
         error = max(error, maxval(abs(reconstruction(:, :, regime) - &
            model%covariance(:, :, regime)))/scale)
      end do
      out%relative_variance = lambda
      out%covariance = reconstruction
      out%reference_regime = reference
      out%reconstruction_error = error
      if (error > selected_tolerance) out%info = 3
   end function gmvarkit_identify_structural

   pure function gmvarkit_reorder_structural(structural, permutation) result(out)
      !! Reorder structural shocks without changing reduced-form covariances.
      type(gmvarkit_structural_t), intent(in) :: structural !! Identified structural parameters.
      integer, intent(in) :: permutation(:) !! New ordering of the structural shocks.
      type(gmvarkit_structural_t) :: out
      logical, allocatable :: seen(:)
      integer :: variables, shock

      if (.not. valid_structure(structural)) then
         out%info = 1
         return
      end if
      variables = size(structural%impact, 1)
      if (size(permutation) /= variables .or. any(permutation < 1) .or. &
         any(permutation > variables)) then
         out%info = 2
         return
      end if
      allocate(seen(variables), source=.false.)
      do shock = 1, variables
         if (seen(permutation(shock))) then
            out%info = 3
            return
         end if
         seen(permutation(shock)) = .true.
      end do
      out = structural
      out%impact = structural%impact(:, permutation)
      out%relative_variance = structural%relative_variance(permutation, :)
   end function gmvarkit_reorder_structural

   pure function gmvarkit_swap_structural_signs(structural, shocks) result(out)
      !! Reverse selected impact-matrix columns without changing covariances.
      type(gmvarkit_structural_t), intent(in) :: structural !! Identified structural parameters.
      integer, intent(in) :: shocks(:) !! Structural shock columns whose signs are reversed.
      type(gmvarkit_structural_t) :: out
      logical, allocatable :: seen(:)
      integer :: variables, index

      if (.not. valid_structure(structural)) then
         out%info = 1
         return
      end if
      variables = size(structural%impact, 1)
      if (any(shocks < 1) .or. any(shocks > variables)) then
         out%info = 2
         return
      end if
      allocate(seen(variables), source=.false.)
      out = structural
      do index = 1, size(shocks)
         if (seen(shocks(index))) cycle
         out%impact(:, shocks(index)) = -out%impact(:, shocks(index))
         seen(shocks(index)) = .true.
      end do
   end function gmvarkit_swap_structural_signs

   pure function gmvarkit_rebase_structural(structural, reference_regime) &
      result(out)
      !! Change the unit-variance reference regime without changing covariances.
      type(gmvarkit_structural_t), intent(in) :: structural !! Identified structural parameters.
      integer, intent(in) :: reference_regime !! Regime assigned an identity Lambda matrix.
      type(gmvarkit_structural_t) :: out
      real(dp), allocatable :: reference_variance(:)
      integer :: regimes, variables

      if (.not. valid_structure(structural)) then
         out%info = 1
         return
      end if
      regimes = size(structural%relative_variance, 2)
      variables = size(structural%impact, 1)
      if (reference_regime < 1 .or. reference_regime > regimes) then
         out%info = 2
         return
      end if
      out = structural
      reference_variance = structural%relative_variance(:, reference_regime)
      out%impact = structural%impact*spread(sqrt(reference_variance), 1, variables)
      out%relative_variance = structural%relative_variance/ &
         spread(reference_variance, 2, regimes)
      out%reference_regime = reference_regime
   end function gmvarkit_rebase_structural

   pure function gmvarkit_companion_eigenvalues(model, tolerance) result(out)
      !! Diagnose regime stability from VAR companion eigenvalues.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      real(dp), intent(in), optional :: tolerance !! Distance from one flagged as near-unit; default 0.0015.
      type(gmvarkit_companion_eigen_t) :: out
      real(dp), allocatable :: companion(:, :)
      complex(dp), allocatable :: roots(:)
      real(dp) :: selected_tolerance
      integer :: variables, lags, regimes, states, regime, lag, info

      selected_tolerance = 0.0015_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (.not. valid_model(model) .or. selected_tolerance < 0.0_dp .or. &
         selected_tolerance >= 1.0_dp) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      states = variables*lags
      allocate(companion(states, states))
      allocate(out%eigenvalue(states, regimes))
      allocate(out%modulus(states, regimes))
      allocate(out%spectral_radius(regimes))
      allocate(out%stationary(regimes))
      allocate(out%near_unit_root(regimes))
      do regime = 1, regimes
         companion = 0.0_dp
         do lag = 1, lags
            companion(:variables, (lag - 1)*variables + 1:lag*variables) = &
               model%ar(:, :, lag, regime)
         end do
         if (lags > 1) companion(variables + 1:, :states - variables) = &
            identity_matrix(states - variables)
         call general_eigenvalues(companion, roots, info)
         if (info /= 0) then
            out%info = 100 + regime
            return
         end if
         out%eigenvalue(:, regime) = roots
         out%modulus(:, regime) = abs(roots)
         out%spectral_radius(regime) = maxval(out%modulus(:, regime))
         out%stationary(regime) = out%spectral_radius(regime) < 1.0_dp
         out%near_unit_root(regime) = any(out%modulus(:, regime) > &
            1.0_dp - selected_tolerance)
      end do
      out%tolerance = selected_tolerance
   end function gmvarkit_companion_eigenvalues

   pure function gmvarkit_covariance_eigenvalues(model, &
      positive_definite_tolerance, identification_tolerance) result(out)
      !! Diagnose covariance singularity and heteroskedastic identification.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      real(dp), intent(in), optional :: positive_definite_tolerance !! Small covariance eigenvalue threshold; default 0.0002.
      real(dp), intent(in), optional :: identification_tolerance !! Relative generalized-eigenvalue separation; default 0.001.
      type(gmvarkit_covariance_eigen_t) :: out
      real(dp), allocatable :: values(:), vectors(:, :), lower(:, :)
      real(dp), allocatable :: inverse(:, :), normalized(:, :)
      real(dp) :: posdef_tolerance, id_tolerance, separation, denominator
      integer :: variables, regimes, pairs, regime, first, second, pair_index
      integer :: left, right, info

      posdef_tolerance = 0.0002_dp
      if (present(positive_definite_tolerance)) posdef_tolerance = &
         positive_definite_tolerance
      id_tolerance = 0.001_dp
      if (present(identification_tolerance)) id_tolerance = &
         identification_tolerance
      if (.not. valid_model(model) .or. posdef_tolerance <= 0.0_dp .or. &
         id_tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      pairs = regimes*(regimes - 1)/2
      allocate(out%eigenvalue(variables, regimes))
      allocate(out%near_singular(regimes))
      do regime = 1, regimes
         call symmetric_eigen(model%covariance(:, :, regime), values, vectors, &
            info)
         if (info /= 0) then
            out%info = 100 + regime
            return
         end if
         out%eigenvalue(:, regime) = values
         out%near_singular(regime) = any(abs(values) < posdef_tolerance)
      end do
      allocate(out%pair(2, pairs))
      allocate(out%ratio_eigenvalue(variables, pairs))
      allocate(out%minimum_separation(pairs), source=huge(1.0_dp))
      allocate(out%weakly_identified(pairs), source=.false.)
      pair_index = 0
      do first = 1, regimes - 1
         call cholesky_lower(model%covariance(:, :, first), lower, info)
         if (info /= 0) then
            out%info = 200 + first
            return
         end if
         call invert_matrix(lower, inverse, info)
         if (info /= 0) then
            out%info = 300 + first
            return
         end if
         do second = first + 1, regimes
            pair_index = pair_index + 1
            out%pair(:, pair_index) = [first, second]
            normalized = matmul(inverse, matmul( &
               model%covariance(:, :, second), transpose(inverse)))
            normalized = 0.5_dp*(normalized + transpose(normalized))
            call symmetric_eigen(normalized, values, vectors, info)
            if (info /= 0) then
               out%info = 400 + pair_index
               return
            end if
            out%ratio_eigenvalue(:, pair_index) = values
            do left = 1, variables - 1
               do right = left + 1, variables
                  denominator = max(1.0_dp, abs(values(left)), &
                     abs(values(right)))
                  separation = abs(values(left) - values(right))/denominator
                  out%minimum_separation(pair_index) = min( &
                     out%minimum_separation(pair_index), separation)
               end do
            end do
            if (variables > 1) out%weakly_identified(pair_index) = &
               out%minimum_separation(pair_index) <= id_tolerance
         end do
      end do
      out%positive_definite_tolerance = posdef_tolerance
      out%identification_tolerance = id_tolerance
   end function gmvarkit_covariance_eigenvalues

   pure function gmvarkit_location_parameters(model) result(out)
      !! Return equivalent intercept and stationary-mean location parameters.
      type(gmvarkit_model_t), intent(in) :: model !! Valid intercept-parameterized model.
      type(gmvarkit_location_parameters_t) :: out
      real(dp), allocatable :: ar_sum(:, :), inverse(:, :)
      integer :: variables, regimes, regime, info

      if (.not. valid_model(model)) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      out%intercept = model%intercept
      allocate(out%mean(variables, regimes))
      allocate(ar_sum(variables, variables))
      do regime = 1, regimes
         ar_sum = sum(model%ar(:, :, :, regime), dim=3)
         call invert_matrix(identity_matrix(variables) - ar_sum, inverse, info)
         if (info /= 0) then
            out%info = 10 + regime
            return
         end if
         out%mean(:, regime) = matmul(inverse, model%intercept(:, regime))
      end do
   end function gmvarkit_location_parameters

   pure function gmvarkit_model_from_regime_means(template, mean) result(out)
      !! Recalculate model intercepts from supplied stationary regime means.
      type(gmvarkit_model_t), intent(in) :: template !! Model supplying AR and distribution parameters.
      real(dp), intent(in) :: mean(:, :) !! Stationary mean by variable and regime.
      type(gmvarkit_model_t) :: out
      real(dp), allocatable :: ar_sum(:, :)
      integer :: variables, regimes, regime

      if (.not. valid_model(template)) then
         out%info = 1
         return
      end if
      variables = size(template%intercept, 1)
      regimes = size(template%intercept, 2)
      if (any(shape(mean) /= [variables, regimes])) then
         out%info = 2
         return
      end if
      out = template
      allocate(ar_sum(variables, variables))
      do regime = 1, regimes
         ar_sum = sum(template%ar(:, :, :, regime), dim=3)
         out%intercept(:, regime) = matmul( &
            identity_matrix(variables) - ar_sum, mean(:, regime))
      end do
      out%info = 0
   end function gmvarkit_model_from_regime_means

   pure elemental function gmvarkit_information_criteria(log_likelihood, &
      parameter_count, observations) result(out)
      !! Calculate likelihood-based model-selection criteria.
      real(dp), intent(in) :: log_likelihood !! Maximized model log likelihood.
      integer, intent(in) :: parameter_count !! Number of freely estimated parameters.
      integer, intent(in) :: observations !! Effective observations excluding starting values.
      type(gmvarkit_information_criteria_t) :: out

      if (parameter_count < 0 .or. observations < 2) then
         out%info = 1
         return
      end if
      out%aic = -2.0_dp*log_likelihood + 2.0_dp*real(parameter_count, dp)
      out%hqic = -2.0_dp*log_likelihood + 2.0_dp*real(parameter_count, dp)* &
         log(log(real(observations, dp)))
      out%bic = -2.0_dp*log_likelihood + real(parameter_count, dp)* &
         log(real(observations, dp))
   end function gmvarkit_information_criteria

   pure function gmvarkit_pearson_residuals(series, model, standardize) &
      result(out)
      !! Calculate raw or symmetrically standardized conditional residuals.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      logical, intent(in), optional :: standardize !! Apply the inverse symmetric covariance square root.
      type(gmvarkit_pearson_residuals_t) :: out
      type(gmvarkit_evaluation_t) :: evaluation
      real(dp), allocatable :: values(:), vectors(:, :), inverse_sqrt(:, :)
      real(dp), allocatable :: covariance(:, :)
      integer :: lags, usable, variables, time, info
      logical :: use_standardization

      use_standardization = .true.
      if (present(standardize)) use_standardization = standardize
      evaluation = gmvarkit_evaluate(series, model)
      if (evaluation%info /= 0) then
         out%info = 100 + evaluation%info
         return
      end if
      lags = size(model%ar, 3)
      usable = size(series, 1) - lags
      variables = size(series, 2)
      allocate(out%residual(usable, variables))
      out%residual = series(lags + 1:, :) - transpose(evaluation%conditional_mean)
      if (.not. use_standardization) return
      do time = 1, usable
         covariance = 0.5_dp*(evaluation%conditional_covariance(:, :, time) + &
            transpose(evaluation%conditional_covariance(:, :, time)))
         call symmetric_eigen(covariance, values, vectors, info)
         if (info /= 0 .or. any(values <= 100.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(values))))) then
            out%info = 200 + time
            return
         end if
         inverse_sqrt = matmul(vectors, matmul( &
            spread(1.0_dp/sqrt(values), 1, variables)*identity_matrix(variables), &
            transpose(vectors)))
         out%residual(time, :) = matmul(inverse_sqrt, out%residual(time, :))
      end do
      out%standardized = .true.
   end function gmvarkit_pearson_residuals

   pure function gmvarkit_unconditional_moments(model) result(out)
      !! Calculate unconditional mixture moments through lag p.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      type(gmvarkit_unconditional_moments_t) :: out
      type(gmvarkit_regime_moments_t) :: moments
      real(dp), allocatable :: difference(:), between(:, :)
      real(dp) :: denominator
      integer :: variables, lags, regimes, regime, lag, ar_lag
      integer :: row_variable, column_variable

      moments = gmvarkit_regime_moments(model)
      if (moments%info /= 0) then
         out%info = 100 + moments%info
         return
      end if
      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      allocate(out%mean(variables), source=0.0_dp)
      allocate(out%autocovariance(variables, variables, lags + 1), source=0.0_dp)
      allocate(out%autocorrelation(variables, variables, lags + 1), source=0.0_dp)
      allocate(out%regime_autocovariance(variables, variables, lags + 1, &
         regimes), source=0.0_dp)
      allocate(difference(variables), between(variables, variables))
      do regime = 1, regimes
         out%mean = out%mean + model%weight(regime)*moments%mean(:, regime)
         do lag = 0, lags - 1
            out%regime_autocovariance(:, :, lag + 1, regime) = &
               moments%lag_covariance(:variables, &
               lag*variables + 1:(lag + 1)*variables, regime)
         end do
         do ar_lag = 1, lags
            out%regime_autocovariance(:, :, lags + 1, regime) = &
               out%regime_autocovariance(:, :, lags + 1, regime) + &
               matmul(model%ar(:, :, ar_lag, regime), &
               out%regime_autocovariance(:, :, lags + 1 - ar_lag, regime))
         end do
      end do
      do regime = 1, regimes
         difference = moments%mean(:, regime) - out%mean
         between = spread(difference, 2, variables)* &
            spread(difference, 1, variables)
         do lag = 0, lags
            out%autocovariance(:, :, lag + 1) = &
               out%autocovariance(:, :, lag + 1) + model%weight(regime)* &
               (out%regime_autocovariance(:, :, lag + 1, regime) + between)
         end do
      end do
      do lag = 0, lags
         do column_variable = 1, variables
            do row_variable = 1, variables
               denominator = sqrt(max(0.0_dp, &
                  out%autocovariance(row_variable, row_variable, 1)* &
                  out%autocovariance(column_variable, column_variable, 1)))
               if (denominator > tiny(1.0_dp)) then
                  out%autocorrelation(row_variable, column_variable, lag + 1) = &
                     out%autocovariance(row_variable, column_variable, lag + 1)/ &
                     denominator
               end if
            end do
         end do
      end do
   end function gmvarkit_unconditional_moments

   pure function gmvarkit_linear_irf(model, horizon, regime, &
      cumulative_variables, scale_shocks, scale_variables, scale_values, &
      structural) &
      result(out)
      !! Calculate a regime VAR's recursively identified linear responses.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      integer, intent(in) :: horizon !! Largest response horizon, including zero impact.
      integer, intent(in), optional :: regime !! Regime whose linear dynamics are used.
      integer, intent(in), optional :: cumulative_variables(:) !! Variables accumulated over horizons.
      integer, intent(in), optional :: scale_shocks(:) !! Shocks assigned impact-response targets.
      integer, intent(in), optional :: scale_variables(:) !! Variables whose impact responses are targeted.
      real(dp), intent(in), optional :: scale_values(:) !! Nonzero target impact responses.
      type(gmvarkit_structural_t), intent(in), optional :: structural !! Optional heteroskedastic structural identification.
      type(gmvarkit_linear_irf_t) :: out
      real(dp), allocatable :: moving_average(:, :, :), factor(:, :)
      real(dp) :: multiplier
      logical, allocatable :: selected_variable(:), selected_shock(:)
      integer :: variables, lags, regimes, selected_regime
      integer :: step, lag, variable, scaling, shock, info
      logical :: has_scaling

      selected_regime = 1
      if (present(regime)) selected_regime = regime
      if (.not. valid_model(model)) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      if (horizon < 0 .or. selected_regime < 1 .or. &
         selected_regime > regimes) then
         out%info = 2
         return
      end if
      if (present(structural)) then
         if (.not. compatible_structure(structural, model)) then
            out%info = 10
            return
         end if
      end if
      has_scaling = present(scale_shocks) .or. present(scale_variables) .or. &
         present(scale_values)
      if (has_scaling .and. .not. (present(scale_shocks) .and. &
         present(scale_variables) .and. present(scale_values))) then
         out%info = 3
         return
      end if
      if (has_scaling) then
         if (size(scale_shocks) /= size(scale_variables) .or. &
            size(scale_shocks) /= size(scale_values)) then
            out%info = 4
            return
         end if
         if (any(scale_shocks < 1) .or. any(scale_shocks > variables) .or. &
            any(scale_variables < 1) .or. any(scale_variables > variables) .or. &
            any(abs(scale_values) <= tiny(1.0_dp))) then
            out%info = 5
            return
         end if
      end if
      allocate(selected_variable(variables), source=.false.)
      if (present(cumulative_variables)) then
         if (any(cumulative_variables < 1) .or. &
            any(cumulative_variables > variables)) then
            out%info = 6
            return
         end if
         do variable = 1, size(cumulative_variables)
            if (selected_variable(cumulative_variables(variable))) then
               out%info = 7
               return
            end if
            selected_variable(cumulative_variables(variable)) = .true.
         end do
         out%cumulative_variable = cumulative_variables
      else
         allocate(out%cumulative_variable(0))
      end if
      if (has_scaling) then
         allocate(selected_shock(variables), source=.false.)
         do scaling = 1, size(scale_shocks)
            if (selected_shock(scale_shocks(scaling))) then
               out%info = 8
               return
            end if
            selected_shock(scale_shocks(scaling)) = .true.
         end do
      end if
      if (present(structural)) then
         factor = structural%impact*spread( &
            sqrt(structural%relative_variance(:, selected_regime)), 1, variables)
      else
         call cholesky_lower(model%covariance(:, :, selected_regime), factor, info)
         if (info /= 0) then
            out%info = 100 + info
            return
         end if
      end if
      allocate(moving_average(variables, variables, horizon + 1), source=0.0_dp)
      moving_average(:, :, 1) = identity_matrix(variables)
      do step = 1, horizon
         do lag = 1, min(lags, step)
            moving_average(:, :, step + 1) = &
               moving_average(:, :, step + 1) + &
               matmul(model%ar(:, :, lag, selected_regime), &
               moving_average(:, :, step - lag + 1))
         end do
      end do
      allocate(out%response(horizon + 1, variables, variables))
      do step = 0, horizon
         out%response(step + 1, :, :) = &
            matmul(moving_average(:, :, step + 1), factor)
      end do
      do variable = 1, variables
         if (.not. selected_variable(variable)) cycle
         do step = 1, horizon
            out%response(step + 1, variable, :) = &
               out%response(step, variable, :) + &
               out%response(step + 1, variable, :)
         end do
      end do
      if (has_scaling) then
         do scaling = 1, size(scale_shocks)
            shock = scale_shocks(scaling)
            variable = scale_variables(scaling)
            if (abs(out%response(1, variable, shock)) <= tiny(1.0_dp)) then
               out%info = 9
               return
            end if
            multiplier = scale_values(scaling)/out%response(1, variable, shock)
            out%response(:, :, shock) = multiplier*out%response(:, :, shock)
         end do
      end if
      out%impact = factor
      out%regime = selected_regime
   end function gmvarkit_linear_irf

   function gmvarkit_linear_irf_bootstrap(series, fit, horizon, regime, &
      replications, confidence_level, cumulative_variables, scale_shocks, &
      scale_variables, scale_values, structural, max_iterations, tolerance) &
      result(out)
      !! Form fixed-design wild-bootstrap intervals for linear responses.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_fit_t), intent(in) :: fit !! Fitted model used for residuals and starting values.
      integer, intent(in) :: horizon !! Largest response horizon, including zero impact.
      integer, intent(in), optional :: regime !! Regime whose linear dynamics are used.
      integer, intent(in), optional :: replications !! Number of wild-bootstrap fits requested.
      real(dp), intent(in), optional :: confidence_level !! Equal-tail interval coverage in (0,1).
      integer, intent(in), optional :: cumulative_variables(:) !! Variables accumulated over horizons.
      integer, intent(in), optional :: scale_shocks(:) !! Shocks assigned impact-response targets.
      integer, intent(in), optional :: scale_variables(:) !! Variables whose impact responses are targeted.
      real(dp), intent(in), optional :: scale_values(:) !! Nonzero target impact responses.
      type(gmvarkit_structural_t), intent(in), optional :: structural !! Optional baseline structural identification.
      integer, intent(in), optional :: max_iterations !! Maximum iterations for each bootstrap fit.
      real(dp), intent(in), optional :: tolerance !! Gradient tolerance for each bootstrap fit.
      type(gmvarkit_linear_irf_bootstrap_t) :: out
      type(gmvarkit_evaluation_t) :: evaluation
      type(gmvarkit_linear_irf_t) :: response
      type(gmvarkit_fit_t) :: bootstrap_fit
      type(gmvarkit_constraints_t) :: constraints
      type(gmvarkit_structural_t) :: bootstrap_structure
      real(dp), allocatable :: residual(:, :), bootstrap_series(:, :)
      real(dp), allocatable :: draws(:, :, :, :), ordered(:)
      real(dp) :: coverage, selected_tolerance, sign_value
      integer :: count, selected_regime, limit, lags, usable, variables
      integer :: replication, success, time, response_variable, shock, step

      count = 100
      if (present(replications)) count = replications
      coverage = 0.95_dp
      if (present(confidence_level)) coverage = confidence_level
      selected_regime = 1
      if (present(regime)) selected_regime = regime
      limit = 100
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      out%requested_replications = count
      out%confidence_level = coverage
      if (.not. valid_model(fit%model) .or. horizon < 0 .or. count < 2 .or. &
         coverage <= 0.0_dp .or. coverage >= 1.0_dp .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. &
         .not. linear_bootstrap_supported(fit%model)) then
         out%info = 1
         return
      end if
      if (present(structural)) then
         if (.not. compatible_structure(structural, fit%model)) then
            out%info = 2
            return
         end if
      end if
      lags = size(fit%model%ar, 3)
      variables = size(fit%model%intercept, 1)
      usable = size(series, 1) - lags
      if (usable < 2 .or. size(series, 2) /= variables) then
         out%info = 3
         return
      end if
      if (present(structural)) then
         response = gmvarkit_linear_irf(fit%model, horizon, selected_regime, &
            cumulative_variables, scale_shocks, scale_variables, scale_values, &
            structural)
      else
         response = gmvarkit_linear_irf(fit%model, horizon, selected_regime, &
            cumulative_variables, scale_shocks, scale_variables, scale_values)
      end if
      if (response%info /= 0) then
         out%info = 100 + response%info
         return
      end if
      out%point = response%response
      evaluation = gmvarkit_evaluate(series, fit%model)
      if (evaluation%info /= 0) then
         out%info = 200 + evaluation%info
         return
      end if
      residual = series(lags + 1:, :) - transpose(evaluation%conditional_mean)
      allocate(bootstrap_series(size(series, 1), variables))
      allocate(draws(horizon + 1, variables, variables, count))
      constraints%fixed_weight = fit%model%weight
      success = 0
      do replication = 1, count
         bootstrap_series(:lags, :) = series(:lags, :)
         do time = 1, usable
            sign_value = merge(1.0_dp, -1.0_dp, random_uniform() >= 0.5_dp)
            bootstrap_series(lags + time, :) = &
               evaluation%conditional_mean(:, time) + sign_value*residual(time, :)
         end do
         bootstrap_fit = gmvarkit_estimate_constrained(bootstrap_series, &
            fit%model, constraints, max_iterations=limit, &
            tolerance=selected_tolerance)
         if (.not. valid_model(bootstrap_fit%model)) cycle
         if (present(structural)) then
            bootstrap_structure = gmvarkit_identify_structural( &
               bootstrap_fit%model, structural%reference_regime)
            if (bootstrap_structure%info /= 0) cycle
            bootstrap_structure = align_structure(bootstrap_structure, structural)
            if (bootstrap_structure%info /= 0) cycle
            response = gmvarkit_linear_irf(bootstrap_fit%model, horizon, &
               selected_regime, cumulative_variables, scale_shocks, &
               scale_variables, scale_values, bootstrap_structure)
         else
            response = gmvarkit_linear_irf(bootstrap_fit%model, horizon, &
               selected_regime, cumulative_variables, scale_shocks, &
               scale_variables, scale_values)
         end if
         if (response%info /= 0) cycle
         success = success + 1
         draws(:, :, :, success) = response%response
      end do
      out%successful_replications = success
      if (success < 2) then
         out%info = 4
         return
      end if
      allocate(out%draw(horizon + 1, variables, variables, success))
      out%draw = draws(:, :, :, :success)
      allocate(out%lower(horizon + 1, variables, variables))
      allocate(out%upper(horizon + 1, variables, variables))
      allocate(ordered(success))
      do step = 1, horizon + 1
         do shock = 1, variables
            do response_variable = 1, variables
               ordered = sorted(out%draw(step, response_variable, shock, :))
               out%lower(step, response_variable, shock) = &
                  quantile(ordered, 0.5_dp*(1.0_dp - coverage))
               out%upper(step, response_variable, shock) = &
                  quantile(ordered, 0.5_dp*(1.0_dp + coverage))
            end do
         end do
      end do
   end function gmvarkit_linear_irf_bootstrap

   function gmvarkit_girf(model, initial_values, horizon, simulations, &
      shock_size, structural) result(out)
      !! Estimate nonlinear generalized impulse responses by paired simulation.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      real(dp), intent(in) :: initial_values(:, :) !! Initial observations by time and variable.
      integer, intent(in) :: horizon !! Largest response horizon, including zero impact.
      integer, intent(in), optional :: simulations !! Paired paths averaged per structural shock.
      real(dp), intent(in), optional :: shock_size !! Replacement structural innovation size.
      type(gmvarkit_structural_t), intent(in), optional :: structural !! Optional heteroskedastic structural identification.
      type(gmvarkit_girf_t) :: out
      type(gmvarkit_regime_moments_t) :: moments
      real(dp), allocatable :: baseline_history(:, :), shocked_history(:, :)
      real(dp), allocatable :: baseline_lag(:), shocked_lag(:)
      real(dp), allocatable :: history_inverse(:, :, :), history_logdet(:)
      real(dp), allocatable :: factor(:, :, :), lower(:, :), impact(:, :)
      real(dp), allocatable :: impact_inverse(:, :), total_covariance(:, :)
      real(dp), allocatable :: baseline_weight(:), shocked_weight(:)
      real(dp), allocatable :: baseline_mean(:, :), shocked_mean(:, :)
      real(dp), allocatable :: baseline_arch(:), shocked_arch(:)
      real(dp), allocatable :: normal(:), baseline_shock(:), shocked_shock(:)
      real(dp), allocatable :: structural_shock(:), chi_square(:)
      real(dp) :: uniform, degrees, replacement
      integer :: variables, lags, regimes, count, shock_index, path, step, lag
      integer :: regime, baseline_regime, shocked_regime, info

      count = 1000
      if (present(simulations)) count = simulations
      replacement = 1.0_dp
      if (present(shock_size)) replacement = shock_size
      if (.not. valid_model(model)) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      if (horizon < 0 .or. count < 1 .or. &
         size(initial_values, 1) < lags .or. &
         size(initial_values, 2) /= variables) then
         out%info = 2
         return
      end if
      if (present(structural)) then
         if (.not. compatible_structure(structural, model)) then
            out%info = 3
            return
         end if
      end if
      moments = gmvarkit_regime_moments(model)
      if (moments%info /= 0) then
         out%info = 100 + moments%info
         return
      end if
      allocate(baseline_history(lags, variables), shocked_history(lags, variables))
      allocate(baseline_lag(lags*variables), shocked_lag(lags*variables))
      allocate(history_inverse(lags*variables, lags*variables, regimes))
      allocate(history_logdet(regimes), factor(variables, variables, regimes))
      allocate(baseline_weight(regimes), shocked_weight(regimes))
      allocate(baseline_mean(variables, regimes), shocked_mean(variables, regimes))
      allocate(baseline_arch(regimes), shocked_arch(regimes), chi_square(regimes))
      allocate(normal(variables), baseline_shock(variables), shocked_shock(variables))
      allocate(structural_shock(variables), total_covariance(variables, variables))
      do regime = 1, regimes
         call inverse_logdet(moments%lag_covariance(:, :, regime), &
            history_inverse(:, :, regime), history_logdet(regime), info, &
            100.0_dp*epsilon(1.0_dp))
         if (info /= 0) then
            out%info = 200 + regime
            return
         end if
         if (present(structural)) then
            factor(:, :, regime) = structural%impact*spread( &
               sqrt(structural%relative_variance(:, regime)), 1, variables)
         else
            call cholesky_lower(model%covariance(:, :, regime), lower, info)
            if (info /= 0) then
               out%info = 300 + regime
               return
            end if
            factor(:, :, regime) = lower
         end if
      end do
      allocate(out%response(horizon + 1, variables, variables), source=0.0_dp)
      allocate(out%mixing_weight_response(horizon + 1, regimes, variables), &
         source=0.0_dp)
      do shock_index = 1, variables
         do path = 1, count
            baseline_history = initial_values(size(initial_values, 1) - lags + 1:, :)
            shocked_history = baseline_history
            do step = 0, horizon
               do lag = 1, lags
                  baseline_lag((lag - 1)*variables + 1:lag*variables) = &
                     baseline_history(lags + 1 - lag, :)
                  shocked_lag((lag - 1)*variables + 1:lag*variables) = &
                     shocked_history(lags + 1 - lag, :)
               end do
               call mixture_state(model, moments, baseline_lag, history_inverse, &
                  history_logdet, baseline_weight, baseline_mean, baseline_arch)
               call mixture_state(model, moments, shocked_lag, history_inverse, &
                  history_logdet, shocked_weight, shocked_mean, shocked_arch)
               uniform = random_uniform()
               baseline_regime = regime_from_uniform(baseline_weight, uniform)
               shocked_regime = regime_from_uniform(shocked_weight, uniform)
               do regime = 1, variables
                  normal(regime) = random_standard_normal()
               end do
               chi_square = 1.0_dp
               do regime = model%gaussian_regimes + 1, regimes
                  degrees = model%degrees_of_freedom(regime) + &
                     real(variables*lags, dp)
                  chi_square(regime) = 2.0_dp*random_gamma(0.5_dp*degrees)
               end do
               baseline_shock = matmul(factor(:, :, baseline_regime), normal)
               if (baseline_regime > model%gaussian_regimes) then
                  degrees = model%degrees_of_freedom(baseline_regime) + &
                     real(variables*lags, dp)
                  baseline_shock = baseline_shock*sqrt( &
                     baseline_arch(baseline_regime)*(degrees - 2.0_dp)/ &
                     max(chi_square(baseline_regime), tiny(1.0_dp)))
               end if
               shocked_shock = matmul(factor(:, :, shocked_regime), normal)
               if (shocked_regime > model%gaussian_regimes) then
                  degrees = model%degrees_of_freedom(shocked_regime) + &
                     real(variables*lags, dp)
                  shocked_shock = shocked_shock*sqrt( &
                     shocked_arch(shocked_regime)*(degrees - 2.0_dp)/ &
                     max(chi_square(shocked_regime), tiny(1.0_dp)))
               end if
               if (step == 0) then
                  if (present(structural)) then
                     structural_shock = 0.0_dp
                     do regime = 1, regimes
                        structural_shock = structural_shock + &
                           baseline_weight(regime)*baseline_arch(regime)* &
                           structural%relative_variance(:, regime)
                     end do
                     impact = structural%impact*spread( &
                        sqrt(structural_shock), 1, variables)
                  else
                     total_covariance = 0.0_dp
                     do regime = 1, regimes
                        total_covariance = total_covariance + &
                           baseline_weight(regime)*baseline_arch(regime)* &
                           model%covariance(:, :, regime)
                     end do
                     call cholesky_lower(total_covariance, impact, info)
                     if (info /= 0) then
                        out%info = 400
                        return
                     end if
                  end if
                  call invert_matrix(impact, impact_inverse, info)
                  if (info /= 0) then
                     out%info = 401
                     return
                  end if
                  structural_shock = matmul(impact_inverse, baseline_shock)
                  structural_shock(shock_index) = replacement
                  shocked_shock = matmul(impact, structural_shock)
               end if
               out%response(step + 1, :, shock_index) = &
                  out%response(step + 1, :, shock_index) + &
                  shocked_mean(:, shocked_regime) + shocked_shock - &
                  baseline_mean(:, baseline_regime) - baseline_shock
               out%mixing_weight_response(step + 1, :, shock_index) = &
                  out%mixing_weight_response(step + 1, :, shock_index) + &
                  shocked_weight - baseline_weight
               if (lags > 1) then
                  baseline_history(:lags - 1, :) = baseline_history(2:, :)
                  shocked_history(:lags - 1, :) = shocked_history(2:, :)
               end if
               baseline_history(lags, :) = &
                  baseline_mean(:, baseline_regime) + baseline_shock
               shocked_history(lags, :) = &
                  shocked_mean(:, shocked_regime) + shocked_shock
            end do
         end do
      end do
      out%response = out%response/real(count, dp)
      out%mixing_weight_response = out%mixing_weight_response/real(count, dp)
      out%simulations = count
   end function gmvarkit_girf

   function gmvarkit_girf_inference(model, series, horizon, &
      inner_simulations, outer_replications, confidence_level, &
      initial_value_mode, initial_regimes, fixed_initial_values, shock_size, &
      cumulative_variables, scale_shocks, scale_variables, scale_values, &
      scale_type, scale_horizon, structural) result(out)
      !! Estimate GIRF distributions across a collection of initial histories.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      real(dp), intent(in) :: series(:, :) !! Data supplying observed or default fixed histories.
      integer, intent(in) :: horizon !! Largest response horizon, including zero impact.
      integer, intent(in), optional :: inner_simulations !! Paired paths averaged for each history.
      integer, intent(in), optional :: outer_replications !! Random histories drawn for inference.
      real(dp), intent(in), optional :: confidence_level !! Equal-tail interval coverage in (0,1).
      character(len=*), intent(in), optional :: initial_value_mode !! One of random, data, or fixed.
      integer, intent(in), optional :: initial_regimes(:) !! Regimes eligible for random histories.
      real(dp), intent(in), optional :: fixed_initial_values(:, :) !! User-supplied fixed history.
      real(dp), intent(in), optional :: shock_size !! Replacement structural innovation size.
      integer, intent(in), optional :: cumulative_variables(:) !! Variables accumulated over horizons.
      integer, intent(in), optional :: scale_shocks(:) !! Shocks assigned response targets.
      integer, intent(in), optional :: scale_variables(:) !! Variables defining response targets.
      real(dp), intent(in), optional :: scale_values(:) !! Nonzero target response magnitudes.
      character(len=*), intent(in), optional :: scale_type !! Instantaneous or peak scaling rule.
      integer, intent(in), optional :: scale_horizon !! Last horizon considered for peak scaling.
      type(gmvarkit_structural_t), intent(in), optional :: structural !! Optional structural identification.
      type(gmvarkit_girf_inference_t) :: out
      type(gmvarkit_regime_moments_t) :: moments
      type(gmvarkit_girf_t) :: one
      real(dp), allocatable :: history(:, :), state(:), state_mean(:), normal(:)
      real(dp), allocatable :: factor(:, :), ordered(:), eligible_weight(:)
      integer, allocatable :: eligible_regime(:)
      logical, allocatable :: selected_variable(:), selected_shock(:)
      real(dp) :: coverage, replacement, uniform, cumulative, denominator
      real(dp) :: magnitude, degrees, chi_square
      character(len=6) :: mode
      character(len=7) :: scaling_rule
      integer :: inner_count, outer_count, requested_outer, peak_limit
      integer :: variables, lags, regimes, outer, regime, selected, lag
      integer :: variable, shock, scaling, step, peak_step, info
      logical :: has_scaling

      inner_count = 250
      if (present(inner_simulations)) inner_count = inner_simulations
      requested_outer = 250
      if (present(outer_replications)) requested_outer = outer_replications
      coverage = 0.95_dp
      if (present(confidence_level)) coverage = confidence_level
      replacement = 1.0_dp
      if (present(shock_size)) replacement = shock_size
      mode = 'random'
      if (present(initial_value_mode)) mode = trim(initial_value_mode)
      scaling_rule = 'instant'
      if (present(scale_type)) scaling_rule = trim(scale_type)
      peak_limit = horizon
      if (present(scale_horizon)) peak_limit = scale_horizon
      if (.not. valid_model(model) .or. horizon < 0 .or. inner_count < 1 .or. &
         requested_outer < 1 .or. coverage <= 0.0_dp .or. coverage >= 1.0_dp .or. &
         replacement == 0.0_dp .or. peak_limit < 0 .or. peak_limit > horizon .or. &
         (mode /= 'random' .and. mode /= 'data' .and. mode /= 'fixed') .or. &
         (scaling_rule /= 'instant' .and. scaling_rule /= 'peak')) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      if (size(series, 2) /= variables .or. size(series, 1) < lags) then
         out%info = 2
         return
      end if
      if (present(structural)) then
         if (.not. compatible_structure(structural, model)) then
            out%info = 3
            return
         end if
      end if
      has_scaling = present(scale_shocks) .or. present(scale_variables) .or. &
         present(scale_values)
      if (has_scaling .and. .not. (present(scale_shocks) .and. &
         present(scale_variables) .and. present(scale_values))) then
         out%info = 4
         return
      end if
      allocate(selected_variable(variables), selected_shock(variables))
      selected_variable = .false.
      selected_shock = .false.
      if (present(cumulative_variables)) then
         if (any(cumulative_variables < 1) .or. &
            any(cumulative_variables > variables)) then
            out%info = 5
            return
         end if
         do variable = 1, size(cumulative_variables)
            if (selected_variable(cumulative_variables(variable))) then
               out%info = 6
               return
            end if
            selected_variable(cumulative_variables(variable)) = .true.
         end do
      end if
      if (has_scaling) then
         if (size(scale_shocks) /= size(scale_variables) .or. &
            size(scale_shocks) /= size(scale_values) .or. &
            any(scale_shocks < 1) .or. any(scale_shocks > variables) .or. &
            any(scale_variables < 1) .or. any(scale_variables > variables) .or. &
            any(abs(scale_values) <= tiny(1.0_dp))) then
            out%info = 7
            return
         end if
         do scaling = 1, size(scale_shocks)
            if (selected_shock(scale_shocks(scaling))) then
               out%info = 8
               return
            end if
            selected_shock(scale_shocks(scaling)) = .true.
         end do
      end if
      if (mode == 'data') then
         outer_count = size(series, 1) - lags + 1
      else if (mode == 'fixed') then
         outer_count = 1
         if (present(fixed_initial_values)) then
            if (size(fixed_initial_values, 1) < lags .or. &
               size(fixed_initial_values, 2) /= variables) then
               out%info = 9
               return
            end if
         end if
      else
         outer_count = requested_outer
      end if
      if (present(initial_regimes)) then
         if (size(initial_regimes) < 1 .or. any(initial_regimes < 1) .or. &
            any(initial_regimes > regimes)) then
            out%info = 10
            return
         end if
         eligible_regime = initial_regimes
      else
         eligible_regime = [(regime, regime=1, regimes)]
      end if
      allocate(eligible_weight(size(eligible_regime)))
      eligible_weight = model%weight(eligible_regime)
      eligible_weight = eligible_weight/sum(eligible_weight)
      moments = gmvarkit_regime_moments(model)
      if (moments%info /= 0) then
         out%info = 100 + moments%info
         return
      end if
      allocate(history(lags, variables), state(lags*variables))
      allocate(state_mean(lags*variables), normal(lags*variables))
      allocate(out%response_draw(horizon + 1, variables, variables, outer_count))
      allocate(out%mixing_weight_draw(horizon + 1, regimes, variables, outer_count))
      do outer = 1, outer_count
         if (mode == 'data') then
            history = series(outer:outer + lags - 1, :)
         else if (mode == 'fixed') then
            if (present(fixed_initial_values)) then
               history = fixed_initial_values( &
                  size(fixed_initial_values, 1) - lags + 1:, :)
            else
               history = series(size(series, 1) - lags + 1:, :)
            end if
         else
            uniform = random_uniform()
            cumulative = 0.0_dp
            selected = eligible_regime(size(eligible_regime))
            do regime = 1, size(eligible_regime)
               cumulative = cumulative + eligible_weight(regime)
               if (uniform <= cumulative) then
                  selected = eligible_regime(regime)
                  exit
               end if
            end do
            do variable = 1, lags*variables
               normal(variable) = random_standard_normal()
            end do
            call cholesky_lower(moments%lag_covariance(:, :, selected), &
               factor, info)
            if (info /= 0) then
               out%info = 200 + selected
               return
            end if
            do lag = 1, lags
               state_mean((lag - 1)*variables + 1:lag*variables) = &
                  moments%mean(:, selected)
            end do
            state = state_mean + matmul(factor, normal)
            if (selected > model%gaussian_regimes) then
               degrees = model%degrees_of_freedom(selected)
               chi_square = 2.0_dp*random_gamma(0.5_dp*degrees)
               state = state_mean + (state - state_mean)* &
                  sqrt((degrees - 2.0_dp)/max(chi_square, tiny(1.0_dp)))
            end if
            do lag = 1, lags
               history(lags + 1 - lag, :) = &
                  state((lag - 1)*variables + 1:lag*variables)
            end do
         end if
         if (present(structural)) then
            one = gmvarkit_girf(model, history, horizon, inner_count, &
               replacement, structural)
         else
            one = gmvarkit_girf(model, history, horizon, inner_count, replacement)
         end if
         if (one%info /= 0) then
            out%info = 300 + one%info
            return
         end if
         out%response_draw(:, :, :, outer) = one%response
         out%mixing_weight_draw(:, :, :, outer) = one%mixing_weight_response
         do variable = 1, variables
            if (.not. selected_variable(variable)) cycle
            do step = 2, horizon + 1
               out%response_draw(step, variable, :, outer) = &
                  out%response_draw(step - 1, variable, :, outer) + &
                  out%response_draw(step, variable, :, outer)
            end do
         end do
         if (has_scaling) then
            do scaling = 1, size(scale_shocks)
               shock = scale_shocks(scaling)
               variable = scale_variables(scaling)
               magnitude = scale_values(scaling)
               if (scaling_rule == 'instant') then
                  peak_step = 1
               else
                  peak_step = maxloc(abs(out%response_draw( &
                     :peak_limit + 1, variable, shock, outer)), dim=1)
               end if
               denominator = out%response_draw(peak_step, variable, shock, outer)
               if (abs(denominator) <= tiny(1.0_dp)) then
                  out%info = 11
                  return
               end if
               out%response_draw(:, :, shock, outer) = magnitude/denominator* &
                  out%response_draw(:, :, shock, outer)
               out%mixing_weight_draw(:, :, shock, outer) = magnitude/denominator* &
                  out%mixing_weight_draw(:, :, shock, outer)
            end do
         end if
      end do
      allocate(out%point_response(horizon + 1, variables, variables))
      allocate(out%lower_response(horizon + 1, variables, variables))
      allocate(out%upper_response(horizon + 1, variables, variables))
      allocate(out%point_mixing_weight(horizon + 1, regimes, variables))
      allocate(out%lower_mixing_weight(horizon + 1, regimes, variables))
      allocate(out%upper_mixing_weight(horizon + 1, regimes, variables))
      out%point_response = sum(out%response_draw, dim=4)/real(outer_count, dp)
      out%point_mixing_weight = sum(out%mixing_weight_draw, dim=4)/ &
         real(outer_count, dp)
      allocate(ordered(outer_count))
      do shock = 1, variables
         do variable = 1, variables
            do step = 1, horizon + 1
               ordered = sorted(out%response_draw(step, variable, shock, :))
               out%lower_response(step, variable, shock) = &
                  quantile(ordered, 0.5_dp*(1.0_dp - coverage))
               out%upper_response(step, variable, shock) = &
                  quantile(ordered, 0.5_dp*(1.0_dp + coverage))
            end do
         end do
         do regime = 1, regimes
            do step = 1, horizon + 1
               ordered = sorted(out%mixing_weight_draw(step, regime, shock, :))
               out%lower_mixing_weight(step, regime, shock) = &
                  quantile(ordered, 0.5_dp*(1.0_dp - coverage))
               out%upper_mixing_weight(step, regime, shock) = &
                  quantile(ordered, 0.5_dp*(1.0_dp + coverage))
            end do
         end do
      end do
      out%confidence_level = coverage
      out%inner_simulations = inner_count
      out%outer_replications = outer_count
   end function gmvarkit_girf_inference

   pure function gmvarkit_gfevd_inference(girf) result(out)
      !! Average history-specific GFEVDs formed from outer GIRF draws.
      type(gmvarkit_girf_inference_t), intent(in) :: girf !! Outer GIRF response distribution.
      type(gmvarkit_gfevd_inference_t) :: out
      real(dp), allocatable :: cumulative(:, :, :)
      real(dp) :: denominator
      integer :: horizons, variables, regimes, shocks, histories
      integer :: history, horizon, variable, regime, shock

      if (girf%info /= 0 .or. .not. allocated(girf%response_draw) .or. &
         .not. allocated(girf%mixing_weight_draw)) then
         out%info = 1
         return
      end if
      horizons = size(girf%response_draw, 1)
      variables = size(girf%response_draw, 2)
      shocks = size(girf%response_draw, 3)
      histories = size(girf%response_draw, 4)
      regimes = size(girf%mixing_weight_draw, 2)
      allocate(out%individual(horizons, variables, shocks, histories), source=0.0_dp)
      allocate(out%mixing_weight_individual(horizons, regimes, shocks, histories), &
         source=0.0_dp)
      allocate(cumulative(horizons, max(variables, regimes), shocks))
      do history = 1, histories
         cumulative = 0.0_dp
         do shock = 1, shocks
            do variable = 1, variables
               cumulative(1, variable, shock) = &
                  girf%response_draw(1, variable, shock, history)**2
               do horizon = 2, horizons
                  cumulative(horizon, variable, shock) = &
                     cumulative(horizon - 1, variable, shock) + &
                     girf%response_draw(horizon, variable, shock, history)**2
               end do
            end do
         end do
         do horizon = 1, horizons
            do variable = 1, variables
               denominator = sum(cumulative(horizon, variable, :))
               if (denominator > tiny(1.0_dp)) out%individual( &
                  horizon, variable, :, history) = &
                  cumulative(horizon, variable, :)/denominator
            end do
         end do
         cumulative = 0.0_dp
         do shock = 1, shocks
            do regime = 1, regimes
               cumulative(1, regime, shock) = &
                  girf%mixing_weight_draw(1, regime, shock, history)**2
               do horizon = 2, horizons
                  cumulative(horizon, regime, shock) = &
                     cumulative(horizon - 1, regime, shock) + &
                     girf%mixing_weight_draw(horizon, regime, shock, history)**2
               end do
            end do
         end do
         do horizon = 1, horizons
            do regime = 1, regimes
               denominator = sum(cumulative(horizon, regime, :))
               if (denominator > tiny(1.0_dp)) out%mixing_weight_individual( &
                  horizon, regime, :, history) = &
                  cumulative(horizon, regime, :)/denominator
            end do
         end do
      end do
      allocate(out%decomposition(horizons, variables, shocks))
      allocate(out%mixing_weight_decomposition(horizons, regimes, shocks))
      out%decomposition = sum(out%individual, dim=4)/real(histories, dp)
      out%mixing_weight_decomposition = &
         sum(out%mixing_weight_individual, dim=4)/real(histories, dp)
      out%outer_replications = histories
   end function gmvarkit_gfevd_inference

   pure function gmvarkit_gfevd(girf) result(out)
      !! Form normalized cumulative-square decompositions from mean GIRFs.
      type(gmvarkit_girf_t), intent(in) :: girf !! Generalized impulse-response estimates.
      type(gmvarkit_gfevd_t) :: out
      real(dp), allocatable :: cumulative(:, :, :)
      real(dp) :: denominator
      integer :: horizons, variables, shocks, horizon, variable, shock

      if (girf%info /= 0 .or. .not. allocated(girf%response)) then
         out%info = 1
         return
      end if
      horizons = size(girf%response, 1)
      variables = size(girf%response, 2)
      shocks = size(girf%response, 3)
      if (horizons < 1 .or. variables < 1 .or. shocks < 1) then
         out%info = 2
         return
      end if
      allocate(cumulative(horizons, variables, shocks), source=0.0_dp)
      allocate(out%decomposition(horizons, variables, shocks), source=0.0_dp)
      do shock = 1, shocks
         do variable = 1, variables
            cumulative(1, variable, shock) = girf%response(1, variable, shock)**2
            do horizon = 2, horizons
               cumulative(horizon, variable, shock) = &
                  cumulative(horizon - 1, variable, shock) + &
                  girf%response(horizon, variable, shock)**2
            end do
         end do
      end do
      do horizon = 1, horizons
         do variable = 1, variables
            denominator = sum(cumulative(horizon, variable, :))
            if (denominator > tiny(1.0_dp)) then
               out%decomposition(horizon, variable, :) = &
                  cumulative(horizon, variable, :)/denominator
            end if
         end do
      end do
   end function gmvarkit_gfevd

   pure function gmvarkit_score_matrix(series, model, conditional, &
      difference_step) result(out)
      !! Calculate observationwise finite-difference likelihood scores and OPG.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: model !! Full model at the null estimate.
      logical, intent(in), optional :: conditional !! Use conditional likelihood terms.
      real(dp), intent(in), optional :: difference_step !! Relative finite-difference step.
      type(gmvarkit_score_t) :: out
      type(gmvarkit_evaluation_t) :: plus_evaluation, minus_evaluation
      real(dp), allocatable :: shifted(:)
      real(dp) :: step, increment
      integer :: parameter, info, observations
      logical :: use_conditional

      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      step = 1.0e-5_dp
      if (present(difference_step)) step = difference_step
      if (.not. valid_model(model) .or. step <= 0.0_dp .or. &
         size(series, 2) /= size(model%intercept, 1) .or. &
         size(series, 1) <= size(model%ar, 3)) then
         out%info = 1
         return
      end if
      call pack_model(model, out%parameter, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      observations = size(series, 1) - size(model%ar, 3)
      allocate(out%observation(observations, size(out%parameter)))
      allocate(out%total(size(out%parameter)), shifted(size(out%parameter)))
      do parameter = 1, size(out%parameter)
         increment = step*max(1.0_dp, abs(out%parameter(parameter)))
         shifted = out%parameter
         shifted(parameter) = shifted(parameter) + increment
         plus_evaluation = gmvarkit_evaluate(series, &
            unpack_model(shifted, model), conditional=use_conditional)
         shifted(parameter) = out%parameter(parameter) - increment
         minus_evaluation = gmvarkit_evaluate(series, &
            unpack_model(shifted, model), conditional=use_conditional)
         if (plus_evaluation%info /= 0 .or. minus_evaluation%info /= 0) then
            out%info = 3
            return
         end if
         out%observation(:, parameter) = &
            (plus_evaluation%log_likelihood_term - &
            minus_evaluation%log_likelihood_term)/(2.0_dp*increment)
      end do
      out%total = sum(out%observation, dim=1)
      out%opg = matmul(transpose(out%observation), out%observation)
   end function gmvarkit_score_matrix

   pure function gmvarkit_rao_test(score, degrees_of_freedom) result(out)
      !! Evaluate the outer-product-of-gradients Rao score statistic.
      type(gmvarkit_score_t), intent(in) :: score !! Observationwise score information.
      integer, intent(in) :: degrees_of_freedom !! Number of restrictions under the null.
      type(gmvarkit_hypothesis_test_t) :: out
      real(dp), allocatable :: inverse(:, :)
      integer :: info

      if (score%info /= 0 .or. .not. allocated(score%total) .or. &
         .not. allocated(score%opg) .or. degrees_of_freedom < 1) then
         out%info = 1
         return
      end if
      call invert_matrix(score%opg, inverse, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      out%statistic = max(0.0_dp, dot_product(score%total, &
         matmul(inverse, score%total)))
      out%degrees_of_freedom = degrees_of_freedom
      out%p_value = regularized_gamma_q(0.5_dp*real(degrees_of_freedom, dp), &
         0.5_dp*out%statistic)
   end function gmvarkit_rao_test

   pure function gmvarkit_quantile_residuals(series, model) result(out)
      !! Calculate sequential multivariate Gaussian quantile residuals.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: model !! Fitted G-StMVAR model.
      type(gmvarkit_quantile_residuals_t) :: out
      type(gmvarkit_evaluation_t) :: evaluation
      real(dp), allocatable :: inverse(:, :), difference(:), beta(:)
      real(dp), allocatable :: log_density(:), unnormalized(:), regime_cdf(:)
      real(dp) :: log_determinant, quadratic, arch, conditional_mean
      real(dp) :: conditional_variance, schur, degrees, maximum_log
      real(dp) :: probability
      integer :: variables, lags, regimes, observations, time, coordinate
      integer :: regime, info

      if (.not. valid_model(model) .or. &
         size(series, 2) /= size(model%intercept, 1) .or. &
         size(series, 1) <= size(model%ar, 3)) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      observations = size(series, 1) - lags
      evaluation = gmvarkit_evaluate(series, model)
      if (evaluation%info /= 0) then
         out%info = 2
         return
      end if
      allocate(out%residual(observations, variables))
      allocate(out%cdf(observations, variables), inverse(variables, variables))
      allocate(difference(variables), beta(regimes), log_density(regimes))
      allocate(unnormalized(regimes), regime_cdf(regimes))
      do time = 1, observations
         do coordinate = 1, variables
            if (coordinate == 1) then
               beta = evaluation%mixing_weight(time, :)
            else
               do regime = 1, regimes
                  call inverse_logdet( &
                     model%covariance(:coordinate - 1, :coordinate - 1, regime), &
                     inverse(:coordinate - 1, :coordinate - 1), &
                     log_determinant, info, 100.0_dp*epsilon(1.0_dp))
                  if (info /= 0) then
                     out%info = 10 + regime
                     return
                  end if
                  difference(:coordinate - 1) = &
                     series(lags + time, :coordinate - 1) - &
                     evaluation%regime_mean(:coordinate - 1, time, regime)
                  quadratic = dot_product(difference(:coordinate - 1), &
                     matmul(inverse(:coordinate - 1, :coordinate - 1), &
                     difference(:coordinate - 1)))
                  if (regime <= model%gaussian_regimes) then
                     log_density(regime) = gaussian_log_density(coordinate - 1, &
                        quadratic, log_determinant)
                  else
                     arch = evaluation%arch_scalar(time, regime)
                     degrees = model%degrees_of_freedom(regime) + &
                        real(variables*lags, dp)
                     log_density(regime) = student_log_density(coordinate - 1, &
                        quadratic/arch, log_determinant + &
                        real(coordinate - 1, dp)*log(arch), degrees)
                  end if
               end do
               maximum_log = maxval(log(evaluation%mixing_weight(time, :)) + &
                  log_density)
               unnormalized = exp(log(evaluation%mixing_weight(time, :)) + &
                  log_density - maximum_log)
               beta = unnormalized/sum(unnormalized)
            end if
            do regime = 1, regimes
               conditional_mean = evaluation%regime_mean(coordinate, time, regime)
               arch = evaluation%arch_scalar(time, regime)
               if (coordinate == 1) then
                  schur = model%covariance(1, 1, regime)
                  conditional_variance = arch*schur
               else
                  call inverse_logdet( &
                     model%covariance(:coordinate - 1, :coordinate - 1, regime), &
                     inverse(:coordinate - 1, :coordinate - 1), &
                     log_determinant, info, 100.0_dp*epsilon(1.0_dp))
                  difference(:coordinate - 1) = &
                     series(lags + time, :coordinate - 1) - &
                     evaluation%regime_mean(:coordinate - 1, time, regime)
                  conditional_mean = conditional_mean + &
                     dot_product(model%covariance(coordinate, &
                     :coordinate - 1, regime), &
                     matmul(inverse(:coordinate - 1, :coordinate - 1), &
                     difference(:coordinate - 1)))
                  schur = model%covariance(coordinate, coordinate, regime) - &
                     dot_product(model%covariance(coordinate, &
                     :coordinate - 1, regime), &
                     matmul(inverse(:coordinate - 1, :coordinate - 1), &
                     model%covariance(:coordinate - 1, coordinate, regime)))
                  if (regime <= model%gaussian_regimes) then
                     conditional_variance = schur
                  else
                     quadratic = dot_product(difference(:coordinate - 1), &
                        matmul(inverse(:coordinate - 1, :coordinate - 1), &
                        difference(:coordinate - 1)))
                     degrees = model%degrees_of_freedom(regime) + &
                        real(variables*lags, dp)
                     conditional_variance = arch*schur*(degrees + &
                        quadratic/arch)/(degrees + real(coordinate, dp) - 3.0_dp)
                  end if
               end if
               if (regime <= model%gaussian_regimes) then
                  regime_cdf(regime) = 0.5_dp*(1.0_dp + erf( &
                     (series(lags + time, coordinate) - conditional_mean)/ &
                     sqrt(2.0_dp*conditional_variance)))
               else
                  degrees = model%degrees_of_freedom(regime) + &
                     real(variables*lags + coordinate - 1, dp)
                  regime_cdf(regime) = student_covariance_cdf( &
                     series(lags + time, coordinate) - conditional_mean, &
                     conditional_variance, degrees)
               end if
            end do
            probability = sum(beta*regime_cdf)
            probability = max(0.5_dp*epsilon(1.0_dp), &
               min(1.0_dp - 0.5_dp*epsilon(1.0_dp), probability))
            out%cdf(time, coordinate) = probability
            out%residual(time, coordinate) = normal_quantile(probability)
         end do
      end do
   end function gmvarkit_quantile_residuals

   pure function gmvarkit_quantile_residual_tests(residuals, &
      autocorrelation_lags, heteroskedasticity_lags, series, model, &
      difference_step, conditional) result(out)
      !! Test quantile-residual normality, serial correlation, and variance dependence.
      real(dp), intent(in) :: residuals(:, :) !! Quantile residuals by time and variable.
      integer, intent(in) :: autocorrelation_lags(:) !! Maximum lags for correlation tests.
      integer, intent(in) :: heteroskedasticity_lags(:) !! Maximum lags for squared-residual tests.
      real(dp), intent(in), optional :: series(:, :) !! Data used for parameter-corrected covariance.
      type(gmvarkit_model_t), intent(in), optional :: model !! Model used for parameter-corrected covariance.
      real(dp), intent(in), optional :: difference_step !! Relative derivative step for the correction.
      logical, intent(in), optional :: conditional !! Use conditional likelihood scores in the correction.
      type(gmvarkit_residual_tests_t) :: out
      real(dp), allocatable :: values(:, :), centered_square(:, :)
      real(dp) :: step
      integer :: observations, variables, columns, time, variable, position
      integer :: lag, delay, other, test
      logical :: correct_parameters
      logical :: use_conditional

      observations = size(residuals, 1)
      variables = size(residuals, 2)
      correct_parameters = present(series) .and. present(model)
      step = 1.0e-5_dp
      if (present(difference_step)) step = difference_step
      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      if (observations < 2 .or. variables < 1 .or. &
         any(autocorrelation_lags < 1) .or. &
         any(heteroskedasticity_lags < 1) .or. &
         any(autocorrelation_lags >= observations) .or. &
         any(heteroskedasticity_lags >= observations) .or. step <= 0.0_dp .or. &
         (present(series) .neqv. present(model))) then
         out%info = 1
         return
      end if
      allocate(values(observations, 3*variables))
      do variable = 1, variables
         values(:, 3*variable - 2) = residuals(:, variable)**2 - 1.0_dp
         values(:, 3*variable - 1) = residuals(:, variable)**3
         values(:, 3*variable) = residuals(:, variable)**4 - 3.0_dp
      end do
      out%normality = residual_moment_test(values)
      if (correct_parameters) out%normality = &
         corrected_residual_moment_test(series, model, values, 1, 0, step, &
         use_conditional)
      out%autocorrelation_lag = autocorrelation_lags
      out%heteroskedasticity_lag = heteroskedasticity_lags
      allocate(out%autocorrelation(size(autocorrelation_lags)))
      allocate(out%heteroskedasticity(size(heteroskedasticity_lags)))
      allocate(centered_square(observations, variables))
      centered_square = residuals**2 - 1.0_dp
      do test = 1, size(autocorrelation_lags)
         lag = autocorrelation_lags(test)
         columns = lag*variables**2
         deallocate(values)
         allocate(values(observations - lag, columns))
         do time = lag + 1, observations
            position = 0
            do delay = 1, lag
               do other = 1, variables
                  do variable = 1, variables
                     position = position + 1
                     values(time - lag, position) = residuals(time, variable)* &
                        residuals(time - delay, other)
                  end do
               end do
            end do
         end do
         out%autocorrelation(test) = residual_moment_test(values)
         if (correct_parameters) out%autocorrelation(test) = &
            corrected_residual_moment_test(series, model, values, 2, lag, step, &
            use_conditional)
      end do
      do test = 1, size(heteroskedasticity_lags)
         lag = heteroskedasticity_lags(test)
         columns = lag*variables**2
         if (allocated(values)) deallocate(values)
         allocate(values(observations - lag, columns))
         do time = lag + 1, observations
            position = 0
            do delay = 1, lag
               do other = 1, variables
                  do variable = 1, variables
                     position = position + 1
                     values(time - lag, position) = &
                        centered_square(time, variable)* &
                        centered_square(time - delay, other)
                  end do
               end do
            end do
         end do
         out%heteroskedasticity(test) = residual_moment_test(values)
         if (correct_parameters) out%heteroskedasticity(test) = &
            corrected_residual_moment_test(series, model, values, 3, lag, step, &
            use_conditional)
      end do
      out%parameter_corrected = correct_parameters
   end function gmvarkit_quantile_residual_tests

   pure function gmvarkit_inference(series, fit, conditional, difference_step) &
      result(out)
      !! Estimate transformed-scale covariance and standard errors by Hessian inversion.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_fit_t), intent(in) :: fit !! Fitted unrestricted model.
      logical, intent(in), optional :: conditional !! Use the conditional likelihood.
      real(dp), intent(in), optional :: difference_step !! Relative finite-difference step.
      type(gmvarkit_inference_t) :: out
      real(dp), allocatable :: inverse(:, :), shifted(:)
      real(dp) :: step, increment, plus, minus
      integer :: parameter, info
      logical :: use_conditional

      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      step = epsilon(1.0_dp)**0.25_dp
      if (present(difference_step)) step = difference_step
      if (.not. valid_model(fit%model) .or. step <= 0.0_dp .or. &
         size(series, 2) /= size(fit%model%intercept, 1) .or. &
         size(series, 1) <= size(fit%model%ar, 3)) then
         out%info = 1
         return
      end if
      call pack_model(fit%model, out%parameter, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      out%hessian = finite_difference_hessian(objective, out%parameter, step)
      call invert_matrix(out%hessian, inverse, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      out%covariance = 0.5_dp*(inverse + transpose(inverse))
      allocate(out%gradient(size(out%parameter)), shifted(size(out%parameter)))
      allocate(out%standard_error(size(out%parameter)))
      do parameter = 1, size(out%parameter)
         increment = step*max(1.0_dp, abs(out%parameter(parameter)))
         shifted = out%parameter
         shifted(parameter) = shifted(parameter) + increment
         plus = objective(shifted)
         shifted(parameter) = out%parameter(parameter) - increment
         minus = objective(shifted)
         out%gradient(parameter) = (plus - minus)/(2.0_dp*increment)
         out%standard_error(parameter) = &
            sqrt(max(out%covariance(parameter, parameter), 0.0_dp))
      end do

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Return the negative log likelihood for numerical differentiation.
         real(dp), intent(in) :: parameters(:) !! Unconstrained transformed parameters.
         type(gmvarkit_model_t) :: candidate
         type(gmvarkit_evaluation_t) :: evaluation

         candidate = unpack_model(parameters, fit%model)
         evaluation = gmvarkit_evaluate(series, candidate, &
            conditional=use_conditional)
         if (evaluation%info /= 0) then
            value = 1.0e100_dp
         else
            value = -evaluation%log_likelihood
         end if
      end function objective

   end function gmvarkit_inference

   pure function gmvarkit_wald_test(inference, restriction, null_value) &
      result(out)
      !! Test linear restrictions on transformed model parameters.
      type(gmvarkit_inference_t), intent(in) :: inference !! Hessian-based parameter inference.
      real(dp), intent(in) :: restriction(:, :) !! Full-row-rank restriction matrix.
      real(dp), intent(in) :: null_value(:) !! Values imposed under the null hypothesis.
      type(gmvarkit_hypothesis_test_t) :: out
      real(dp), allocatable :: difference(:), middle(:, :), inverse(:, :)
      integer :: restrictions, info

      restrictions = size(restriction, 1)
      if (inference%info /= 0 .or. .not. allocated(inference%parameter) .or. &
         .not. allocated(inference%covariance) .or. restrictions < 1 .or. &
         size(restriction, 2) /= size(inference%parameter) .or. &
         size(null_value) /= restrictions .or. &
         matrix_rank(restriction, 1.0e-10_dp) /= restrictions) then
         out%info = 1
         return
      end if
      difference = matmul(restriction, inference%parameter) - null_value
      middle = matmul(matmul(restriction, inference%covariance), &
         transpose(restriction))
      call invert_matrix(middle, inverse, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      out%statistic = max(0.0_dp, dot_product(difference, &
         matmul(inverse, difference)))
      out%degrees_of_freedom = restrictions
      out%p_value = regularized_gamma_q(0.5_dp*real(restrictions, dp), &
         0.5_dp*out%statistic)
   end function gmvarkit_wald_test

   pure function gmvarkit_likelihood_ratio(unrestricted, restricted) result(out)
      !! Compare nested unrestricted and restricted maximum-likelihood fits.
      type(gmvarkit_fit_t), intent(in) :: unrestricted !! Freely estimated fit.
      type(gmvarkit_fit_t), intent(in) :: restricted !! Nested restricted fit.
      type(gmvarkit_hypothesis_test_t) :: out

      out%degrees_of_freedom = unrestricted%parameter_count - &
         restricted%parameter_count
      if (out%degrees_of_freedom < 1 .or. &
         unrestricted%log_likelihood < restricted%log_likelihood .or. &
         unrestricted%log_likelihood <= -huge(1.0_dp) .or. &
         restricted%log_likelihood <= -huge(1.0_dp)) then
         out%info = 1
         return
      end if
      out%statistic = 2.0_dp*(unrestricted%log_likelihood - &
         restricted%log_likelihood)
      out%p_value = regularized_gamma_q( &
         0.5_dp*real(out%degrees_of_freedom, dp), 0.5_dp*out%statistic)
   end function gmvarkit_likelihood_ratio

   pure function gmvarkit_convert_student_regimes(model, &
      maximum_degrees_of_freedom, series, structural, estimate, conditional, &
      max_iterations, tolerance) result(out)
      !! Replace large-df Student regimes by Gaussian regimes and reorder them.
      type(gmvarkit_model_t), intent(in) :: model !! Source G-StMVAR or StMVAR model.
      real(dp), intent(in), optional :: maximum_degrees_of_freedom !! Strict conversion threshold; default 100.
      real(dp), intent(in), optional :: series(:, :) !! Data used for likelihood evaluation or re-estimation.
      type(gmvarkit_structural_t), intent(in), optional :: structural !! Compatible structural representation.
      logical, intent(in), optional :: estimate !! Re-estimate the converted model; default false.
      logical, intent(in), optional :: conditional !! Use the conditional likelihood when data are supplied.
      integer, intent(in), optional :: max_iterations !! Maximum iterations during optional re-estimation.
      real(dp), intent(in), optional :: tolerance !! Gradient tolerance during optional re-estimation.
      type(gmvarkit_regime_conversion_t) :: out
      type(gmvarkit_evaluation_t) :: evaluation
      type(gmvarkit_structural_fit_t) :: structural_fit
      type(gmvarkit_fit_t) :: reduced_fit
      logical, allocatable :: used(:), gaussian(:)
      real(dp) :: threshold, best_weight
      integer :: regimes, old_gaussian, new_gaussian, old_regime, new_regime
      integer :: selected, group
      logical :: reestimate, use_conditional

      threshold = 100.0_dp
      if (present(maximum_degrees_of_freedom)) threshold = &
         maximum_degrees_of_freedom
      reestimate = .false.
      if (present(estimate)) reestimate = estimate
      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      if (.not. valid_model(model) .or. threshold <= 2.0_dp .or. &
         (reestimate .and. .not. present(series))) then
         out%info = 1
         return
      end if
      if (present(series)) then
         if (size(series, 2) /= size(model%intercept, 1) .or. &
            size(series, 1) <= size(model%ar, 3)) then
            out%info = 2
            return
         end if
      end if
      if (present(structural)) then
         if (.not. compatible_structure(structural, model)) then
            out%info = 3
            return
         end if
      end if
      regimes = size(model%intercept, 2)
      old_gaussian = model%gaussian_regimes
      allocate(out%converted(regimes), source=.false.)
      do old_regime = old_gaussian + 1, regimes
         out%converted(old_regime) = &
            model%degrees_of_freedom(old_regime) > threshold
      end do
      allocate(gaussian(regimes), source=.false.)
      if (old_gaussian > 0) gaussian(:old_gaussian) = .true.
      gaussian = gaussian .or. out%converted
      new_gaussian = count(gaussian)
      allocate(out%new_to_old(regimes), out%old_to_new(regimes))
      if (any(out%converted)) then
         allocate(used(regimes), source=.false.)
         new_regime = 0
         do group = 1, 2
            do while (count(used .and. &
               (gaussian .eqv. (group == 1))) < &
               count(gaussian .eqv. (group == 1)))
               selected = 0
               best_weight = -huge(1.0_dp)
               do old_regime = 1, regimes
                  if (used(old_regime) .or. &
                     (gaussian(old_regime) .neqv. (group == 1))) cycle
                  if (model%weight(old_regime) > best_weight) then
                     selected = old_regime
                     best_weight = model%weight(old_regime)
                  end if
               end do
               if (selected == 0) exit
               new_regime = new_regime + 1
               out%new_to_old(new_regime) = selected
               out%old_to_new(selected) = new_regime
               used(selected) = .true.
            end do
         end do
      else
         out%new_to_old = [(old_regime, old_regime=1, regimes)]
         out%old_to_new = out%new_to_old
      end if
      out%fit%model = model
      out%fit%model%intercept = model%intercept(:, out%new_to_old)
      out%fit%model%ar = model%ar(:, :, :, out%new_to_old)
      out%fit%model%covariance = model%covariance(:, :, out%new_to_old)
      out%fit%model%weight = model%weight(out%new_to_old)
      out%fit%model%degrees_of_freedom = &
         model%degrees_of_freedom(out%new_to_old)
      out%fit%model%gaussian_regimes = new_gaussian
      if (new_gaussian > 0) &
         out%fit%model%degrees_of_freedom(:new_gaussian) = 0.0_dp
      out%fit%model%info = 0
      if (present(structural)) then
         out%structural = structural
         out%structural%relative_variance = &
            structural%relative_variance(:, out%new_to_old)
         out%structural%covariance = structural%covariance(:, :, out%new_to_old)
         out%structural%reference_regime = &
            out%old_to_new(structural%reference_regime)
         out%has_structural = .true.
      end if
      if (present(series)) then
         evaluation = gmvarkit_evaluate(series, out%fit%model, &
            conditional=use_conditional)
         if (evaluation%info /= 0) then
            out%info = 100 + evaluation%info
            return
         end if
         out%fit%log_likelihood = evaluation%log_likelihood
      end if
      if (.not. reestimate) return
      out%reestimated = .true.
      if (present(structural)) then
         structural_fit = gmvarkit_estimate_structural(series, out%fit%model, &
            out%structural, conditional=use_conditional, &
            max_iterations=max_iterations, tolerance=tolerance, &
            calculate_inference=.false.)
         out%fit = structural_fit%fit
         out%structural = structural_fit%structural
         out%info = structural_fit%info
      else
         reduced_fit = gmvarkit_estimate(series, out%fit%model, &
            conditional=use_conditional, max_iterations=max_iterations, &
            tolerance=tolerance)
         out%fit = reduced_fit
         out%info = reduced_fit%info
      end if
   end function gmvarkit_convert_student_regimes

   pure function gmvarkit_profile_likelihood_reduced(series, model, parameter, &
      scale, points, conditional) result(out)
      !! Evaluate reduced-form likelihood slices on transformed coordinates.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: model !! Model at the center of every profile.
      integer, intent(in), optional :: parameter(:) !! One-based transformed coordinates to profile.
      real(dp), intent(in), optional :: scale !! Relative half-width with a unit magnitude floor.
      integer, intent(in), optional :: points !! Number of equally spaced values per profile.
      logical, intent(in), optional :: conditional !! Evaluate the conditional likelihood.
      type(gmvarkit_likelihood_profile_t) :: out
      type(gmvarkit_model_t) :: candidate
      type(gmvarkit_evaluation_t) :: evaluation
      real(dp), allocatable :: parameters(:), trial(:)
      real(dp) :: selected_scale
      integer :: selected_points, info, profile, grid
      logical :: use_conditional

      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      selected_scale = 0.02_dp
      if (present(scale)) selected_scale = scale
      selected_points = 101
      if (present(points)) selected_points = points
      if (.not. valid_model(model) .or. &
         size(series, 2) /= size(model%intercept, 1) .or. &
         size(series, 1) <= size(model%ar, 3)) then
         out%info = 1
         return
      end if
      call pack_model(model, parameters, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      if (present(parameter)) then
         call initialize_likelihood_profile(parameters, selected_scale, &
            selected_points, .false., out, parameter)
      else
         call initialize_likelihood_profile(parameters, selected_scale, &
            selected_points, .false., out)
      end if
      if (out%info /= 0) return
      allocate(trial(size(parameters)))
      do profile = 1, size(out%parameter)
         do grid = 1, selected_points
            trial = parameters
            trial(out%parameter(profile)) = out%value(grid, profile)
            candidate = unpack_model(trial, model)
            evaluation = gmvarkit_evaluate(series, candidate, &
               conditional=use_conditional)
            if (evaluation%info == 0) then
               out%log_likelihood(grid, profile) = evaluation%log_likelihood
               out%valid(grid, profile) = .true.
            end if
         end do
      end do
   end function gmvarkit_profile_likelihood_reduced

   pure function gmvarkit_profile_likelihood_structural(series, model, &
      structural, constraints, parameter, scale, points, conditional) result(out)
      !! Evaluate structural likelihood slices on transformed coordinates.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: model !! Reduced-form model at the profile center.
      type(gmvarkit_structural_t), intent(in) :: structural !! Structural parameters at the profile center.
      type(gmvarkit_structural_constraints_t), intent(in), optional :: &
         constraints !! Structural restrictions defining free coordinates.
      integer, intent(in), optional :: parameter(:) !! One-based transformed coordinates to profile.
      real(dp), intent(in), optional :: scale !! Relative half-width with a unit magnitude floor.
      integer, intent(in), optional :: points !! Number of equally spaced values per profile.
      logical, intent(in), optional :: conditional !! Evaluate the conditional likelihood.
      type(gmvarkit_likelihood_profile_t) :: out
      type(gmvarkit_structural_constraints_t) :: restrictions
      type(gmvarkit_structural_t) :: candidate_structural
      type(gmvarkit_model_t) :: candidate_model
      type(gmvarkit_evaluation_t) :: evaluation
      real(dp), allocatable :: parameters(:), trial(:)
      real(dp) :: selected_scale
      integer :: selected_points, info, profile, grid
      logical :: use_conditional

      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      selected_scale = 0.02_dp
      if (present(scale)) selected_scale = scale
      selected_points = 101
      if (present(points)) selected_points = points
      if (present(constraints)) restrictions = constraints
      if (.not. compatible_structure(structural, model) .or. &
         .not. valid_structural_constraints(model, structural, restrictions) .or. &
         size(series, 2) /= size(model%intercept, 1) .or. &
         size(series, 1) <= size(model%ar, 3)) then
         out%info = 1
         return
      end if
      call pack_structural_model(model, structural, restrictions, parameters, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      if (present(parameter)) then
         call initialize_likelihood_profile(parameters, selected_scale, &
            selected_points, .true., out, parameter)
      else
         call initialize_likelihood_profile(parameters, selected_scale, &
            selected_points, .true., out)
      end if
      if (out%info /= 0) return
      allocate(trial(size(parameters)))
      do profile = 1, size(out%parameter)
         do grid = 1, selected_points
            trial = parameters
            trial(out%parameter(profile)) = out%value(grid, profile)
            call unpack_structural_model(trial, model, structural, restrictions, &
               candidate_model, candidate_structural, info)
            if (info /= 0) cycle
            evaluation = gmvarkit_evaluate(series, candidate_model, &
               conditional=use_conditional)
            if (evaluation%info == 0) then
               out%log_likelihood(grid, profile) = evaluation%log_likelihood
               out%valid(grid, profile) = .true.
            end if
         end do
      end do
   end function gmvarkit_profile_likelihood_structural

   pure function gmvarkit_estimate_structural(series, initial_model, &
      initial_structural, constraints, conditional, max_iterations, tolerance, &
      calculate_inference) result(out)
      !! Estimate W and lambda parameters directly by structural likelihood.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: initial_model !! Valid reduced-form starting model.
      type(gmvarkit_structural_t), intent(in) :: initial_structural !! Compatible structural starting parameters.
      type(gmvarkit_structural_constraints_t), intent(in), optional :: constraints !! Fixed, sign, or lambda restrictions.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      logical, intent(in), optional :: calculate_inference !! Calculate a numerical Hessian covariance.
      type(gmvarkit_structural_fit_t) :: out
      type(gmvarkit_structural_constraints_t) :: restrictions
      type(optimization_result_t) :: optimized
      type(gmvarkit_evaluation_t) :: evaluation
      real(dp), allocatable :: initial(:), inverse(:, :), shifted(:)
      real(dp) :: selected_tolerance, step, increment, plus, minus
      integer :: limit, info, parameter
      logical :: use_conditional, infer

      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      infer = .true.
      if (present(calculate_inference)) infer = calculate_inference
      if (present(constraints)) restrictions = constraints
      if (.not. compatible_structure(initial_structural, initial_model) .or. &
         .not. valid_structural_constraints(initial_model, initial_structural, &
         restrictions) .or. size(series, 2) /= size(initial_model%intercept, 1) .or. &
         size(series, 1) <= size(initial_model%ar, 3) .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      call pack_structural_model(initial_model, initial_structural, restrictions, &
         initial, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      optimized = bfgs_minimize_fd(objective, initial, max_iterations=limit, &
         gradient_tolerance=selected_tolerance)
      call unpack_structural_model(optimized%parameters, initial_model, &
         initial_structural, restrictions, out%fit%model, out%structural, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      evaluation = gmvarkit_evaluate(series, out%fit%model, &
         conditional=use_conditional)
      if (evaluation%info /= 0) then
         out%info = 100 + evaluation%info
         return
      end if
      out%fit%log_likelihood = evaluation%log_likelihood
      out%fit%iterations = optimized%iterations
      out%fit%parameter_count = size(initial)
      out%fit%converged = optimized%converged
      out%fit%info = optimized%info
      out%info = optimized%info
      if (.not. infer) return
      out%inference%parameter = optimized%parameters
      step = epsilon(1.0_dp)**0.25_dp
      out%inference%hessian = finite_difference_hessian(objective, &
         optimized%parameters, step)
      call invert_matrix(out%inference%hessian, inverse, info)
      if (info /= 0) then
         out%inference%info = 3
         return
      end if
      out%inference%covariance = 0.5_dp*(inverse + transpose(inverse))
      allocate(out%inference%gradient(size(initial)))
      allocate(out%inference%standard_error(size(initial)), shifted(size(initial)))
      do parameter = 1, size(initial)
         increment = step*max(1.0_dp, abs(optimized%parameters(parameter)))
         shifted = optimized%parameters
         shifted(parameter) = shifted(parameter) + increment
         plus = objective(shifted)
         shifted(parameter) = optimized%parameters(parameter) - increment
         minus = objective(shifted)
         out%inference%gradient(parameter) = (plus - minus)/(2.0_dp*increment)
         out%inference%standard_error(parameter) = sqrt(max(0.0_dp, &
            out%inference%covariance(parameter, parameter)))
      end do

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Return a penalized structural negative log likelihood.
         real(dp), intent(in) :: parameters(:) !! Free transformed structural coordinates.
         type(gmvarkit_model_t) :: candidate_model
         type(gmvarkit_structural_t) :: candidate_structural
         type(gmvarkit_evaluation_t) :: candidate_evaluation
         integer :: status

         call unpack_structural_model(parameters, initial_model, &
            initial_structural, restrictions, candidate_model, &
            candidate_structural, status)
         if (status /= 0) then
            value = 1.0e100_dp
            return
         end if
         candidate_evaluation = gmvarkit_evaluate(series, candidate_model, &
            conditional=use_conditional)
         if (candidate_evaluation%info /= 0) then
            value = 1.0e100_dp
         else
            value = -candidate_evaluation%log_likelihood
         end if
      end function objective

   end function gmvarkit_estimate_structural

   pure function gmvarkit_structural_multistart_estimate(series, initial_model, &
      initial_structural, constraints, conditional, max_iterations, tolerance, &
      likelihood_tolerance, parameter_tolerance, calculate_inference) result(out)
      !! Fit and deduplicate multiple direct structural starting values.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: initial_model(:) !! Reduced-form starting models.
      type(gmvarkit_structural_t), intent(in) :: initial_structural(:) !! Structural starting values paired with the models.
      type(gmvarkit_structural_constraints_t), intent(in), optional :: constraints !! Restrictions applied to every start.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations for each start.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      real(dp), intent(in), optional :: likelihood_tolerance !! Relative log-likelihood duplicate tolerance.
      real(dp), intent(in), optional :: parameter_tolerance !! Relative reduced-parameter duplicate tolerance.
      logical, intent(in), optional :: calculate_inference !! Calculate Hessian inference for the best fit only.
      type(gmvarkit_structural_multistart_fit_t) :: out
      type(gmvarkit_structural_constraints_t) :: restrictions
      real(dp), allocatable :: candidate_parameter(:), previous_parameter(:)
      real(dp) :: log_tolerance, par_tolerance, likelihood_scale
      real(dp) :: parameter_scale, parameter_distance
      integer :: starts, start, rank, previous_rank, candidate, previous
      integer :: selected, swap_index, info
      logical :: use_conditional, infer

      starts = size(initial_model)
      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      log_tolerance = 1.0e-6_dp
      if (present(likelihood_tolerance)) log_tolerance = likelihood_tolerance
      par_tolerance = 1.0e-4_dp
      if (present(parameter_tolerance)) par_tolerance = parameter_tolerance
      infer = .true.
      if (present(calculate_inference)) infer = calculate_inference
      if (present(constraints)) restrictions = constraints
      if (starts < 1 .or. size(initial_structural) /= starts .or. &
         log_tolerance < 0.0_dp .or. par_tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      do start = 1, starts
         if (.not. valid_model(initial_model(start))) then
            out%info = 2
            return
         end if
         if (.not. compatible_structure(initial_structural(start), &
            initial_model(start))) then
            out%info = 2
            return
         end if
         if (.not. same_model_specification(initial_model(1), &
            initial_model(start))) then
            out%info = 3
            return
         end if
         if (.not. valid_structural_constraints(initial_model(start), &
            initial_structural(start), restrictions)) then
            out%info = 4
            return
         end if
      end do
      if (size(series, 2) /= size(initial_model(1)%intercept, 1) .or. &
         size(series, 1) <= size(initial_model(1)%ar, 3)) then
         out%info = 5
         return
      end if
      allocate(out%fit(starts), out%order(starts))
      allocate(out%successful(starts), source=.false.)
      allocate(out%distinct(starts), source=.false.)
      allocate(out%duplicate_of(starts), source=0)
      do start = 1, starts
         out%fit(start) = gmvarkit_estimate_structural(series, &
            initial_model(start), initial_structural(start), restrictions, &
            conditional=use_conditional, max_iterations=max_iterations, &
            tolerance=tolerance, calculate_inference=.false.)
         out%successful(start) = out%fit(start)%fit%log_likelihood > &
            -0.5_dp*huge(1.0_dp)
         out%order(start) = start
      end do
      out%successful_count = count(out%successful)
      do rank = 1, starts - 1
         selected = rank
         do candidate = rank + 1, starts
            if (out%fit(out%order(candidate))%fit%log_likelihood > &
               out%fit(out%order(selected))%fit%log_likelihood) &
               selected = candidate
         end do
         if (selected /= rank) then
            swap_index = out%order(rank)
            out%order(rank) = out%order(selected)
            out%order(selected) = swap_index
         end if
      end do
      do rank = 1, starts
         candidate = out%order(rank)
         if (.not. out%successful(candidate)) cycle
         call canonical_model_parameters(out%fit(candidate)%fit%model, &
            candidate_parameter, info)
         if (info /= 0) then
            out%info = 6
            return
         end if
         out%distinct(candidate) = .true.
         do previous_rank = 1, rank - 1
            previous = out%order(previous_rank)
            if (.not. out%distinct(previous)) cycle
            likelihood_scale = max(1.0_dp, &
               abs(out%fit(candidate)%fit%log_likelihood), &
               abs(out%fit(previous)%fit%log_likelihood))
            if (abs(out%fit(candidate)%fit%log_likelihood - &
               out%fit(previous)%fit%log_likelihood) > &
               log_tolerance*likelihood_scale) cycle
            call canonical_model_parameters(out%fit(previous)%fit%model, &
               previous_parameter, info)
            if (info /= 0) then
               out%info = 7
               return
            end if
            parameter_scale = max(1.0_dp, maxval(abs(candidate_parameter)), &
               maxval(abs(previous_parameter)))
            parameter_distance = maxval(abs(candidate_parameter - &
               previous_parameter))/parameter_scale
            if (parameter_distance <= par_tolerance) then
               out%distinct(candidate) = .false.
               out%duplicate_of(candidate) = previous
               exit
            end if
         end do
         if (out%distinct(candidate)) then
            out%distinct_count = out%distinct_count + 1
            if (out%best_index == 0) out%best_index = candidate
         end if
      end do
      if (out%successful_count == 0) then
         out%info = 8
         return
      end if
      if (infer) out%fit(out%best_index) = gmvarkit_estimate_structural(series, &
         initial_model(out%best_index), initial_structural(out%best_index), &
         restrictions, conditional=use_conditional, &
         max_iterations=max_iterations, tolerance=tolerance, &
         calculate_inference=.true.)
   end function gmvarkit_structural_multistart_estimate

   pure function gmvarkit_estimate_constrained(series, initial_model, &
      constraints, conditional, max_iterations, tolerance) result(out)
      !! Estimate a model subject to linear AR, shared-mean, and weight constraints.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: initial_model !! Valid starting model.
      type(gmvarkit_constraints_t), intent(in) :: constraints !! Parameter restrictions.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(gmvarkit_fit_t) :: out
      type(optimization_result_t) :: optimized
      type(gmvarkit_evaluation_t) :: evaluation
      real(dp), allocatable :: initial(:)
      real(dp) :: selected_tolerance
      integer :: limit, info
      logical :: use_conditional

      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (.not. valid_model(initial_model) .or. &
         .not. valid_constraints(initial_model, constraints) .or. &
         limit < 1 .or. selected_tolerance <= 0.0_dp .or. &
         size(series, 2) /= size(initial_model%intercept, 1) .or. &
         size(series, 1) <= size(initial_model%ar, 3)) then
         out%info = 1
         return
      end if
      call pack_constrained_model(initial_model, constraints, initial, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      optimized = bfgs_minimize_fd(objective, initial, &
         max_iterations=limit, gradient_tolerance=selected_tolerance)
      out%model = unpack_constrained_model(optimized%parameters, initial_model, &
         constraints)
      evaluation = gmvarkit_evaluate(series, out%model, &
         conditional=use_conditional)
      if (evaluation%info /= 0) then
         out%info = 3
         return
      end if
      out%log_likelihood = evaluation%log_likelihood
      out%iterations = optimized%iterations
      out%parameter_count = size(initial)
      out%converged = optimized%converged
      out%info = optimized%info

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Return the constrained negative log likelihood for BFGS.
         real(dp), intent(in) :: parameters(:) !! Free transformed parameters.
         type(gmvarkit_model_t) :: candidate
         type(gmvarkit_evaluation_t) :: candidate_evaluation

         candidate = unpack_constrained_model(parameters, initial_model, constraints)
         candidate_evaluation = gmvarkit_evaluate(series, candidate, &
            conditional=use_conditional)
         if (candidate_evaluation%info /= 0) then
            value = 1.0e100_dp
         else
            value = -candidate_evaluation%log_likelihood
         end if
      end function objective

   end function gmvarkit_estimate_constrained

   function gmvarkit_genetic_estimate(series, initial_model, conditional, &
      population_size, generations, mutation_scale, local_iterations) result(out)
      !! Search globally with an elitist genetic algorithm and refine by BFGS.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: initial_model !! Center of the initial population.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: population_size !! Number of candidates per generation.
      integer, intent(in), optional :: generations !! Number of genetic generations.
      real(dp), intent(in), optional :: mutation_scale !! Initial Gaussian mutation scale.
      integer, intent(in), optional :: local_iterations !! Final BFGS iteration limit.
      type(gmvarkit_genetic_fit_t) :: out
      real(dp), allocatable :: center(:), population(:, :), next_population(:, :)
      real(dp), allocatable :: objective_value(:)
      integer, allocatable :: order(:)
      type(gmvarkit_model_t) :: best_model
      type(gmvarkit_evaluation_t) :: evaluation
      real(dp) :: scale, generation_scale
      integer :: count, generation_count, refinement, elite_count
      integer :: parameter_count, candidate, parameter, parent1, parent2
      integer :: generation, info
      logical :: use_conditional

      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      count = 40
      if (present(population_size)) count = population_size
      generation_count = 50
      if (present(generations)) generation_count = generations
      scale = 0.25_dp
      if (present(mutation_scale)) scale = mutation_scale
      refinement = 100
      if (present(local_iterations)) refinement = local_iterations
      if (.not. valid_model(initial_model) .or. count < 4 .or. &
         generation_count < 1 .or. scale <= 0.0_dp .or. refinement < 1 .or. &
         size(series, 2) /= size(initial_model%intercept, 1) .or. &
         size(series, 1) <= size(initial_model%ar, 3)) then
         out%info = 1
         return
      end if
      call pack_model(initial_model, center, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      parameter_count = size(center)
      elite_count = max(2, count/5)
      allocate(population(parameter_count, count))
      allocate(next_population(parameter_count, count), objective_value(count))
      allocate(order(count), out%best_objective(generation_count))
      population(:, 1) = center
      do candidate = 2, count
         do parameter = 1, parameter_count
            population(parameter, candidate) = center(parameter) + &
               scale*max(1.0_dp, abs(center(parameter)))*random_standard_normal()
         end do
      end do
      do generation = 1, generation_count
         do candidate = 1, count
            objective_value(candidate) = genetic_objective(population(:, candidate))
         end do
         out%evaluations = out%evaluations + count
         call order_values(objective_value, order)
         out%best_objective(generation) = objective_value(order(1))
         next_population(:, :elite_count) = population(:, order(:elite_count))
         generation_scale = scale*sqrt(1.0_dp - &
            0.9_dp*real(generation - 1, dp)/real(max(1, generation_count - 1), dp))
         do candidate = elite_count + 1, count
            parent1 = 1 + int(random_uniform()*real(elite_count, dp))
            parent2 = 1 + int(random_uniform()*real(elite_count, dp))
            parent1 = min(parent1, elite_count)
            parent2 = min(parent2, elite_count)
            do parameter = 1, parameter_count
               if (random_uniform() < 0.5_dp) then
                  next_population(parameter, candidate) = &
                     population(parameter, order(parent1))
               else
                  next_population(parameter, candidate) = &
                     population(parameter, order(parent2))
               end if
               next_population(parameter, candidate) = &
                  next_population(parameter, candidate) + generation_scale* &
                  max(1.0_dp, abs(next_population(parameter, candidate)))* &
                  random_standard_normal()
            end do
         end do
         population = next_population
      end do
      do candidate = 1, count
         objective_value(candidate) = genetic_objective(population(:, candidate))
      end do
      out%evaluations = out%evaluations + count
      call order_values(objective_value, order)
      if (objective_value(order(1)) < out%best_objective(generation_count)) then
         out%best_objective(generation_count) = objective_value(order(1))
      end if
      best_model = unpack_model(population(:, order(1)), initial_model)
      evaluation = gmvarkit_evaluate(series, best_model, conditional=use_conditional)
      if (evaluation%info /= 0) then
         out%info = 3
         return
      end if
      out%fit = gmvarkit_estimate(series, best_model, &
         conditional=use_conditional, max_iterations=refinement)
      out%population_size = count
      out%generations = generation_count

   contains

      pure real(dp) function genetic_objective(parameters) result(value)
         !! Return a finite penalized objective for one genetic candidate.
         real(dp), intent(in) :: parameters(:) !! Unconstrained model parameters.
         type(gmvarkit_model_t) :: candidate_model
         type(gmvarkit_evaluation_t) :: candidate_evaluation

         candidate_model = unpack_model(parameters, initial_model)
         candidate_evaluation = gmvarkit_evaluate(series, candidate_model, &
            conditional=use_conditional)
         if (candidate_evaluation%info /= 0) then
            value = 1.0e100_dp
         else
            value = -candidate_evaluation%log_likelihood
         end if
      end function genetic_objective

   end function gmvarkit_genetic_estimate

   pure function gmvarkit_estimate(series, initial_model, conditional, &
      max_iterations, tolerance) result(out)
      !! Refine an unconstrained reduced-form G-StMVAR model by BFGS.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: initial_model !! Valid starting model.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(gmvarkit_fit_t) :: out
      type(optimization_result_t) :: optimized
      type(gmvarkit_evaluation_t) :: evaluation
      real(dp), allocatable :: initial(:)
      real(dp) :: selected_tolerance
      integer :: limit, info
      logical :: use_conditional

      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (.not. valid_model(initial_model) .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. &
         size(series, 2) /= size(initial_model%intercept, 1) .or. &
         size(series, 1) <= size(initial_model%ar, 3)) then
         out%info = 1
         return
      end if
      call pack_model(initial_model, initial, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      optimized = bfgs_minimize_fd(objective, initial, &
         max_iterations=limit, gradient_tolerance=selected_tolerance)
      out%model = unpack_model(optimized%parameters, initial_model)
      evaluation = gmvarkit_evaluate(series, out%model, &
         conditional=use_conditional)
      if (evaluation%info /= 0) then
         out%info = 3
         return
      end if
      out%log_likelihood = evaluation%log_likelihood
      out%iterations = optimized%iterations
      out%parameter_count = size(initial)
      out%converged = optimized%converged
      out%info = optimized%info

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Return the penalized negative log likelihood for BFGS.
         real(dp), intent(in) :: parameters(:) !! Unconstrained model parameters.
         type(gmvarkit_model_t) :: candidate
         type(gmvarkit_evaluation_t) :: candidate_evaluation

         candidate = unpack_model(parameters, initial_model)
         candidate_evaluation = gmvarkit_evaluate(series, candidate, &
            conditional=use_conditional)
         if (candidate_evaluation%info /= 0) then
            value = 1.0e100_dp
         else
            value = -candidate_evaluation%log_likelihood
         end if
      end function objective

   end function gmvarkit_estimate

   pure function gmvarkit_multistart_estimate(series, initial_model, &
      conditional, max_iterations, tolerance, likelihood_tolerance, &
      parameter_tolerance) result(out)
      !! Fit multiple starts and retain ranked, deduplicated local solutions.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: initial_model(:) !! Starting models with a common specification.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations for each start.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      real(dp), intent(in), optional :: likelihood_tolerance !! Relative log-likelihood duplicate tolerance.
      real(dp), intent(in), optional :: parameter_tolerance !! Relative canonical-parameter duplicate tolerance.
      type(gmvarkit_multistart_fit_t) :: out
      real(dp), allocatable :: candidate_parameter(:), previous_parameter(:)
      real(dp) :: log_tolerance, par_tolerance, likelihood_scale
      real(dp) :: parameter_scale, parameter_distance
      integer :: starts, start, rank, previous_rank, candidate, previous
      integer :: selected, swap_index, info
      logical :: use_conditional

      starts = size(initial_model)
      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      log_tolerance = 1.0e-6_dp
      if (present(likelihood_tolerance)) log_tolerance = likelihood_tolerance
      par_tolerance = 1.0e-4_dp
      if (present(parameter_tolerance)) par_tolerance = parameter_tolerance
      if (starts < 1 .or. log_tolerance < 0.0_dp .or. &
         par_tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      do start = 1, starts
         if (.not. valid_model(initial_model(start))) then
            out%info = 2
            return
         end if
         if (.not. same_model_specification(initial_model(1), &
            initial_model(start))) then
            out%info = 3
            return
         end if
      end do
      if (size(series, 2) /= size(initial_model(1)%intercept, 1) .or. &
         size(series, 1) <= size(initial_model(1)%ar, 3)) then
         out%info = 4
         return
      end if
      allocate(out%fit(starts), out%order(starts))
      allocate(out%successful(starts), source=.false.)
      allocate(out%distinct(starts), source=.false.)
      allocate(out%duplicate_of(starts), source=0)
      do start = 1, starts
         out%fit(start) = gmvarkit_estimate(series, initial_model(start), &
            conditional=use_conditional, max_iterations=max_iterations, &
            tolerance=tolerance)
         out%successful(start) = out%fit(start)%log_likelihood > &
            -0.5_dp*huge(1.0_dp)
         out%order(start) = start
      end do
      out%successful_count = count(out%successful)
      do rank = 1, starts - 1
         selected = rank
         do candidate = rank + 1, starts
            if (out%fit(out%order(candidate))%log_likelihood > &
               out%fit(out%order(selected))%log_likelihood) selected = candidate
         end do
         if (selected /= rank) then
            swap_index = out%order(rank)
            out%order(rank) = out%order(selected)
            out%order(selected) = swap_index
         end if
      end do
      do rank = 1, starts
         candidate = out%order(rank)
         if (.not. out%successful(candidate)) cycle
         call canonical_model_parameters(out%fit(candidate)%model, &
            candidate_parameter, info)
         if (info /= 0) then
            out%info = 5
            return
         end if
         out%distinct(candidate) = .true.
         do previous_rank = 1, rank - 1
            previous = out%order(previous_rank)
            if (.not. out%distinct(previous)) cycle
            likelihood_scale = max(1.0_dp, &
               abs(out%fit(candidate)%log_likelihood), &
               abs(out%fit(previous)%log_likelihood))
            if (abs(out%fit(candidate)%log_likelihood - &
               out%fit(previous)%log_likelihood) > &
               log_tolerance*likelihood_scale) cycle
            call canonical_model_parameters(out%fit(previous)%model, &
               previous_parameter, info)
            if (info /= 0) then
               out%info = 6
               return
            end if
            parameter_scale = max(1.0_dp, maxval(abs(candidate_parameter)), &
               maxval(abs(previous_parameter)))
            parameter_distance = maxval(abs(candidate_parameter - &
               previous_parameter))/parameter_scale
            if (parameter_distance <= par_tolerance) then
               out%distinct(candidate) = .false.
               out%duplicate_of(candidate) = previous
               exit
            end if
         end do
         if (out%distinct(candidate)) then
            out%distinct_count = out%distinct_count + 1
            if (out%best_index == 0) out%best_index = candidate
         end if
      end do
      if (out%successful_count == 0) out%info = 7
   end function gmvarkit_multistart_estimate

   function gmvarkit_simulate(model, initial_values, observations, paths) result(out)
      !! Simulate recursive paths from a Gaussian and Student-t mixture VAR.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      real(dp), intent(in) :: initial_values(:, :) !! Initial observations by time and variable.
      integer, intent(in) :: observations !! Number of observations generated per path.
      integer, intent(in), optional :: paths !! Number of independent simulation paths.
      type(gmvarkit_simulation_t) :: out
      type(gmvarkit_regime_moments_t) :: moments
      real(dp), allocatable :: history(:, :), lag_vector(:), history_inverse(:, :, :)
      real(dp), allocatable :: history_logdet(:), factor(:, :, :), weight(:)
      real(dp), allocatable :: regime_mean(:, :), arch_scalar(:), normal(:), shock(:)
      real(dp), allocatable :: lower(:, :)
      real(dp) :: uniform, cumulative, degrees, chi_square
      integer :: variables, lags, regimes, simulations, path, time, lag
      integer :: regime, selected, info

      simulations = 1
      if (present(paths)) simulations = paths
      if (.not. valid_model(model)) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      if (observations < 1 .or. simulations < 1 .or. &
         size(initial_values, 1) < lags .or. &
         size(initial_values, 2) /= variables) then
         out%info = 2
         return
      end if
      moments = gmvarkit_regime_moments(model)
      if (moments%info /= 0) then
         out%info = 100 + moments%info
         return
      end if
      allocate(history(lags, variables), lag_vector(lags*variables))
      allocate(history_inverse(lags*variables, lags*variables, regimes))
      allocate(history_logdet(regimes), factor(variables, variables, regimes))
      allocate(weight(regimes), regime_mean(variables, regimes))
      allocate(arch_scalar(regimes), normal(variables), shock(variables))
      do regime = 1, regimes
         call cholesky_lower(moments%lag_covariance(:, :, regime), lower, info)
         if (info /= 0) then
            out%info = 200 + regime
            return
         end if
         call inverse_logdet(moments%lag_covariance(:, :, regime), &
            history_inverse(:, :, regime), history_logdet(regime), info, &
            100.0_dp*epsilon(1.0_dp))
         if (info /= 0) then
            out%info = 200 + regime
            return
         end if
         call cholesky_lower(model%covariance(:, :, regime), &
            lower=lower, info=info)
         if (info /= 0) then
            out%info = 300 + regime
            return
         end if
         factor(:, :, regime) = lower
      end do
      allocate(out%series(observations, variables, simulations))
      allocate(out%regime(observations, simulations))
      allocate(out%mixing_weight(observations, regimes, simulations))
      do path = 1, simulations
         history = initial_values(size(initial_values, 1) - lags + 1:, :)
         do time = 1, observations
            do lag = 1, lags
               lag_vector((lag - 1)*variables + 1:lag*variables) = &
                  history(lags + 1 - lag, :)
            end do
            call mixture_state(model, moments, lag_vector, history_inverse, &
               history_logdet, weight, regime_mean, arch_scalar)
            uniform = random_uniform()
            cumulative = 0.0_dp
            selected = regimes
            do regime = 1, regimes
               cumulative = cumulative + weight(regime)
               if (uniform <= cumulative) then
                  selected = regime
                  exit
               end if
            end do
            do regime = 1, variables
               normal(regime) = random_standard_normal()
            end do
            shock = matmul(factor(:, :, selected), normal)
            if (selected > model%gaussian_regimes) then
               degrees = model%degrees_of_freedom(selected) + &
                  real(variables*lags, dp)
               chi_square = 2.0_dp*random_gamma(0.5_dp*degrees)
               shock = shock*sqrt(arch_scalar(selected)*(degrees - 2.0_dp)/ &
                  max(chi_square, tiny(1.0_dp)))
            end if
            out%series(time, :, path) = regime_mean(:, selected) + shock
            out%regime(time, path) = selected
            out%mixing_weight(time, :, path) = weight
            if (lags > 1) history(:lags - 1, :) = history(2:, :)
            history(lags, :) = out%series(time, :, path)
         end do
      end do
   end function gmvarkit_simulate

   function gmvarkit_forecast(model, initial_values, horizon, simulations, &
      probabilities) result(out)
      !! Summarize independent recursive simulations as predictive forecasts.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      real(dp), intent(in) :: initial_values(:, :) !! Initial observations by time and variable.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      integer, intent(in), optional :: simulations !! Number of Monte Carlo paths.
      real(dp), intent(in), optional :: probabilities(:) !! Requested empirical quantile probabilities.
      type(gmvarkit_forecast_t) :: out
      type(gmvarkit_simulation_t) :: draws
      real(dp), allocatable :: ordered(:)
      integer :: count, variables, regimes, probability_count
      integer :: time, variable, probability, regime

      count = 2000
      if (present(simulations)) count = simulations
      if (horizon < 1 .or. count < 1) then
         out%info = 1
         return
      end if
      if (present(probabilities)) then
         if (any(probabilities < 0.0_dp) .or. any(probabilities > 1.0_dp)) then
            out%info = 2
            return
         end if
         out%probability = probabilities
      else
         out%probability = [0.025_dp, 0.1_dp, 0.9_dp, 0.975_dp]
      end if
      draws = gmvarkit_simulate(model, initial_values, horizon, paths=count)
      if (draws%info /= 0) then
         out%info = 100 + draws%info
         return
      end if
      variables = size(draws%series, 2)
      regimes = size(draws%mixing_weight, 2)
      probability_count = size(out%probability)
      allocate(out%mean(horizon, variables), out%median(horizon, variables))
      allocate(out%quantile(horizon, probability_count, variables))
      allocate(out%mixing_weight_mean(horizon, regimes), ordered(count))
      do time = 1, horizon
         do variable = 1, variables
            out%mean(time, variable) = &
               sum(draws%series(time, variable, :))/real(count, dp)
            out%median(time, variable) = median(draws%series(time, variable, :))
            ordered = sorted(draws%series(time, variable, :))
            do probability = 1, probability_count
               out%quantile(time, probability, variable) = &
                  quantile(ordered, out%probability(probability))
            end do
         end do
         do regime = 1, regimes
            out%mixing_weight_mean(time, regime) = &
               sum(draws%mixing_weight(time, regime, :))/real(count, dp)
         end do
      end do
      out%simulations = count
   end function gmvarkit_forecast

   pure function gmvarkit_regime_moments(model) result(out)
      !! Calculate regime means and stationary covariance matrices of p lags.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      type(gmvarkit_regime_moments_t) :: out
      real(dp), allocatable :: companion(:, :), innovation(:, :), system(:, :)
      real(dp), allocatable :: inverse(:, :), vectorized(:), covariance(:, :)
      real(dp), allocatable :: ar_sum(:, :), mean_inverse(:, :)
      integer :: variables, lags, regimes, state_size, regime, lag, info

      if (.not. valid_model(model)) then
         out%info = 1
         return
      end if
      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      lags = size(model%ar, 3)
      state_size = variables*lags
      allocate(out%mean(variables, regimes))
      allocate(out%lag_covariance(state_size, state_size, regimes))
      allocate(companion(state_size, state_size), innovation(state_size, state_size))
      allocate(system(state_size**2, state_size**2), vectorized(state_size**2))
      allocate(covariance(state_size, state_size), ar_sum(variables, variables))
      do regime = 1, regimes
         companion = 0.0_dp
         do lag = 1, lags
            companion(:variables, (lag - 1)*variables + 1:lag*variables) = &
               model%ar(:, :, lag, regime)
         end do
         if (lags > 1) companion(variables + 1:, :state_size - variables) = &
            identity_matrix(state_size - variables)
         innovation = 0.0_dp
         innovation(:variables, :variables) = model%covariance(:, :, regime)
         system = identity_matrix(state_size**2) - &
            kronecker_product(companion, companion)
         call invert_matrix(system, inverse, info)
         if (info /= 0) then
            out%info = 10 + regime
            return
         end if
         vectorized = matmul(inverse, reshape(innovation, [state_size**2]))
         covariance = reshape(vectorized, [state_size, state_size])
         out%lag_covariance(:, :, regime) = &
            0.5_dp*(covariance + transpose(covariance))
         ar_sum = sum(model%ar(:, :, :, regime), dim=3)
         call invert_matrix(identity_matrix(variables) - ar_sum, &
            mean_inverse, info)
         if (info /= 0) then
            out%info = 20 + regime
            return
         end if
         out%mean(:, regime) = matmul(mean_inverse, model%intercept(:, regime))
      end do
   end function gmvarkit_regime_moments

   pure function gmvarkit_evaluate(series, model, conditional) result(out)
      !! Evaluate G-StMVAR mixing weights, moments, and log likelihood.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      logical, intent(in), optional :: conditional !! Omit the initial stationary density term.
      type(gmvarkit_evaluation_t) :: out
      type(gmvarkit_regime_moments_t) :: moments
      real(dp), allocatable :: lag_vector(:), repeated_mean(:), log_history(:)
      real(dp), allocatable :: history_inverse(:, :, :), history_logdet(:)
      real(dp), allocatable :: innovation_inverse(:, :, :), innovation_logdet(:)
      real(dp), allocatable :: factor(:, :)
      real(dp), allocatable :: difference(:), regime_density(:), total_mean(:)
      real(dp), allocatable :: mean_difference(:)
      real(dp) :: quadratic, density, maximum_log, denominator, scalar
      real(dp) :: initial_log_likelihood
      integer :: observations, variables, lags, regimes, usable
      integer :: time, regime, lag, info
      logical :: use_conditional

      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      if (.not. valid_model(model)) then
         out%info = 1
         return
      end if
      observations = size(series, 1)
      variables = size(series, 2)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      usable = observations - lags
      if (variables /= size(model%intercept, 1) .or. usable < 1) then
         out%info = 2
         return
      end if
      moments = gmvarkit_regime_moments(model)
      if (moments%info /= 0) then
         out%info = 100 + moments%info
         return
      end if
      allocate(out%mixing_weight(usable, regimes))
      allocate(out%regime_mean(variables, usable, regimes))
      allocate(out%arch_scalar(usable, regimes))
      allocate(out%conditional_mean(variables, usable))
      allocate(out%conditional_covariance(variables, variables, usable))
      allocate(out%log_likelihood_term(usable))
      allocate(lag_vector(variables*lags), repeated_mean(variables*lags))
      allocate(log_history(regimes), history_inverse(variables*lags, &
         variables*lags, regimes), history_logdet(regimes))
      allocate(innovation_inverse(variables, variables, regimes))
      allocate(innovation_logdet(regimes), difference(max(variables, variables*lags)))
      allocate(regime_density(regimes), total_mean(variables))
      allocate(mean_difference(variables))
      initial_log_likelihood = 0.0_dp
      do regime = 1, regimes
         call cholesky_lower(moments%lag_covariance(:, :, regime), factor, info)
         if (info /= 0) then
            out%info = 200 + regime
            return
         end if
         call inverse_logdet(moments%lag_covariance(:, :, regime), &
            history_inverse(:, :, regime), history_logdet(regime), info, &
            100.0_dp*epsilon(1.0_dp))
         if (info /= 0) then
            out%info = 200 + regime
            return
         end if
         call cholesky_lower(model%covariance(:, :, regime), factor, info)
         if (info /= 0) then
            out%info = 300 + regime
            return
         end if
         call inverse_logdet(model%covariance(:, :, regime), &
            innovation_inverse(:, :, regime), innovation_logdet(regime), info, &
            100.0_dp*epsilon(1.0_dp))
         if (info /= 0) then
            out%info = 300 + regime
            return
         end if
         do lag = 1, lags
            repeated_mean((lag - 1)*variables + 1:lag*variables) = &
               moments%mean(:, regime)
         end do
      end do
      do time = 1, usable
         do lag = 1, lags
            lag_vector((lag - 1)*variables + 1:lag*variables) = &
               series(lags + time - lag, :)
         end do
         do regime = 1, regimes
            do lag = 1, lags
               repeated_mean((lag - 1)*variables + 1:lag*variables) = &
                  moments%mean(:, regime)
            end do
            difference(:variables*lags) = lag_vector - repeated_mean
            quadratic = dot_product(difference(:variables*lags), &
               matmul(history_inverse(:, :, regime), &
               difference(:variables*lags)))
            if (regime <= model%gaussian_regimes) then
               log_history(regime) = gaussian_log_density(variables*lags, &
                  quadratic, history_logdet(regime))
               out%arch_scalar(time, regime) = 1.0_dp
            else
               log_history(regime) = student_log_density(variables*lags, &
                  quadratic, history_logdet(regime), &
                  model%degrees_of_freedom(regime))
               out%arch_scalar(time, regime) = &
                  (model%degrees_of_freedom(regime) - 2.0_dp + quadratic)/ &
                  (model%degrees_of_freedom(regime) - 2.0_dp + &
                  real(variables*lags, dp))
            end if
         end do
         maximum_log = maxval(log(model%weight) + log_history)
         regime_density = exp(log(model%weight) + log_history - maximum_log)
         denominator = sum(regime_density)
         out%mixing_weight(time, :) = regime_density/denominator
         if (time == 1) initial_log_likelihood = maximum_log + log(denominator)
         do regime = 1, regimes
            out%regime_mean(:, time, regime) = model%intercept(:, regime)
            do lag = 1, lags
               out%regime_mean(:, time, regime) = &
                  out%regime_mean(:, time, regime) + &
                  matmul(model%ar(:, :, lag, regime), &
                  series(lags + time - lag, :))
            end do
            difference(:variables) = series(lags + time, :) - &
               out%regime_mean(:, time, regime)
            quadratic = dot_product(difference(:variables), &
               matmul(innovation_inverse(:, :, regime), difference(:variables)))
            scalar = out%arch_scalar(time, regime)
            if (regime <= model%gaussian_regimes) then
               density = gaussian_log_density(variables, quadratic, &
                  innovation_logdet(regime))
            else
               density = student_log_density(variables, quadratic/scalar, &
                  innovation_logdet(regime) + real(variables, dp)*log(scalar), &
                  model%degrees_of_freedom(regime) + real(variables*lags, dp))
            end if
            regime_density(regime) = log(out%mixing_weight(time, regime)) + density
         end do
         maximum_log = maxval(regime_density)
         out%log_likelihood_term(time) = maximum_log + &
            log(sum(exp(regime_density - maximum_log)))
         total_mean = 0.0_dp
         do regime = 1, regimes
            total_mean = total_mean + out%mixing_weight(time, regime)* &
               out%regime_mean(:, time, regime)
         end do
         out%conditional_mean(:, time) = total_mean
         out%conditional_covariance(:, :, time) = 0.0_dp
         do regime = 1, regimes
            mean_difference = out%regime_mean(:, time, regime) - total_mean
            out%conditional_covariance(:, :, time) = &
               out%conditional_covariance(:, :, time) + &
               out%mixing_weight(time, regime)*( &
               out%arch_scalar(time, regime)*model%covariance(:, :, regime) + &
               spread(mean_difference, 2, variables)* &
               spread(mean_difference, 1, variables))
         end do
      end do
      out%log_likelihood = sum(out%log_likelihood_term)
      if (.not. use_conditional) out%log_likelihood = &
         out%log_likelihood + initial_log_likelihood
   end function gmvarkit_evaluate

   pure logical function linear_bootstrap_supported(model) result(supported)
      !! Test the common-dynamics condition required by the wild bootstrap.
      type(gmvarkit_model_t), intent(in) :: model !! Model whose regime dynamics are compared.
      type(gmvarkit_regime_moments_t) :: moments
      real(dp) :: scale
      integer :: regimes, regime

      supported = .false.
      if (.not. valid_model(model)) return
      regimes = size(model%intercept, 2)
      if (regimes == 1) then
         supported = .true.
         return
      end if
      moments = gmvarkit_regime_moments(model)
      if (moments%info /= 0) return
      scale = max(1.0_dp, maxval(abs(model%ar)))
      do regime = 2, regimes
         if (maxval(abs(model%ar(:, :, :, regime) - model%ar(:, :, :, 1))) > &
            1.0e-10_dp*scale) return
         if (maxval(abs(moments%mean(:, regime) - moments%mean(:, 1))) > &
            1.0e-10_dp*max(1.0_dp, maxval(abs(moments%mean)))) return
      end do
      supported = .true.
   end function linear_bootstrap_supported

   pure function align_structure(candidate, target) result(out)
      !! Align candidate structural columns and signs to a target ordering.
      type(gmvarkit_structural_t), intent(in) :: candidate !! Newly identified structural parameters.
      type(gmvarkit_structural_t), intent(in) :: target !! Baseline shock ordering and signs.
      type(gmvarkit_structural_t) :: out
      integer, allocatable :: permutation(:)
      logical, allocatable :: used(:)
      real(dp) :: distance, best_distance, cosine_denominator
      integer :: variables, target_shock, candidate_shock, selected

      if (.not. valid_structure(candidate) .or. .not. valid_structure(target) .or. &
         any(shape(candidate%relative_variance) /= &
         shape(target%relative_variance))) then
         out%info = 1
         return
      end if
      variables = size(target%impact, 1)
      allocate(permutation(variables), used(variables))
      used = .false.
      do target_shock = 1, variables
         selected = 0
         best_distance = huge(1.0_dp)
         do candidate_shock = 1, variables
            if (used(candidate_shock)) cycle
            distance = sum((log(candidate%relative_variance(candidate_shock, :)) - &
               log(target%relative_variance(target_shock, :)))**2)
            cosine_denominator = sqrt(sum(candidate%impact(:, candidate_shock)**2)* &
               sum(target%impact(:, target_shock)**2))
            if (cosine_denominator > tiny(1.0_dp)) then
               distance = distance + 1.0_dp - abs(dot_product( &
                  candidate%impact(:, candidate_shock), &
                  target%impact(:, target_shock)))/cosine_denominator
            end if
            if (distance < best_distance) then
               best_distance = distance
               selected = candidate_shock
            end if
         end do
         if (selected == 0) then
            out%info = 2
            return
         end if
         permutation(target_shock) = selected
         used(selected) = .true.
      end do
      out = gmvarkit_reorder_structural(candidate, permutation)
      if (out%info /= 0) return
      do target_shock = 1, variables
         if (dot_product(out%impact(:, target_shock), &
            target%impact(:, target_shock)) < 0.0_dp) then
            out%impact(:, target_shock) = -out%impact(:, target_shock)
         end if
      end do
   end function align_structure

   pure logical function valid_structure(structural) result(valid)
      !! Test internal dimensions and positivity of structural parameters.
      type(gmvarkit_structural_t), intent(in) :: structural !! Structural parameters to validate.
      integer :: variables, regimes

      valid = .false.
      if (structural%info /= 0 .or. .not. allocated(structural%impact) .or. &
         .not. allocated(structural%relative_variance) .or. &
         .not. allocated(structural%covariance)) return
      variables = size(structural%impact, 1)
      regimes = size(structural%relative_variance, 2)
      if (variables < 1 .or. regimes < 1) return
      if (size(structural%impact, 2) /= variables .or. &
         size(structural%relative_variance, 1) /= variables .or. &
         any(shape(structural%covariance) /= [variables, variables, regimes]) .or. &
         any(structural%relative_variance <= 0.0_dp)) return
      valid = structural%reference_regime >= 1 .and. &
         structural%reference_regime <= regimes
   end function valid_structure

   pure logical function compatible_structure(structural, model) result(valid)
      !! Test whether structural parameters reproduce a model's covariances.
      type(gmvarkit_structural_t), intent(in) :: structural !! Structural parameters to validate.
      type(gmvarkit_model_t), intent(in) :: model !! Reduced-form model to reproduce.
      real(dp) :: scale

      valid = .false.
      if (.not. valid_structure(structural) .or. .not. valid_model(model)) return
      if (any(shape(structural%covariance) /= shape(model%covariance))) return
      scale = max(1.0_dp, maxval(abs(model%covariance)))
      valid = maxval(abs(structural%covariance - model%covariance))/scale <= 1.0e-7_dp
   end function compatible_structure

   pure integer function regime_from_uniform(weight, uniform) result(selected)
      !! Map a common uniform draw through a regime probability vector.
      real(dp), intent(in) :: weight(:) !! Regime probabilities summing to one.
      real(dp), intent(in) :: uniform !! Uniform draw on the unit interval.
      real(dp) :: cumulative
      integer :: regime

      selected = size(weight)
      cumulative = 0.0_dp
      do regime = 1, size(weight)
         cumulative = cumulative + weight(regime)
         if (uniform <= cumulative) then
            selected = regime
            return
         end if
      end do
   end function regime_from_uniform

   pure subroutine mixture_state(model, moments, lag_vector, history_inverse, &
      history_logdet, weight, regime_mean, arch_scalar)
      !! Compute regime probabilities and moments conditional on one lag vector.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      type(gmvarkit_regime_moments_t), intent(in) :: moments !! Stationary regime moments.
      real(dp), intent(in) :: lag_vector(:) !! Most-recent-first stacked lag vector.
      real(dp), intent(in) :: history_inverse(:, :, :) !! Inverse stationary lag covariances.
      real(dp), intent(in) :: history_logdet(:) !! Stationary lag-covariance log determinants.
      real(dp), intent(out) :: weight(:) !! Conditional regime probabilities.
      real(dp), intent(out) :: regime_mean(:, :) !! Regimewise conditional means.
      real(dp), intent(out) :: arch_scalar(:) !! Regimewise conditional covariance scalars.
      real(dp) :: repeated_mean(size(lag_vector)), difference(size(lag_vector))
      real(dp) :: log_density(size(weight)), unnormalized(size(weight))
      real(dp) :: quadratic, maximum_log
      integer :: variables, lags, regimes, regime, lag

      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      do regime = 1, regimes
         do lag = 1, lags
            repeated_mean((lag - 1)*variables + 1:lag*variables) = &
               moments%mean(:, regime)
         end do
         difference = lag_vector - repeated_mean
         quadratic = dot_product(difference, &
            matmul(history_inverse(:, :, regime), difference))
         if (regime <= model%gaussian_regimes) then
            log_density(regime) = gaussian_log_density(variables*lags, &
               quadratic, history_logdet(regime))
            arch_scalar(regime) = 1.0_dp
         else
            log_density(regime) = student_log_density(variables*lags, &
               quadratic, history_logdet(regime), &
               model%degrees_of_freedom(regime))
            arch_scalar(regime) = &
               (model%degrees_of_freedom(regime) - 2.0_dp + quadratic)/ &
               (model%degrees_of_freedom(regime) - 2.0_dp + &
               real(variables*lags, dp))
         end if
         regime_mean(:, regime) = model%intercept(:, regime)
         do lag = 1, lags
            regime_mean(:, regime) = regime_mean(:, regime) + &
               matmul(model%ar(:, :, lag, regime), &
               lag_vector((lag - 1)*variables + 1:lag*variables))
         end do
      end do
      maximum_log = maxval(log(model%weight) + log_density)
      unnormalized = exp(log(model%weight) + log_density - maximum_log)
      weight = unnormalized/sum(unnormalized)
   end subroutine mixture_state

   pure logical function same_model_specification(first, second) result(same)
      !! Check whether two valid models share dimensions and regime families.
      type(gmvarkit_model_t), intent(in) :: first !! First valid model.
      type(gmvarkit_model_t), intent(in) :: second !! Second valid model.

      same = all(shape(first%intercept) == shape(second%intercept)) .and. &
         all(shape(first%ar) == shape(second%ar)) .and. &
         all(shape(first%covariance) == shape(second%covariance)) .and. &
         size(first%weight) == size(second%weight) .and. &
         first%gaussian_regimes == second%gaussian_regimes
   end function same_model_specification

   pure subroutine canonical_model_parameters(model, parameters, info)
      !! Pack a model after deterministic within-family regime ordering.
      type(gmvarkit_model_t), intent(in) :: model !! Valid model to canonicalize.
      real(dp), allocatable, intent(out) :: parameters(:) !! Canonically ordered transformed parameters.
      integer, intent(out) :: info !! Status code; zero indicates success.
      type(gmvarkit_model_t) :: canonical
      integer, allocatable :: permutation(:)
      logical, allocatable :: used(:)
      integer :: regimes, gaussian, group, lower, upper, position, candidate
      integer :: selected

      regimes = size(model%intercept, 2)
      gaussian = model%gaussian_regimes
      allocate(permutation(regimes), source=0)
      allocate(used(regimes), source=.false.)
      position = 0
      do group = 1, 2
         if (group == 1) then
            lower = 1
            upper = gaussian
         else
            lower = gaussian + 1
            upper = regimes
         end if
         do while (position < merge(gaussian, regimes, group == 1))
            selected = 0
            do candidate = lower, upper
               if (used(candidate)) cycle
               if (selected == 0) then
                  selected = candidate
               else if (regime_precedes(model, candidate, selected)) then
                  selected = candidate
               end if
            end do
            if (selected == 0) exit
            position = position + 1
            permutation(position) = selected
            used(selected) = .true.
         end do
      end do
      canonical = model
      canonical%intercept = model%intercept(:, permutation)
      canonical%ar = model%ar(:, :, :, permutation)
      canonical%covariance = model%covariance(:, :, permutation)
      canonical%weight = model%weight(permutation)
      canonical%degrees_of_freedom = model%degrees_of_freedom(permutation)
      call pack_model(canonical, parameters, info)
   end subroutine canonical_model_parameters

   pure logical function regime_precedes(model, first, second) result(precedes)
      !! Compare two regimes by weight and then lexicographic parameter values.
      type(gmvarkit_model_t), intent(in) :: model !! Model containing both regimes.
      integer, intent(in) :: first !! Candidate regime index.
      integer, intent(in) :: second !! Current regime index.
      real(dp), allocatable :: first_value(:), second_value(:)
      integer :: variables, lags, covariance_size, position, index

      if (model%weight(first) /= model%weight(second)) then
         precedes = model%weight(first) > model%weight(second)
         return
      end if
      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      covariance_size = variables**2
      allocate(first_value(variables + variables**2*lags + covariance_size + 1))
      allocate(second_value(size(first_value)))
      position = variables
      first_value(:position) = model%intercept(:, first)
      second_value(:position) = model%intercept(:, second)
      first_value(position + 1:position + variables**2*lags) = &
         reshape(model%ar(:, :, :, first), [variables**2*lags])
      second_value(position + 1:position + variables**2*lags) = &
         reshape(model%ar(:, :, :, second), [variables**2*lags])
      position = position + variables**2*lags
      first_value(position + 1:position + covariance_size) = &
         reshape(model%covariance(:, :, first), [covariance_size])
      second_value(position + 1:position + covariance_size) = &
         reshape(model%covariance(:, :, second), [covariance_size])
      position = position + covariance_size
      first_value(position + 1) = model%degrees_of_freedom(first)
      second_value(position + 1) = model%degrees_of_freedom(second)
      precedes = .false.
      do index = 1, size(first_value)
         if (first_value(index) < second_value(index)) then
            precedes = .true.
            return
         else if (first_value(index) > second_value(index)) then
            return
         end if
      end do
   end function regime_precedes

   pure subroutine initialize_likelihood_profile(parameters, scale, points, &
      structural, profile, selection)
      !! Allocate equally spaced transformed-coordinate likelihood grids.
      real(dp), intent(in) :: parameters(:) !! Complete transformed parameter vector.
      real(dp), intent(in) :: scale !! Relative half-width with a unit magnitude floor.
      integer, intent(in) :: points !! Number of grid values per coordinate.
      logical, intent(in) :: structural !! Whether structural coordinates are used.
      type(gmvarkit_likelihood_profile_t), intent(out) :: profile !! Initialized profile result.
      integer, intent(in), optional :: selection(:) !! Coordinates included in the result.
      real(dp) :: half_width
      integer :: coordinate, grid, other

      if (points < 2 .or. scale <= 0.0_dp .or. size(parameters) < 1) then
         profile%info = 1
         return
      end if
      if (present(selection)) then
         if (size(selection) < 1 .or. any(selection < 1) .or. &
            any(selection > size(parameters))) then
            profile%info = 2
            return
         end if
         do coordinate = 1, size(selection)
            do other = coordinate + 1, size(selection)
               if (selection(coordinate) == selection(other)) then
                  profile%info = 3
                  return
               end if
            end do
         end do
         profile%parameter = selection
      else
         profile%parameter = [(coordinate, coordinate=1, size(parameters))]
      end if
      profile%center = parameters(profile%parameter)
      allocate(profile%value(points, size(profile%parameter)))
      allocate(profile%log_likelihood(points, size(profile%parameter)), &
         source=-huge(1.0_dp))
      allocate(profile%valid(points, size(profile%parameter)), source=.false.)
      do coordinate = 1, size(profile%parameter)
         half_width = scale*max(1.0_dp, abs(profile%center(coordinate)))
         do grid = 1, points
            profile%value(grid, coordinate) = profile%center(coordinate) - &
               half_width + 2.0_dp*half_width*real(grid - 1, dp)/ &
               real(points - 1, dp)
         end do
      end do
      profile%structural = structural
      profile%info = 0
   end subroutine initialize_likelihood_profile

   pure logical function valid_structural_constraints(model, structural, &
      constraints) result(valid)
      !! Check dimensions and consistency of direct structural restrictions.
      type(gmvarkit_model_t), intent(in) :: model !! Reduced-form model dimensions.
      type(gmvarkit_structural_t), intent(in) :: structural !! Structural starting values.
      type(gmvarkit_structural_constraints_t), intent(in) :: constraints !! Restrictions to validate.
      integer :: variables, regimes, reference, row, column

      valid = .false.
      if (constraints%info /= 0 .or. .not. compatible_structure(structural, model)) &
         return
      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      reference = structural%reference_regime
      if (allocated(constraints%impact_fixed) .neqv. &
         allocated(constraints%impact_value)) return
      if (allocated(constraints%impact_fixed)) then
         if (any(shape(constraints%impact_fixed) /= [variables, variables]) .or. &
            any(shape(constraints%impact_value) /= [variables, variables])) return
         if (any(constraints%impact_fixed)) then
            if (maxval(abs(pack(structural%impact - constraints%impact_value, &
               constraints%impact_fixed))) > 1.0e-10_dp* &
               max(1.0_dp, maxval(abs(structural%impact)))) return
         end if
      end if
      if (allocated(constraints%impact_sign)) then
         if (any(shape(constraints%impact_sign) /= [variables, variables]) .or. &
            any(abs(constraints%impact_sign) > 1)) return
         do column = 1, variables
            do row = 1, variables
               if (allocated(constraints%impact_fixed)) then
                  if (constraints%impact_fixed(row, column)) cycle
               end if
               if (constraints%impact_sign(row, column) /= 0 .and. &
                  sign(1.0_dp, structural%impact(row, column)) /= &
                  real(constraints%impact_sign(row, column), dp)) return
            end do
         end do
      end if
      if (allocated(constraints%lambda_mapping) .and. &
         allocated(constraints%fixed_relative_variance)) return
      if (allocated(constraints%fixed_relative_variance)) then
         if (any(shape(constraints%fixed_relative_variance) /= &
            [variables, regimes]) .or. &
            any(constraints%fixed_relative_variance <= 0.0_dp) .or. &
            maxval(abs(constraints%fixed_relative_variance(:, reference) - &
            1.0_dp)) > 1.0e-12_dp .or. &
            maxval(abs(constraints%fixed_relative_variance - &
            structural%relative_variance)) > 1.0e-10_dp*max(1.0_dp, &
            maxval(structural%relative_variance))) return
      end if
      if (allocated(constraints%lambda_mapping)) then
         if (regimes < 2 .or. size(constraints%lambda_mapping, 1) /= &
            variables*(regimes - 1) .or. &
            size(constraints%lambda_mapping, 2) < 1 .or. &
            any(constraints%lambda_mapping < 0.0_dp) .or. &
            any(sum(constraints%lambda_mapping, dim=2) <= 0.0_dp) .or. &
            matrix_rank(constraints%lambda_mapping, 1.0e-10_dp) /= &
            size(constraints%lambda_mapping, 2)) return
      end if
      valid = .true.
   end function valid_structural_constraints

   pure subroutine pack_structural_model(model, structural, constraints, &
      parameters, info)
      !! Map a structural model to free unconstrained coordinates.
      type(gmvarkit_model_t), intent(in) :: model !! Reduced-form model values.
      type(gmvarkit_structural_t), intent(in) :: structural !! Structural model values.
      type(gmvarkit_structural_constraints_t), intent(in) :: constraints !! Structural restrictions.
      real(dp), allocatable, intent(out) :: parameters(:) !! Free transformed coordinates.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: lambda_vector(:), cross(:, :), inverse(:, :), gamma(:)
      integer :: variables, regimes, lags, reference, impact_count, lambda_count
      integer :: total, position, regime, shock, row, column, sign_constraint
      logical :: fixed

      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      lags = size(model%ar, 3)
      reference = structural%reference_regime
      impact_count = variables**2
      if (allocated(constraints%impact_fixed)) impact_count = &
         count(.not. constraints%impact_fixed)
      lambda_count = variables*(regimes - 1)
      if (allocated(constraints%fixed_relative_variance)) lambda_count = 0
      if (allocated(constraints%lambda_mapping)) lambda_count = &
         size(constraints%lambda_mapping, 2)
      total = variables*regimes + variables**2*lags*regimes + impact_count + &
         lambda_count + regimes - 1 + regimes - model%gaussian_regimes
      allocate(parameters(total))
      position = 0
      parameters(1:variables*regimes) = reshape(model%intercept, &
         [variables*regimes])
      position = variables*regimes
      parameters(position + 1:position + variables**2*lags*regimes) = &
         reshape(model%ar, [variables**2*lags*regimes])
      position = position + variables**2*lags*regimes
      do column = 1, variables
         do row = 1, variables
            fixed = .false.
            if (allocated(constraints%impact_fixed)) fixed = &
               constraints%impact_fixed(row, column)
            if (fixed) cycle
            sign_constraint = 0
            if (allocated(constraints%impact_sign)) sign_constraint = &
               constraints%impact_sign(row, column)
            position = position + 1
            if (sign_constraint == 0) then
               parameters(position) = structural%impact(row, column)
            else
               parameters(position) = log(max(abs(structural%impact(row, column)), &
                  tiny(1.0_dp)))
            end if
         end do
      end do
      if (lambda_count > 0) then
         allocate(lambda_vector(variables*(regimes - 1)))
         shock = 0
         do regime = 1, regimes
            if (regime == reference) cycle
            lambda_vector(shock + 1:shock + variables) = &
               structural%relative_variance(:, regime)
            shock = shock + variables
         end do
         if (allocated(constraints%lambda_mapping)) then
            cross = matmul(transpose(constraints%lambda_mapping), &
               constraints%lambda_mapping)
            call invert_matrix(cross, inverse, info)
            if (info /= 0) return
            gamma = matmul(inverse, matmul(transpose( &
               constraints%lambda_mapping), lambda_vector))
            if (any(gamma <= 0.0_dp) .or. maxval(abs(matmul( &
               constraints%lambda_mapping, gamma) - lambda_vector)) > &
               1.0e-8_dp*max(1.0_dp, maxval(lambda_vector))) then
               info = 2
               return
            end if
            parameters(position + 1:position + lambda_count) = log(gamma)
         else
            parameters(position + 1:position + lambda_count) = log(lambda_vector)
         end if
         position = position + lambda_count
      end if
      do regime = 1, regimes - 1
         position = position + 1
         parameters(position) = log(model%weight(regime)/model%weight(regimes))
      end do
      do regime = model%gaussian_regimes + 1, regimes
         position = position + 1
         parameters(position) = log(model%degrees_of_freedom(regime) - 2.0_dp)
      end do
      info = 0
   end subroutine pack_structural_model

   pure subroutine unpack_structural_model(parameters, template_model, &
      template_structural, constraints, model, structural, info)
      !! Rebuild reduced-form covariances from free structural coordinates.
      real(dp), intent(in) :: parameters(:) !! Free transformed coordinates.
      type(gmvarkit_model_t), intent(in) :: template_model !! Model dimensions and regime types.
      type(gmvarkit_structural_t), intent(in) :: template_structural !! Structural dimensions and reference.
      type(gmvarkit_structural_constraints_t), intent(in) :: constraints !! Structural restrictions.
      type(gmvarkit_model_t), intent(out) :: model !! Reconstructed reduced-form model.
      type(gmvarkit_structural_t), intent(out) :: structural !! Reconstructed structural model.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: impact(:, :), lambda(:, :), lambda_vector(:)
      real(dp), allocatable :: logits(:), gamma(:)
      real(dp) :: maximum_logit
      integer :: variables, regimes, lags, reference, position
      integer :: regime, row, column, sign_constraint, lambda_count, location
      logical :: fixed

      model = template_model
      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      lags = size(model%ar, 3)
      reference = template_structural%reference_regime
      position = variables*regimes
      model%intercept = reshape(parameters(:position), [variables, regimes])
      model%ar = reshape(parameters(position + 1: &
         position + variables**2*lags*regimes), &
         [variables, variables, lags, regimes])
      position = position + variables**2*lags*regimes
      allocate(impact(variables, variables))
      do column = 1, variables
         do row = 1, variables
            fixed = .false.
            if (allocated(constraints%impact_fixed)) fixed = &
               constraints%impact_fixed(row, column)
            if (fixed) then
               impact(row, column) = constraints%impact_value(row, column)
            else
               position = position + 1
               sign_constraint = 0
               if (allocated(constraints%impact_sign)) sign_constraint = &
                  constraints%impact_sign(row, column)
               if (sign_constraint == 0) then
                  impact(row, column) = parameters(position)
               else
                  impact(row, column) = real(sign_constraint, dp)* &
                     exp(max(-30.0_dp, min(30.0_dp, parameters(position))))
               end if
            end if
         end do
      end do
      allocate(lambda(variables, regimes), source=1.0_dp)
      if (allocated(constraints%fixed_relative_variance)) then
         lambda = constraints%fixed_relative_variance
      else if (regimes > 1) then
         lambda_count = variables*(regimes - 1)
         if (allocated(constraints%lambda_mapping)) then
            lambda_count = size(constraints%lambda_mapping, 2)
            gamma = exp(max(-30.0_dp, min(30.0_dp, &
               parameters(position + 1:position + lambda_count))))
            lambda_vector = matmul(constraints%lambda_mapping, gamma)
         else
            lambda_vector = exp(max(-30.0_dp, min(30.0_dp, &
               parameters(position + 1:position + lambda_count))))
         end if
         position = position + lambda_count
         location = 0
         do regime = 1, regimes
            if (regime == reference) cycle
            lambda(:, regime) = lambda_vector(location + 1:location + variables)
            location = location + variables
         end do
      end if
      structural = gmvarkit_structural_covariances(impact, lambda, reference)
      if (structural%info /= 0) then
         info = 1
         return
      end if
      model%covariance = structural%covariance
      allocate(logits(regimes), source=0.0_dp)
      do regime = 1, regimes - 1
         position = position + 1
         logits(regime) = parameters(position)
      end do
      maximum_logit = maxval(logits)
      model%weight = exp(logits - maximum_logit)
      model%weight = model%weight/sum(model%weight)
      model%degrees_of_freedom(:model%gaussian_regimes) = 0.0_dp
      do regime = model%gaussian_regimes + 1, regimes
         position = position + 1
         model%degrees_of_freedom(regime) = 2.0_dp + &
            exp(max(-30.0_dp, min(30.0_dp, parameters(position))))
      end do
      model%info = 0
      info = 0
   end subroutine unpack_structural_model

   pure subroutine pack_model(model, parameters, info)
      !! Map a valid model to the unconstrained optimization parameter vector.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      real(dp), allocatable, intent(out) :: parameters(:) !! Unconstrained parameter vector.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: lower(:, :)
      integer :: variables, lags, regimes, covariance_size, total
      integer :: regime, row, column, position, block_size

      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      lags = size(model%ar, 3)
      covariance_size = variables*(variables + 1)/2
      block_size = variables + lags*variables**2 + covariance_size
      total = regimes*block_size + regimes - 1 + regimes - model%gaussian_regimes
      allocate(parameters(total))
      position = 0
      info = 0
      do regime = 1, regimes
         parameters(position + 1:position + variables) = model%intercept(:, regime)
         position = position + variables
         parameters(position + 1:position + lags*variables**2) = &
            reshape(model%ar(:, :, :, regime), [lags*variables**2])
         position = position + lags*variables**2
         call cholesky_lower(model%covariance(:, :, regime), lower, info)
         if (info /= 0) return
         do row = 1, variables
            do column = 1, row
               position = position + 1
               if (row == column) then
                  parameters(position) = log(lower(row, column))
               else
                  parameters(position) = lower(row, column)
               end if
            end do
         end do
      end do
      do regime = 1, regimes - 1
         position = position + 1
         parameters(position) = log(model%weight(regime)/model%weight(regimes))
      end do
      do regime = model%gaussian_regimes + 1, regimes
         position = position + 1
         parameters(position) = log(model%degrees_of_freedom(regime) - 2.0_dp)
      end do
   end subroutine pack_model

   pure function unpack_model(parameters, template) result(model)
      !! Transform unconstrained optimization parameters into a valid model form.
      real(dp), intent(in) :: parameters(:) !! Unconstrained parameter vector.
      type(gmvarkit_model_t), intent(in) :: template !! Model dimensions and regime types.
      type(gmvarkit_model_t) :: model
      real(dp), allocatable :: lower(:, :), logits(:)
      real(dp) :: maximum_logit
      integer :: variables, lags, regimes, regime, row, column, position

      model = template
      variables = size(template%intercept, 1)
      regimes = size(template%intercept, 2)
      lags = size(template%ar, 3)
      allocate(lower(variables, variables), logits(regimes))
      position = 0
      do regime = 1, regimes
         model%intercept(:, regime) = parameters(position + 1:position + variables)
         position = position + variables
         model%ar(:, :, :, regime) = reshape( &
            parameters(position + 1:position + lags*variables**2), &
            [variables, variables, lags])
         position = position + lags*variables**2
         lower = 0.0_dp
         do row = 1, variables
            do column = 1, row
               position = position + 1
               if (row == column) then
                  lower(row, column) = exp(max(-30.0_dp, &
                     min(30.0_dp, parameters(position))))
               else
                  lower(row, column) = parameters(position)
               end if
            end do
         end do
         model%covariance(:, :, regime) = matmul(lower, transpose(lower))
      end do
      logits = 0.0_dp
      do regime = 1, regimes - 1
         position = position + 1
         logits(regime) = parameters(position)
      end do
      maximum_logit = maxval(logits)
      model%weight = exp(logits - maximum_logit)
      model%weight = model%weight/sum(model%weight)
      model%degrees_of_freedom(:model%gaussian_regimes) = 0.0_dp
      do regime = model%gaussian_regimes + 1, regimes
         position = position + 1
         model%degrees_of_freedom(regime) = 2.0_dp + &
            exp(max(-30.0_dp, min(30.0_dp, parameters(position))))
      end do
      model%info = 0
   end function unpack_model

   pure subroutine pack_constrained_model(model, constraints, parameters, info)
      !! Map a model to free coordinates under the supplied restrictions.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      type(gmvarkit_constraints_t), intent(in) :: constraints !! Parameter restrictions.
      real(dp), allocatable, intent(out) :: parameters(:) !! Free transformed parameters.
      integer, intent(out) :: info !! Status code; zero indicates success.
      type(gmvarkit_regime_moments_t) :: moments
      real(dp), allocatable :: full_ar(:), free_ar(:), cross(:, :), inverse(:, :)
      real(dp), allocatable :: lower(:, :)
      integer :: variables, lags, regimes, groups, location_count, ar_count
      integer :: covariance_count, weight_count, df_count, total, position
      integer :: regime, group, member_count, row, column

      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      groups = regimes
      if (allocated(constraints%mean_group)) groups = maxval(constraints%mean_group)
      location_count = variables*regimes
      if (allocated(constraints%mean_group)) location_count = variables*groups
      ar_count = variables**2*lags*regimes
      if (allocated(constraints%ar_mapping)) ar_count = size(constraints%ar_mapping, 2)
      covariance_count = regimes*variables*(variables + 1)/2
      weight_count = regimes - 1
      if (allocated(constraints%fixed_weight)) weight_count = 0
      df_count = regimes - model%gaussian_regimes
      total = location_count + ar_count + covariance_count + weight_count + df_count
      allocate(parameters(total), full_ar(variables**2*lags*regimes))
      position = 0
      info = 0
      if (allocated(constraints%mean_group)) then
         moments = gmvarkit_regime_moments(model)
         if (moments%info /= 0) then
            info = 1
            return
         end if
         do group = 1, groups
            member_count = count(constraints%mean_group == group)
            parameters(position + 1:position + variables) = &
               sum(moments%mean(:, pack([(regime, regime=1, regimes)], &
               constraints%mean_group == group)), dim=2)/real(member_count, dp)
            position = position + variables
         end do
      else
         parameters(position + 1:position + variables*regimes) = &
            reshape(model%intercept, [variables*regimes])
         position = position + variables*regimes
      end if
      do regime = 1, regimes
         full_ar((regime - 1)*variables**2*lags + 1: &
            regime*variables**2*lags) = &
            reshape(model%ar(:, :, :, regime), [variables**2*lags])
      end do
      if (allocated(constraints%ar_mapping)) then
         cross = matmul(transpose(constraints%ar_mapping), &
            constraints%ar_mapping)
         call invert_matrix(cross, inverse, info)
         if (info /= 0) return
         free_ar = matmul(inverse, &
            matmul(transpose(constraints%ar_mapping), full_ar))
         parameters(position + 1:position + ar_count) = free_ar
      else
         parameters(position + 1:position + ar_count) = full_ar
      end if
      position = position + ar_count
      do regime = 1, regimes
         call cholesky_lower(model%covariance(:, :, regime), lower, info)
         if (info /= 0) return
         do row = 1, variables
            do column = 1, row
               position = position + 1
               if (row == column) then
                  parameters(position) = log(lower(row, column))
               else
                  parameters(position) = lower(row, column)
               end if
            end do
         end do
      end do
      if (.not. allocated(constraints%fixed_weight)) then
         do regime = 1, regimes - 1
            position = position + 1
            parameters(position) = log(model%weight(regime)/model%weight(regimes))
         end do
      end if
      do regime = model%gaussian_regimes + 1, regimes
         position = position + 1
         parameters(position) = log(model%degrees_of_freedom(regime) - 2.0_dp)
      end do
   end subroutine pack_constrained_model

   pure function unpack_constrained_model(parameters, template, constraints) &
      result(model)
      !! Reconstruct a full model from restricted free coordinates.
      real(dp), intent(in) :: parameters(:) !! Free transformed parameters.
      type(gmvarkit_model_t), intent(in) :: template !! Model dimensions and regime types.
      type(gmvarkit_constraints_t), intent(in) :: constraints !! Parameter restrictions.
      type(gmvarkit_model_t) :: model
      real(dp), allocatable :: group_mean(:, :), free_ar(:), full_ar(:)
      real(dp), allocatable :: lower(:, :), logits(:), ar_sum(:, :)
      real(dp) :: maximum_logit
      integer :: variables, lags, regimes, groups, ar_count, position
      integer :: regime, row, column

      model = template
      variables = size(template%intercept, 1)
      lags = size(template%ar, 3)
      regimes = size(template%intercept, 2)
      groups = regimes
      if (allocated(constraints%mean_group)) groups = maxval(constraints%mean_group)
      allocate(lower(variables, variables), logits(regimes), &
         full_ar(variables**2*lags*regimes), ar_sum(variables, variables))
      position = 0
      if (allocated(constraints%mean_group)) then
         allocate(group_mean(variables, groups))
         group_mean = reshape(parameters(position + 1:position + variables*groups), &
            [variables, groups])
         position = position + variables*groups
      else
         model%intercept = reshape( &
            parameters(position + 1:position + variables*regimes), &
            [variables, regimes])
         position = position + variables*regimes
      end if
      ar_count = variables**2*lags*regimes
      if (allocated(constraints%ar_mapping)) then
         ar_count = size(constraints%ar_mapping, 2)
         free_ar = parameters(position + 1:position + ar_count)
         full_ar = matmul(constraints%ar_mapping, free_ar)
      else
         full_ar = parameters(position + 1:position + ar_count)
      end if
      position = position + ar_count
      do regime = 1, regimes
         model%ar(:, :, :, regime) = reshape( &
            full_ar((regime - 1)*variables**2*lags + 1: &
            regime*variables**2*lags), [variables, variables, lags])
      end do
      if (allocated(constraints%mean_group)) then
         do regime = 1, regimes
            ar_sum = sum(model%ar(:, :, :, regime), dim=3)
            model%intercept(:, regime) = matmul(identity_matrix(variables) - &
               ar_sum, group_mean(:, constraints%mean_group(regime)))
         end do
      end if
      do regime = 1, regimes
         lower = 0.0_dp
         do row = 1, variables
            do column = 1, row
               position = position + 1
               if (row == column) then
                  lower(row, column) = exp(max(-30.0_dp, &
                     min(30.0_dp, parameters(position))))
               else
                  lower(row, column) = parameters(position)
               end if
            end do
         end do
         model%covariance(:, :, regime) = matmul(lower, transpose(lower))
      end do
      if (allocated(constraints%fixed_weight)) then
         model%weight = constraints%fixed_weight
      else
         logits = 0.0_dp
         do regime = 1, regimes - 1
            position = position + 1
            logits(regime) = parameters(position)
         end do
         maximum_logit = maxval(logits)
         model%weight = exp(logits - maximum_logit)
         model%weight = model%weight/sum(model%weight)
      end if
      model%degrees_of_freedom(:model%gaussian_regimes) = 0.0_dp
      do regime = model%gaussian_regimes + 1, regimes
         position = position + 1
         model%degrees_of_freedom(regime) = 2.0_dp + &
            exp(max(-30.0_dp, min(30.0_dp, parameters(position))))
      end do
      model%info = 0
   end function unpack_constrained_model

   pure logical function valid_constraints(model, constraints) result(valid)
      !! Check constraint dimensions, group labels, and fixed weights.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      type(gmvarkit_constraints_t), intent(in) :: constraints !! Parameter restrictions.
      integer :: variables, lags, regimes, groups, group

      valid = .false.
      if (constraints%info /= 0) return
      variables = size(model%intercept, 1)
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      if (allocated(constraints%ar_mapping)) then
         if (size(constraints%ar_mapping, 1) /= variables**2*lags*regimes .or. &
            size(constraints%ar_mapping, 2) < 1) return
      end if
      if (allocated(constraints%mean_group)) then
         if (size(constraints%mean_group) /= regimes .or. &
            any(constraints%mean_group < 1)) return
         groups = maxval(constraints%mean_group)
         do group = 1, groups
            if (count(constraints%mean_group == group) < 1) return
         end do
      end if
      if (allocated(constraints%fixed_weight)) then
         if (size(constraints%fixed_weight) /= regimes .or. &
            any(constraints%fixed_weight <= 0.0_dp) .or. &
            abs(sum(constraints%fixed_weight) - 1.0_dp) > 1.0e-10_dp) return
      end if
      valid = .true.
   end function valid_constraints

   pure subroutine order_values(values, order)
      !! Return indices that arrange objective values in ascending order.
      real(dp), intent(in) :: values(:) !! Objective values.
      integer, intent(out) :: order(:) !! Ascending index permutation.
      integer :: item, position, held

      do item = 1, size(values)
         order(item) = item
      end do
      do item = 2, size(values)
         held = order(item)
         position = item - 1
         do while (position >= 1)
            if (values(order(position)) <= values(held)) exit
            order(position + 1) = order(position)
            position = position - 1
         end do
         order(position + 1) = held
      end do
   end subroutine order_values

   pure logical function valid_model(model) result(valid)
      !! Test the dimensions and scalar restrictions of a mixture VAR model.
      type(gmvarkit_model_t), intent(in) :: model !! G-StMVAR model parameters.
      integer :: variables, regimes

      valid = .false.
      if (model%info /= 0 .or. .not. allocated(model%intercept) .or. &
         .not. allocated(model%ar) .or. .not. allocated(model%covariance) .or. &
         .not. allocated(model%weight) .or. &
         .not. allocated(model%degrees_of_freedom)) return
      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      if (variables < 1 .or. regimes < 1 .or. size(model%ar, 1) /= variables .or. &
         size(model%ar, 2) /= variables .or. size(model%ar, 3) < 1 .or. &
         size(model%ar, 4) /= regimes .or. &
         any(shape(model%covariance) /= [variables, variables, regimes]) .or. &
         size(model%weight) /= regimes .or. &
         size(model%degrees_of_freedom) /= regimes .or. &
         model%gaussian_regimes < 0 .or. model%gaussian_regimes > regimes .or. &
         any(model%weight <= 0.0_dp) .or. &
         abs(sum(model%weight) - 1.0_dp) > 1.0e-10_dp) return
      if (model%gaussian_regimes < regimes) then
         if (any(model%degrees_of_freedom(model%gaussian_regimes + 1:) <= &
            2.0_dp)) return
      end if
      valid = .true.
   end function valid_model

   pure function residual_moment_test(values) result(out)
      !! Form a chi-square test from sample moment vectors and their second moment.
      real(dp), intent(in) :: values(:, :) !! Moment vectors by usable observation.
      type(gmvarkit_hypothesis_test_t) :: out
      real(dp), allocatable :: covariance(:, :), inverse(:, :), total(:)
      integer :: observations, moments, info

      observations = size(values, 1)
      moments = size(values, 2)
      if (observations < 1 .or. moments < 1) then
         out%info = 1
         return
      end if
      covariance = matmul(transpose(values), values)/real(observations, dp)
      call invert_matrix(covariance, inverse, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      total = sum(values, dim=1)
      out%statistic = max(0.0_dp, dot_product(total, matmul(inverse, total))/ &
         real(observations, dp))
      out%degrees_of_freedom = moments
      out%p_value = regularized_gamma_q(0.5_dp*real(moments, dp), &
         0.5_dp*out%statistic)
   end function residual_moment_test

   pure function corrected_residual_moment_test(series, model, base_values, &
      transformation, lag, difference_step, conditional) result(out)
      !! Apply the parameter-estimation correction to a residual moment test.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(gmvarkit_model_t), intent(in) :: model !! Fitted full model.
      real(dp), intent(in) :: base_values(:, :) !! Moment vectors at the estimate.
      integer, intent(in) :: transformation !! One normality, two correlation, or three variance.
      integer, intent(in) :: lag !! Maximum lag used by serial moment transformations.
      real(dp), intent(in) :: difference_step !! Relative finite-difference step.
      logical, intent(in) :: conditional !! Use conditional likelihood scores.
      type(gmvarkit_hypothesis_test_t) :: out
      type(gmvarkit_score_t) :: score
      type(gmvarkit_quantile_residuals_t) :: plus_residuals, minus_residuals
      real(dp), allocatable :: parameters(:), shifted(:), plus_values(:, :)
      real(dp), allocatable :: minus_values(:, :), derivative(:, :), fisher(:, :)
      real(dp), allocatable :: inverse_fisher(:, :), psi(:, :), h(:, :), fg(:, :)
      real(dp), allocatable :: omega(:, :), inverse_omega(:, :), total(:)
      real(dp) :: increment
      integer :: parameter, parameter_count, observations, moments, info

      observations = size(base_values, 1)
      moments = size(base_values, 2)
      call pack_model(model, parameters, info)
      if (info /= 0 .or. observations < 1 .or. moments < 1) then
         out%info = 1
         return
      end if
      parameter_count = size(parameters)
      allocate(shifted(parameter_count), derivative(moments, parameter_count))
      do parameter = 1, parameter_count
         increment = difference_step*max(1.0_dp, abs(parameters(parameter)))
         shifted = parameters
         shifted(parameter) = shifted(parameter) + increment
         plus_residuals = gmvarkit_quantile_residuals(series, &
            unpack_model(shifted, model))
         shifted(parameter) = parameters(parameter) - increment
         minus_residuals = gmvarkit_quantile_residuals(series, &
            unpack_model(shifted, model))
         if (plus_residuals%info /= 0 .or. minus_residuals%info /= 0) then
            out%info = 2
            return
         end if
         call residual_moment_values(plus_residuals%residual, transformation, &
            lag, plus_values)
         call residual_moment_values(minus_residuals%residual, transformation, &
            lag, minus_values)
         if (any(shape(plus_values) /= shape(base_values)) .or. &
            any(shape(minus_values) /= shape(base_values))) then
            out%info = 3
            return
         end if
         derivative(:, parameter) = sum(plus_values - minus_values, dim=1)/ &
            (2.0_dp*increment*real(observations, dp))
      end do
      score = gmvarkit_score_matrix(series, model, conditional=conditional, &
         difference_step=difference_step)
      if (score%info /= 0) then
         out%info = 4
         return
      end if
      fisher = score%opg/real(size(score%observation, 1), dp)
      call invert_matrix(fisher, inverse_fisher, info)
      if (info /= 0) then
         out%info = 5
         return
      end if
      if (transformation == 1) then
         psi = matmul(transpose(base_values), score%observation)/ &
            real(observations, dp)
      else
         psi = matmul(transpose(base_values), score%observation(lag + 1:, :))/ &
            real(observations, dp)
      end if
      h = matmul(transpose(base_values), base_values)/real(observations, dp)
      fg = matmul(inverse_fisher, transpose(derivative))
      omega = h + matmul(derivative, fg) + matmul(psi, fg) + &
         matmul(matmul(derivative, inverse_fisher), transpose(psi))
      omega = 0.5_dp*(omega + transpose(omega))
      call invert_matrix(omega, inverse_omega, info)
      if (info /= 0) then
         out%info = 6
         return
      end if
      total = sum(base_values, dim=1)
      out%statistic = max(0.0_dp, dot_product(total, &
         matmul(inverse_omega, total))/real(observations, dp))
      out%degrees_of_freedom = moments
      out%p_value = regularized_gamma_q(0.5_dp*real(moments, dp), &
         0.5_dp*out%statistic)
   end function corrected_residual_moment_test

   pure subroutine residual_moment_values(residuals, transformation, lag, values)
      !! Construct normality, correlation, or squared-correlation moment vectors.
      real(dp), intent(in) :: residuals(:, :) !! Quantile residuals by time and variable.
      integer, intent(in) :: transformation !! One normality, two correlation, or three variance.
      integer, intent(in) :: lag !! Maximum lag for serial transformations.
      real(dp), allocatable, intent(out) :: values(:, :) !! Moment vectors by usable observation.
      real(dp), allocatable :: centered_square(:, :)
      integer :: observations, variables, time, delay, variable, other, position

      observations = size(residuals, 1)
      variables = size(residuals, 2)
      if (transformation == 1) then
         allocate(values(observations, 3*variables))
         do variable = 1, variables
            values(:, 3*variable - 2) = residuals(:, variable)**2 - 1.0_dp
            values(:, 3*variable - 1) = residuals(:, variable)**3
            values(:, 3*variable) = residuals(:, variable)**4 - 3.0_dp
         end do
         return
      end if
      allocate(values(observations - lag, lag*variables**2))
      if (transformation == 3) then
         allocate(centered_square(observations, variables))
         centered_square = residuals**2 - 1.0_dp
      end if
      do time = lag + 1, observations
         position = 0
         do delay = 1, lag
            do other = 1, variables
               do variable = 1, variables
                  position = position + 1
                  if (transformation == 2) then
                     values(time - lag, position) = residuals(time, variable)* &
                        residuals(time - delay, other)
                  else
                     values(time - lag, position) = &
                        centered_square(time, variable)* &
                        centered_square(time - delay, other)
                  end if
               end do
            end do
         end do
      end do
   end subroutine residual_moment_values

   pure real(dp) function student_covariance_cdf(value, variance, &
      degrees_of_freedom) result(probability)
      !! Evaluate a univariate Student-t CDF parameterized by its variance.
      real(dp), intent(in) :: value !! Centered evaluation value.
      real(dp), intent(in) :: variance !! Positive distribution variance.
      real(dp), intent(in) :: degrees_of_freedom !! Degrees of freedom above two.
      real(dp) :: standardized, beta_argument, beta_value

      if (variance <= 0.0_dp .or. degrees_of_freedom <= 2.0_dp) then
         probability = 0.5_dp
         return
      end if
      standardized = value/sqrt(variance*(degrees_of_freedom - 2.0_dp)/ &
         degrees_of_freedom)
      beta_argument = degrees_of_freedom/(degrees_of_freedom + standardized**2)
      beta_value = regularized_beta(beta_argument, &
         0.5_dp*degrees_of_freedom, 0.5_dp)
      if (standardized >= 0.0_dp) then
         probability = 1.0_dp - 0.5_dp*beta_value
      else
         probability = 0.5_dp*beta_value
      end if
   end function student_covariance_cdf

   pure real(dp) function gaussian_log_density(dimension, quadratic, &
      log_determinant) result(value)
      !! Evaluate a multivariate Gaussian log density from sufficient scalars.
      integer, intent(in) :: dimension !! Distribution dimension.
      real(dp), intent(in) :: quadratic !! Mahalanobis quadratic form.
      real(dp), intent(in) :: log_determinant !! Log determinant of covariance.

      value = -0.5_dp*(real(dimension, dp)*log(2.0_dp*acos(-1.0_dp)) + &
         log_determinant + quadratic)
   end function gaussian_log_density

   pure real(dp) function student_log_density(dimension, quadratic, &
      log_determinant, degrees_of_freedom) result(value)
      !! Evaluate a covariance-parameterized multivariate Student-t log density.
      integer, intent(in) :: dimension !! Distribution dimension.
      real(dp), intent(in) :: quadratic !! Mahalanobis quadratic form.
      real(dp), intent(in) :: log_determinant !! Log determinant of covariance.
      real(dp), intent(in) :: degrees_of_freedom !! Degrees of freedom above two.

      value = log_gamma(0.5_dp*(real(dimension, dp) + degrees_of_freedom)) - &
         log_gamma(0.5_dp*degrees_of_freedom) - &
         0.5_dp*real(dimension, dp)*log(acos(-1.0_dp)) - &
         0.5_dp*real(dimension, dp)*log(degrees_of_freedom - 2.0_dp) - &
         0.5_dp*log_determinant - &
         0.5_dp*(real(dimension, dp) + degrees_of_freedom)* &
         log(1.0_dp + quadratic/(degrees_of_freedom - 2.0_dp))
   end function student_log_density

end module gmvarkit_mod
