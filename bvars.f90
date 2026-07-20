! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Distinct algorithms translated from the R bvars package.
module bvars_mod
   !! Bayesian large-VAR estimation with common error scaling.
   use kind_mod, only: dp
   use linalg_mod, only: invert_matrix, cholesky_lower, identity_matrix
   use random_mod, only: random_standard_normal, random_standard_normal_matrix, &
      random_gamma, random_uniform, multivariate_normal_from_standard
   use mixture_mod, only: gaussian_mixture_t, fit_gaussian_mixture
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use time_series_var_utils_mod, only: build_var_data
   use bvartools_mod, only: bvartools_bvar_draws_t, bvartools_fevd_t, &
      bvartools_bvar_fevd
   implicit none
   private

   type, public :: bvars_data_t
      !! Aligned BVAR response, regressors, and source series.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: series(:, :)
      integer :: variables = 0
      integer :: observations = 0
      integer :: regressors = 0
      integer :: lags = 0
      integer :: exogenous_count = 0
      integer :: info = 0
   end type bvars_data_t

   type, public :: bvars_prior_t
      !! Matrix-normal and inverse-Wishart BVAR prior.
      real(dp), allocatable :: coefficient_mean(:, :)
      real(dp), allocatable :: coefficient_covariance(:, :)
      real(dp), allocatable :: covariance_scale(:, :)
      real(dp) :: covariance_degrees_of_freedom = 0.0_dp
      integer :: info = 0
   end type bvars_prior_t

   type, public :: bvars_draws_t
      !! Retained coefficient and covariance draws.
      real(dp), allocatable :: coefficient(:, :, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: coefficient_covariance(:, :, :)
      real(dp), allocatable :: common_variance(:, :)
      real(dp), allocatable :: scale_mixture(:, :)
      real(dp), allocatable :: degrees_of_freedom(:)
      logical, allocatable :: degrees_accepted(:)
      real(dp), allocatable :: log_variance(:, :)
      integer, allocatable :: mixture_component(:, :)
      real(dp), allocatable :: persistence(:)
      real(dp), allocatable :: volatility_loading(:)
      real(dp), allocatable :: volatility_innovation_variance(:)
      real(dp), allocatable :: loading_prior_variance(:)
      real(dp), allocatable :: volatility_hyper_scale(:)
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvars_draws_t

   type, public :: bvars_predictive_t
      !! Draw-wise forecast paths, means, and covariances.
      real(dp), allocatable :: path(:, :, :)
      real(dp), allocatable :: mean(:, :, :)
      real(dp), allocatable :: covariance(:, :, :, :)
      integer :: info = 0
   end type bvars_predictive_t

   type, public :: bvars_common_sv_t
      !! Common stochastic-volatility state and parameter draw.
      real(dp), allocatable :: log_variance(:)
      integer, allocatable :: mixture_component(:)
      real(dp) :: persistence = 0.0_dp
      real(dp) :: loading = 0.0_dp
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: loading_prior_variance = 0.0_dp
      real(dp) :: hyper_scale = 0.0_dp
      integer :: info = 0
   end type bvars_common_sv_t

   type, public :: bvars_fitted_t
      !! Posterior fitted-density draws.
      real(dp), allocatable :: value(:, :, :)
      integer :: info = 0
   end type bvars_fitted_t

   public :: bvars_prepare
   public :: bvars_default_prior
   public :: bvars_matrix_normal_from_standard, bvars_matrix_normal
   public :: bvars_conjugate_draws
   public :: bvars_student_log_kernel, bvars_student_t_draws
   public :: bvars_student_scale_forecast
   public :: bvars_sv_auxiliary_mixture
   public :: bvars_centered_sv_update, bvars_centered_sv_draws
   public :: bvars_noncentered_sv_update, bvars_noncentered_sv_draws
   public :: bvars_student_sv_draws
   public :: bvars_common_variance_forecast
   public :: bvars_shocks
   public :: bvars_fitted_from_random, bvars_fitted
   public :: bvars_forecast_from_random, bvars_forecast
   public :: bvars_conditional_forecast_from_random, bvars_conditional_forecast
   public :: bvars_fevd

contains

   pure function bvars_prepare(series, lags, exogenous) result(out)
      !! Construct the lag-major BVAR response and design matrices.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      integer, intent(in) :: lags !! Positive autoregressive lag order.
      real(dp), intent(in), optional :: exogenous(:, :) !! Contemporaneous exogenous observations.
      type(bvars_data_t) :: out
      real(dp), allocatable :: response(:, :), lag_design(:, :)
      integer :: rows, variables, exogenous_count

      rows = size(series, 1)
      variables = size(series, 2)
      exogenous_count = 0
      if (present(exogenous)) exogenous_count = size(exogenous, 2)
      if (variables < 1 .or. lags < 1 .or. rows <= lags .or. &
         (present(exogenous) .and. size(exogenous, 1) /= rows)) then
         out%info = 1
         return
      end if
      call build_var_data(series, lags, response, lag_design)
      out%observations = rows - lags
      out%variables = variables
      out%lags = lags
      out%exogenous_count = exogenous_count
      out%regressors = variables*lags + 1 + exogenous_count
      allocate(out%series(variables, rows))
      allocate(out%y(variables, out%observations))
      allocate(out%x(out%regressors, out%observations))
      out%series = transpose(series)
      out%y = transpose(response)
      out%x = 0.0_dp
      out%x(1:variables*lags, :) = transpose(lag_design)
      out%x(variables*lags + 1, :) = 1.0_dp
      if (present(exogenous)) then
         out%x(variables*lags + 2:, :) = transpose(exogenous(lags + 1:, :))
      end if
   end function bvars_prepare

   pure function bvars_default_prior(data, stationary, innovation_variance, &
      tightness, deterministic_variance) result(out)
      !! Construct the default bvars matrix-normal and inverse-Wishart prior.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      logical, intent(in), optional :: stationary(:) !! Equations centred on white noise instead of random walks.
      real(dp), intent(in), optional :: innovation_variance(:) !! Marginal innovation variance estimates.
      real(dp), intent(in), optional :: tightness !! Overall Minnesota variance scale.
      real(dp), intent(in), optional :: deterministic_variance !! Constant and exogenous variance scale.
      type(bvars_prior_t) :: out
      real(dp), allocatable :: scale(:)
      real(dp) :: lambda, deterministic_scale
      integer :: variable, lag, position

      lambda = 0.2_dp**2
      deterministic_scale = 10.0_dp**2
      if (present(tightness)) lambda = tightness
      if (present(deterministic_variance)) deterministic_scale = &
         deterministic_variance
      if (data%info /= 0 .or. data%variables < 1 .or. lambda <= 0.0_dp .or. &
         deterministic_scale <= 0.0_dp .or. &
         (present(stationary) .and. size(stationary) /= data%variables) .or. &
         (present(innovation_variance) .and. &
         (size(innovation_variance) /= data%variables .or. &
         any(innovation_variance <= 0.0_dp)))) then
         out%info = 1
         return
      end if
      allocate(scale(data%variables))
      if (present(innovation_variance)) then
         scale = innovation_variance
      else
         do variable = 1, data%variables
            scale(variable) = sum((data%series(variable, 2:) - &
               data%series(variable, :size(data%series, 2) - 1))**2)/ &
               real(size(data%series, 2) - 1, dp)
         end do
         scale = max(scale, sqrt(epsilon(1.0_dp)))
      end if
      allocate(out%coefficient_mean(data%variables, data%regressors))
      allocate(out%coefficient_covariance(data%regressors, data%regressors))
      allocate(out%covariance_scale(data%variables, data%variables))
      out%coefficient_mean = 0.0_dp
      do variable = 1, data%variables
         if (.not. present(stationary)) then
            out%coefficient_mean(variable, variable) = 1.0_dp
         else if (.not. stationary(variable)) then
            out%coefficient_mean(variable, variable) = 1.0_dp
         end if
      end do
      out%coefficient_covariance = 0.0_dp
      do lag = 1, data%lags
         do variable = 1, data%variables
            position = (lag - 1)*data%variables + variable
            out%coefficient_covariance(position, position) = &
               lambda/(scale(variable)*real(lag*lag, dp))
         end do
      end do
      do position = data%variables*data%lags + 1, data%regressors
         out%coefficient_covariance(position, position) = deterministic_scale
      end do
      out%covariance_scale = identity_matrix(data%variables)
      out%covariance_degrees_of_freedom = real(data%variables + 3, dp)
   end function bvars_default_prior

   pure function bvars_matrix_normal_from_standard(mean, row_covariance, &
      column_covariance, standard_normal) result(draw)
      !! Transform standard normals to a matrix-normal draw.
      real(dp), intent(in) :: mean(:, :) !! Matrix-normal mean.
      real(dp), intent(in) :: row_covariance(:, :) !! Row covariance.
      real(dp), intent(in) :: column_covariance(:, :) !! Column covariance.
      real(dp), intent(in) :: standard_normal(:, :) !! Independent standard normals.
      real(dp), allocatable :: draw(:, :)
      real(dp), allocatable :: row_factor(:, :), column_factor(:, :)
      integer :: info

      if (any(shape(row_covariance) /= [size(mean, 1), size(mean, 1)]) .or. &
         any(shape(column_covariance) /= [size(mean, 2), size(mean, 2)]) .or. &
         any(shape(standard_normal) /= shape(mean))) then
         allocate(draw(0, 0))
         return
      end if
      call cholesky_lower(row_covariance, row_factor, info)
      if (info /= 0) then
         allocate(draw(0, 0))
         return
      end if
      call cholesky_lower(column_covariance, column_factor, info)
      if (info /= 0) then
         allocate(draw(0, 0))
         return
      end if
      draw = mean + matmul(row_factor, &
         matmul(standard_normal, transpose(column_factor)))
   end function bvars_matrix_normal_from_standard

   function bvars_matrix_normal(mean, row_covariance, column_covariance) &
      result(draw)
      !! Generate a matrix-normal draw using the shared RNG.
      real(dp), intent(in) :: mean(:, :) !! Matrix-normal mean.
      real(dp), intent(in) :: row_covariance(:, :) !! Row covariance.
      real(dp), intent(in) :: column_covariance(:, :) !! Column covariance.
      real(dp), allocatable :: draw(:, :)
      real(dp), allocatable :: standard_normal(:, :)

      allocate(standard_normal(size(mean, 1), size(mean, 2)))
      call random_standard_normal_matrix(standard_normal)
      draw = bvars_matrix_normal_from_standard(mean, row_covariance, &
         column_covariance, standard_normal)
   end function bvars_matrix_normal

   function bvars_conjugate_draws(data, prior, draws, burnin, observation_scale) &
      result(out)
      !! Sample the conditional matrix-normal/inverse-Wishart posterior.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      type(bvars_prior_t), intent(in) :: prior !! BVAR prior specification.
      integer, intent(in) :: draws !! Number of retained posterior draws.
      integer, intent(in), optional :: burnin !! Number of discarded posterior draws.
      real(dp), intent(in), optional :: observation_scale(:) !! Common variance multiplier by observation.
      type(bvars_draws_t) :: out
      real(dp), allocatable :: weight(:), weighted_y(:, :), weighted_x(:, :)
      real(dp), allocatable :: prior_precision(:, :), posterior_precision(:, :)
      real(dp), allocatable :: posterior_covariance(:, :), posterior_mean(:, :)
      real(dp), allocatable :: posterior_scale(:, :), covariance(:, :), precision(:, :)
      integer :: total, discard, iteration, retained, info

      discard = 0
      if (present(burnin)) discard = burnin
      if (data%info /= 0 .or. prior%info /= 0 .or. draws < 1 .or. discard < 0 .or. &
         any(shape(prior%coefficient_mean) /= [data%variables, data%regressors]) .or. &
         any(shape(prior%coefficient_covariance) /= [data%regressors, data%regressors]) .or. &
         any(shape(prior%covariance_scale) /= [data%variables, data%variables]) .or. &
         prior%covariance_degrees_of_freedom <= real(data%variables - 1, dp) .or. &
         (present(observation_scale) .and. &
         (size(observation_scale) /= data%observations .or. &
         any(observation_scale <= 0.0_dp)))) then
         out%info = 1
         return
      end if
      allocate(weight(data%observations))
      weight = 1.0_dp
      if (present(observation_scale)) weight = 1.0_dp/sqrt(observation_scale)
      weighted_y = data%y*spread(weight, 1, data%variables)
      weighted_x = data%x*spread(weight, 1, data%regressors)
      call invert_matrix(prior%coefficient_covariance, prior_precision, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      posterior_precision = matmul(weighted_x, transpose(weighted_x)) + &
         prior_precision
      call invert_matrix(posterior_precision, posterior_covariance, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      posterior_mean = matmul(matmul(weighted_y, transpose(weighted_x)) + &
         matmul(prior%coefficient_mean, prior_precision), posterior_covariance)
      posterior_scale = prior%covariance_scale + &
         matmul(weighted_y, transpose(weighted_y)) + &
         matmul(matmul(prior%coefficient_mean, prior_precision), &
         transpose(prior%coefficient_mean)) - &
         matmul(matmul(posterior_mean, posterior_precision), &
         transpose(posterior_mean))
      posterior_scale = 0.5_dp*(posterior_scale + transpose(posterior_scale))
      allocate(out%coefficient(data%variables, data%regressors, draws))
      allocate(out%covariance(data%variables, data%variables, draws))
      allocate(out%coefficient_covariance(data%regressors, data%regressors, draws))
      allocate(out%common_variance(data%observations, draws))
      allocate(out%scale_mixture(data%observations, draws))
      out%common_variance = 1.0_dp
      out%scale_mixture = 1.0_dp
      total = draws + discard
      retained = 0
      do iteration = 1, total
         call inverse_wishart_draw(posterior_scale, &
            prior%covariance_degrees_of_freedom + real(data%observations, dp), &
            covariance, precision, info)
         if (info /= 0) then
            out%info = 10 + info
            return
         end if
         if (iteration > discard) then
            retained = retained + 1
            out%coefficient(:, :, retained) = bvars_matrix_normal( &
               posterior_mean, covariance, posterior_covariance)
            out%covariance(:, :, retained) = covariance
            out%coefficient_covariance(:, :, retained) = &
               prior%coefficient_covariance
         end if
      end do
      out%retained_draws = retained
   end function bvars_conjugate_draws

   pure real(dp) function bvars_student_log_kernel(degrees_of_freedom, &
      scale_mixture) result(log_kernel)
      !! Evaluate the bvars Student-t degrees-of-freedom posterior kernel.
      real(dp), intent(in) :: degrees_of_freedom !! Degrees of freedom greater than two.
      real(dp), intent(in) :: scale_mixture(:) !! Positive latent variance multipliers.
      integer :: observations

      if (degrees_of_freedom <= 2.0_dp .or. size(scale_mixture) < 1 .or. &
         any(scale_mixture <= 0.0_dp)) then
         log_kernel = -huge(1.0_dp)
         return
      end if
      observations = size(scale_mixture)
      log_kernel = -real(observations, dp)* &
         log_gamma(0.5_dp*degrees_of_freedom) + &
         0.5_dp*real(observations, dp)*degrees_of_freedom* &
         log(0.5_dp*(degrees_of_freedom - 2.0_dp)) - &
         0.5_dp*(degrees_of_freedom + 2.0_dp)*sum(log(scale_mixture)) - &
         0.5_dp*(degrees_of_freedom - 2.0_dp)*sum(1.0_dp/scale_mixture) - &
         2.0_dp*log(degrees_of_freedom - 1.0_dp)
   end function bvars_student_log_kernel

   function bvars_student_t_draws(data, prior, draws, burnin, initial_df, &
      initial_proposal_scale) result(out)
      !! Sample the homoskedastic BVAR with multivariate Student-t errors.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      type(bvars_prior_t), intent(in) :: prior !! BVAR prior specification.
      integer, intent(in) :: draws !! Number of retained posterior draws.
      integer, intent(in), optional :: burnin !! Number of discarded posterior draws.
      real(dp), intent(in), optional :: initial_df !! Initial degrees of freedom greater than two.
      real(dp), intent(in), optional :: initial_proposal_scale !! Initial transformed random-walk variance.
      type(bvars_draws_t) :: out
      real(dp), allocatable :: coefficient(:, :), covariance(:, :), precision(:, :)
      real(dp), allocatable :: residual(:, :), quadratic(:), mixture(:), weight(:)
      real(dp), allocatable :: weighted_y(:, :), weighted_x(:, :)
      real(dp), allocatable :: prior_precision(:, :), posterior_precision(:, :)
      real(dp), allocatable :: posterior_covariance(:, :), posterior_mean(:, :)
      real(dp), allocatable :: posterior_scale(:, :)
      real(dp) :: degrees, candidate, eta, candidate_eta, proposal_scale
      real(dp) :: old_kernel, new_kernel, log_ratio, acceptance
      integer :: discard, total, iteration, retained, time, info
      logical :: accepted

      discard = 0
      degrees = 30.0_dp
      proposal_scale = 0.05_dp
      if (present(burnin)) discard = burnin
      if (present(initial_df)) degrees = initial_df
      if (present(initial_proposal_scale)) proposal_scale = initial_proposal_scale
      if (data%info /= 0 .or. prior%info /= 0 .or. draws < 1 .or. discard < 0 .or. &
         degrees <= 2.0_dp .or. proposal_scale <= 0.0_dp .or. &
         any(shape(prior%coefficient_mean) /= [data%variables, data%regressors]) .or. &
         any(shape(prior%coefficient_covariance) /= &
         [data%regressors, data%regressors]) .or. &
         any(shape(prior%covariance_scale) /= [data%variables, data%variables]) .or. &
         prior%covariance_degrees_of_freedom <= real(data%variables + 1, dp)) then
         out%info = 1
         return
      end if
      call invert_matrix(prior%coefficient_covariance, prior_precision, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      coefficient = prior%coefficient_mean
      covariance = prior%covariance_scale/ &
         (prior%covariance_degrees_of_freedom - real(data%variables + 1, dp))
      call invert_matrix(covariance, precision, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      allocate(mixture(data%observations), quadratic(data%observations))
      allocate(weight(data%observations))
      mixture = 1.0_dp
      allocate(out%coefficient(data%variables, data%regressors, draws))
      allocate(out%covariance(data%variables, data%variables, draws))
      allocate(out%coefficient_covariance(data%regressors, data%regressors, draws))
      allocate(out%common_variance(data%observations, draws))
      allocate(out%scale_mixture(data%observations, draws))
      allocate(out%degrees_of_freedom(draws), out%degrees_accepted(draws))
      out%common_variance = 1.0_dp
      total = draws + discard
      retained = 0
      do iteration = 1, total
         eta = log(degrees - 2.0_dp)
         candidate_eta = eta + sqrt(proposal_scale)*random_standard_normal()
         candidate = 2.0_dp + exp(candidate_eta)
         old_kernel = bvars_student_log_kernel(degrees, mixture) + eta
         new_kernel = bvars_student_log_kernel(candidate, mixture) + candidate_eta
         log_ratio = min(0.0_dp, new_kernel - old_kernel)
         acceptance = exp(log_ratio)
         accepted = log(max(random_uniform(), tiny(1.0_dp))) < log_ratio
         if (accepted) degrees = candidate
         if (iteration > 1) proposal_scale = exp(log(proposal_scale) + &
            0.5_dp*log(1.0_dp + real(iteration, dp)**(-0.6_dp)* &
            (acceptance - 0.44_dp)))

         residual = data%y - matmul(coefficient, data%x)
         do time = 1, data%observations
            quadratic(time) = dot_product(residual(:, time), &
               matmul(precision, residual(:, time)))/real(data%variables, dp)
         end do
         if (sum(quadratic) <= tiny(1.0_dp)) then
            quadratic = 1.0_dp
         else
            quadratic = quadratic/(sum(quadratic)/real(data%observations, dp))
         end if
         do time = 1, data%observations
            mixture(time) = (quadratic(time) + degrees - 2.0_dp)/ &
               (2.0_dp*random_gamma(0.5_dp*(degrees + real(data%variables, dp))))
         end do
         weight = 1.0_dp/sqrt(mixture)
         weighted_y = data%y*spread(weight, 1, data%variables)
         weighted_x = data%x*spread(weight, 1, data%regressors)
         posterior_precision = matmul(weighted_x, transpose(weighted_x)) + &
            prior_precision
         call invert_matrix(posterior_precision, posterior_covariance, info)
         if (info /= 0) then
            out%info = 10
            return
         end if
         posterior_mean = matmul(matmul(weighted_y, transpose(weighted_x)) + &
            matmul(prior%coefficient_mean, prior_precision), posterior_covariance)
         posterior_scale = prior%covariance_scale + &
            matmul(weighted_y, transpose(weighted_y)) + &
            matmul(matmul(prior%coefficient_mean, prior_precision), &
            transpose(prior%coefficient_mean)) - &
            matmul(matmul(posterior_mean, posterior_precision), &
            transpose(posterior_mean))
         posterior_scale = 0.5_dp*(posterior_scale + transpose(posterior_scale))
         call inverse_wishart_draw(posterior_scale, &
            prior%covariance_degrees_of_freedom + real(data%observations, dp), &
            covariance, precision, info)
         if (info /= 0) then
            out%info = 11
            return
         end if
         coefficient = bvars_matrix_normal(posterior_mean, covariance, &
            posterior_covariance)
         if (iteration > discard) then
            retained = retained + 1
            out%coefficient(:, :, retained) = coefficient
            out%covariance(:, :, retained) = covariance
            out%coefficient_covariance(:, :, retained) = &
               prior%coefficient_covariance
            out%scale_mixture(:, retained) = mixture
            out%degrees_of_freedom(retained) = degrees
            out%degrees_accepted(retained) = accepted
         end if
      end do
      out%retained_draws = retained
   end function bvars_student_t_draws

   function bvars_student_scale_forecast(degrees_of_freedom, horizon) &
      result(scale_mixture)
      !! Simulate future standardized Student-t variance multipliers.
      real(dp), intent(in) :: degrees_of_freedom(:) !! Draw-specific degrees of freedom.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), allocatable :: scale_mixture(:, :)
      integer :: draw, step

      if (horizon < 1 .or. size(degrees_of_freedom) < 1 .or. &
         any(degrees_of_freedom <= 2.0_dp)) then
         allocate(scale_mixture(0, 0))
         return
      end if
      allocate(scale_mixture(horizon, size(degrees_of_freedom)))
      do draw = 1, size(degrees_of_freedom)
         do step = 1, horizon
            scale_mixture(step, draw) = (degrees_of_freedom(draw) - 2.0_dp)/ &
               (2.0_dp*random_gamma(0.5_dp*degrees_of_freedom(draw)))
         end do
      end do
   end function bvars_student_scale_forecast

   function bvars_sv_auxiliary_mixture(variables, sample_count, components, &
      max_iterations, tolerance) result(mixture)
      !! Construct the common-SV auxiliary Gaussian-mixture approximation.
      integer, intent(in) :: variables !! Positive number of jointly standardized variables.
      integer, intent(in), optional :: sample_count !! Number of simulated log-chi-square averages.
      integer, intent(in), optional :: components !! Positive Gaussian-mixture component count.
      integer, intent(in), optional :: max_iterations !! Maximum mixture EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative mixture EM convergence tolerance.
      real(dp), allocatable :: mixture(:, :)
      type(gaussian_mixture_t) :: fit
      real(dp), allocatable :: simulated(:)
      real(dp) :: tolerance_value, standard_normal
      integer :: count_value, components_value, iterations_value
      integer :: sample, variable

      count_value = 100000
      components_value = 10
      iterations_value = 200
      tolerance_value = 1.0e-8_dp
      if (present(sample_count)) count_value = sample_count
      if (present(components)) components_value = components
      if (present(max_iterations)) iterations_value = max_iterations
      if (present(tolerance)) tolerance_value = tolerance
      if (variables < 1 .or. count_value < 2 .or. components_value < 1 .or. &
         components_value > count_value .or. iterations_value < 1 .or. &
         tolerance_value <= 0.0_dp) then
         allocate(mixture(0, 0))
         return
      end if
      allocate(simulated(count_value))
      simulated = 0.0_dp
      do variable = 1, variables
         do sample = 1, count_value
            standard_normal = random_standard_normal()
            simulated(sample) = simulated(sample) + &
               log(max(standard_normal**2, tiny(1.0_dp)))/real(variables, dp)
         end do
      end do
      fit = fit_gaussian_mixture(simulated, components_value, iterations_value, &
         tolerance_value)
      if (fit%info /= 0) then
         allocate(mixture(0, 0))
         return
      end if
      allocate(mixture(3, components_value))
      mixture(1, :) = fit%weight
      mixture(2, :) = fit%mean
      mixture(3, :) = fit%variance
   end function bvars_sv_auxiliary_mixture

   function bvars_centered_sv_update(log_observation, current_log_variance, &
      current_persistence, current_innovation_variance, &
      current_loading_prior_variance, current_hyper_scale, mixture, &
      prior_shape, prior_scale) result(out)
      !! Perform one centred common-SV Gibbs update from bvars.
      real(dp), intent(in) :: log_observation(:) !! Aggregated log-squared standardized residuals.
      real(dp), intent(in) :: current_log_variance(:) !! Current centred log-variance states.
      real(dp), intent(in) :: current_persistence !! Current AR(1) persistence.
      real(dp), intent(in) :: current_innovation_variance !! Current state innovation variance.
      real(dp), intent(in) :: current_loading_prior_variance !! Current hierarchical loading variance.
      real(dp), intent(in) :: current_hyper_scale !! Current hierarchical scale.
      real(dp), intent(in) :: mixture(:, :) !! Mixture probability scores, means, and variances.
      real(dp), intent(in) :: prior_shape !! Positive hierarchical prior shape.
      real(dp), intent(in) :: prior_scale !! Positive hierarchical prior scale.
      type(bvars_common_sv_t) :: out
      real(dp), allocatable :: mixture_mean(:), mixture_precision(:)
      real(dp), allocatable :: transition_difference(:, :), state_precision(:, :)
      real(dp), allocatable :: state_covariance(:, :), state_mean(:), normal(:)
      real(dp) :: log_weight, maximum_weight, total_weight, uniform
      real(dp) :: persistence_variance, persistence_mean, sum_lagged
      integer :: observations, components, time, component, selected, info

      observations = size(log_observation)
      components = size(mixture, 2)
      if (observations < 2 .or. size(current_log_variance) /= observations .or. &
         any(shape(mixture) /= [3, components]) .or. components < 1 .or. &
         any(mixture(3, :) <= 0.0_dp) .or. abs(current_persistence) >= 1.0_dp .or. &
         current_innovation_variance <= 0.0_dp .or. &
         current_loading_prior_variance <= 0.0_dp .or. &
         current_hyper_scale <= 0.0_dp .or. prior_shape <= 0.0_dp .or. &
         prior_scale <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%mixture_component(observations))
      allocate(mixture_mean(observations), mixture_precision(observations))
      do time = 1, observations
         maximum_weight = -huge(1.0_dp)
         do component = 1, components
            log_weight = mixture(1, component) - &
               (log_observation(time) - current_log_variance(time) - &
               mixture(2, component))**2/mixture(3, component)
            maximum_weight = max(maximum_weight, log_weight)
         end do
         total_weight = 0.0_dp
         do component = 1, components
            total_weight = total_weight + exp(mixture(1, component) - &
               (log_observation(time) - current_log_variance(time) - &
               mixture(2, component))**2/mixture(3, component) - maximum_weight)
         end do
         uniform = random_uniform()*total_weight
         total_weight = 0.0_dp
         selected = components
         do component = 1, components
            total_weight = total_weight + exp(mixture(1, component) - &
               (log_observation(time) - current_log_variance(time) - &
               mixture(2, component))**2/mixture(3, component) - maximum_weight)
            if (uniform <= total_weight) then
               selected = component
               exit
            end if
         end do
         out%mixture_component(time) = selected
         mixture_mean(time) = mixture(2, selected)
         mixture_precision(time) = 1.0_dp/mixture(3, selected)
      end do

      out%hyper_scale = (1.0_dp + 2.0_dp*current_loading_prior_variance)/ &
         (2.0_dp*random_gamma(0.5_dp*(3.0_dp + 2.0_dp*prior_shape)))
      out%loading_prior_variance = random_gamma(1.0_dp + 0.5_dp*prior_shape, &
         1.0_dp/(1.0_dp/prior_scale + 1.0_dp/ &
         (2.0_dp*current_innovation_variance)))
      sum_lagged = dot_product(current_log_variance(:observations - 1), &
         current_log_variance(:observations - 1))
      if (sum_lagged <= tiny(1.0_dp)) then
         persistence_mean = 0.0_dp
         persistence_variance = 1.0_dp
      else
         persistence_variance = current_innovation_variance/sum_lagged
         persistence_mean = dot_product(current_log_variance(:observations - 1), &
            current_log_variance(2:))/sum_lagged
      end if
      out%persistence = truncated_normal_random(persistence_mean, &
         sqrt(persistence_variance), -1.0_dp, 1.0_dp)
      allocate(transition_difference(observations, observations))
      transition_difference = identity_matrix(observations)
      do time = 2, observations
         transition_difference(time, time - 1) = -out%persistence
      end do
      out%innovation_variance = (out%loading_prior_variance + &
         sum(matmul(transition_difference, current_log_variance)**2))/ &
         (2.0_dp*random_gamma(0.5_dp*(3.0_dp + real(observations, dp))))
      out%loading = sqrt(out%innovation_variance)
      state_precision = matmul(transpose(transition_difference), &
         transition_difference)/out%innovation_variance
      do time = 1, observations
         state_precision(time, time) = state_precision(time, time) + &
            mixture_precision(time)
      end do
      call invert_matrix(state_precision, state_covariance, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      state_mean = matmul(state_covariance, mixture_precision* &
         (log_observation - mixture_mean))
      allocate(normal(observations))
      do time = 1, observations
         normal(time) = random_standard_normal()
      end do
      allocate(out%log_variance(observations))
      call multivariate_normal_from_standard(state_mean, state_covariance, &
         normal, out%log_variance, info)
      if (info /= 0) out%info = 3
   end function bvars_centered_sv_update

   function bvars_noncentered_sv_update(log_observation, current_state, &
      current_persistence, current_loading, current_innovation_variance, &
      current_loading_prior_variance, current_hyper_scale, mixture, &
      prior_shape, prior_scale) result(out)
      !! Perform one non-centred common-SV update with ASIS interweaving.
      real(dp), intent(in) :: log_observation(:) !! Aggregated log-squared standardized residuals.
      real(dp), intent(in) :: current_state(:) !! Current unit-innovation non-centred states.
      real(dp), intent(in) :: current_persistence !! Current AR(1) persistence.
      real(dp), intent(in) :: current_loading !! Current signed volatility loading.
      real(dp), intent(in) :: current_innovation_variance !! Current centred innovation variance.
      real(dp), intent(in) :: current_loading_prior_variance !! Current loading prior variance.
      real(dp), intent(in) :: current_hyper_scale !! Current hierarchical scale.
      real(dp), intent(in) :: mixture(:, :) !! Mixture probability scores, means, and variances.
      real(dp), intent(in) :: prior_shape !! Positive hierarchical prior shape.
      real(dp), intent(in) :: prior_scale !! Positive hierarchical prior scale.
      type(bvars_common_sv_t) :: out
      real(dp), allocatable :: mixture_mean(:), mixture_precision(:)
      real(dp), allocatable :: difference(:, :), precision(:, :), covariance(:, :)
      real(dp), allocatable :: mean(:), normal(:), state(:), centred_state(:)
      real(dp) :: maximum_weight, total_weight, log_weight, uniform
      real(dp) :: persistence_variance, persistence_mean, bound
      real(dp) :: loading_variance, loading_mean, quadratic
      integer :: observations, components, time, component, selected, info, sign_value

      observations = size(log_observation)
      components = size(mixture, 2)
      if (observations < 2 .or. size(current_state) /= observations .or. &
         size(mixture, 1) /= 3 .or. components < 1 .or. &
         any(mixture(3, :) <= 0.0_dp) .or. abs(current_persistence) >= 1.0_dp .or. &
         current_innovation_variance <= 0.0_dp .or. &
         current_loading_prior_variance <= 0.0_dp .or. current_hyper_scale <= 0.0_dp .or. &
         prior_shape <= 0.0_dp .or. prior_scale <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%mixture_component(observations), mixture_mean(observations))
      allocate(mixture_precision(observations))
      do time = 1, observations
         maximum_weight = -huge(1.0_dp)
         do component = 1, components
            log_weight = mixture(1, component) - (log_observation(time) - &
               current_loading*current_state(time) - mixture(2, component))**2/ &
               mixture(3, component)
            maximum_weight = max(maximum_weight, log_weight)
         end do
         total_weight = 0.0_dp
         do component = 1, components
            total_weight = total_weight + exp(mixture(1, component) - &
               (log_observation(time) - current_loading*current_state(time) - &
               mixture(2, component))**2/mixture(3, component) - maximum_weight)
         end do
         uniform = random_uniform()*total_weight
         total_weight = 0.0_dp
         selected = components
         do component = 1, components
            total_weight = total_weight + exp(mixture(1, component) - &
               (log_observation(time) - current_loading*current_state(time) - &
               mixture(2, component))**2/mixture(3, component) - maximum_weight)
            if (uniform <= total_weight) then
               selected = component
               exit
            end if
         end do
         out%mixture_component(time) = selected
         mixture_mean(time) = mixture(2, selected)
         mixture_precision(time) = 1.0_dp/mixture(3, selected)
      end do
      out%hyper_scale = (prior_scale + 2.0_dp*current_loading_prior_variance)/ &
         (2.0_dp*random_gamma(0.5_dp*(3.0_dp + 2.0_dp*prior_shape)))
      out%loading_prior_variance = gig_slice_draw(prior_shape - 0.5_dp, &
         current_loading**2, 2.0_dp/out%hyper_scale, &
         current_loading_prior_variance)
      bound = sqrt(max(epsilon(1.0_dp), 1.0_dp - out%loading_prior_variance))
      quadratic = dot_product(current_state(:observations - 1), &
         current_state(:observations - 1))
      persistence_variance = 1.0_dp/max(quadratic, tiny(1.0_dp))
      persistence_mean = persistence_variance* &
         dot_product(current_state(:observations - 1), current_state(2:))
      out%persistence = truncated_normal_random(persistence_mean, &
         sqrt(persistence_variance), -bound, bound)
      allocate(difference(observations, observations))
      difference = identity_matrix(observations)
      do time = 2, observations
         difference(time, time - 1) = -out%persistence
      end do
      loading_variance = 1.0_dp/(sum(mixture_precision*current_state**2) + &
         1.0_dp/out%loading_prior_variance)
      loading_mean = loading_variance*sum(current_state*mixture_precision* &
         (log_observation - mixture_mean))
      out%loading = loading_mean + sqrt(loading_variance)*random_standard_normal()
      precision = out%loading**2*diagonal_from_vector(mixture_precision) + &
         matmul(transpose(difference), difference)
      call invert_matrix(precision, covariance, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      mean = matmul(covariance, out%loading*mixture_precision* &
         (log_observation - mixture_mean))
      allocate(normal(observations), state(observations))
      do time = 1, observations
         normal(time) = random_standard_normal()
      end do
      call multivariate_normal_from_standard(mean, covariance, normal, state, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      centred_state = out%loading*state
      quadratic = sum(matmul(difference, centred_state)**2)
      out%innovation_variance = gig_slice_draw(-0.5_dp*real(observations - 1, dp), &
         quadratic, 1.0_dp/out%loading_prior_variance, current_innovation_variance)
      sign_value = 1
      if (random_uniform() < 0.5_dp) sign_value = -1
      out%loading = real(sign_value, dp)*sqrt(out%innovation_variance)
      state = centred_state/out%loading
      quadratic = dot_product(state(:observations - 1), state(:observations - 1))
      persistence_variance = 1.0_dp/max(quadratic, tiny(1.0_dp))
      persistence_mean = persistence_variance* &
         dot_product(state(:observations - 1), state(2:))
      bound = sqrt(max(epsilon(1.0_dp), 1.0_dp - out%loading_prior_variance))
      out%persistence = truncated_normal_random(persistence_mean, &
         sqrt(persistence_variance), -bound, bound)
      out%log_variance = out%loading*state
   end function bvars_noncentered_sv_update

   function bvars_centered_sv_draws(data, prior, mixture, draws, burnin, &
      prior_shape, prior_scale) result(out)
      !! Sample the BVAR with centred common stochastic volatility.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      type(bvars_prior_t), intent(in) :: prior !! BVAR prior specification.
      real(dp), intent(in) :: mixture(:, :) !! Mixture probability scores, means, and variances.
      integer, intent(in) :: draws !! Number of retained posterior draws.
      integer, intent(in), optional :: burnin !! Number of discarded posterior draws.
      real(dp), intent(in), optional :: prior_shape !! Hierarchical volatility prior shape.
      real(dp), intent(in), optional :: prior_scale !! Hierarchical volatility prior scale.
      type(bvars_draws_t) :: out
      type(bvars_common_sv_t) :: volatility
      real(dp), allocatable :: coefficient(:, :), covariance(:, :), precision(:, :)
      real(dp), allocatable :: residual(:, :), lower(:, :), lower_inverse(:, :)
      real(dp), allocatable :: standardized(:, :), log_observation(:), weight(:)
      real(dp), allocatable :: weighted_y(:, :), weighted_x(:, :)
      real(dp), allocatable :: prior_precision(:, :), posterior_precision(:, :)
      real(dp), allocatable :: posterior_covariance(:, :), posterior_mean(:, :)
      real(dp), allocatable :: posterior_scale(:, :), log_variance(:)
      real(dp) :: persistence, innovation_variance, loading_prior_variance
      real(dp) :: hyper_scale, shape_value, scale_value
      integer :: discard, total, iteration, retained, time, info

      discard = 0
      shape_value = 1.0_dp
      scale_value = 0.1_dp
      if (present(burnin)) discard = burnin
      if (present(prior_shape)) shape_value = prior_shape
      if (present(prior_scale)) scale_value = prior_scale
      if (data%info /= 0 .or. prior%info /= 0 .or. draws < 1 .or. discard < 0 .or. &
         size(mixture, 1) /= 3 .or. size(mixture, 2) < 1 .or. &
         any(mixture(3, :) <= 0.0_dp) .or. shape_value <= 0.0_dp .or. &
         scale_value <= 0.0_dp .or. &
         prior%covariance_degrees_of_freedom <= real(data%variables + 1, dp)) then
         out%info = 1
         return
      end if
      call invert_matrix(prior%coefficient_covariance, prior_precision, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      coefficient = prior%coefficient_mean
      covariance = prior%covariance_scale/ &
         (prior%covariance_degrees_of_freedom - real(data%variables + 1, dp))
      allocate(log_variance(data%observations), log_observation(data%observations))
      allocate(weight(data%observations))
      log_variance = 0.0_dp
      persistence = 0.5_dp
      innovation_variance = 0.01_dp
      loading_prior_variance = 1.0_dp
      hyper_scale = 0.05_dp
      allocate(out%coefficient(data%variables, data%regressors, draws))
      allocate(out%covariance(data%variables, data%variables, draws))
      allocate(out%coefficient_covariance(data%regressors, data%regressors, draws))
      allocate(out%common_variance(data%observations, draws))
      allocate(out%scale_mixture(data%observations, draws))
      allocate(out%log_variance(data%observations, draws))
      allocate(out%mixture_component(data%observations, draws))
      allocate(out%persistence(draws), out%volatility_loading(draws))
      allocate(out%volatility_innovation_variance(draws))
      allocate(out%loading_prior_variance(draws), out%volatility_hyper_scale(draws))
      out%scale_mixture = 1.0_dp
      total = draws + discard
      retained = 0
      do iteration = 1, total
         residual = data%y - matmul(coefficient, data%x)
         call cholesky_lower(covariance, lower, info)
         if (info /= 0) then
            out%info = 3
            return
         end if
         call invert_matrix(lower, lower_inverse, info)
         if (info /= 0) then
            out%info = 4
            return
         end if
         standardized = matmul(lower_inverse, residual)
         do time = 1, data%observations
            log_observation(time) = sum(log(standardized(:, time)**2 + &
               1.0e-9_dp))/real(data%variables, dp)
         end do
         volatility = bvars_centered_sv_update(log_observation, log_variance, &
            persistence, innovation_variance, loading_prior_variance, &
            hyper_scale, mixture, shape_value, scale_value)
         if (volatility%info /= 0) then
            out%info = 10 + volatility%info
            return
         end if
         log_variance = volatility%log_variance
         persistence = volatility%persistence
         innovation_variance = volatility%innovation_variance
         loading_prior_variance = volatility%loading_prior_variance
         hyper_scale = volatility%hyper_scale
         weight = exp(-0.5_dp*log_variance)
         weighted_y = data%y*spread(weight, 1, data%variables)
         weighted_x = data%x*spread(weight, 1, data%regressors)
         posterior_precision = matmul(weighted_x, transpose(weighted_x)) + &
            prior_precision
         call invert_matrix(posterior_precision, posterior_covariance, info)
         if (info /= 0) then
            out%info = 20
            return
         end if
         posterior_mean = matmul(matmul(weighted_y, transpose(weighted_x)) + &
            matmul(prior%coefficient_mean, prior_precision), posterior_covariance)
         posterior_scale = prior%covariance_scale + &
            matmul(weighted_y, transpose(weighted_y)) + &
            matmul(matmul(prior%coefficient_mean, prior_precision), &
            transpose(prior%coefficient_mean)) - &
            matmul(matmul(posterior_mean, posterior_precision), &
            transpose(posterior_mean))
         posterior_scale = 0.5_dp*(posterior_scale + transpose(posterior_scale))
         call inverse_wishart_draw(posterior_scale, &
            prior%covariance_degrees_of_freedom + real(data%observations, dp), &
            covariance, precision, info)
         if (info /= 0) then
            out%info = 21
            return
         end if
         coefficient = bvars_matrix_normal(posterior_mean, covariance, &
            posterior_covariance)
         if (iteration > discard) then
            retained = retained + 1
            out%coefficient(:, :, retained) = coefficient
            out%covariance(:, :, retained) = covariance
            out%coefficient_covariance(:, :, retained) = &
               prior%coefficient_covariance
            out%log_variance(:, retained) = log_variance
            out%common_variance(:, retained) = exp(log_variance)
            out%mixture_component(:, retained) = volatility%mixture_component
            out%persistence(retained) = persistence
            out%volatility_loading(retained) = volatility%loading
            out%volatility_innovation_variance(retained) = innovation_variance
            out%loading_prior_variance(retained) = loading_prior_variance
            out%volatility_hyper_scale(retained) = hyper_scale
         end if
      end do
      out%retained_draws = retained
   end function bvars_centered_sv_draws

   function bvars_noncentered_sv_draws(data, prior, mixture, draws, burnin, &
      prior_shape, prior_scale) result(out)
      !! Sample the BVAR with non-centred common stochastic volatility and ASIS.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      type(bvars_prior_t), intent(in) :: prior !! BVAR prior specification.
      real(dp), intent(in) :: mixture(:, :) !! Mixture probability scores, means, and variances.
      integer, intent(in) :: draws !! Number of retained posterior draws.
      integer, intent(in), optional :: burnin !! Number of discarded posterior draws.
      real(dp), intent(in), optional :: prior_shape !! Hierarchical volatility prior shape.
      real(dp), intent(in), optional :: prior_scale !! Hierarchical volatility prior scale.
      type(bvars_draws_t) :: out
      type(bvars_common_sv_t) :: volatility
      real(dp), allocatable :: coefficient(:, :), covariance(:, :), precision(:, :)
      real(dp), allocatable :: residual(:, :), lower(:, :), lower_inverse(:, :)
      real(dp), allocatable :: standardized(:, :), log_observation(:), weight(:)
      real(dp), allocatable :: weighted_y(:, :), weighted_x(:, :)
      real(dp), allocatable :: prior_precision(:, :), posterior_precision(:, :)
      real(dp), allocatable :: posterior_covariance(:, :), posterior_mean(:, :)
      real(dp), allocatable :: posterior_scale(:, :), state(:), log_variance(:)
      real(dp) :: persistence, loading, innovation_variance
      real(dp) :: loading_prior_variance, hyper_scale, shape_value, scale_value
      integer :: discard, total, iteration, retained, time, info

      discard = 0
      shape_value = 1.0_dp
      scale_value = 0.1_dp
      if (present(burnin)) discard = burnin
      if (present(prior_shape)) shape_value = prior_shape
      if (present(prior_scale)) scale_value = prior_scale
      if (data%info /= 0 .or. prior%info /= 0 .or. draws < 1 .or. discard < 0 .or. &
         size(mixture, 1) /= 3 .or. size(mixture, 2) < 1 .or. &
         any(mixture(3, :) <= 0.0_dp) .or. shape_value <= 0.0_dp .or. &
         scale_value <= 0.0_dp .or. &
         prior%covariance_degrees_of_freedom <= real(data%variables + 1, dp)) then
         out%info = 1
         return
      end if
      call invert_matrix(prior%coefficient_covariance, prior_precision, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      coefficient = prior%coefficient_mean
      covariance = prior%covariance_scale/ &
         (prior%covariance_degrees_of_freedom - real(data%variables + 1, dp))
      allocate(state(data%observations), log_variance(data%observations))
      allocate(log_observation(data%observations), weight(data%observations))
      state = 0.0_dp
      log_variance = 0.0_dp
      persistence = 0.5_dp
      loading = 0.1_dp
      innovation_variance = loading**2
      loading_prior_variance = 0.5_dp
      hyper_scale = 0.05_dp
      allocate(out%coefficient(data%variables, data%regressors, draws))
      allocate(out%covariance(data%variables, data%variables, draws))
      allocate(out%coefficient_covariance(data%regressors, data%regressors, draws))
      allocate(out%common_variance(data%observations, draws))
      allocate(out%scale_mixture(data%observations, draws))
      allocate(out%log_variance(data%observations, draws))
      allocate(out%mixture_component(data%observations, draws))
      allocate(out%persistence(draws), out%volatility_loading(draws))
      allocate(out%volatility_innovation_variance(draws))
      allocate(out%loading_prior_variance(draws), out%volatility_hyper_scale(draws))
      out%scale_mixture = 1.0_dp
      total = draws + discard
      retained = 0
      do iteration = 1, total
         residual = data%y - matmul(coefficient, data%x)
         call cholesky_lower(covariance, lower, info)
         if (info /= 0) then
            out%info = 3
            return
         end if
         call invert_matrix(lower, lower_inverse, info)
         if (info /= 0) then
            out%info = 4
            return
         end if
         standardized = matmul(lower_inverse, residual)
         do time = 1, data%observations
            log_observation(time) = sum(log(standardized(:, time)**2 + &
               1.0e-9_dp))/real(data%variables, dp)
         end do
         volatility = bvars_noncentered_sv_update(log_observation, state, &
            persistence, loading, innovation_variance, &
            loading_prior_variance, hyper_scale, mixture, shape_value, scale_value)
         if (volatility%info /= 0) then
            out%info = 10 + volatility%info
            return
         end if
         log_variance = volatility%log_variance
         loading = volatility%loading
         state = log_variance/loading
         persistence = volatility%persistence
         innovation_variance = volatility%innovation_variance
         loading_prior_variance = volatility%loading_prior_variance
         hyper_scale = volatility%hyper_scale
         weight = exp(-0.5_dp*log_variance)
         weighted_y = data%y*spread(weight, 1, data%variables)
         weighted_x = data%x*spread(weight, 1, data%regressors)
         posterior_precision = matmul(weighted_x, transpose(weighted_x)) + &
            prior_precision
         call invert_matrix(posterior_precision, posterior_covariance, info)
         if (info /= 0) then
            out%info = 20
            return
         end if
         posterior_mean = matmul(matmul(weighted_y, transpose(weighted_x)) + &
            matmul(prior%coefficient_mean, prior_precision), posterior_covariance)
         posterior_scale = prior%covariance_scale + &
            matmul(weighted_y, transpose(weighted_y)) + &
            matmul(matmul(prior%coefficient_mean, prior_precision), &
            transpose(prior%coefficient_mean)) - &
            matmul(matmul(posterior_mean, posterior_precision), &
            transpose(posterior_mean))
         posterior_scale = 0.5_dp*(posterior_scale + transpose(posterior_scale))
         call inverse_wishart_draw(posterior_scale, &
            prior%covariance_degrees_of_freedom + real(data%observations, dp), &
            covariance, precision, info)
         if (info /= 0) then
            out%info = 21
            return
         end if
         coefficient = bvars_matrix_normal(posterior_mean, covariance, &
            posterior_covariance)
         if (iteration > discard) then
            retained = retained + 1
            out%coefficient(:, :, retained) = coefficient
            out%covariance(:, :, retained) = covariance
            out%coefficient_covariance(:, :, retained) = &
               prior%coefficient_covariance
            out%log_variance(:, retained) = log_variance
            out%common_variance(:, retained) = exp(log_variance)
            out%mixture_component(:, retained) = volatility%mixture_component
            out%persistence(retained) = persistence
            out%volatility_loading(retained) = loading
            out%volatility_innovation_variance(retained) = innovation_variance
            out%loading_prior_variance(retained) = loading_prior_variance
            out%volatility_hyper_scale(retained) = hyper_scale
         end if
      end do
      out%retained_draws = retained
   end function bvars_noncentered_sv_draws

   function bvars_student_sv_draws(data, prior, sv_mixture, centered, draws, &
      burnin, initial_df, initial_proposal_scale, prior_shape, prior_scale) &
      result(out)
      !! Sample a Student-t BVAR with centred or non-centred common SV.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      type(bvars_prior_t), intent(in) :: prior !! BVAR prior specification.
      real(dp), intent(in) :: sv_mixture(:, :) !! SV mixture probability scores, means, and variances.
      logical, intent(in) :: centered !! Use centred rather than non-centred SV.
      integer, intent(in) :: draws !! Number of retained posterior draws.
      integer, intent(in), optional :: burnin !! Number of discarded posterior draws.
      real(dp), intent(in), optional :: initial_df !! Initial degrees of freedom greater than two.
      real(dp), intent(in), optional :: initial_proposal_scale !! Initial transformed proposal variance.
      real(dp), intent(in), optional :: prior_shape !! Hierarchical volatility prior shape.
      real(dp), intent(in), optional :: prior_scale !! Hierarchical volatility prior scale.
      type(bvars_draws_t) :: out
      type(bvars_common_sv_t) :: volatility
      real(dp), allocatable :: coefficient(:, :), covariance(:, :), precision(:, :)
      real(dp), allocatable :: residual(:, :), lower(:, :), lower_inverse(:, :)
      real(dp), allocatable :: standardized(:, :), quadratic(:), scale_mixture(:)
      real(dp), allocatable :: log_observation(:), observation_weight(:)
      real(dp), allocatable :: weighted_y(:, :), weighted_x(:, :)
      real(dp), allocatable :: prior_precision(:, :), posterior_precision(:, :)
      real(dp), allocatable :: posterior_covariance(:, :), posterior_mean(:, :)
      real(dp), allocatable :: posterior_scale(:, :), state(:), log_variance(:)
      real(dp) :: degrees, candidate, eta, candidate_eta, proposal_variance
      real(dp) :: old_kernel, new_kernel, log_ratio, acceptance
      real(dp) :: persistence, loading, innovation_variance
      real(dp) :: loading_prior_variance, hyper_scale, shape_value, scale_value
      integer :: discard, total, iteration, retained, time, info
      logical :: accepted

      discard = 0
      degrees = 30.0_dp
      proposal_variance = 0.05_dp
      shape_value = 1.0_dp
      scale_value = 0.1_dp
      if (present(burnin)) discard = burnin
      if (present(initial_df)) degrees = initial_df
      if (present(initial_proposal_scale)) proposal_variance = &
         initial_proposal_scale
      if (present(prior_shape)) shape_value = prior_shape
      if (present(prior_scale)) scale_value = prior_scale
      if (data%info /= 0 .or. prior%info /= 0 .or. draws < 1 .or. discard < 0 .or. &
         degrees <= 2.0_dp .or. proposal_variance <= 0.0_dp .or. &
         size(sv_mixture, 1) /= 3 .or. size(sv_mixture, 2) < 1 .or. &
         any(sv_mixture(3, :) <= 0.0_dp) .or. shape_value <= 0.0_dp .or. &
         scale_value <= 0.0_dp .or. &
         prior%covariance_degrees_of_freedom <= real(data%variables + 1, dp)) then
         out%info = 1
         return
      end if
      call invert_matrix(prior%coefficient_covariance, prior_precision, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      coefficient = prior%coefficient_mean
      covariance = prior%covariance_scale/ &
         (prior%covariance_degrees_of_freedom - real(data%variables + 1, dp))
      allocate(state(data%observations), log_variance(data%observations))
      allocate(quadratic(data%observations), scale_mixture(data%observations))
      allocate(log_observation(data%observations))
      allocate(observation_weight(data%observations))
      state = 0.0_dp
      log_variance = 0.0_dp
      scale_mixture = 1.0_dp
      persistence = 0.5_dp
      loading = 0.1_dp
      innovation_variance = loading**2
      loading_prior_variance = 0.5_dp
      hyper_scale = 0.05_dp
      allocate(out%coefficient(data%variables, data%regressors, draws))
      allocate(out%covariance(data%variables, data%variables, draws))
      allocate(out%coefficient_covariance(data%regressors, data%regressors, draws))
      allocate(out%common_variance(data%observations, draws))
      allocate(out%scale_mixture(data%observations, draws))
      allocate(out%degrees_of_freedom(draws), out%degrees_accepted(draws))
      allocate(out%log_variance(data%observations, draws))
      allocate(out%mixture_component(data%observations, draws))
      allocate(out%persistence(draws), out%volatility_loading(draws))
      allocate(out%volatility_innovation_variance(draws))
      allocate(out%loading_prior_variance(draws), out%volatility_hyper_scale(draws))
      total = draws + discard
      retained = 0
      do iteration = 1, total
         eta = log(degrees - 2.0_dp)
         candidate_eta = eta + sqrt(proposal_variance)*random_standard_normal()
         candidate = 2.0_dp + exp(candidate_eta)
         old_kernel = bvars_student_log_kernel(degrees, scale_mixture) + eta
         new_kernel = bvars_student_log_kernel(candidate, scale_mixture) + &
            candidate_eta
         log_ratio = min(0.0_dp, new_kernel - old_kernel)
         acceptance = exp(log_ratio)
         accepted = log(max(random_uniform(), tiny(1.0_dp))) < log_ratio
         if (accepted) degrees = candidate
         if (iteration > 1) proposal_variance = exp(log(proposal_variance) + &
            0.5_dp*log(1.0_dp + real(iteration, dp)**(-0.6_dp)* &
            (acceptance - 0.44_dp)))

         residual = data%y - matmul(coefficient, data%x)
         call cholesky_lower(covariance, lower, info)
         if (info /= 0) then
            out%info = 3
            return
         end if
         call invert_matrix(lower, lower_inverse, info)
         if (info /= 0) then
            out%info = 4
            return
         end if
         standardized = matmul(lower_inverse, residual)
         do time = 1, data%observations
            quadratic(time) = sum(standardized(:, time)**2)/ &
               (real(data%variables, dp)*exp(log_variance(time)))
         end do
         if (sum(quadratic) <= tiny(1.0_dp)) then
            quadratic = 1.0_dp
         else
            quadratic = quadratic/(sum(quadratic)/real(data%observations, dp))
         end if
         do time = 1, data%observations
            scale_mixture(time) = (quadratic(time) + degrees - 2.0_dp)/ &
               (2.0_dp*random_gamma(0.5_dp*(degrees + &
               real(data%variables, dp))))
         end do

         standardized = standardized/ &
            spread(sqrt(scale_mixture), 1, data%variables)
         do time = 1, data%observations
            log_observation(time) = sum(log(standardized(:, time)**2 + &
               1.0e-9_dp))/real(data%variables, dp)
         end do
         if (centered) then
            volatility = bvars_centered_sv_update(log_observation, &
               log_variance, persistence, innovation_variance, &
               loading_prior_variance, hyper_scale, sv_mixture, &
               shape_value, scale_value)
         else
            volatility = bvars_noncentered_sv_update(log_observation, state, &
               persistence, loading, innovation_variance, &
               loading_prior_variance, hyper_scale, sv_mixture, &
               shape_value, scale_value)
         end if
         if (volatility%info /= 0) then
            out%info = 10 + volatility%info
            return
         end if
         log_variance = volatility%log_variance
         loading = volatility%loading
         if (.not. centered) state = log_variance/loading
         persistence = volatility%persistence
         innovation_variance = volatility%innovation_variance
         loading_prior_variance = volatility%loading_prior_variance
         hyper_scale = volatility%hyper_scale

         observation_weight = 1.0_dp/ &
            sqrt(exp(log_variance)*scale_mixture)
         weighted_y = data%y*spread(observation_weight, 1, data%variables)
         weighted_x = data%x*spread(observation_weight, 1, data%regressors)
         posterior_precision = matmul(weighted_x, transpose(weighted_x)) + &
            prior_precision
         call invert_matrix(posterior_precision, posterior_covariance, info)
         if (info /= 0) then
            out%info = 20
            return
         end if
         posterior_mean = matmul(matmul(weighted_y, transpose(weighted_x)) + &
            matmul(prior%coefficient_mean, prior_precision), posterior_covariance)
         posterior_scale = prior%covariance_scale + &
            matmul(weighted_y, transpose(weighted_y)) + &
            matmul(matmul(prior%coefficient_mean, prior_precision), &
            transpose(prior%coefficient_mean)) - &
            matmul(matmul(posterior_mean, posterior_precision), &
            transpose(posterior_mean))
         posterior_scale = 0.5_dp*(posterior_scale + transpose(posterior_scale))
         call inverse_wishart_draw(posterior_scale, &
            prior%covariance_degrees_of_freedom + real(data%observations, dp), &
            covariance, precision, info)
         if (info /= 0) then
            out%info = 21
            return
         end if
         coefficient = bvars_matrix_normal(posterior_mean, covariance, &
            posterior_covariance)
         if (iteration > discard) then
            retained = retained + 1
            out%coefficient(:, :, retained) = coefficient
            out%covariance(:, :, retained) = covariance
            out%coefficient_covariance(:, :, retained) = &
               prior%coefficient_covariance
            out%common_variance(:, retained) = exp(log_variance)
            out%scale_mixture(:, retained) = scale_mixture
            out%degrees_of_freedom(retained) = degrees
            out%degrees_accepted(retained) = accepted
            out%log_variance(:, retained) = log_variance
            out%mixture_component(:, retained) = volatility%mixture_component
            out%persistence(retained) = persistence
            out%volatility_loading(retained) = loading
            out%volatility_innovation_variance(retained) = innovation_variance
            out%loading_prior_variance(retained) = loading_prior_variance
            out%volatility_hyper_scale(retained) = hyper_scale
         end if
      end do
      out%retained_draws = retained
   end function bvars_student_sv_draws

   function bvars_common_variance_forecast(last_log_variance, persistence, &
      innovation_standard_deviation, horizon) result(common_variance)
      !! Simulate future common stochastic-variance paths.
      real(dp), intent(in) :: last_log_variance(:) !! Draw-specific terminal log variance.
      real(dp), intent(in) :: persistence(:) !! Draw-specific AR(1) persistence.
      real(dp), intent(in) :: innovation_standard_deviation(:) !! Draw-specific innovation deviation.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), allocatable :: common_variance(:, :)
      real(dp), allocatable :: state(:)
      integer :: draw, step

      if (horizon < 1 .or. size(last_log_variance) < 1 .or. &
         size(persistence) /= size(last_log_variance) .or. &
         size(innovation_standard_deviation) /= size(last_log_variance) .or. &
         any(abs(persistence) >= 1.0_dp) .or. &
         any(innovation_standard_deviation <= 0.0_dp)) then
         allocate(common_variance(0, 0))
         return
      end if
      allocate(common_variance(horizon, size(last_log_variance)))
      state = last_log_variance
      do step = 1, horizon
         do draw = 1, size(state)
            state(draw) = persistence(draw)*state(draw) + &
               innovation_standard_deviation(draw)*random_standard_normal()
            common_variance(step, draw) = exp(state(draw))
         end do
      end do
   end function bvars_common_variance_forecast

   pure function bvars_shocks(draws, data) result(shocks)
      !! Compute reduced-form residual draws.
      type(bvars_draws_t), intent(in) :: draws !! Posterior BVAR draws.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      real(dp), allocatable :: shocks(:, :, :)
      integer :: draw

      if (.not. allocated(draws%coefficient) .or. data%info /= 0) then
         allocate(shocks(0, 0, 0))
         return
      end if
      allocate(shocks(data%variables, data%observations, draws%retained_draws))
      do draw = 1, draws%retained_draws
         shocks(:, :, draw) = data%y - &
            matmul(draws%coefficient(:, :, draw), data%x)
      end do
   end function bvars_shocks

   pure function bvars_fitted_from_random(draws, data, standard_normal) result(out)
      !! Draw fitted observations from supplied standard normals.
      type(bvars_draws_t), intent(in) :: draws !! Posterior BVAR draws.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      real(dp), intent(in) :: standard_normal(:, :, :) !! Normals by variable, observation, and draw.
      type(bvars_fitted_t) :: out
      real(dp), allocatable :: factor(:, :)
      integer :: draw, info

      if (.not. allocated(draws%coefficient) .or. &
         any(shape(standard_normal) /= &
         [data%variables, data%observations, draws%retained_draws])) then
         out%info = 1
         return
      end if
      allocate(out%value(data%variables, data%observations, draws%retained_draws))
      do draw = 1, draws%retained_draws
         call cholesky_lower(draws%covariance(:, :, draw), factor, info)
         if (info /= 0) then
            out%info = 2
            return
         end if
         out%value(:, :, draw) = matmul(draws%coefficient(:, :, draw), data%x) + &
            matmul(factor, standard_normal(:, :, draw))*spread(sqrt( &
            draws%common_variance(:, draw)*draws%scale_mixture(:, draw)), &
            1, data%variables)
      end do
   end function bvars_fitted_from_random

   function bvars_fitted(draws, data) result(out)
      !! Draw fitted observations using the shared RNG.
      type(bvars_draws_t), intent(in) :: draws !! Posterior BVAR draws.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      type(bvars_fitted_t) :: out
      real(dp), allocatable :: standard_normal(:, :, :)
      integer :: draw

      allocate(standard_normal(data%variables, data%observations, &
         draws%retained_draws))
      do draw = 1, draws%retained_draws
         call random_standard_normal_matrix(standard_normal(:, :, draw))
      end do
      out = bvars_fitted_from_random(draws, data, standard_normal)
   end function bvars_fitted

   pure function bvars_forecast_from_random(draws, data, horizon, &
      standard_normal, future_exogenous, future_common_variance, &
      future_scale_mixture, conditional_path) result(out)
      !! Generate recursive BVAR forecasts from supplied normal draws.
      type(bvars_draws_t), intent(in) :: draws !! Posterior BVAR draws.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in) :: standard_normal(:, :, :) !! Normals by variable, horizon, and draw.
      real(dp), intent(in), optional :: future_exogenous(:, :) !! Exogenous values by horizon and variable.
      real(dp), intent(in), optional :: future_common_variance(:, :) !! Common variance by horizon and draw.
      real(dp), intent(in), optional :: future_scale_mixture(:, :) !! Scale mixture by horizon and draw.
      real(dp), intent(in), optional :: conditional_path(:, :) !! Fixed values by horizon and variable, with NaN for free values.
      type(bvars_predictive_t) :: out
      real(dp), allocatable :: regressor(:), history(:, :), factor(:, :)
      real(dp) :: scale
      integer :: draw, step, lag, first, info

      if (horizon < 1 .or. draws%retained_draws < 1 .or. &
         any(shape(standard_normal) /= &
         [data%variables, horizon, draws%retained_draws]) .or. &
         (data%exogenous_count > 0 .and. .not. present(future_exogenous)) .or. &
         (present(future_exogenous) .and. any(shape(future_exogenous) /= &
         [horizon, data%exogenous_count])) .or. &
         (present(future_common_variance) .and. &
         any(shape(future_common_variance) /= [horizon, draws%retained_draws])) .or. &
         (present(future_scale_mixture) .and. &
         any(shape(future_scale_mixture) /= [horizon, draws%retained_draws])) .or. &
         (present(conditional_path) .and. any(shape(conditional_path) /= &
         [horizon, data%variables]))) then
         out%info = 1
         return
      end if
      allocate(out%path(data%variables, horizon, draws%retained_draws))
      allocate(out%mean(data%variables, horizon, draws%retained_draws))
      allocate(out%covariance(data%variables, data%variables, horizon, &
         draws%retained_draws))
      allocate(regressor(data%regressors), history(data%variables, data%lags))
      do draw = 1, draws%retained_draws
         history = data%series(:, size(data%series, 2) - data%lags + 1:)
         call cholesky_lower(draws%covariance(:, :, draw), factor, info)
         if (info /= 0) then
            out%info = 2
            return
         end if
         do step = 1, horizon
            do lag = 1, data%lags
               first = (lag - 1)*data%variables + 1
               regressor(first:first + data%variables - 1) = &
                  history(:, data%lags + 1 - lag)
            end do
            regressor(data%variables*data%lags + 1) = 1.0_dp
            if (data%exogenous_count > 0) regressor(data%variables*data%lags + 2:) = &
               future_exogenous(step, :)
            scale = 1.0_dp
            if (present(future_common_variance)) scale = scale* &
               future_common_variance(step, draw)
            if (present(future_scale_mixture)) scale = scale* &
               future_scale_mixture(step, draw)
            if (scale <= 0.0_dp) then
               out%info = 3
               return
            end if
            out%mean(:, step, draw) = &
               matmul(draws%coefficient(:, :, draw), regressor)
            out%covariance(:, :, step, draw) = &
               scale*draws%covariance(:, :, draw)
            if (present(conditional_path)) then
               call conditional_normal_from_standard(out%mean(:, step, draw), &
                  out%covariance(:, :, step, draw), conditional_path(step, :), &
                  standard_normal(:, step, draw), out%path(:, step, draw), info)
               if (info /= 0) then
                  out%info = 4
                  return
               end if
            else
               out%path(:, step, draw) = out%mean(:, step, draw) + &
                  sqrt(scale)*matmul(factor, standard_normal(:, step, draw))
            end if
            if (data%lags > 1) history(:, :data%lags - 1) = history(:, 2:)
            history(:, data%lags) = out%path(:, step, draw)
         end do
      end do
   end function bvars_forecast_from_random

   pure function bvars_conditional_forecast_from_random(draws, data, horizon, &
      standard_normal, conditional_path, future_exogenous, &
      future_common_variance, future_scale_mixture) result(out)
      !! Generate BVAR forecasts conditionally on supplied future values.
      type(bvars_draws_t), intent(in) :: draws !! Posterior BVAR draws.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in) :: standard_normal(:, :, :) !! Normals by variable, horizon, and draw.
      real(dp), intent(in) :: conditional_path(:, :) !! Fixed values by horizon and variable, with NaN for free values.
      real(dp), intent(in), optional :: future_exogenous(:, :) !! Exogenous values by horizon and variable.
      real(dp), intent(in), optional :: future_common_variance(:, :) !! Common variance by horizon and draw.
      real(dp), intent(in), optional :: future_scale_mixture(:, :) !! Scale mixture by horizon and draw.
      type(bvars_predictive_t) :: out

      out = bvars_forecast_from_random(draws, data, horizon, standard_normal, &
         future_exogenous, future_common_variance, future_scale_mixture, &
         conditional_path)
   end function bvars_conditional_forecast_from_random

   function bvars_forecast(draws, data, horizon, future_exogenous, &
      future_common_variance, future_scale_mixture, conditional_path) result(out)
      !! Generate recursive BVAR forecasts using the shared RNG.
      type(bvars_draws_t), intent(in) :: draws !! Posterior BVAR draws.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: future_exogenous(:, :) !! Exogenous values by horizon and variable.
      real(dp), intent(in), optional :: future_common_variance(:, :) !! Common variance by horizon and draw.
      real(dp), intent(in), optional :: future_scale_mixture(:, :) !! Scale mixture by horizon and draw.
      real(dp), intent(in), optional :: conditional_path(:, :) !! Fixed values by horizon and variable, with NaN for free values.
      type(bvars_predictive_t) :: out
      real(dp), allocatable :: standard_normal(:, :, :)
      integer :: draw

      allocate(standard_normal(data%variables, horizon, draws%retained_draws))
      do draw = 1, draws%retained_draws
         call random_standard_normal_matrix(standard_normal(:, :, draw))
      end do
      if (present(conditional_path)) then
         out = bvars_forecast_from_random(draws, data, horizon, standard_normal, &
            future_exogenous, future_common_variance, future_scale_mixture, &
            conditional_path)
      else
         out = bvars_forecast_from_random(draws, data, horizon, standard_normal, &
            future_exogenous, future_common_variance, future_scale_mixture)
      end if
   end function bvars_forecast

   function bvars_conditional_forecast(draws, data, horizon, conditional_path, &
      future_exogenous, future_common_variance, future_scale_mixture) result(out)
      !! Generate BVAR forecasts conditionally using the shared RNG.
      type(bvars_draws_t), intent(in) :: draws !! Posterior BVAR draws.
      type(bvars_data_t), intent(in) :: data !! Prepared BVAR data.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in) :: conditional_path(:, :) !! Fixed values by horizon and variable, with NaN for free values.
      real(dp), intent(in), optional :: future_exogenous(:, :) !! Exogenous values by horizon and variable.
      real(dp), intent(in), optional :: future_common_variance(:, :) !! Common variance by horizon and draw.
      real(dp), intent(in), optional :: future_scale_mixture(:, :) !! Scale mixture by horizon and draw.
      type(bvars_predictive_t) :: out
      real(dp), allocatable :: standard_normal(:, :, :)
      integer :: draw

      allocate(standard_normal(data%variables, horizon, draws%retained_draws))
      do draw = 1, draws%retained_draws
         call random_standard_normal_matrix(standard_normal(:, :, draw))
      end do
      out = bvars_conditional_forecast_from_random(draws, data, horizon, &
         standard_normal, conditional_path, future_exogenous, &
         future_common_variance, future_scale_mixture)
   end function bvars_conditional_forecast

   pure subroutine conditional_normal_from_standard(mean, covariance, condition, &
      standard_normal, draw, info)
      !! Draw unconstrained entries from a conditional multivariate normal law.
      real(dp), intent(in) :: mean(:) !! Unconditional multivariate normal mean.
      real(dp), intent(in) :: covariance(:, :) !! Unconditional positive-definite covariance.
      real(dp), intent(in) :: condition(:) !! Fixed values, with NaN entries left unconstrained.
      real(dp), intent(in) :: standard_normal(:) !! Independent standard-normal variates.
      real(dp), intent(out) :: draw(:) !! Conditional multivariate normal draw.
      integer, intent(out) :: info !! Status code; zero indicates success.
      integer, allocatable :: constrained(:), unconstrained(:)
      real(dp), allocatable :: inverse_fixed(:, :), conditional_mean(:)
      real(dp), allocatable :: conditional_covariance(:, :), factor(:, :)
      integer :: variables, constrained_count, unconstrained_count, variable

      info = 0
      variables = size(mean)
      if (variables < 1 .or. size(condition) /= variables .or. &
         size(standard_normal) < variables .or. size(draw) /= variables .or. &
         any(shape(covariance) /= [variables, variables])) then
         info = 1
         return
      end if
      allocate(constrained(variables), unconstrained(variables))
      constrained_count = 0
      unconstrained_count = 0
      do variable = 1, variables
         if (ieee_is_nan(condition(variable))) then
            unconstrained_count = unconstrained_count + 1
            unconstrained(unconstrained_count) = variable
         else if (ieee_is_finite(condition(variable))) then
            constrained_count = constrained_count + 1
            constrained(constrained_count) = variable
         else
            info = 2
            return
         end if
      end do
      if (constrained_count == 0) then
         call cholesky_lower(covariance, factor, info)
         if (info /= 0) return
         draw = mean + matmul(factor, standard_normal)
         return
      end if
      draw(constrained(:constrained_count)) = condition(constrained(:constrained_count))
      if (unconstrained_count == 0) return
      call invert_matrix(covariance(constrained(:constrained_count), &
         constrained(:constrained_count)), inverse_fixed, info)
      if (info /= 0) then
         info = 3
         return
      end if
      conditional_mean = mean(unconstrained(:unconstrained_count)) + matmul( &
         matmul(covariance(unconstrained(:unconstrained_count), &
         constrained(:constrained_count)), inverse_fixed), &
         condition(constrained(:constrained_count)) - &
         mean(constrained(:constrained_count)))
      conditional_covariance = covariance(unconstrained(:unconstrained_count), &
         unconstrained(:unconstrained_count)) - matmul(matmul( &
         covariance(unconstrained(:unconstrained_count), &
         constrained(:constrained_count)), inverse_fixed), covariance( &
         constrained(:constrained_count), unconstrained(:unconstrained_count)))
      conditional_covariance = 0.5_dp*(conditional_covariance + &
         transpose(conditional_covariance))
      call cholesky_lower(conditional_covariance, factor, info)
      if (info /= 0) then
         info = 4
         return
      end if
      draw(unconstrained(:unconstrained_count)) = conditional_mean + matmul( &
         factor, standard_normal(:unconstrained_count))
   end subroutine conditional_normal_from_standard

   pure function bvars_fevd(draws, lags, horizon) result(out)
      !! Compute orthogonalized forecast-error variance decompositions.
      type(bvars_draws_t), intent(in) :: draws !! Posterior BVAR draws.
      integer, intent(in) :: lags !! Autoregressive lag order.
      integer, intent(in) :: horizon !! Largest forecast-error horizon.
      type(bvartools_fevd_t) :: out
      type(bvartools_bvar_draws_t) :: compatible
      integer :: draw, coefficient_count

      if (.not. allocated(draws%coefficient) .or. lags < 1 .or. horizon < 0) then
         out%info = 1
         return
      end if
      coefficient_count = size(draws%coefficient, 1)*size(draws%coefficient, 2)
      allocate(compatible%coefficients(coefficient_count, draws%retained_draws))
      compatible%covariance = draws%covariance
      compatible%retained_draws = draws%retained_draws
      do draw = 1, draws%retained_draws
         compatible%coefficients(:, draw) = reshape( &
            draws%coefficient(:, :, draw), [coefficient_count])
      end do
      out = bvartools_bvar_fevd(compatible, lags, horizon, 'oir')
   end function bvars_fevd

   subroutine inverse_wishart_draw(scale, degrees_of_freedom, covariance, &
      precision, info)
      !! Draw an inverse-Wishart matrix by Bartlett decomposition.
      real(dp), intent(in) :: scale(:, :) !! Inverse-Wishart scale matrix.
      real(dp), intent(in) :: degrees_of_freedom !! Degrees of freedom.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Drawn covariance.
      real(dp), allocatable, intent(out) :: precision(:, :) !! Drawn precision.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: inverse_scale(:, :), lower(:, :)
      real(dp), allocatable :: bartlett(:, :), factor(:, :)
      integer :: variables, row, column

      variables = size(scale, 1)
      call invert_matrix(scale, inverse_scale, info)
      if (info /= 0) return
      call cholesky_lower(inverse_scale, lower, info)
      if (info /= 0) return
      allocate(bartlett(variables, variables))
      bartlett = 0.0_dp
      do row = 1, variables
         bartlett(row, row) = sqrt(2.0_dp*random_gamma(0.5_dp* &
            (degrees_of_freedom - real(row, dp) + 1.0_dp)))
         do column = 1, row - 1
            bartlett(row, column) = random_standard_normal()
         end do
      end do
      factor = matmul(lower, bartlett)
      precision = matmul(factor, transpose(factor))
      call invert_matrix(precision, covariance, info)
   end subroutine inverse_wishart_draw

   real(dp) function truncated_normal_random(mean, standard_deviation, lower, &
      upper) result(draw)
      !! Draw a normal variate restricted to a finite interval.
      real(dp), intent(in) :: mean !! Untruncated normal mean.
      real(dp), intent(in) :: standard_deviation !! Positive normal deviation.
      real(dp), intent(in) :: lower !! Inclusive lower bound.
      real(dp), intent(in) :: upper !! Inclusive upper bound.
      integer :: attempt

      if (standard_deviation <= 0.0_dp .or. lower >= upper) then
         draw = max(lower, min(upper, mean))
         return
      end if
      do attempt = 1, 100000
         draw = mean + standard_deviation*random_standard_normal()
         if (draw >= lower .and. draw <= upper) return
      end do
      draw = max(lower, min(upper, mean))
   end function truncated_normal_random

   real(dp) function gig_slice_draw(shape, inverse_scale, scale, current) &
      result(draw)
      !! Draw a generalized inverse-Gaussian variate by log-scale slice sampling.
      real(dp), intent(in) :: shape !! GIG power parameter.
      real(dp), intent(in) :: inverse_scale !! Coefficient of the reciprocal variate.
      real(dp), intent(in) :: scale !! Coefficient of the variate.
      real(dp), intent(in) :: current !! Positive current state.
      real(dp) :: log_current, log_level, left, right, proposal
      integer :: step

      if (inverse_scale < 0.0_dp .or. scale <= 0.0_dp .or. current <= 0.0_dp) then
         draw = max(current, tiny(1.0_dp))
         return
      end if
      log_current = log(current)
      log_level = gig_log_density(log_current, shape, inverse_scale, scale) + &
         log(max(random_uniform(), tiny(1.0_dp)))
      left = log_current - random_uniform()
      right = left + 1.0_dp
      do step = 1, 100
         if (gig_log_density(left, shape, inverse_scale, scale) <= log_level) exit
         left = left - 1.0_dp
      end do
      do step = 1, 100
         if (gig_log_density(right, shape, inverse_scale, scale) <= log_level) exit
         right = right + 1.0_dp
      end do
      do step = 1, 10000
         proposal = left + (right - left)*random_uniform()
         if (gig_log_density(proposal, shape, inverse_scale, scale) >= log_level) then
            draw = exp(proposal)
            return
         end if
         if (proposal < log_current) then
            left = proposal
         else
            right = proposal
         end if
      end do
      draw = current
   end function gig_slice_draw

   pure real(dp) function gig_log_density(log_value, shape, inverse_scale, &
      scale) result(log_density)
      !! Evaluate an unnormalized log-GIG density after a logarithmic transform.
      real(dp), intent(in) :: log_value !! Logarithm of the positive variate.
      real(dp), intent(in) :: shape !! GIG power parameter.
      real(dp), intent(in) :: inverse_scale !! Reciprocal coefficient.
      real(dp), intent(in) :: scale !! Linear coefficient.

      log_density = shape*log_value - 0.5_dp*( &
         inverse_scale*exp(-log_value) + scale*exp(log_value))
   end function gig_log_density

   pure function diagonal_from_vector(values) result(matrix)
      !! Construct a diagonal matrix from a vector.
      real(dp), intent(in) :: values(:) !! Diagonal values.
      real(dp) :: matrix(size(values), size(values))
      integer :: item

      matrix = 0.0_dp
      do item = 1, size(values)
         matrix(item, item) = values(item)
      end do
   end function diagonal_from_vector

end module bvars_mod
