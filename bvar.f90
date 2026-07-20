! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Distinct algorithms translated from the R BVAR package.
module bvar_mod
   !! Hierarchical conjugate priors and evidence calculations for Bayesian VARs.
   use kind_mod, only: dp
   use linalg_mod, only: identity_matrix, inverse_logdet, solve_matrix
   use random_mod, only: random_standard_normal_matrix, random_uniform
   implicit none
   private

   type, public :: bvar_dummy_observations_t
      !! Dummy response and regressor observations for a conjugate VAR prior.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: x(:, :)
      integer :: info = 0
   end type bvar_dummy_observations_t

   type, public :: bvar_evidence_t
      !! Conjugate posterior sufficient statistics and log marginal likelihood.
      real(dp), allocatable :: beta_hat(:, :)
      real(dp), allocatable :: sse(:, :)
      real(dp), allocatable :: posterior_scale(:, :)
      real(dp), allocatable :: coefficient_precision(:, :)
      real(dp) :: omega_log_determinant = 0.0_dp
      real(dp) :: psi_log_determinant = 0.0_dp
      real(dp) :: log_marginal_likelihood = -huge(1.0_dp)
      integer :: info = 0
   end type bvar_evidence_t

   type, public :: bvar_hyper_evaluation_t
      !! Marginal posterior evaluation at one vector of BVAR hyperparameters.
      type(bvar_evidence_t) :: evidence
      real(dp) :: log_posterior = -huge(1.0_dp)
      integer :: info = 0
   end type bvar_hyper_evaluation_t

   type, public :: bvar_metropolis_result_t
      !! Retained hierarchical-prior draws and Metropolis acceptance diagnostics.
      real(dp), allocatable :: hyperparameters(:, :)
      real(dp), allocatable :: log_posterior(:)
      logical, allocatable :: accepted(:)
      real(dp) :: acceptance_rate = 0.0_dp
      integer :: info = 0
   end type bvar_metropolis_result_t

   public :: bvar_minnesota_variance
   public :: bvar_soc_dummy
   public :: bvar_sur_dummy
   public :: bvar_conjugate_evidence
   public :: bvar_gamma_log_density
   public :: bvar_inverse_gamma_log_density
   public :: bvar_hyper_log_posterior
   public :: bvar_hierarchical_metropolis_from_random
   public :: bvar_hierarchical_metropolis

contains

   pure function bvar_minnesota_variance(innovation_scale, lags, lambda, alpha, &
      constant_variance) result(variance)
      !! Construct BVAR's diagonal Minnesota row variances.
      real(dp), intent(in) :: innovation_scale(:) !! Variable-specific innovation scales.
      integer, intent(in) :: lags !! VAR lag order.
      real(dp), intent(in) :: lambda !! Overall Minnesota tightness.
      real(dp), intent(in) :: alpha !! Lag-decay exponent.
      real(dp), intent(in) :: constant_variance !! Prior variance of the intercept row.
      real(dp), allocatable :: variance(:)
      integer :: variables, lag, first, last

      variables = size(innovation_scale)
      if (variables < 1 .or. lags < 0 .or. lambda <= 0.0_dp .or. &
         alpha < 0.0_dp .or. constant_variance <= 0.0_dp .or. &
         any(innovation_scale <= 0.0_dp)) then
         allocate(variance(0))
         return
      end if
      allocate(variance(1 + variables*lags))
      variance(1) = constant_variance
      do lag = 1, lags
         first = 2 + variables*(lag - 1)
         last = 1 + variables*lag
         variance(first:last) = lambda**2/ &
            (real(lag, dp)**alpha*innovation_scale)
      end do
   end function bvar_minnesota_variance

   pure function bvar_soc_dummy(y, lags, tightness) result(out)
      !! Construct sum-of-coefficients dummy observations.
      real(dp), intent(in) :: y(:, :) !! Variables by observation before lag trimming.
      integer, intent(in) :: lags !! VAR lag order.
      real(dp), intent(in) :: tightness !! Positive SOC prior tightness.
      type(bvar_dummy_observations_t) :: out
      real(dp), allocatable :: center(:)
      integer :: variables, lag, item, first

      variables = size(y, 1)
      if (variables < 1 .or. lags < 1 .or. size(y, 2) < lags .or. &
         tightness <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (lags == 1) then
         center = y(:, 1)/tightness
      else
         center = sum(y(:, :lags), dim=2)/(real(lags, dp)*tightness)
      end if
      allocate(out%y(variables, variables), &
         out%x(1 + variables*lags, variables))
      out%y = 0.0_dp
      out%x = 0.0_dp
      do item = 1, variables
         out%y(item, item) = center(item)
      end do
      do lag = 1, lags
         first = 2 + variables*(lag - 1)
         do item = 1, variables
            out%x(first + item - 1, item) = center(item)
         end do
      end do
   end function bvar_soc_dummy

   pure function bvar_sur_dummy(y, lags, tightness) result(out)
      !! Construct single-unit-root dummy observations.
      real(dp), intent(in) :: y(:, :) !! Variables by observation before lag trimming.
      integer, intent(in) :: lags !! VAR lag order.
      real(dp), intent(in) :: tightness !! Positive SUR prior tightness.
      type(bvar_dummy_observations_t) :: out
      real(dp), allocatable :: center(:)
      integer :: variables, lag, first, last

      variables = size(y, 1)
      if (variables < 1 .or. lags < 1 .or. size(y, 2) < lags .or. &
         tightness <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (lags == 1) then
         center = y(:, 1)/tightness
      else
         center = sum(y(:, :lags), dim=2)/(real(lags, dp)*tightness)
      end if
      allocate(out%y(variables, 1), out%x(1 + variables*lags, 1))
      out%y(:, 1) = center
      out%x(1, 1) = 1.0_dp/tightness
      do lag = 1, lags
         first = 2 + variables*(lag - 1)
         last = 1 + variables*lag
         out%x(first:last, 1) = center
      end do
   end function bvar_sur_dummy

   pure function bvar_conjugate_evidence(y, x, prior_mean, omega_variance, &
      psi_scale, estimate_beta) result(out)
      !! Compute conjugate BVAR sufficient statistics and log marginal likelihood.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: x(:, :) !! Regressors by observation.
      real(dp), intent(in) :: prior_mean(:, :) !! Minnesota coefficient prior mean.
      real(dp), intent(in) :: omega_variance(:) !! Diagonal coefficient row variances.
      real(dp), intent(in) :: psi_scale(:) !! Diagonal inverse-Wishart scale entries.
      logical, intent(in), optional :: estimate_beta !! Estimate coefficients instead of using their prior mean.
      type(bvar_evidence_t) :: out
      real(dp), allocatable :: cross(:, :), right(:, :), omega_inverse(:, :)
      real(dp), allocatable :: omega_sqrt(:, :), psi_inverse_sqrt(:, :)
      real(dp), allocatable :: omega_matrix(:, :), psi_matrix(:, :)
      real(dp), allocatable :: mostly_harmless(:, :), difference(:, :)
      real(dp), allocatable :: residual(:, :), unused_inverse(:, :)
      real(dp) :: omega_logdet, psi_logdet
      integer :: variables, observations, regressors, item, status
      logical :: estimate

      variables = size(y, 1)
      observations = size(y, 2)
      regressors = size(x, 1)
      estimate = .true.
      if (present(estimate_beta)) estimate = estimate_beta
      if (variables < 1 .or. observations < 1 .or. regressors < 1 .or. &
         size(x, 2) /= observations .or. &
         any(shape(prior_mean) /= [regressors, variables]) .or. &
         size(omega_variance) /= regressors .or. size(psi_scale) /= variables .or. &
         any(omega_variance <= 0.0_dp) .or. any(psi_scale <= 0.0_dp)) then
         out%info = 1
         return
      end if
      cross = matmul(x, transpose(x))
      allocate(omega_inverse(regressors, regressors), &
         omega_sqrt(regressors, regressors), &
         psi_inverse_sqrt(variables, variables), psi_matrix(variables, variables))
      omega_inverse = 0.0_dp
      omega_sqrt = 0.0_dp
      psi_inverse_sqrt = 0.0_dp
      psi_matrix = 0.0_dp
      do item = 1, regressors
         omega_inverse(item, item) = 1.0_dp/omega_variance(item)
         omega_sqrt(item, item) = sqrt(omega_variance(item))
      end do
      do item = 1, variables
         psi_inverse_sqrt(item, item) = 1.0_dp/sqrt(psi_scale(item))
         psi_matrix(item, item) = psi_scale(item)
      end do
      out%coefficient_precision = cross + omega_inverse
      right = matmul(x, transpose(y)) + matmul(omega_inverse, prior_mean)
      allocate(out%beta_hat(regressors, variables))
      if (estimate) then
         call solve_matrix(out%coefficient_precision, right, out%beta_hat, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
      else
         out%beta_hat = prior_mean
      end if
      residual = y - matmul(transpose(out%beta_hat), x)
      out%sse = matmul(residual, transpose(residual))
      difference = out%beta_hat - prior_mean
      mostly_harmless = out%sse + matmul(transpose(difference), &
         matmul(omega_inverse, difference))
      out%posterior_scale = psi_matrix + mostly_harmless
      omega_matrix = identity_matrix(regressors) + &
         matmul(omega_sqrt, matmul(cross, omega_sqrt))
      allocate(unused_inverse(regressors, regressors))
      call inverse_logdet(omega_matrix, unused_inverse, omega_logdet, status, &
         100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         out%info = 3
         return
      end if
      psi_matrix = identity_matrix(variables) + &
         matmul(psi_inverse_sqrt, matmul(mostly_harmless, psi_inverse_sqrt))
      deallocate(unused_inverse)
      allocate(unused_inverse(variables, variables))
      call inverse_logdet(psi_matrix, unused_inverse, psi_logdet, status, &
         100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         out%info = 4
         return
      end if
      out%omega_log_determinant = omega_logdet
      out%psi_log_determinant = psi_logdet
      out%log_marginal_likelihood = &
         -0.5_dp*real(variables*observations, dp)*log(acos(-1.0_dp)) + &
         sum([(log_gamma(0.5_dp*real(observations + variables + 3 - item, dp)) - &
         log_gamma(0.5_dp*real(variables + 3 - item, dp)), &
         item=1, variables)]) - &
         0.5_dp*real(observations, dp)*sum(log(psi_scale)) - &
         0.5_dp*real(variables, dp)*omega_logdet - &
         0.5_dp*real(observations + variables + 2, dp)*psi_logdet
   end function bvar_conjugate_evidence

   pure elemental function bvar_gamma_log_density(value, shape, scale) result(log_density)
      !! Evaluate a Gamma log density using its shape and scale parameterization.
      real(dp), intent(in) :: value !! Positive evaluation point.
      real(dp), intent(in) :: shape !! Positive Gamma shape.
      real(dp), intent(in) :: scale !! Positive Gamma scale.
      real(dp) :: log_density

      if (value <= 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
         log_density = -huge(1.0_dp)
      else
         log_density = (shape - 1.0_dp)*log(value) - value/scale - &
            log_gamma(shape) - shape*log(scale)
      end if
   end function bvar_gamma_log_density

   pure elemental function bvar_inverse_gamma_log_density(value, shape, scale) &
      result(log_density)
      !! Evaluate the inverse-Gamma log density used for BVAR innovation scales.
      real(dp), intent(in) :: value !! Positive evaluation point.
      real(dp), intent(in) :: shape !! Positive inverse-Gamma shape.
      real(dp), intent(in) :: scale !! Positive inverse-Gamma scale.
      real(dp) :: log_density

      if (value <= 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
         log_density = -huge(1.0_dp)
      else
         log_density = shape*log(scale) - (shape + 1.0_dp)*log(value) - &
            scale/value - log_gamma(shape)
      end if
   end function bvar_inverse_gamma_log_density

   pure function bvar_hyper_log_posterior(y, x, dummy_data, prior_mean, lags, &
      constant_variance, hyperparameters, lower, upper, prior_shape, prior_scale, &
      use_soc, use_sur) result(out)
      !! Evaluate BVAR's hierarchical marginal posterior, including dummy priors.
      real(dp), intent(in) :: y(:, :) !! Lag-trimmed endogenous variables by observation.
      real(dp), intent(in) :: x(:, :) !! VAR regressors by observation.
      real(dp), intent(in) :: dummy_data(:, :) !! Untrimmed data used to form dummy observations.
      real(dp), intent(in) :: prior_mean(:, :) !! Minnesota coefficient prior mean.
      integer, intent(in) :: lags !! VAR lag order.
      real(dp), intent(in) :: constant_variance !! Prior variance of the intercept row.
      real(dp), intent(in) :: hyperparameters(:) !! Lambda, alpha, psi values, then optional SOC and SUR tightnesses.
      real(dp), intent(in) :: lower(:) !! Inclusive lower bounds for the hyperparameters.
      real(dp), intent(in) :: upper(:) !! Inclusive upper bounds for the hyperparameters.
      real(dp), intent(in) :: prior_shape(:) !! Hyperprior shapes in hyperparameter order.
      real(dp), intent(in) :: prior_scale(:) !! Hyperprior scales in hyperparameter order.
      logical, intent(in) :: use_soc !! Include a hierarchical sum-of-coefficients dummy prior.
      logical, intent(in) :: use_sur !! Include a hierarchical single-unit-root dummy prior.
      type(bvar_hyper_evaluation_t) :: out
      type(bvar_dummy_observations_t) :: soc, sur
      type(bvar_evidence_t) :: dummy_evidence
      real(dp), allocatable :: variance(:), augmented_y(:, :), augmented_x(:, :)
      real(dp), allocatable :: dummy_y(:, :), dummy_x(:, :)
      integer :: variables, regressors, count, expected, first, item, soc_index, sur_index

      variables = size(y, 1)
      regressors = size(x, 1)
      expected = 2 + variables
      if (use_soc) expected = expected + 1
      if (use_sur) expected = expected + 1
      if (variables < 1 .or. lags < 1 .or. size(x, 2) /= size(y, 2) .or. &
         any(shape(prior_mean) /= [regressors, variables]) .or. &
         size(dummy_data, 1) /= variables .or. size(dummy_data, 2) < lags .or. &
         size(hyperparameters) /= expected .or. size(lower) /= expected .or. &
         size(upper) /= expected .or. size(prior_shape) /= expected .or. &
         size(prior_scale) /= expected .or. constant_variance <= 0.0_dp .or. &
         any(hyperparameters < lower) .or. any(hyperparameters > upper) .or. &
         any(prior_shape <= 0.0_dp) .or. any(prior_scale <= 0.0_dp)) then
         out%info = 1
         return
      end if

      variance = bvar_minnesota_variance(hyperparameters(3:2 + variables), &
         lags, hyperparameters(1), hyperparameters(2), constant_variance)
      if (size(variance) /= regressors) then
         out%info = 2
         return
      end if
      count = 0
      if (use_soc) count = count + variables
      if (use_sur) count = count + 1
      allocate(dummy_y(variables, count), dummy_x(regressors, count))
      first = 1
      soc_index = 3 + variables
      sur_index = soc_index
      if (use_soc) then
         soc = bvar_soc_dummy(dummy_data, lags, hyperparameters(soc_index))
         if (soc%info /= 0) then
            out%info = 3
            return
         end if
         dummy_y(:, first:first + variables - 1) = soc%y
         dummy_x(:, first:first + variables - 1) = soc%x
         first = first + variables
         sur_index = sur_index + 1
      end if
      if (use_sur) then
         sur = bvar_sur_dummy(dummy_data, lags, hyperparameters(sur_index))
         if (sur%info /= 0) then
            out%info = 4
            return
         end if
         dummy_y(:, first) = sur%y(:, 1)
         dummy_x(:, first) = sur%x(:, 1)
      end if

      allocate(augmented_y(variables, count + size(y, 2)))
      allocate(augmented_x(regressors, count + size(x, 2)))
      if (count > 0) then
         augmented_y(:, :count) = dummy_y
         augmented_x(:, :count) = dummy_x
      end if
      augmented_y(:, count + 1:) = y
      augmented_x(:, count + 1:) = x
      out%evidence = bvar_conjugate_evidence(augmented_y, augmented_x, prior_mean, &
         variance, hyperparameters(3:2 + variables))
      if (out%evidence%info /= 0) then
         out%info = 5
         return
      end if
      out%log_posterior = out%evidence%log_marginal_likelihood
      out%log_posterior = out%log_posterior + &
         bvar_gamma_log_density(hyperparameters(1), prior_shape(1), prior_scale(1)) + &
         bvar_gamma_log_density(hyperparameters(2), prior_shape(2), prior_scale(2))
      do item = 1, variables
         out%log_posterior = out%log_posterior + &
            bvar_inverse_gamma_log_density(hyperparameters(2 + item), &
            prior_shape(2 + item), prior_scale(2 + item))
      end do
      do item = 3 + variables, expected
         out%log_posterior = out%log_posterior + &
            bvar_gamma_log_density(hyperparameters(item), prior_shape(item), &
            prior_scale(item))
      end do
      if (count > 0) then
         dummy_evidence = bvar_conjugate_evidence(dummy_y, dummy_x, prior_mean, &
            variance, hyperparameters(3:2 + variables), .false.)
         if (dummy_evidence%info /= 0) then
            out%info = 6
            return
         end if
         out%log_posterior = out%log_posterior - &
            dummy_evidence%log_marginal_likelihood
      end if
   end function bvar_hyper_log_posterior

   pure function bvar_hierarchical_metropolis_from_random(y, x, dummy_data, &
      prior_mean, lags, constant_variance, initial_hyperparameters, lower, upper, &
      prior_shape, prior_scale, use_soc, use_sur, burnin, thin, proposal_factor, &
      standard_normals, acceptance_uniforms, adjust_acceptance, adjust_iterations, &
      acceptance_lower, acceptance_upper, acceptance_change) result(out)
      !! Draw bounded hierarchical BVAR parameters from supplied random variates.
      real(dp), intent(in) :: y(:, :) !! Lag-trimmed endogenous variables by observation.
      real(dp), intent(in) :: x(:, :) !! VAR regressors by observation.
      real(dp), intent(in) :: dummy_data(:, :) !! Untrimmed data used to form dummy observations.
      real(dp), intent(in) :: prior_mean(:, :) !! Minnesota coefficient prior mean.
      integer, intent(in) :: lags !! VAR lag order.
      real(dp), intent(in) :: constant_variance !! Prior variance of the intercept row.
      real(dp), intent(in) :: initial_hyperparameters(:) !! Initial bounded hyperparameter vector.
      real(dp), intent(in) :: lower(:) !! Inclusive lower bounds for the hyperparameters.
      real(dp), intent(in) :: upper(:) !! Inclusive upper bounds for the hyperparameters.
      real(dp), intent(in) :: prior_shape(:) !! Hyperprior shapes in hyperparameter order.
      real(dp), intent(in) :: prior_scale(:) !! Hyperprior scales in hyperparameter order.
      logical, intent(in) :: use_soc !! Include a hierarchical SOC dummy tightness.
      logical, intent(in) :: use_sur !! Include a hierarchical SUR dummy tightness.
      integer, intent(in) :: burnin !! Number of initial proposals to discard.
      integer, intent(in) :: thin !! Number of proposals between retained draws.
      real(dp), intent(in) :: proposal_factor(:, :) !! Factor mapping standard normals to proposal increments.
      real(dp), intent(in) :: standard_normals(:, :) !! Standard-normal proposal innovations.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Uniform variates used in acceptance tests.
      logical, intent(in), optional :: adjust_acceptance !! Tune proposal scale during burn-in.
      integer, intent(in), optional :: adjust_iterations !! Burn-in proposals eligible for tuning.
      real(dp), intent(in), optional :: acceptance_lower !! Lower target acceptance rate.
      real(dp), intent(in), optional :: acceptance_upper !! Upper target acceptance rate.
      real(dp), intent(in), optional :: acceptance_change !! Fractional proposal-variance adjustment.
      type(bvar_metropolis_result_t) :: out
      type(bvar_hyper_evaluation_t) :: current_evaluation, candidate_evaluation
      real(dp), allocatable :: current(:), candidate(:), factor(:, :)
      real(dp) :: lower_rate, upper_rate, change, rate
      integer :: hyper_count, total, saved, iteration, save_index, accepted_total
      integer :: accepted_adjustment, adjustment_limit
      logical :: accepted, adjust

      hyper_count = size(initial_hyperparameters)
      total = size(acceptance_uniforms)
      if (hyper_count < 1 .or. burnin < 0 .or. thin < 1 .or. total <= burnin .or. &
         mod(total - burnin, thin) /= 0 .or. &
         any(shape(proposal_factor) /= [hyper_count, hyper_count]) .or. &
         any(shape(standard_normals) /= [hyper_count, total]) .or. &
         any(acceptance_uniforms < 0.0_dp) .or. &
         any(acceptance_uniforms > 1.0_dp)) then
         out%info = 1
         return
      end if
      adjust = .false.
      if (present(adjust_acceptance)) adjust = adjust_acceptance
      adjustment_limit = (3*burnin)/4
      if (present(adjust_iterations)) adjustment_limit = adjust_iterations
      lower_rate = 0.25_dp
      if (present(acceptance_lower)) lower_rate = acceptance_lower
      upper_rate = 0.45_dp
      if (present(acceptance_upper)) upper_rate = acceptance_upper
      change = 0.01_dp
      if (present(acceptance_change)) change = acceptance_change
      if (adjustment_limit < 0 .or. adjustment_limit > burnin .or. &
         lower_rate < 0.0_dp .or. upper_rate > 1.0_dp .or. &
         lower_rate > upper_rate .or. change <= 0.0_dp .or. change >= 1.0_dp) then
         out%info = 2
         return
      end if
      current = initial_hyperparameters
      factor = proposal_factor
      current_evaluation = bvar_hyper_log_posterior(y, x, dummy_data, prior_mean, &
         lags, constant_variance, current, lower, upper, prior_shape, prior_scale, &
         use_soc, use_sur)
      if (current_evaluation%info /= 0) then
         out%info = 3
         return
      end if
      saved = (total - burnin)/thin
      allocate(out%hyperparameters(hyper_count, saved), out%log_posterior(saved))
      allocate(out%accepted(total))
      accepted_total = 0
      accepted_adjustment = 0
      save_index = 0
      do iteration = 1, total
         candidate = current + matmul(factor, standard_normals(:, iteration))
         candidate_evaluation = bvar_hyper_log_posterior(y, x, dummy_data, &
            prior_mean, lags, constant_variance, candidate, lower, upper, &
            prior_shape, prior_scale, use_soc, use_sur)
         accepted = .false.
         if (candidate_evaluation%info == 0) then
            accepted = log(max(acceptance_uniforms(iteration), tiny(1.0_dp))) < &
               candidate_evaluation%log_posterior - current_evaluation%log_posterior
         end if
         if (accepted) then
            current = candidate
            current_evaluation = candidate_evaluation
            if (iteration <= burnin) accepted_adjustment = accepted_adjustment + 1
            if (iteration > burnin) accepted_total = accepted_total + 1
         end if
         out%accepted(iteration) = accepted
         if (adjust .and. iteration <= adjustment_limit .and. mod(iteration, 10) == 0) then
            rate = real(accepted_adjustment, dp)/real(iteration, dp)
            if (rate < lower_rate) factor = sqrt(1.0_dp - change)*factor
            if (rate > upper_rate) factor = sqrt(1.0_dp + change)*factor
         end if
         if (iteration > burnin .and. mod(iteration - burnin, thin) == 0) then
            save_index = save_index + 1
            out%hyperparameters(:, save_index) = current
            out%log_posterior(save_index) = current_evaluation%log_posterior
         end if
      end do
      out%acceptance_rate = real(accepted_total, dp)/real(total - burnin, dp)
   end function bvar_hierarchical_metropolis_from_random

   function bvar_hierarchical_metropolis(y, x, dummy_data, prior_mean, lags, &
      constant_variance, initial_hyperparameters, lower, upper, prior_shape, &
      prior_scale, use_soc, use_sur, burnin, draws, thin, proposal_factor, &
      adjust_acceptance, adjust_iterations, acceptance_lower, acceptance_upper, &
      acceptance_change) result(out)
      !! Draw hierarchical BVAR parameters using the library's shared RNG.
      real(dp), intent(in) :: y(:, :) !! Lag-trimmed endogenous variables by observation.
      real(dp), intent(in) :: x(:, :) !! VAR regressors by observation.
      real(dp), intent(in) :: dummy_data(:, :) !! Untrimmed data used to form dummy observations.
      real(dp), intent(in) :: prior_mean(:, :) !! Minnesota coefficient prior mean.
      integer, intent(in) :: lags !! VAR lag order.
      real(dp), intent(in) :: constant_variance !! Prior variance of the intercept row.
      real(dp), intent(in) :: initial_hyperparameters(:) !! Initial bounded hyperparameter vector.
      real(dp), intent(in) :: lower(:) !! Inclusive lower bounds for the hyperparameters.
      real(dp), intent(in) :: upper(:) !! Inclusive upper bounds for the hyperparameters.
      real(dp), intent(in) :: prior_shape(:) !! Hyperprior shapes in hyperparameter order.
      real(dp), intent(in) :: prior_scale(:) !! Hyperprior scales in hyperparameter order.
      logical, intent(in) :: use_soc !! Include a hierarchical SOC dummy tightness.
      logical, intent(in) :: use_sur !! Include a hierarchical SUR dummy tightness.
      integer, intent(in) :: burnin !! Number of initial proposals to discard.
      integer, intent(in) :: draws !! Number of retained draws.
      integer, intent(in) :: thin !! Number of proposals between retained draws.
      real(dp), intent(in) :: proposal_factor(:, :) !! Factor mapping standard normals to proposal increments.
      logical, intent(in), optional :: adjust_acceptance !! Tune proposal scale during burn-in.
      integer, intent(in), optional :: adjust_iterations !! Burn-in proposals eligible for tuning.
      real(dp), intent(in), optional :: acceptance_lower !! Lower target acceptance rate.
      real(dp), intent(in), optional :: acceptance_upper !! Upper target acceptance rate.
      real(dp), intent(in), optional :: acceptance_change !! Fractional proposal-variance adjustment.
      type(bvar_metropolis_result_t) :: out
      real(dp), allocatable :: normals(:, :), uniforms(:)
      integer :: total, item

      if (burnin < 0 .or. draws < 1 .or. thin < 1) then
         out%info = 1
         return
      end if
      total = burnin + draws*thin
      allocate(normals(size(initial_hyperparameters), total), uniforms(total))
      call random_standard_normal_matrix(normals)
      do item = 1, total
         uniforms(item) = random_uniform()
      end do
      out = bvar_hierarchical_metropolis_from_random(y, x, dummy_data, prior_mean, &
         lags, constant_variance, initial_hyperparameters, lower, upper, prior_shape, &
         prior_scale, use_soc, use_sur, burnin, thin, proposal_factor, normals, &
         uniforms, adjust_acceptance, adjust_iterations, acceptance_lower, &
         acceptance_upper, acceptance_change)
   end function bvar_hierarchical_metropolis

end module bvar_mod
