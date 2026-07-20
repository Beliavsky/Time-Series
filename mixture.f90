! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original finite-mixture infrastructure for this Fortran library.
module mixture_mod
   !! Finite-mixture estimation algorithms reusable across statistical domains.
   use kind_mod, only: dp
   use stats_mod, only: sorted, quantile, variance
   implicit none
   private

   real(dp), parameter :: log_two_pi = 1.83787706640934548356_dp

   type, public :: gaussian_mixture_t
      !! Fitted univariate Gaussian-mixture parameters and diagnostics.
      real(dp), allocatable :: weight(:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: variance(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: info = 0
   end type gaussian_mixture_t

   public :: fit_gaussian_mixture

contains

   pure function fit_gaussian_mixture(values, components, max_iterations, &
      tolerance, minimum_variance) result(out)
      !! Fit a univariate Gaussian mixture by deterministic EM iterations.
      real(dp), intent(in) :: values(:) !! Observed scalar values.
      integer, intent(in) :: components !! Positive number of mixture components.
      integer, intent(in), optional :: max_iterations !! Maximum number of EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative log-likelihood convergence tolerance.
      real(dp), intent(in), optional :: minimum_variance !! Positive component variance floor.
      type(gaussian_mixture_t) :: out
      real(dp), allocatable :: ordered(:), effective_count(:), weighted_sum(:)
      real(dp), allocatable :: weighted_square_sum(:), log_density(:)
      real(dp), allocatable :: new_mean(:), new_variance(:), new_weight(:)
      real(dp) :: tolerance_value, variance_floor, total_variance
      real(dp) :: maximum_log_density, density_sum, responsibility
      real(dp) :: log_likelihood, previous_log_likelihood, probability
      integer :: iterations_limit, observation, component, iteration

      iterations_limit = 200
      tolerance_value = 1.0e-8_dp
      if (present(max_iterations)) iterations_limit = max_iterations
      if (present(tolerance)) tolerance_value = tolerance
      if (size(values) < 2 .or. components < 1 .or. &
         components > size(values) .or. iterations_limit < 1 .or. &
         tolerance_value <= 0.0_dp) then
         out%info = 1
         return
      end if
      ordered = sorted(values)
      total_variance = max(variance(values), epsilon(1.0_dp))
      variance_floor = max(1.0e-8_dp*total_variance, epsilon(1.0_dp))
      if (present(minimum_variance)) variance_floor = minimum_variance
      if (variance_floor <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%weight(components), out%mean(components))
      allocate(out%variance(components))
      do component = 1, components
         probability = (real(component, dp) - 0.5_dp)/real(components, dp)
         out%mean(component) = quantile(ordered, probability)
      end do
      out%weight = 1.0_dp/real(components, dp)
      out%variance = total_variance
      allocate(effective_count(components), weighted_sum(components))
      allocate(weighted_square_sum(components), log_density(components))
      allocate(new_mean(components), new_variance(components), new_weight(components))
      previous_log_likelihood = -huge(1.0_dp)
      do iteration = 1, iterations_limit
         effective_count = 0.0_dp
         weighted_sum = 0.0_dp
         weighted_square_sum = 0.0_dp
         log_likelihood = 0.0_dp
         do observation = 1, size(values)
            do component = 1, components
               log_density(component) = log(max(out%weight(component), &
                  tiny(1.0_dp))) - 0.5_dp*(log_two_pi + &
                  log(out%variance(component)) + (values(observation) - &
                  out%mean(component))**2/out%variance(component))
            end do
            maximum_log_density = maxval(log_density)
            density_sum = sum(exp(log_density - maximum_log_density))
            log_likelihood = log_likelihood + maximum_log_density + &
               log(density_sum)
            do component = 1, components
               responsibility = exp(log_density(component) - &
                  maximum_log_density)/density_sum
               effective_count(component) = effective_count(component) + &
                  responsibility
               weighted_sum(component) = weighted_sum(component) + &
                  responsibility*values(observation)
               weighted_square_sum(component) = weighted_square_sum(component) + &
                  responsibility*values(observation)**2
            end do
         end do
         do component = 1, components
            if (effective_count(component) <= tiny(1.0_dp)) then
               new_weight(component) = 1.0_dp/real(size(values), dp)
               new_mean(component) = out%mean(component)
               new_variance(component) = total_variance
            else
               new_weight(component) = effective_count(component)/ &
                  real(size(values), dp)
               new_mean(component) = weighted_sum(component)/ &
                  effective_count(component)
               new_variance(component) = max(variance_floor, &
                  weighted_square_sum(component)/effective_count(component) - &
                  new_mean(component)**2)
            end if
         end do
         new_weight = new_weight/sum(new_weight)
         out%weight = new_weight
         out%mean = new_mean
         out%variance = new_variance
         out%iterations = iteration
         out%log_likelihood = log_likelihood
         if (iteration > 1) then
            if (abs(log_likelihood - previous_log_likelihood) <= &
               tolerance_value*(1.0_dp + abs(previous_log_likelihood))) exit
         end if
         previous_log_likelihood = log_likelihood
      end do
      call order_gaussian_components(out%weight, out%mean, out%variance)
   end function fit_gaussian_mixture

   pure subroutine order_gaussian_components(weight, mean, component_variance)
      !! Order mixture components by ascending mean for stable output.
      real(dp), intent(inout) :: weight(:) !! Component weights, reordered in place.
      real(dp), intent(inout) :: mean(:) !! Component means, reordered in place.
      real(dp), intent(inout) :: component_variance(:) !! Component variances, reordered in place.
      real(dp) :: held_weight, held_mean, held_variance
      integer :: component, position

      do component = 2, size(mean)
         held_weight = weight(component)
         held_mean = mean(component)
         held_variance = component_variance(component)
         position = component - 1
         do while (position >= 1)
            if (mean(position) <= held_mean) exit
            weight(position + 1) = weight(position)
            mean(position + 1) = mean(position)
            component_variance(position + 1) = component_variance(position)
            position = position - 1
         end do
         weight(position + 1) = held_weight
         mean(position + 1) = held_mean
         component_variance(position + 1) = held_variance
      end do
   end subroutine order_gaussian_components

end module mixture_mod
