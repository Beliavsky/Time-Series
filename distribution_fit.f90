! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Distribution fitting translated from the R rugarch package.
module distribution_fit_mod
   !! Maximum-likelihood fitting for standardized innovation distributions.
   use kind_mod, only: dp
   use distribution_mod, only: standardized_log_density, &
      distribution_has_skew, distribution_has_shape, distribution_has_lambda, &
      distribution_skew_normal, distribution_skew_student, &
      distribution_skew_ged, distribution_johnson_su, distribution_nig, &
      distribution_ghyp, distribution_gh_skew_student, distribution_student, &
      distribution_ged
   use linalg_mod, only: symmetric_pseudoinverse
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd, &
      finite_difference_hessian
   use stats_mod, only: variance
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   type, public :: distribution_fit_t
      !! Maximum-likelihood location, scale, skew, shape, and lambda estimates.
      integer :: distribution = 0
      real(dp) :: location = 0.0_dp
      real(dp) :: scale = 1.0_dp
      real(dp) :: skew = 1.0_dp
      real(dp) :: shape = 0.0_dp
      real(dp) :: lambda = -0.5_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: standard_errors(:)
      real(dp), allocatable :: covariance(:, :)
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: converged = .false.
   end type distribution_fit_t

   interface display
      module procedure display_distribution_fit
   end interface display

   public :: fit_distribution, display

contains

   pure function fit_distribution(data, distribution, initial_shape, &
      initial_skew, initial_lambda, max_iterations, tolerance) result(out)
      !! Estimate a location-scale innovation distribution by maximum likelihood.
      real(dp), intent(in) :: data(:) !! Finite sample observations.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in), optional :: initial_shape !! Optional starting shape.
      real(dp), intent(in), optional :: initial_skew !! Optional starting skew.
      real(dp), intent(in), optional :: initial_lambda !! Optional starting GH lambda.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(distribution_fit_t) :: out
      type(optimization_result_t) :: optimum
      real(dp), allocatable :: raw(:), hessian(:, :), inverse(:, :)
      real(dp), allocatable :: jacobian(:, :), physical_covariance(:, :)
      real(dp) :: selected_tolerance
      integer :: count, limit, inverse_info, parameter

      out%distribution = distribution
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (size(data) < 5 .or. .not. all(ieee_is_finite(data)) .or. &
         distribution < 1 .or. distribution > 10 .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      raw = initial_raw(data, distribution, initial_shape, initial_skew, &
         initial_lambda)
      optimum = bfgs_minimize_fd(objective, raw, limit, selected_tolerance)
      if (.not. allocated(optimum%parameters)) then
         out%info = 2
         return
      end if
      call unpack_parameters(distribution, optimum%parameters, out%location, &
         out%scale, out%shape, out%skew, out%lambda)
      out%coefficients = physical_coefficients(out)
      out%log_likelihood = -objective(optimum%parameters)
      out%iterations = optimum%iterations
      out%optimizer_info = optimum%info
      out%converged = optimum%converged
      if (.not. ieee_is_finite(out%log_likelihood)) then
         out%info = 3
         return
      end if
      hessian = finite_difference_hessian(objective, optimum%parameters)
      allocate(inverse(size(hessian, 1), size(hessian, 2)))
      call symmetric_pseudoinverse(hessian, inverse, inverse_info)
      count = size(out%coefficients)
      if (inverse_info == 0) then
         jacobian = physical_jacobian(distribution, optimum%parameters)
         physical_covariance = matmul(matmul(jacobian, inverse), &
            transpose(jacobian))
         out%covariance = physical_covariance
         allocate(out%standard_errors(count))
         do parameter = 1, count
            out%standard_errors(parameter) = sqrt(max(physical_covariance( &
               parameter, parameter), 0.0_dp))
         end do
      else
         allocate(out%covariance(0, 0), out%standard_errors(0))
      end if

   contains

      pure real(dp) function objective(candidate) result(value)
         !! Evaluate the negative distribution log likelihood.
         real(dp), intent(in) :: candidate(:) !! Unconstrained parameters.
         real(dp) :: location, scale, shape, skew, lambda

         call unpack_parameters(distribution, candidate, location, scale, &
            shape, skew, lambda)
         value = -sum(standardized_log_density((data - location)/scale, &
            distribution, shape, skew, lambda) - log(scale))
         if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
      end function objective

   end function fit_distribution

   pure function initial_raw(data, distribution, initial_shape, initial_skew, &
      initial_lambda) result(raw)
      !! Construct unconstrained starting coordinates for distribution fitting.
      real(dp), intent(in) :: data(:) !! Sample observations.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in), optional :: initial_shape !! Optional starting shape.
      real(dp), intent(in), optional :: initial_skew !! Optional starting skew.
      real(dp), intent(in), optional :: initial_lambda !! Optional starting lambda.
      real(dp), allocatable :: raw(:)
      real(dp) :: shape, skew, lambda
      integer :: index

      shape = 8.0_dp
      if (distribution == distribution_ged .or. &
         distribution == distribution_skew_ged) shape = 1.5_dp
      if (distribution == distribution_johnson_su) shape = 1.5_dp
      if (distribution == distribution_nig .or. &
         distribution == distribution_ghyp) shape = 3.0_dp
      if (present(initial_shape)) shape = initial_shape
      skew = 1.0_dp
      if (distribution == distribution_johnson_su .or. &
         distribution == distribution_nig .or. &
         distribution == distribution_ghyp .or. &
         distribution == distribution_gh_skew_student) skew = 0.0_dp
      if (present(initial_skew)) skew = initial_skew
      lambda = -0.5_dp
      if (present(initial_lambda)) lambda = initial_lambda
      allocate(raw(2 + merge(1, 0, distribution_has_skew(distribution)) + &
         merge(1, 0, distribution_has_shape(distribution)) + &
         merge(1, 0, distribution_has_lambda(distribution))))
      raw(1) = sum(data)/real(size(data), dp)
      raw(2) = log(sqrt(max(variance(data), tiny(1.0_dp))))
      index = 3
      if (distribution_has_skew(distribution)) then
         select case (distribution)
         case (distribution_skew_normal, distribution_skew_student, &
            distribution_skew_ged)
            raw(index) = log(max(skew, 1.0e-6_dp))
         case (distribution_nig, distribution_ghyp)
            raw(index) = atanh(max(-0.999_dp, min(0.999_dp, skew/0.999_dp)))
         case default
            raw(index) = skew
         end select
         index = index + 1
      end if
      if (distribution_has_shape(distribution)) then
         select case (distribution)
         case (distribution_student, distribution_skew_student)
            raw(index) = log(max(shape - 2.01_dp, 1.0e-6_dp))
         case (distribution_ged, distribution_skew_ged)
            raw(index) = log(max(shape - 0.20_dp, 1.0e-6_dp))
         case (distribution_johnson_su)
            raw(index) = logit((shape - 0.20_dp)/19.80_dp)
         case (distribution_gh_skew_student)
            raw(index) = log(max(shape - 4.01_dp, 1.0e-6_dp))
         case default
            raw(index) = log(max(shape - 0.10_dp, 1.0e-6_dp))
         end select
         index = index + 1
      end if
      if (distribution_has_lambda(distribution)) raw(index) = lambda
   end function initial_raw

   pure subroutine unpack_parameters(distribution, raw, location, scale, shape, &
      skew, lambda)
      !! Transform unconstrained coordinates to physical distribution parameters.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: raw(:) !! Unconstrained parameters.
      real(dp), intent(out) :: location !! Location parameter.
      real(dp), intent(out) :: scale !! Positive scale parameter.
      real(dp), intent(out) :: shape !! Tail or shape parameter.
      real(dp), intent(out) :: skew !! Skew parameter.
      real(dp), intent(out) :: lambda !! Generalized-hyperbolic lambda.
      integer :: index

      location = raw(1)
      scale = exp(min(raw(2), 50.0_dp))
      shape = 0.0_dp
      skew = 1.0_dp
      lambda = -0.5_dp
      index = 3
      if (distribution_has_skew(distribution)) then
         select case (distribution)
         case (distribution_skew_normal, distribution_skew_student, &
            distribution_skew_ged)
            skew = exp(min(raw(index), 20.0_dp))
         case (distribution_nig, distribution_ghyp)
            skew = 0.999_dp*tanh(raw(index))
         case default
            skew = raw(index)
         end select
         index = index + 1
      end if
      if (distribution_has_shape(distribution)) then
         select case (distribution)
         case (distribution_student, distribution_skew_student)
            shape = 2.01_dp + exp(min(raw(index), 20.0_dp))
         case (distribution_ged, distribution_skew_ged)
            shape = 0.20_dp + exp(min(raw(index), 20.0_dp))
         case (distribution_johnson_su)
            shape = 0.20_dp + 19.80_dp*logistic(raw(index))
         case (distribution_gh_skew_student)
            shape = 4.01_dp + exp(min(raw(index), 20.0_dp))
         case default
            shape = 0.10_dp + exp(min(raw(index), 20.0_dp))
         end select
         index = index + 1
      end if
      if (distribution_has_lambda(distribution)) lambda = raw(index)
   end subroutine unpack_parameters

   pure function physical_coefficients(fit) result(coefficients)
      !! Pack the included physical parameters in fitting order.
      type(distribution_fit_t), intent(in) :: fit !! Fitted distribution.
      real(dp), allocatable :: coefficients(:)

      coefficients = [fit%location, fit%scale]
      if (distribution_has_skew(fit%distribution)) &
         coefficients = [coefficients, fit%skew]
      if (distribution_has_shape(fit%distribution)) &
         coefficients = [coefficients, fit%shape]
      if (distribution_has_lambda(fit%distribution)) &
         coefficients = [coefficients, fit%lambda]
   end function physical_coefficients

   pure function physical_jacobian(distribution, raw) result(jacobian)
      !! Numerically differentiate physical parameters by raw coordinates.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: raw(:) !! Unconstrained parameters.
      real(dp), allocatable :: jacobian(:, :)
      real(dp) :: shifted(size(raw)), step
      real(dp), allocatable :: upper(:), lower(:)
      integer :: parameter

      upper = raw_physical_coefficients(distribution, raw)
      allocate(jacobian(size(upper), size(raw)))
      do parameter = 1, size(raw)
         step = epsilon(1.0_dp)**(1.0_dp/3.0_dp)* &
            max(1.0_dp, abs(raw(parameter)))
         shifted = raw
         shifted(parameter) = raw(parameter) + step
         upper = raw_physical_coefficients(distribution, shifted)
         shifted(parameter) = raw(parameter) - step
         lower = raw_physical_coefficients(distribution, shifted)
         jacobian(:, parameter) = (upper - lower)/(2.0_dp*step)
      end do
   end function physical_jacobian

   pure function raw_physical_coefficients(distribution, raw) result(values)
      !! Transform and pack one raw distribution parameter vector.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: raw(:) !! Unconstrained parameters.
      real(dp), allocatable :: values(:)
      real(dp) :: location, scale, shape, skew, lambda

      call unpack_parameters(distribution, raw, location, scale, shape, skew, &
         lambda)
      values = [location, scale]
      if (distribution_has_skew(distribution)) values = [values, skew]
      if (distribution_has_shape(distribution)) values = [values, shape]
      if (distribution_has_lambda(distribution)) values = [values, lambda]
   end function raw_physical_coefficients

   pure elemental real(dp) function logistic(value) result(transformed)
      !! Apply a stable logistic transform.
      real(dp), intent(in) :: value !! Unconstrained value.

      transformed = 1.0_dp/(1.0_dp + exp(max(-50.0_dp, min(50.0_dp, -value))))
   end function logistic

   pure elemental real(dp) function logit(probability) result(value)
      !! Apply a clipped inverse logistic transform.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp) :: clipped

      clipped = max(1.0e-8_dp, min(1.0_dp - 1.0e-8_dp, probability))
      value = log(clipped/(1.0_dp - clipped))
   end function logit

   subroutine display_distribution_fit(value, unit)
      !! Display a fitted innovation distribution.
      type(distribution_fit_t), intent(in) :: value !! Distribution fit.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'distribution fit'
      write(output, '(a,i0)') 'distribution: ', value%distribution
      write(output, '(a,es14.6)') 'location: ', value%location
      write(output, '(a,es14.6)') 'scale: ', value%scale
      write(output, '(a,es14.6)') 'skew: ', value%skew
      write(output, '(a,es14.6)') 'shape: ', value%shape
      write(output, '(a,es14.6)') 'lambda: ', value%lambda
      write(output, '(a,es14.6)') 'log likelihood: ', value%log_likelihood
      write(output, '(a,l1)') 'converged: ', value%converged
      write(output, '(a,i0)') 'info: ', value%info
   end subroutine display_distribution_fit

end module distribution_fit_mod
