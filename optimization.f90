! SPDX-License-Identifier: MIT
! SPDX-FileComment: Optimization routines adapted from GARCH-BFGS.
! Shared numerical optimization algorithms for time-series estimation.
! Adapted from GARCH-BFGS/bfgs.f90, Copyright (c) 2026 Beliavsky, MIT licensed.
module optimization_mod
   use kind_mod, only: dp
   use linalg_mod, only: outer_product
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   type, public :: optimization_result_t
      ! Optimizer parameters, objective value, convergence, and status.
      real(dp), allocatable :: parameters(:)
      real(dp) :: objective = huge(1.0_dp)
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type

   abstract interface
      pure function objective_function_t(parameters) result(value)
         !! Evaluate a scalar optimization objective.
         import :: dp
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: value
      end function objective_function_t

      pure function gradient_function_t(parameters) result(gradient)
         !! Evaluate an optimization-objective gradient.
         import :: dp
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: gradient(size(parameters))
      end function gradient_function_t
   end interface

   public :: bfgs_minimize, bfgs_minimize_fd
   public :: nelder_mead_minimize, finite_difference_hessian

contains

   pure function bfgs_minimize(objective, gradient_function, initial, &
      max_iterations, gradient_tolerance) result(out)
      !! Minimize a pure objective using BFGS and an analytic gradient.
      procedure(objective_function_t) :: objective !! Objective callback procedure.
      procedure(gradient_function_t) :: gradient_function !! Gradient callback procedure.
      real(dp), intent(in) :: initial(:) !! Initial parameter values.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: gradient_tolerance !! Gradient convergence tolerance.
      type(optimization_result_t) :: out
      real(dp), allocatable :: x(:), candidate(:), gradient(:), new_gradient(:)
      real(dp), allocatable :: direction(:), inverse_hessian(:, :), identity(:, :)
      real(dp), allocatable :: s(:), y(:), update_left(:, :)
      real(dp) :: value, new_value, tolerance, sy, rho
      integer :: limit, iteration, n, i
      logical :: accepted

      n = size(initial)
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      tolerance = 1.0e-6_dp
      if (present(gradient_tolerance)) tolerance = gradient_tolerance
      if (n < 1 .or. limit < 1 .or. tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(x(n), candidate(n), gradient(n), new_gradient(n), direction(n))
      allocate(inverse_hessian(n, n), identity(n, n), s(n), y(n))
      allocate(update_left(n, n))
      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
      inverse_hessian = identity
      x = initial
      value = objective(x)
      gradient = gradient_function(x)
      if (.not. ieee_is_finite(value) .or. &
         .not. all(ieee_is_finite(gradient))) then
         out%info = 2
         out%parameters = x
         return
      end if
      do iteration = 1, limit
         out%iterations = iteration
         if (norm2(gradient) <= tolerance) then
            out%converged = .true.
            exit
         end if
         direction = -matmul(inverse_hessian, gradient)
         if (dot_product(gradient, direction) >= 0.0_dp) then
            inverse_hessian = identity
            direction = -gradient
         end if
         call armijo_search(objective, x, value, gradient, direction, &
            candidate, new_value, accepted)
         if (.not. accepted) then
            out%info = 3
            exit
         end if
         new_gradient = gradient_function(candidate)
         if (.not. all(ieee_is_finite(new_gradient))) then
            out%info = 2
            exit
         end if
         s = candidate - x
         y = new_gradient - gradient
         sy = dot_product(s, y)
         x = candidate
         value = new_value
         gradient = new_gradient
         if (sy > 1.0e-10_dp*norm2(s)*norm2(y)) then
            rho = 1.0_dp/sy
            update_left = identity - rho*outer_product(s, y)
            inverse_hessian = matmul(matmul(update_left, inverse_hessian), &
               transpose(update_left)) + rho*outer_product(s, s)
         end if
      end do
      if (.not. out%converged .and. out%info == 0) out%info = 4
      out%parameters = x
      out%objective = value
   end function bfgs_minimize

   pure function nelder_mead_minimize(objective, initial, max_iterations, &
      tolerance, initial_step) result(out)
      !! Minimize a pure scalar objective with the derivative-free Nelder-Mead simplex.
      procedure(objective_function_t) :: objective !! Objective callback procedure.
      real(dp), intent(in) :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: initial_step !! Initial step.
      type(optimization_result_t) :: out
      real(dp), allocatable :: simplex(:, :), values(:), centroid(:), trial(:)
      real(dp) :: tol, step, reflected_value, expanded_value, contracted_value
      integer :: n, limit, iteration, i, best, worst, second_worst

      n = size(initial)
      limit = 1000
      if (present(max_iterations)) limit = max_iterations
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = tolerance
      step = 0.1_dp
      if (present(initial_step)) step = initial_step
      if (n < 1 .or. limit < 1 .or. tol <= 0.0_dp .or. step <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(simplex(n, n + 1), values(n + 1), centroid(n), trial(n))
      simplex = spread(initial, 2, n + 1)
      do i = 1, n
         simplex(i, i + 1) = simplex(i, i + 1) + step*max(1.0_dp, abs(initial(i)))
      end do
      do i = 1, n + 1
         values(i) = objective(simplex(:, i))
      end do
      do iteration = 1, limit
         call simplex_extrema(values, best, worst, second_worst)
         out%iterations = iteration
         if (maxval(abs(values - values(best))) <= tol*max(1.0_dp, abs(values(best))) .and. &
            maxval(abs(simplex - spread(simplex(:, best), 2, n + 1))) <= &
            sqrt(tol)*max(1.0_dp, maxval(abs(simplex(:, best))))) then
            out%converged = .true.
            exit
         end if
         centroid = (sum(simplex, dim=2) - simplex(:, worst))/real(n, dp)
         trial = centroid + (centroid - simplex(:, worst))
         reflected_value = objective(trial)
         if (reflected_value < values(best)) then
            trial = centroid + 2.0_dp*(centroid - simplex(:, worst))
            expanded_value = objective(trial)
            if (expanded_value < reflected_value) then
               simplex(:, worst) = trial
               values(worst) = expanded_value
            else
               simplex(:, worst) = centroid + (centroid - simplex(:, worst))
               values(worst) = reflected_value
            end if
         else if (reflected_value < values(second_worst)) then
            simplex(:, worst) = trial
            values(worst) = reflected_value
         else
            if (reflected_value < values(worst)) then
               trial = centroid + 0.5_dp*(centroid - simplex(:, worst))
            else
               trial = centroid + 0.5_dp*(simplex(:, worst) - centroid)
            end if
            contracted_value = objective(trial)
            if (contracted_value < min(reflected_value, values(worst))) then
               simplex(:, worst) = trial
               values(worst) = contracted_value
            else
               do i = 1, n + 1
                  if (i == best) cycle
                  simplex(:, i) = simplex(:, best) + 0.5_dp*(simplex(:, i) - simplex(:, best))
                  values(i) = objective(simplex(:, i))
               end do
            end if
         end if
      end do
      call simplex_extrema(values, best, worst, second_worst)
      out%parameters = simplex(:, best)
      out%objective = values(best)
      if (.not. out%converged) out%info = 4
   end function nelder_mead_minimize

   pure subroutine simplex_extrema(values, best, worst, second_worst)
      !! Locate the best and two worst simplex vertices.
      real(dp), intent(in) :: values(:) !! Input values.
      integer, intent(out) :: best !! Best.
      integer, intent(out) :: worst !! Worst.
      integer, intent(out) :: second_worst !! Second worst.
      integer :: i

      best = minloc(values, dim=1)
      worst = maxloc(values, dim=1)
      second_worst = best
      do i = 1, size(values)
         if (i == worst) cycle
         if (second_worst == best .or. values(i) > values(second_worst)) second_worst = i
      end do
   end subroutine simplex_extrema

   pure function bfgs_minimize_fd(objective, initial, max_iterations, gradient_tolerance, &
      difference_step) result(out)
      !! Minimize a pure scalar objective using finite-difference BFGS.
      procedure(objective_function_t) :: objective !! Objective callback procedure.
      real(dp), intent(in) :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: gradient_tolerance !! Gradient tolerance.
      real(dp), intent(in), optional :: difference_step !! Difference step.
      type(optimization_result_t) :: out
      real(dp), allocatable :: x(:), candidate(:), gradient(:), new_gradient(:), direction(:)
      real(dp), allocatable :: inverse_hessian(:, :), identity(:, :), s(:), y(:), a(:, :)
      real(dp) :: value, new_value, tolerance, step, sy, rho
      integer :: limit, iteration, n, i
      logical :: accepted

      n = size(initial)
      if (n < 1) then
         out%info = 1
         return
      end if
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      tolerance = 1.0e-6_dp
      if (present(gradient_tolerance)) tolerance = gradient_tolerance
      step = 1.0e-5_dp
      if (present(difference_step)) step = difference_step
      if (limit < 1 .or. tolerance <= 0.0_dp .or. step <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(x(n), candidate(n), gradient(n), new_gradient(n), direction(n))
      allocate(inverse_hessian(n, n), identity(n, n), s(n), y(n), a(n, n))
      x = initial
      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
      inverse_hessian = identity
      value = objective(x)
      if (.not. ieee_is_finite(value)) then
         out%info = 2
         out%parameters = x
         return
      end if
      gradient = central_gradient(objective, x, step)
      do iteration = 1, limit
         out%iterations = iteration
         if (norm2(gradient) <= tolerance) then
            out%converged = .true.
            exit
         end if
         direction = -matmul(inverse_hessian, gradient)
         if (dot_product(gradient, direction) >= 0.0_dp) then
            inverse_hessian = identity
            direction = -gradient
         end if
         call armijo_search(objective, x, value, gradient, direction, candidate, new_value, accepted)
         if (.not. accepted) then
            out%info = 3
            exit
         end if
         new_gradient = central_gradient(objective, candidate, step)
         s = candidate - x
         y = new_gradient - gradient
         sy = dot_product(s, y)
         x = candidate
         value = new_value
         gradient = new_gradient
         if (sy > 1.0e-10_dp*norm2(s)*norm2(y)) then
            rho = 1.0_dp/sy
            a = identity - rho*outer_product(s, y)
            inverse_hessian = matmul(matmul(a, inverse_hessian), transpose(a)) + &
               rho*outer_product(s, s)
         end if
      end do
      if (.not. out%converged .and. out%info == 0) out%info = 4
      out%parameters = x
      out%objective = value
   end function bfgs_minimize_fd

   pure function central_gradient(objective, parameters, step) result(gradient)
      !! Compute a scale-aware central finite-difference gradient.
      procedure(objective_function_t) :: objective !! Objective callback procedure.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: step !! Step.
      real(dp) :: gradient(size(parameters)), shifted(size(parameters)), h, plus, minus
      integer :: i

      do i = 1, size(parameters)
         h = step*max(1.0_dp, abs(parameters(i)))
         shifted = parameters
         shifted(i) = parameters(i) + h
         plus = objective(shifted)
         shifted(i) = parameters(i) - h
         minus = objective(shifted)
         if (ieee_is_finite(plus) .and. ieee_is_finite(minus)) then
            gradient(i) = (plus - minus)/(2.0_dp*h)
         else
            gradient(i) = 0.0_dp
         end if
      end do
   end function central_gradient

   pure function finite_difference_hessian(objective, parameters, difference_step) result(hessian)
      !! Compute a symmetric finite-difference Hessian of a scalar objective.
      procedure(objective_function_t) :: objective !! Objective callback procedure.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in), optional :: difference_step !! Difference step.
      real(dp) :: hessian(size(parameters), size(parameters))
      real(dp) :: shifted(size(parameters)), hi, hj, center, fpp, fpm, fmp, fmm, step
      integer :: i, j

      step = epsilon(1.0_dp)**0.25_dp
      if (present(difference_step)) step = difference_step
      center = objective(parameters)
      hessian = 0.0_dp
      do i = 1, size(parameters)
         hi = step*max(1.0_dp, abs(parameters(i)))
         shifted = parameters
         shifted(i) = parameters(i) + hi
         fpp = objective(shifted)
         shifted(i) = parameters(i) - hi
         fmm = objective(shifted)
         hessian(i, i) = (fpp - 2.0_dp*center + fmm)/(hi*hi)
         do j = i + 1, size(parameters)
            hj = step*max(1.0_dp, abs(parameters(j)))
            shifted = parameters
            shifted(i) = parameters(i) + hi
            shifted(j) = parameters(j) + hj
            fpp = objective(shifted)
            shifted(j) = parameters(j) - hj
            fpm = objective(shifted)
            shifted(i) = parameters(i) - hi
            shifted(j) = parameters(j) + hj
            fmp = objective(shifted)
            shifted(j) = parameters(j) - hj
            fmm = objective(shifted)
            hessian(i, j) = (fpp - fpm - fmp + fmm)/(4.0_dp*hi*hj)
            hessian(j, i) = hessian(i, j)
         end do
      end do
   end function finite_difference_hessian

   pure subroutine armijo_search(objective, x, value, gradient, direction, candidate, &
      candidate_value, accepted)
      !! Backtrack until the Armijo decrease condition is satisfied.
      procedure(objective_function_t) :: objective !! Objective callback procedure.
      real(dp), intent(in) :: x(:) !! Input data or predictor values.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: gradient(:) !! Gradient.
      real(dp), intent(in) :: direction(:) !! Direction.
      real(dp), intent(out) :: candidate(:) !! Candidate.
      real(dp), intent(out) :: candidate_value !! Candidate value.
      logical, intent(out) :: accepted !! Flag controlling accepted.
      real(dp) :: alpha, slope
      integer :: iteration

      alpha = 1.0_dp
      slope = dot_product(gradient, direction)
      accepted = .false.
      do iteration = 1, 60
         candidate = x + alpha*direction
         candidate_value = objective(candidate)
         if (ieee_is_finite(candidate_value)) then
            if (candidate_value <= value + 1.0e-4_dp*alpha*slope) then
               accepted = .true.
               return
            end if
         end if
         alpha = 0.5_dp*alpha
         if (alpha < 1.0e-12_dp) return
      end do
   end subroutine armijo_search

end module optimization_mod
