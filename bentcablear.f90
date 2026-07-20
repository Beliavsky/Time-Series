! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Algorithms translated from the R bentcableAR package.
module bentcablear_mod
   !! Bent-cable regression with independent or autoregressive errors.
   use kind_mod, only: dp
   use stats_mod, only: ols_fit, normal_quantile
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd
   use polynomial_mod, only: polynomial_roots_t, polynomial_roots
   use linalg_mod, only: invert_matrix, outer_product
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   type, public :: bentcable_residuals_t
      !! Regression residuals and conditional AR innovations.
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: innovations(:)
      real(dp) :: rss = huge(1.0_dp)
      integer :: first_innovation = 0
      integer :: info = 0
   end type bentcable_residuals_t

   type, public :: bentcable_fit_t
      !! Bent-cable or broken-stick CSS fit with autoregressive errors.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: innovations(:)
      real(dp) :: tau = 0.0_dp
      real(dp) :: gamma = 0.0_dp
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: innovation_variance = huge(1.0_dp)
      integer :: ar_order = 0
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: stick = .false.
      logical :: stationary = .false.
      logical :: converged = .false.
   end type bentcable_fit_t

   type, public :: bentcable_profile_t
      !! Profile-deviance surface and its best fixed transition point.
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: deviance(:, :)
      real(dp), allocatable :: initial_cable(:)
      real(dp), allocatable :: initial_ar(:)
      integer :: selected_tau = 0
      integer :: selected_gamma = 0
      integer :: info = 0
   end type bentcable_profile_t

   type, public :: bentcable_change_t
      !! Critical transition point estimate, variance, and normal interval.
      real(dp) :: estimate = 0.0_dp
      real(dp) :: variance = huge(1.0_dp)
      real(dp) :: lower = 0.0_dp
      real(dp) :: upper = 0.0_dp
      real(dp) :: level = 0.95_dp
      integer :: info = 0
   end type bentcable_change_t

   type, public :: bentcable_ar_covariance_t
      !! Yule-Walker AR coefficients, correlations, and Toeplitz covariance.
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: correlation(:)
      real(dp), allocatable :: covariance(:, :)
      integer :: info = 0
   end type bentcable_ar_covariance_t

   public :: bentcable_basis, bentcable_value, bentcable_design_matrix
   public :: bentcable_residuals, bentcable_sse, bentcable_stationary
   public :: bentcable_profile, bentcable_fit
   public :: bentcable_fisher_information, bentcable_change_interval
   public :: bentcable_ar_covariance, bentcable_fit_iterative_yw

contains

   pure elemental real(dp) function bentcable_basis(offset, gamma) result(value)
      !! Evaluate the linear-quadratic-linear cable basis.
      real(dp), intent(in) :: offset !! Covariate minus transition center.
      real(dp), intent(in) :: gamma !! Half-width of the quadratic transition.

      if (gamma > 0.0_dp) then
         if (offset < -gamma) then
            value = 0.0_dp
         else if (offset <= gamma) then
            value = (offset + gamma)**2/(4.0_dp*gamma)
         else
            value = offset
         end if
      else
         value = merge(offset, 0.0_dp, offset > 0.0_dp)
      end if
   end function bentcable_basis

   pure elemental real(dp) function bentcable_value(time, beta0, beta1, beta2, &
      tau, gamma) result(value)
      !! Evaluate a bent cable or broken stick at one covariate value.
      real(dp), intent(in) :: time !! Covariate value.
      real(dp), intent(in) :: beta0 !! Regression intercept.
      real(dp), intent(in) :: beta1 !! Initial linear slope.
      real(dp), intent(in) :: beta2 !! Slope change after the transition.
      real(dp), intent(in) :: tau !! Transition center.
      real(dp), intent(in) :: gamma !! Transition half-width, zero for a stick.

      value = beta0 + beta1*time + beta2*bentcable_basis(time - tau, gamma)
   end function bentcable_value

   pure function bentcable_design_matrix(time, tau, gamma) result(design)
      !! Construct intercept, time, and cable-basis regression columns.
      real(dp), intent(in) :: time(:) !! Covariate values.
      real(dp), intent(in) :: tau !! Transition center.
      real(dp), intent(in) :: gamma !! Transition half-width.
      real(dp) :: design(size(time), 3)

      design(:, 1) = 1.0_dp
      design(:, 2) = time
      design(:, 3) = bentcable_basis(time - tau, gamma)
   end function bentcable_design_matrix

   pure function bentcable_residuals(data, time, beta, tau, gamma, ar) result(out)
      !! Compute fitted-curve residuals and conditional AR innovations.
      real(dp), intent(in) :: data(:) !! Response observations.
      real(dp), intent(in) :: time(:) !! Covariate values.
      real(dp), intent(in) :: beta(:) !! Intercept, initial slope, and slope change.
      real(dp), intent(in) :: tau !! Transition center.
      real(dp), intent(in) :: gamma !! Transition half-width.
      real(dp), intent(in) :: ar(:) !! Autoregressive error coefficients.
      type(bentcable_residuals_t) :: out
      integer :: observations, order, index, lag, time_index

      observations = size(data)
      order = size(ar)
      if (size(time) /= observations .or. size(beta) /= 3 .or. &
         observations <= order .or. gamma < 0.0_dp .or. &
         .not. all(ieee_is_finite(data)) .or. &
         .not. all(ieee_is_finite(time)) .or. &
         .not. all(ieee_is_finite(beta)) .or. &
         .not. all(ieee_is_finite(ar))) then
         out%info = 1
         return
      end if
      allocate(out%residuals(observations), out%innovations(observations - order))
      out%residuals = data - bentcable_value(time, beta(1), beta(2), beta(3), &
         tau, gamma)
      do index = 1, observations - order
         time_index = order + index
         out%innovations(index) = out%residuals(time_index)
         do lag = 1, order
            out%innovations(index) = out%innovations(index) - &
               ar(lag)*out%residuals(time_index - lag)
         end do
      end do
      out%rss = sum(out%innovations**2)
      out%first_innovation = order + 1
   end function bentcable_residuals

   pure real(dp) function bentcable_sse(parameters, data, time, ar_order, &
      stick) result(sse)
      !! Evaluate conditional SSE for a complete cable and AR parameter vector.
      real(dp), intent(in) :: parameters(:) !! Cable parameters followed by AR coefficients.
      real(dp), intent(in) :: data(:) !! Response observations.
      real(dp), intent(in) :: time(:) !! Covariate values.
      integer, intent(in) :: ar_order !! Autoregressive error order.
      logical, intent(in), optional :: stick !! Fit a broken stick with gamma fixed at zero.
      type(bentcable_residuals_t) :: residual
      real(dp) :: gamma
      integer :: cable_count
      logical :: use_stick

      use_stick = .false.
      if (present(stick)) use_stick = stick
      cable_count = merge(4, 5, use_stick)
      sse = huge(1.0_dp)
      if (ar_order < 0 .or. size(parameters) /= cable_count + ar_order) return
      gamma = 0.0_dp
      if (.not. use_stick) gamma = parameters(5)
      if (gamma < 0.0_dp) return
      residual = bentcable_residuals(data, time, parameters(:3), parameters(4), &
         gamma, parameters(cable_count + 1:))
      if (residual%info == 0 .and. ieee_is_finite(residual%rss)) sse = residual%rss
   end function bentcable_sse

   pure logical function bentcable_stationary(ar) result(stationary)
      !! Test whether every root of the AR polynomial lies outside the unit circle.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      type(polynomial_roots_t) :: roots
      real(dp), allocatable :: polynomial(:)

      if (size(ar) == 0) then
         stationary = .true.
         return
      end if
      allocate(polynomial(size(ar) + 1))
      polynomial(1) = 1.0_dp
      polynomial(2:) = -ar
      roots = polynomial_roots(polynomial)
      stationary = roots%info == 0
      if (stationary) stationary = all(abs(roots%roots) > 1.0_dp + 1.0e-8_dp)
   end function bentcable_stationary

   pure function bentcable_profile(data, tau_values, gamma_values, ar_order, &
      time, max_iterations, tolerance) result(out)
      !! Profile conditional deviance over a supplied tau-gamma grid.
      real(dp), intent(in) :: data(:) !! Response observations.
      real(dp), intent(in) :: tau_values(:) !! Candidate transition centers.
      real(dp), intent(in) :: gamma_values(:) !! Candidate nonnegative half-widths.
      integer, intent(in), optional :: ar_order !! Autoregressive error order.
      real(dp), intent(in), optional :: time(:) !! Covariate values.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations per grid point.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(bentcable_profile_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: selected_time(:), design(:, :), coefficient(:)
      real(dp), allocatable :: standard_error(:), residual(:), initial(:), candidate(:)
      real(dp) :: rss, value, best_value, selected_tolerance
      integer :: order, limit, observations, i, j, status, parameter_count

      observations = size(data)
      order = 0
      if (present(ar_order)) order = ar_order
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (present(time)) then
         selected_time = time
      else
         allocate(selected_time(observations))
         selected_time = [(real(i - 1, dp), i=1, observations)]
      end if
      if (observations < 5 .or. size(selected_time) /= observations .or. &
         size(tau_values) < 1 .or. size(gamma_values) < 1 .or. order < 0 .or. &
         any(gamma_values < 0.0_dp) .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. &
         .not. all(ieee_is_finite(data)) .or. &
         .not. all(ieee_is_finite(selected_time))) then
         out%info = 1
         return
      end if
      out%tau = tau_values
      out%gamma = gamma_values
      allocate(out%deviance(size(tau_values), size(gamma_values)))
      out%deviance = -huge(1.0_dp)
      best_value = -huge(1.0_dp)
      parameter_count = 3 + order
      do j = 1, size(gamma_values)
         do i = 1, size(tau_values)
            design = bentcable_design_matrix(selected_time, tau_values(i), &
               gamma_values(j))
            call ols_fit(design, data, coefficient, standard_error, residual, &
               rss, status)
            if (status /= 0) cycle
            allocate(initial(parameter_count))
            initial = 0.0_dp
            initial(:3) = coefficient
            if (order > 0) initial(4:) = alternating_ar_start(order)
            if (order > 0) then
               optimization = bfgs_minimize_fd(objective, initial, limit, &
                  selected_tolerance)
               if (.not. allocated(optimization%parameters)) then
                  deallocate(initial)
                  cycle
               end if
               candidate = optimization%parameters
               rss = optimization%objective
            else
               candidate = initial
            end if
            deallocate(initial)
            if (.not. ieee_is_finite(rss) .or. rss <= 0.0_dp) cycle
            value = -real(observations - order, dp)* &
               log(rss/real(observations - order, dp))
            out%deviance(i, j) = value
            if (value > best_value) then
               best_value = value
               out%selected_tau = i
               out%selected_gamma = j
               out%initial_cable = [candidate(:3), tau_values(i), gamma_values(j)]
               out%initial_ar = candidate(4:)
            end if
         end do
      end do
      if (out%selected_tau == 0) then
         out%info = 2
         return
      end if
      where (out%deviance > -0.5_dp*huge(1.0_dp))
         out%deviance = out%deviance - best_value
      end where

   contains

      pure real(dp) function objective(parameters) result(result_value)
         !! Evaluate fixed-transition CSS for profile optimization.
         real(dp), intent(in) :: parameters(:) !! Beta and AR coefficients.
         type(bentcable_residuals_t) :: evaluated

         evaluated = bentcable_residuals(data, selected_time, parameters(:3), &
            tau_values(i), gamma_values(j), parameters(4:))
         result_value = huge(1.0_dp)
         if (evaluated%info == 0) result_value = evaluated%rss
      end function objective

   end function bentcable_profile

   pure function bentcable_fit(data, ar_order, initial_cable, initial_ar, &
      time, stick, max_iterations, tolerance) result(out)
      !! Fit a bent cable or broken stick by conditional RSS minimization.
      real(dp), intent(in) :: data(:) !! Response observations.
      integer, intent(in), optional :: ar_order !! Autoregressive error order.
      real(dp), intent(in), optional :: initial_cable(:) !! Initial beta, tau, and optional gamma.
      real(dp), intent(in), optional :: initial_ar(:) !! Initial AR coefficients.
      real(dp), intent(in), optional :: time(:) !! Covariate values.
      logical, intent(in), optional :: stick !! Fit a broken stick rather than a bent cable.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(bentcable_fit_t) :: out
      type(optimization_result_t) :: optimization
      type(bentcable_profile_t) :: profile
      type(bentcable_residuals_t) :: residual
      real(dp), allocatable :: selected_time(:), starting(:), tau_grid(:), gamma_grid(:)
      real(dp) :: selected_tolerance, time_range
      integer :: observations, order, cable_count, limit, i
      logical :: use_stick

      observations = size(data)
      order = 0
      if (present(ar_order)) order = ar_order
      use_stick = .false.
      if (present(stick)) use_stick = stick
      cable_count = merge(4, 5, use_stick)
      limit = 1000
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (present(time)) then
         selected_time = time
      else
         allocate(selected_time(observations))
         selected_time = [(real(i - 1, dp), i=1, observations)]
      end if
      if (observations <= cable_count + order .or. order < 0 .or. &
         size(selected_time) /= observations .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. &
         .not. all(ieee_is_finite(data)) .or. &
         .not. all(ieee_is_finite(selected_time))) then
         out%info = 1
         return
      end if
      allocate(starting(cable_count + order))
      if (present(initial_cable)) then
         if (size(initial_cable) /= cable_count) then
            out%info = 2
            return
         end if
         starting(:cable_count) = initial_cable
      else
         time_range = maxval(selected_time) - minval(selected_time)
         allocate(tau_grid(15))
         tau_grid = [(minval(selected_time) + &
            real(i, dp)*time_range/16.0_dp, i=1, 15)]
         if (use_stick) then
            gamma_grid = [0.0_dp]
         else
            allocate(gamma_grid(8))
            gamma_grid = [(max(time_range/160.0_dp, sqrt(epsilon(1.0_dp))) + &
               real(i - 1, dp)*time_range/32.0_dp, i=1, 8)]
         end if
         profile = bentcable_profile(data, tau_grid, gamma_grid, order, &
            selected_time, min(limit, 250), selected_tolerance)
         if (profile%info /= 0) then
            out%info = 3
            return
         end if
         starting(:3) = profile%initial_cable(:3)
         starting(4) = profile%initial_cable(4)
         if (.not. use_stick) starting(5) = profile%initial_cable(5)
      end if
      if (order > 0) then
         if (present(initial_ar)) then
            if (size(initial_ar) /= order) then
               out%info = 2
               return
            end if
            starting(cable_count + 1:) = initial_ar
         else if (allocated(profile%initial_ar)) then
            starting(cable_count + 1:) = profile%initial_ar
         else
            starting(cable_count + 1:) = alternating_ar_start(order)
         end if
      end if
      optimization = bfgs_minimize_fd(objective, starting, limit, &
         selected_tolerance)
      if (.not. allocated(optimization%parameters) .or. &
         .not. ieee_is_finite(optimization%objective)) then
         out%info = 4
         return
      end if
      out%beta = optimization%parameters(:3)
      out%tau = optimization%parameters(4)
      out%gamma = 0.0_dp
      if (.not. use_stick) out%gamma = optimization%parameters(5)
      out%ar = optimization%parameters(cable_count + 1:)
      residual = bentcable_residuals(data, selected_time, out%beta, out%tau, &
         out%gamma, out%ar)
      if (residual%info /= 0) then
         out%info = 5
         return
      end if
      out%data = data
      out%time = selected_time
      out%fitted = data - residual%residuals
      out%residuals = residual%residuals
      out%innovations = residual%innovations
      out%rss = residual%rss
      out%innovation_variance = residual%rss/real(observations - order, dp)
      out%ar_order = order
      out%iterations = optimization%iterations
      out%optimizer_info = optimization%info
      out%converged = optimization%converged
      out%stationary = bentcable_stationary(out%ar)
      out%stick = use_stick

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Evaluate full bent-cable CSS for optimization.
         real(dp), intent(in) :: parameters(:) !! Cable and AR parameters.

         value = bentcable_sse(parameters, data, selected_time, order, use_stick)
      end function objective

   end function bentcable_fit

   pure function bentcable_fisher_information(fit) result(fisher)
      !! Compute the package's conditional Fisher information for cable parameters.
      type(bentcable_fit_t), intent(in) :: fit !! Fitted cable model.
      real(dp), allocatable :: fisher(:, :)
      real(dp), allocatable :: derivative(:), lagged_derivative(:)
      integer :: dimensions, observation, parameter, lag

      if (fit%info /= 0 .or. fit%innovation_variance <= 0.0_dp .or. &
         .not. allocated(fit%time)) then
         allocate(fisher(0, 0))
         return
      end if
      dimensions = merge(4, 5, fit%stick)
      allocate(fisher(dimensions, dimensions), derivative(dimensions), &
         lagged_derivative(dimensions))
      fisher = 0.0_dp
      do observation = 1, size(fit%time)
         derivative = cable_mean_gradient(fit%time(observation), fit%beta, &
            fit%tau, fit%gamma, fit%stick)
         derivative = -derivative
         do lag = 1, min(fit%ar_order, observation - 1)
            lagged_derivative = cable_mean_gradient(fit%time(observation - lag), &
               fit%beta, fit%tau, fit%gamma, fit%stick)
            derivative = derivative + fit%ar(lag)*lagged_derivative
         end do
         fisher = fisher + outer_product(derivative, derivative)
      end do
      fisher = fisher/fit%innovation_variance
   end function bentcable_fisher_information

   pure function bentcable_change_interval(fit, level) result(out)
      !! Estimate the critical transition point and its normal confidence interval.
      type(bentcable_fit_t), intent(in) :: fit !! Fitted cable model.
      real(dp), intent(in), optional :: level !! Confidence level between zero and one.
      type(bentcable_change_t) :: out
      real(dp), allocatable :: fisher(:, :), covariance(:, :), gradient(:)
      real(dp) :: selected_level, quantile_value
      integer :: status

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      if (fit%info /= 0 .or. selected_level <= 0.0_dp .or. &
         selected_level >= 1.0_dp .or. abs(fit%beta(3)) <= sqrt(epsilon(1.0_dp))) then
         out%info = 1
         return
      end if
      if (.not. fit%stick .and. fit%gamma <= sqrt(epsilon(1.0_dp))) then
         out%info = 1
         return
      end if
      fisher = bentcable_fisher_information(fit)
      if (size(fisher, 1) == 0) then
         out%info = 2
         return
      end if
      call invert_matrix(fisher, covariance, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      if (fit%stick) then
         allocate(gradient(4))
         gradient = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
         out%estimate = fit%tau
      else
         allocate(gradient(5))
         gradient = [0.0_dp, -2.0_dp*fit%gamma/fit%beta(3), &
            2.0_dp*fit%beta(2)*fit%gamma/fit%beta(3)**2, 1.0_dp, &
            -(2.0_dp*fit%beta(2) + fit%beta(3))/fit%beta(3)]
         out%estimate = fit%tau - fit%gamma - &
            2.0_dp*fit%beta(2)*fit%gamma/fit%beta(3)
      end if
      out%variance = dot_product(gradient, matmul(covariance, gradient))
      if (out%variance <= 0.0_dp .or. .not. ieee_is_finite(out%variance)) then
         out%info = 4
         return
      end if
      quantile_value = normal_quantile(0.5_dp*(1.0_dp + selected_level))
      out%lower = out%estimate - quantile_value*sqrt(out%variance)
      out%upper = out%estimate + quantile_value*sqrt(out%variance)
      out%level = selected_level
   end function bentcable_change_interval

   pure function bentcable_ar_covariance(residuals, order) result(out)
      !! Estimate Yule-Walker AR dependence and construct its correlation matrix.
      real(dp), intent(in) :: residuals(:) !! Cable-regression residuals.
      integer, intent(in) :: order !! Autoregressive order.
      type(bentcable_ar_covariance_t) :: out
      real(dp), allocatable :: centered(:), gram(:, :), rhs(:), inverse(:, :)
      real(dp) :: denominator
      integer :: observations, lag, row, column, status

      observations = size(residuals)
      if (order < 1 .or. observations <= 2*order + 1 .or. &
         .not. all(ieee_is_finite(residuals))) then
         out%info = 1
         return
      end if
      centered = residuals - sum(residuals)/real(observations, dp)
      denominator = sum(centered**2)
      if (denominator <= 0.0_dp) then
         out%info = 2
         return
      end if
      allocate(out%correlation(0:observations - 1))
      out%correlation(0) = 1.0_dp
      do lag = 1, order
         out%correlation(lag) = dot_product(centered(lag + 1:), &
            centered(:observations - lag))/denominator
      end do
      allocate(gram(order, order), rhs(order))
      do row = 1, order
         rhs(row) = out%correlation(row)
         do column = 1, order
            gram(row, column) = out%correlation(abs(row - column))
         end do
      end do
      call invert_matrix(gram, inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      out%ar = matmul(inverse, rhs)
      do lag = order + 1, observations - 1
         out%correlation(lag) = dot_product(out%ar, &
            out%correlation(lag - 1:lag - order:-1))
      end do
      allocate(out%covariance(observations, observations))
      do row = 1, observations
         do column = 1, observations
            out%covariance(row, column) = &
               out%correlation(abs(row - column))
         end do
      end do
   end function bentcable_ar_covariance

   pure function bentcable_fit_iterative_yw(data, ar_order, initial_cable, &
      initial_ar, time, stick, max_iterations, tolerance) result(out)
      !! Fit a cable by alternating Yule-Walker dependence and cable updates.
      real(dp), intent(in) :: data(:) !! Response observations.
      integer, intent(in) :: ar_order !! Autoregressive error order.
      real(dp), intent(in) :: initial_cable(:) !! Initial beta, tau, and optional gamma.
      real(dp), intent(in), optional :: initial_ar(:) !! Initial AR coefficients.
      real(dp), intent(in), optional :: time(:) !! Covariate values.
      logical, intent(in), optional :: stick !! Fit a broken stick rather than a bent cable.
      integer, intent(in), optional :: max_iterations !! Maximum alternating iterations.
      real(dp), intent(in), optional :: tolerance !! Parameter convergence tolerance.
      type(bentcable_fit_t) :: out
      type(bentcable_fit_t) :: current
      type(bentcable_ar_covariance_t) :: dependence
      type(bentcable_residuals_t) :: residual
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: selected_time(:), cable(:), previous(:), ar(:)
      real(dp) :: selected_tolerance, selected_gamma
      integer :: observations, cable_count, limit, iteration, i
      logical :: use_stick

      observations = size(data)
      use_stick = .false.
      if (present(stick)) use_stick = stick
      cable_count = merge(4, 5, use_stick)
      limit = 100
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (present(time)) then
         selected_time = time
      else
         allocate(selected_time(observations))
         selected_time = [(real(i - 1, dp), i=1, observations)]
      end if
      if (ar_order < 1 .or. size(initial_cable) /= cable_count .or. &
         size(selected_time) /= observations .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      cable = initial_cable
      if (present(initial_ar)) then
         if (size(initial_ar) /= ar_order) then
            out%info = 1
            return
         end if
         ar = initial_ar
      else
         ar = alternating_ar_start(ar_order)
      end if
      do iteration = 1, limit
         previous = [cable, ar]
         selected_gamma = 0.0_dp
         if (.not. use_stick) selected_gamma = cable(5)
         residual = bentcable_residuals(data, selected_time, cable(:3), &
            cable(4), selected_gamma, ar)
         if (residual%info /= 0) then
            out%info = 2
            return
         end if
         dependence = bentcable_ar_covariance(residual%residuals, ar_order)
         if (dependence%info /= 0) then
            out%info = 3
            return
         end if
         ar = dependence%ar
         optimization = bfgs_minimize_fd(objective, cable, 300, &
            0.1_dp*selected_tolerance)
         if (.not. allocated(optimization%parameters)) then
            out%info = 4
            return
         end if
         cable = optimization%parameters
         if (sum(abs([cable, ar] - previous)) <= selected_tolerance) exit
      end do
      current = bentcable_fit(data, ar_order, cable, ar, selected_time, &
         use_stick, 1, selected_tolerance)
      if (current%info /= 0) then
         selected_gamma = 0.0_dp
         if (.not. use_stick) selected_gamma = cable(5)
         residual = bentcable_residuals(data, selected_time, cable(:3), &
            cable(4), selected_gamma, ar)
         if (residual%info /= 0) then
            out%info = 5
            return
         end if
         out%beta = cable(:3)
         out%tau = cable(4)
         if (.not. use_stick) out%gamma = cable(5)
         out%ar = ar
         out%data = data
         out%time = selected_time
         out%fitted = data - residual%residuals
         out%residuals = residual%residuals
         out%innovations = residual%innovations
         out%rss = residual%rss
         out%innovation_variance = residual%rss/real(observations - ar_order, dp)
         out%ar_order = ar_order
         out%stick = use_stick
         out%stationary = bentcable_stationary(ar)
      else
         out = current
      end if
      out%iterations = iteration
      out%converged = iteration <= limit

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Evaluate cable CSS conditional on the current Yule-Walker AR fit.
         real(dp), intent(in) :: parameters(:) !! Candidate cable parameters.
         type(bentcable_residuals_t) :: evaluated
         real(dp) :: selected_gamma

         value = huge(1.0_dp)
         selected_gamma = 0.0_dp
         if (.not. use_stick) selected_gamma = parameters(5)
         if (selected_gamma < 0.0_dp) return
         evaluated = bentcable_residuals(data, selected_time, parameters(:3), &
            parameters(4), selected_gamma, ar)
         if (evaluated%info == 0) value = evaluated%rss
      end function objective

   end function bentcable_fit_iterative_yw

   pure function cable_mean_gradient(time, beta, tau, gamma, stick) &
      result(gradient)
      !! Differentiate the cable mean with respect to its regression parameters.
      real(dp), intent(in) :: time !! Covariate value.
      real(dp), intent(in) :: beta(:) !! Cable regression coefficients.
      real(dp), intent(in) :: tau !! Transition center.
      real(dp), intent(in) :: gamma !! Transition half-width.
      logical, intent(in) :: stick !! Broken-stick indicator.
      real(dp), allocatable :: gradient(:)
      real(dp) :: offset, derivative_tau, derivative_gamma

      offset = time - tau
      if (stick) then
         allocate(gradient(4))
         gradient = [1.0_dp, time, bentcable_basis(offset, 0.0_dp), &
            -beta(3)*merge(1.0_dp, 0.0_dp, offset > 0.0_dp)]
         return
      end if
      allocate(gradient(5))
      derivative_tau = 0.0_dp
      derivative_gamma = 0.0_dp
      if (offset > gamma) then
         derivative_tau = -1.0_dp
      else if (abs(offset) <= gamma) then
         derivative_tau = -(offset + gamma)/(2.0_dp*gamma)
         derivative_gamma = 0.25_dp*(1.0_dp - (offset/gamma)**2)
      end if
      gradient = [1.0_dp, time, bentcable_basis(offset, gamma), &
         beta(3)*derivative_tau, beta(3)*derivative_gamma]
   end function cable_mean_gradient

   pure function alternating_ar_start(order) result(ar)
      !! Construct the package's alternating 0.5 and -0.5 AR start.
      integer, intent(in) :: order !! Autoregressive order.
      real(dp) :: ar(order)
      integer :: lag

      do lag = 1, order
         ar(lag) = merge(0.5_dp, -0.5_dp, mod(lag, 2) == 1)
      end do
   end function alternating_ar_start

end module bentcablear_mod
