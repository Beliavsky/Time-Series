! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Algorithms translated from the R starvars package.
module starvars_mod
   !! Vector logistic smooth-transition autoregression and related utilities.
   use kind_mod, only: dp
   use calendar_mod, only: date_t, date_valid, date_day_number
   use linalg_mod, only: cholesky_lower, inverse_logdet, &
      symmetric_pseudoinverse
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd
   use random_mod, only: random_multivariate_normal, random_uniform, &
      set_random_seed
   use special_functions_mod, only: regularized_gamma_q
   use stats_mod, only: normal_quantile, quantile, sort
   use time_series_stats_mod, only: acf_values
   use utils_mod, only: quiet_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: starvars_method_nls = 1
   integer, parameter, public :: starvars_method_ml = 2
   integer, parameter, public :: starvars_forecast_naive = 1
   integer, parameter, public :: starvars_forecast_monte_carlo = 2
   integer, parameter, public :: starvars_forecast_bootstrap = 3
   integer, parameter, public :: starvars_frequency_daily = 1
   integer, parameter, public :: starvars_frequency_monthly = 2
   integer, parameter, public :: starvars_frequency_quarterly = 3
   integer, parameter, public :: starvars_frequency_yearly = 4

   type, public :: starvars_start_t
      !! Starting logistic slopes and locations for every transition and equation.
      real(dp), allocatable :: gamma(:, :)
      real(dp), allocatable :: location(:, :)
      real(dp) :: objective = huge(1.0_dp)
      integer :: info = 0
   end type starvars_start_t

   type, public :: starvars_fit_t
      !! Profile-estimated vector logistic smooth-transition autoregression.
      real(dp), allocatable :: coefficients(:, :, :)
      real(dp), allocatable :: coefficient_standard_errors(:, :, :)
      real(dp), allocatable :: gamma(:, :)
      real(dp), allocatable :: location(:, :)
      real(dp), allocatable :: transition_weights(:, :, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: series(:, :)
      real(dp), allocatable :: exogenous(:, :)
      real(dp), allocatable :: transition(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: rss = huge(1.0_dp)
      real(dp), allocatable :: aic(:)
      real(dp), allocatable :: bic(:)
      integer :: order = 0
      integer :: regimes = 0
      integer :: regressor_count = 0
      integer :: method = starvars_method_nls
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: constant = .true.
      logical :: common_transition = .false.
      logical :: converged = .false.
   end type starvars_fit_t

   type, public :: starvars_forecast_t
      !! VLSTAR point forecasts, intervals, and optional simulated paths.
      real(dp), allocatable :: point(:, :)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
      real(dp), allocatable :: paths(:, :, :)
      real(dp) :: confidence = 0.95_dp
      integer :: method = starvars_forecast_naive
      integer :: simulations = 0
      integer :: info = 0
   end type starvars_forecast_t

   type, public :: starvars_joint_test_t
      !! Third-order Taylor joint linearity tests for transition candidates.
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: p_value(:)
      real(dp) :: critical_value = 0.0_dp
      integer :: degrees_of_freedom = 0
      integer :: selected = 0
      integer :: info = 0
   end type starvars_joint_test_t

   type, public :: starvars_long_run_variance_t
      !! Bartlett-kernel long-run variance and selected bandwidth.
      real(dp) :: variance = 0.0_dp
      integer :: bandwidth = 0
      integer :: info = 0
   end type starvars_long_run_variance_t

   type, public :: starvars_cumsum_t
      !! Binary-segmentation covariance CUMSUM statistics and break indices.
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: omega(:)
      integer, allocatable :: break_location(:)
      real(dp) :: lambda_critical = 0.0_dp
      real(dp) :: omega_critical = 0.0_dp
      integer :: info = 0
   end type starvars_cumsum_t

   type, public :: starvars_realized_covariance_t
      !! Grouped realized covariance vectors and optional Cholesky factors.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: cholesky(:, :)
      real(dp), allocatable :: returns(:, :)
      integer, allocatable :: groups(:)
      integer :: info = 0
   end type starvars_realized_covariance_t

   type :: starvars_evaluation_t
      !! Internal profile objective and fitted coefficient state.
      real(dp), allocatable :: coefficients(:, :, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: weights(:, :, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: rss = huge(1.0_dp)
      integer :: info = 0
   end type starvars_evaluation_t

   interface display
      module procedure display_starvars_start
      module procedure display_starvars_fit
      module procedure display_starvars_forecast
      module procedure display_starvars_joint_test
      module procedure display_starvars_long_run_variance
      module procedure display_starvars_cumsum
      module procedure display_starvars_realized_covariance
   end interface display

   interface starvars_realized_covariance
      module procedure starvars_realized_covariance_groups
      module procedure starvars_realized_covariance_dates
   end interface starvars_realized_covariance

   public :: starvars_logistic, starvars_starting, starvars_fit
   public :: starvars_forecast, starvars_joint_linearity_test
   public :: starvars_long_run_variance, starvars_multiple_cumsum
   public :: starvars_realized_covariance, display

contains

   pure elemental real(dp) function starvars_logistic(transition, gamma, &
      location) result(weight)
      !! Evaluate a numerically stable logistic transition function.
      real(dp), intent(in) :: transition !! Transition-variable value.
      real(dp), intent(in) :: gamma !! Positive transition slope.
      real(dp), intent(in) :: location !! Transition location.
      real(dp) :: argument

      argument = gamma*(transition - location)
      if (argument >= 0.0_dp) then
         weight = 1.0_dp/(1.0_dp + exp(-min(argument, 700.0_dp)))
      else
         weight = exp(max(argument, -700.0_dp))/ &
            (1.0_dp + exp(max(argument, -700.0_dp)))
      end if
   end function starvars_logistic

   pure function starvars_starting(series, transition, order, regimes, &
      exogenous, constant, combinations, common_transition) result(start)
      !! Select logistic starts using the package's startingVLSTAR grid search.
      real(dp), intent(in) :: series(:, :) !! Observation-by-variable response matrix.
      real(dp), intent(in) :: transition(:) !! Common transition variable.
      integer, intent(in), optional :: order !! Positive VAR lag order.
      integer, intent(in), optional :: regimes !! Number of regimes; defaults to two.
      real(dp), intent(in), optional :: exogenous(:, :) !! Contemporaneous exogenous regressors.
      logical, intent(in), optional :: constant !! Include an intercept.
      integer, intent(in), optional :: combinations !! Grid points per dimension.
      logical, intent(in), optional :: common_transition !! Share gamma and location across equations.
      type(starvars_start_t) :: start
      real(dp), allocatable :: response(:, :), design(:, :), aligned_transition(:)
      real(dp), allocatable :: gamma(:, :), location(:, :), parameters(:)
      type(starvars_evaluation_t) :: evaluated
      real(dp) :: lower, upper, candidate_gamma, candidate_location, score, best
      integer :: selected_order, selected_regimes, grid, variables, transitions
      integer :: width, transition_index, equation, gamma_index, location_index
      logical :: include_constant, common, found

      selected_order = 1
      if (present(order)) selected_order = order
      selected_regimes = 2
      if (present(regimes)) selected_regimes = regimes
      grid = 7
      if (present(combinations)) grid = combinations
      include_constant = .true.
      if (present(constant)) include_constant = constant
      common = .false.
      if (present(common_transition)) common = common_transition
      call prepare_vlstar_data(series, transition, selected_order, exogenous, &
         include_constant, response, design, aligned_transition)
      variables = size(series, 2)
      transitions = selected_regimes - 1
      width = merge(1, variables, common)
      if (size(response, 1) < 2 .or. selected_regimes < 2 .or. grid < 2) then
         start%info = 1
         return
      end if
      allocate(gamma(width, transitions), source=1.0_dp)
      allocate(location(width, transitions), &
         source=sum(aligned_transition)/real(size(aligned_transition), dp))
      lower = minval(aligned_transition)
      upper = maxval(aligned_transition)
      do transition_index = 1, transitions
         do equation = 1, width
            best = huge(1.0_dp)
            found = .false.
            do gamma_index = 1, grid
               candidate_gamma = 100.0_dp*real(gamma_index - 1, dp)/ &
                  real(grid - 1, dp)
               candidate_gamma = max(candidate_gamma, 1.0e-4_dp)
               do location_index = 1, grid
                  candidate_location = lower + (upper - lower)* &
                     real(location_index - 1, dp)/real(grid - 1, dp)
                  gamma(equation, transition_index) = candidate_gamma
                  location(equation, transition_index) = candidate_location
                  parameters = pack_transition_parameters(gamma, location)
                  evaluated = evaluate_vlstar(response, design, &
                     aligned_transition, selected_regimes, common, &
                     starvars_method_nls, parameters)
                  if (evaluated%info /= 0) cycle
                  score = evaluated%rss
                  if (.not. common .and. evaluated%info == 0) &
                     score = sum(evaluated%residuals(:, equation)**2)
                  if (score < best) then
                     best = score
                     start%gamma = gamma
                     start%location = location
                     start%objective = evaluated%rss
                     found = .true.
                  end if
               end do
            end do
            if (.not. found) then
               start%info = 2
               return
            end if
            gamma = start%gamma
            location = start%location
         end do
      end do
   end function starvars_starting

   pure function starvars_fit(series, transition, order, regimes, exogenous, &
      constant, starting, method, common_transition, max_iterations, &
      tolerance) result(out)
      !! Estimate the package's VLSTAR model by profile NLS or Gaussian ML.
      real(dp), intent(in) :: series(:, :) !! Observation-by-variable response matrix.
      real(dp), intent(in) :: transition(:) !! Common transition variable.
      integer, intent(in), optional :: order !! Positive VAR lag order.
      integer, intent(in), optional :: regimes !! Number of regimes; defaults to two.
      real(dp), intent(in), optional :: exogenous(:, :) !! Contemporaneous exogenous regressors.
      logical, intent(in), optional :: constant !! Include an intercept.
      type(starvars_start_t), intent(in), optional :: starting !! Logistic starting values.
      integer, intent(in), optional :: method !! NLS or concentrated-ML method code.
      logical, intent(in), optional :: common_transition !! Share transition parameters.
      integer, intent(in), optional :: max_iterations !! Maximum finite-difference BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Optimizer gradient tolerance.
      type(starvars_fit_t) :: out
      type(starvars_start_t) :: selected_start
      type(optimization_result_t) :: optimum
      type(starvars_evaluation_t) :: evaluated
      real(dp), allocatable :: response(:, :), design(:, :), aligned_transition(:)
      real(dp), allocatable :: initial(:), gamma(:, :), location(:, :)
      integer :: selected_order, selected_regimes, selected_method, limit
      integer :: variables, width, equation, status
      real(dp) :: selected_tolerance
      logical :: include_constant, common

      selected_order = 1
      if (present(order)) selected_order = order
      selected_regimes = 2
      if (present(regimes)) selected_regimes = regimes
      selected_method = starvars_method_nls
      if (present(method)) selected_method = method
      include_constant = .true.
      if (present(constant)) include_constant = constant
      common = .false.
      if (present(common_transition)) common = common_transition
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      variables = size(series, 2)
      width = merge(1, variables, common)
      if (size(series, 1) <= selected_order + 2 .or. variables < 2 .or. &
         size(transition) /= size(series, 1) .or. selected_order < 1 .or. &
         selected_regimes < 2 .or. &
         (selected_method /= starvars_method_nls .and. &
         selected_method /= starvars_method_ml) .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(transition))) then
         out%info = 1
         return
      end if
      call prepare_vlstar_data(series, transition, selected_order, exogenous, &
         include_constant, response, design, aligned_transition)
      if (size(response, 1) < size(design, 2)*selected_regimes + 1) then
         out%info = 2
         return
      end if
      if (present(starting)) then
         selected_start = starting
      else
         selected_start = starvars_starting(series, transition, selected_order, &
            selected_regimes, exogenous, include_constant, 7, common)
      end if
      if (selected_start%info /= 0 .or. &
         .not. allocated(selected_start%gamma) .or. &
         .not. allocated(selected_start%location)) then
         out%info = 3
         return
      end if
      if (size(selected_start%gamma, 1) /= width .or. &
         size(selected_start%gamma, 2) /= selected_regimes - 1 .or. &
         size(selected_start%location, 1) /= width .or. &
         size(selected_start%location, 2) /= selected_regimes - 1 .or. &
         any(selected_start%gamma <= 0.0_dp)) then
         out%info = 3
         return
      end if
      initial = pack_transition_parameters(selected_start%gamma, &
         selected_start%location)
      optimum = bfgs_minimize_fd(objective, initial, limit, selected_tolerance)
      if (.not. allocated(optimum%parameters)) then
         out%info = 4
         return
      end if
      evaluated = evaluate_vlstar(response, design, aligned_transition, &
         selected_regimes, common, selected_method, optimum%parameters)
      if (evaluated%info /= 0) then
         optimum%parameters = initial
         optimum%converged = .false.
         evaluated = evaluate_vlstar(response, design, aligned_transition, &
            selected_regimes, common, selected_method, optimum%parameters)
         if (evaluated%info /= 0) then
            out%info = 5
            return
         end if
      end if
      call unpack_transition_parameters(optimum%parameters, width, &
         selected_regimes - 1, gamma, location)
      out%coefficients = evaluated%coefficients
      out%gamma = gamma
      out%location = location
      out%transition_weights = evaluated%weights
      out%covariance = evaluated%covariance
      out%series = series
      out%transition = transition
      if (present(exogenous)) then
         out%exogenous = exogenous
      else
         allocate(out%exogenous(size(series, 1), 0))
      end if
      allocate(out%fitted(size(series, 1), variables), source=quiet_nan())
      allocate(out%residuals(size(series, 1), variables), source=quiet_nan())
      out%fitted(selected_order + 1:, :) = evaluated%fitted
      out%residuals(selected_order + 1:, :) = evaluated%residuals
      out%rss = evaluated%rss
      out%log_likelihood = evaluated%log_likelihood
      out%order = selected_order
      out%regimes = selected_regimes
      out%regressor_count = size(design, 2)
      out%method = selected_method
      out%constant = include_constant
      out%common_transition = common
      out%iterations = optimum%iterations
      out%optimizer_info = optimum%info
      out%converged = optimum%converged
      allocate(out%aic(variables), out%bic(variables))
      do equation = 1, variables
         out%aic(equation) = real(2*size(design, 2)*selected_regimes, dp) + &
            real(size(response, 1), dp)* &
            (log(max(evaluated%covariance(equation, equation), &
            tiny(1.0_dp))) + 1.0_dp)
         out%bic(equation) = real(size(design, 2)*selected_regimes, dp)* &
            log(real(size(response, 1), dp)) + &
            real(size(response, 1), dp)* &
            (log(max(evaluated%covariance(equation, equation), &
            tiny(1.0_dp))) + 1.0_dp)
      end do
      call coefficient_standard_errors(design, aligned_transition, out, &
         status)
      if (status /= 0) then
         if (allocated(out%coefficient_standard_errors)) &
            deallocate(out%coefficient_standard_errors)
         allocate(out%coefficient_standard_errors(0, 0, 0))
      end if

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Evaluate the profile VLSTAR objective for the shared optimizer.
         real(dp), intent(in) :: parameters(:) !! Packed log-slope and location values.
         type(starvars_evaluation_t) :: candidate

         candidate = evaluate_vlstar(response, design, aligned_transition, &
            selected_regimes, common, selected_method, parameters)
         value = candidate%objective
      end function objective

   end function starvars_fit

   pure subroutine prepare_vlstar_data(series, transition, order, exogenous, &
      constant, response, design, aligned_transition)
      !! Align responses, VAR lags, exogenous inputs, and transition values.
      real(dp), intent(in) :: series(:, :) !! Observation-by-variable response matrix.
      real(dp), intent(in) :: transition(:) !! Transition variable.
      integer, intent(in) :: order !! Positive VAR lag order.
      real(dp), intent(in), optional :: exogenous(:, :) !! Contemporaneous exogenous regressors.
      logical, intent(in) :: constant !! Include an intercept.
      real(dp), allocatable, intent(out) :: response(:, :) !! Aligned responses.
      real(dp), allocatable, intent(out) :: design(:, :) !! Aligned base-regime design.
      real(dp), allocatable, intent(out) :: aligned_transition(:) !! Aligned transition values.
      integer :: observations, variables, exogenous_count, columns
      integer :: row, lag, first_column, time

      observations = size(series, 1)
      variables = size(series, 2)
      exogenous_count = 0
      if (present(exogenous)) exogenous_count = size(exogenous, 2)
      if (order < 1 .or. observations <= order .or. &
         size(transition) /= observations) then
         allocate(response(0, variables), design(0, 0), aligned_transition(0))
         return
      end if
      if (present(exogenous)) then
         if (size(exogenous, 1) /= observations) then
            allocate(response(0, variables), design(0, 0), &
               aligned_transition(0))
            return
         end if
      end if
      columns = merge(1, 0, constant) + order*variables + exogenous_count
      allocate(response(observations - order, variables))
      allocate(design(observations - order, columns))
      aligned_transition = transition(order + 1:)
      response = series(order + 1:, :)
      design = 0.0_dp
      first_column = 0
      if (constant) then
         first_column = 1
         design(:, 1) = 1.0_dp
      end if
      do row = 1, size(response, 1)
         time = order + row
         do lag = 1, order
            design(row, first_column + (lag - 1)*variables + 1: &
               first_column + lag*variables) = series(time - lag, :)
         end do
         if (exogenous_count > 0) &
            design(row, columns - exogenous_count + 1:) = exogenous(time, :)
      end do
   end subroutine prepare_vlstar_data

   pure function pack_transition_parameters(gamma, location) result(parameters)
      !! Pack positive slopes as logarithms followed by transition locations.
      real(dp), intent(in) :: gamma(:, :) !! Positive slopes by equation and transition.
      real(dp), intent(in) :: location(:, :) !! Locations by equation and transition.
      real(dp), allocatable :: parameters(:)
      integer :: count

      count = size(gamma)
      allocate(parameters(2*count))
      parameters(:count) = reshape(log(max(gamma, 1.0e-8_dp)), [count])
      parameters(count + 1:) = reshape(location, [count])
   end function pack_transition_parameters

   pure subroutine unpack_transition_parameters(parameters, width, transitions, &
      gamma, location)
      !! Unpack logged slopes and locations into transition matrices.
      real(dp), intent(in) :: parameters(:) !! Packed log-slopes and locations.
      integer, intent(in) :: width !! One or the number of response equations.
      integer, intent(in) :: transitions !! Number of logistic transitions.
      real(dp), allocatable, intent(out) :: gamma(:, :) !! Positive slopes.
      real(dp), allocatable, intent(out) :: location(:, :) !! Locations.
      integer :: count

      count = width*transitions
      allocate(gamma(width, transitions), location(width, transitions))
      gamma = reshape(exp(min(parameters(:count), log(1.0e6_dp))), &
         [width, transitions])
      location = reshape(parameters(count + 1:2*count), [width, transitions])
   end subroutine unpack_transition_parameters

   pure function evaluate_vlstar(response, design, transition, regimes, &
      common_transition, method, parameters) result(out)
      !! Profile regime coefficients and evaluate NLS or Gaussian ML objective.
      real(dp), intent(in) :: response(:, :) !! Aligned response matrix.
      real(dp), intent(in) :: design(:, :) !! Base-regime design matrix.
      real(dp), intent(in) :: transition(:) !! Aligned transition variable.
      integer, intent(in) :: regimes !! Number of regimes.
      logical, intent(in) :: common_transition !! Share transition parameters.
      integer, intent(in) :: method !! NLS or concentrated-ML method code.
      real(dp), intent(in) :: parameters(:) !! Packed transition parameters.
      type(starvars_evaluation_t) :: out
      real(dp), allocatable :: gamma(:, :), location(:, :)
      real(dp), allocatable :: equation_design(:, :), cross_product(:, :)
      real(dp), allocatable :: inverse(:, :), coefficient(:), covariance_inverse(:, :)
      real(dp) :: log_determinant
      integer :: observations, variables, regressors, transitions, width
      integer :: equation, regime, status

      observations = size(response, 1)
      variables = size(response, 2)
      regressors = size(design, 2)
      transitions = regimes - 1
      width = merge(1, variables, common_transition)
      if (observations /= size(design, 1) .or. &
         observations /= size(transition) .or. regimes < 2 .or. &
         size(parameters) /= 2*width*transitions) then
         out%info = 1
         return
      end if
      call unpack_transition_parameters(parameters, width, transitions, &
         gamma, location)
      allocate(out%weights(observations, variables, transitions))
      do regime = 1, transitions
         do equation = 1, variables
            out%weights(:, equation, regime) = starvars_logistic(transition, &
               gamma(merge(1, equation, common_transition), regime), &
               location(merge(1, equation, common_transition), regime))
         end do
      end do
      allocate(out%coefficients(regressors, variables, regimes))
      allocate(out%fitted(observations, variables))
      allocate(out%residuals(observations, variables))
      allocate(equation_design(observations, regressors*regimes))
      do equation = 1, variables
         equation_design(:, :regressors) = design
         do regime = 2, regimes
            equation_design(:, (regime - 1)*regressors + 1: &
               regime*regressors) = design* &
               spread(out%weights(:, equation, regime - 1), 2, regressors)
         end do
         cross_product = matmul(transpose(equation_design), equation_design)
         if (.not. allocated(inverse)) &
            allocate(inverse(size(cross_product, 1), size(cross_product, 2)))
         call symmetric_pseudoinverse(cross_product, inverse, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         coefficient = matmul(inverse, &
            matmul(transpose(equation_design), response(:, equation)))
         out%coefficients(:, equation, :) = reshape(coefficient, &
            [regressors, regimes])
         out%fitted(:, equation) = matmul(equation_design, coefficient)
      end do
      out%residuals = response - out%fitted
      out%rss = sum(out%residuals**2)
      out%covariance = matmul(transpose(out%residuals), out%residuals)/ &
         real(observations, dp)
      allocate(covariance_inverse(variables, variables))
      call inverse_logdet(out%covariance, covariance_inverse, log_determinant, &
         status, 100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         out%info = 3
         return
      end if
      out%log_likelihood = -0.5_dp*real(observations, dp)* &
         (real(variables, dp)*(1.0_dp + log(2.0_dp*acos(-1.0_dp))) + &
         log_determinant)
      if (method == starvars_method_ml) then
         out%objective = -out%log_likelihood
      else
         out%objective = out%rss
      end if
   end function evaluate_vlstar

   pure subroutine coefficient_standard_errors(design, transition, model, info)
      !! Approximate conditional coefficient standard errors equation by equation.
      real(dp), intent(in) :: design(:, :) !! Base-regime design matrix.
      real(dp), intent(in) :: transition(:) !! Aligned transition variable.
      type(starvars_fit_t), intent(inout) :: model !! Model receiving standard errors.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: equation_design(:, :), cross_product(:, :)
      real(dp), allocatable :: inverse(:, :), diagonal(:)
      integer :: observations, variables, regressors, equation, regime, status

      observations = size(design, 1)
      variables = size(model%coefficients, 2)
      regressors = size(design, 2)
      allocate(model%coefficient_standard_errors(regressors, variables, &
         model%regimes))
      allocate(equation_design(observations, regressors*model%regimes))
      info = 0
      do equation = 1, variables
         equation_design(:, :regressors) = design
         do regime = 2, model%regimes
            equation_design(:, (regime - 1)*regressors + 1: &
               regime*regressors) = design*spread(starvars_logistic(transition, &
               model%gamma(merge(1, equation, model%common_transition), &
               regime - 1), model%location(merge(1, equation, &
               model%common_transition), regime - 1)), 2, regressors)
         end do
         cross_product = matmul(transpose(equation_design), equation_design)
         if (.not. allocated(inverse)) &
            allocate(inverse(size(cross_product, 1), size(cross_product, 2)))
         call symmetric_pseudoinverse(cross_product, inverse, status)
         if (status /= 0) then
            info = equation
            return
         end if
         allocate(diagonal(size(inverse, 1)))
         do regime = 1, size(inverse, 1)
            diagonal(regime) = sqrt(max(inverse(regime, regime)* &
               model%covariance(equation, equation), 0.0_dp))
         end do
         model%coefficient_standard_errors(:, equation, :) = reshape(diagonal, &
            [regressors, model%regimes])
         deallocate(diagonal)
      end do
   end subroutine coefficient_standard_errors

   function starvars_forecast(model, horizon, method, simulations, &
      confidence, exogenous_new, transition_new, transition_variable, seed, &
      keep_paths) result(out)
      !! Translate predict.VLSTAR recursive and simulation forecasting.
      type(starvars_fit_t), intent(in) :: model !! Fitted VLSTAR model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      integer, intent(in), optional :: method !! Forecast-method code.
      integer, intent(in), optional :: simulations !! Number of simulated paths.
      real(dp), intent(in), optional :: confidence !! Central interval probability.
      real(dp), intent(in), optional :: exogenous_new(:, :) !! Future exogenous regressors.
      real(dp), intent(in), optional :: transition_new(:) !! Future transition values.
      integer, intent(in), optional :: transition_variable !! Endogenous transition series column.
      integer, intent(in), optional :: seed !! Shared random-number seed.
      logical, intent(in), optional :: keep_paths !! Retain simulated paths.
      type(starvars_forecast_t) :: out
      real(dp), allocatable :: history(:, :), base(:), innovation(:), zero(:)
      real(dp), allocatable :: ordered(:)
      real(dp) :: probability, z, transition_value
      integer :: selected_method, draws, variables, exogenous_count
      integer :: simulation, step, equation, status, sampled_row
      logical :: retain

      selected_method = starvars_forecast_naive
      if (present(method)) selected_method = method
      draws = 1000
      if (present(simulations)) draws = simulations
      probability = 0.95_dp
      if (present(confidence)) probability = confidence
      retain = .false.
      if (present(keep_paths)) retain = keep_paths
      variables = size(model%series, 2)
      exogenous_count = size(model%exogenous, 2)
      out%method = selected_method
      out%confidence = probability
      if (model%info /= 0 .or. horizon < 1 .or. variables < 1 .or. &
         probability <= 0.0_dp .or. probability >= 1.0_dp .or. &
         selected_method < starvars_forecast_naive .or. &
         selected_method > starvars_forecast_bootstrap) then
         out%info = 1
         return
      end if
      if (exogenous_count > 0 .and. .not. present(exogenous_new)) then
         out%info = 1
         return
      end if
      if (present(exogenous_new)) then
         if (size(exogenous_new, 1) < horizon .or. &
            size(exogenous_new, 2) /= exogenous_count) then
            out%info = 1
            return
         end if
      end if
      if (present(transition_new)) then
         if (size(transition_new) < horizon) then
            out%info = 1
            return
         end if
      end if
      if (present(transition_variable)) then
         if (transition_variable < 1 .or. transition_variable > variables) then
            out%info = 1
            return
         end if
      end if
      if (selected_method /= starvars_forecast_naive .and. draws < 2) then
         out%info = 2
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      allocate(out%point(horizon, variables), out%lower(horizon, variables), &
         out%upper(horizon, variables))
      allocate(history(size(model%series, 1) + horizon, variables))
      history(:size(model%series, 1), :) = model%series
      allocate(base(model%regressor_count), zero(variables), innovation(variables))
      zero = 0.0_dp
      do step = 1, horizon
         call forecast_design(model, history, step, exogenous_new, base)
         transition_value = forecast_transition(model, history, step, &
            transition_new, transition_variable)
         out%point(step, :) = forecast_mean(model, base, transition_value)
         history(size(model%series, 1) + step, :) = out%point(step, :)
      end do
      if (selected_method == starvars_forecast_naive) then
         z = normal_quantile(0.5_dp + 0.5_dp*probability)
         do equation = 1, variables
            out%lower(:, equation) = out%point(:, equation) - z* &
               sqrt(max(model%covariance(equation, equation), 0.0_dp))
            out%upper(:, equation) = out%point(:, equation) + z* &
               sqrt(max(model%covariance(equation, equation), 0.0_dp))
         end do
         out%info = 0
         return
      end if
      allocate(out%paths(horizon, variables, draws))
      do simulation = 1, draws
         history(:size(model%series, 1), :) = model%series
         do step = 1, horizon
            call forecast_design(model, history, step, exogenous_new, base)
            transition_value = forecast_transition(model, history, step, &
               transition_new, transition_variable)
            if (selected_method == starvars_forecast_monte_carlo) then
               call random_multivariate_normal(zero, model%covariance, &
                  innovation, status)
               if (status /= 0) then
                  out%info = 3
                  return
               end if
            else
               sampled_row = 1 + int(random_uniform()* &
                  real(size(model%series, 1) - model%order, dp))
               sampled_row = min(sampled_row, size(model%series, 1) - model%order)
               innovation = model%residuals(model%order + sampled_row, :)
            end if
            out%paths(step, :, simulation) = &
               forecast_mean(model, base, transition_value) + innovation
            history(size(model%series, 1) + step, :) = &
               out%paths(step, :, simulation)
         end do
      end do
      out%point = sum(out%paths, dim=3)/real(draws, dp)
      allocate(ordered(draws))
      do equation = 1, variables
         do step = 1, horizon
            ordered = out%paths(step, equation, :)
            call sort(ordered)
            out%lower(step, equation) = quantile(ordered, &
               0.5_dp*(1.0_dp - probability))
            out%upper(step, equation) = quantile(ordered, &
               0.5_dp*(1.0_dp + probability))
         end do
      end do
      out%simulations = draws
      if (.not. retain) deallocate(out%paths)
   end function starvars_forecast

   pure subroutine forecast_design(model, history, step, exogenous_new, base)
      !! Construct one future base-regime design vector.
      type(starvars_fit_t), intent(in) :: model !! Fitted VLSTAR model.
      real(dp), intent(in) :: history(:, :) !! Observed and forecast history.
      integer, intent(in) :: step !! Forecast step.
      real(dp), intent(in), optional :: exogenous_new(:, :) !! Future exogenous regressors.
      real(dp), intent(out) :: base(:) !! Constructed design vector.
      integer :: offset, lag, variables, time, exogenous_count

      variables = size(model%series, 2)
      exogenous_count = size(model%exogenous, 2)
      time = size(model%series, 1) + step
      base = 0.0_dp
      offset = 0
      if (model%constant) then
         base(1) = 1.0_dp
         offset = 1
      end if
      do lag = 1, model%order
         base(offset + (lag - 1)*variables + 1:offset + lag*variables) = &
            history(time - lag, :)
      end do
      if (exogenous_count > 0 .and. present(exogenous_new)) &
         base(size(base) - exogenous_count + 1:) = exogenous_new(step, :)
   end subroutine forecast_design

   pure real(dp) function forecast_transition(model, history, step, &
      transition_new, transition_variable) result(value)
      !! Select an supplied or endogenous future transition value.
      type(starvars_fit_t), intent(in) :: model !! Fitted VLSTAR model.
      real(dp), intent(in) :: history(:, :) !! Observed and forecast history.
      integer, intent(in) :: step !! Forecast step.
      real(dp), intent(in), optional :: transition_new(:) !! Future transition values.
      integer, intent(in), optional :: transition_variable !! Endogenous transition series column.
      integer :: column, time

      if (present(transition_new)) then
         value = transition_new(step)
      else if (present(transition_variable)) then
         column = transition_variable
         time = size(model%series, 1) + step
         value = history(time - 1, column)
      else
         value = model%transition(size(model%transition))
      end if
   end function forecast_transition

   pure function forecast_mean(model, base, transition) result(value)
      !! Evaluate the conditional mean for one future transition value.
      type(starvars_fit_t), intent(in) :: model !! Fitted VLSTAR model.
      real(dp), intent(in) :: base(:) !! Base-regime design vector.
      real(dp), intent(in) :: transition !! Transition value.
      real(dp) :: value(size(model%coefficients, 2))
      integer :: equation, regime, transition_row

      do equation = 1, size(value)
         value(equation) = dot_product(base, model%coefficients(:, equation, 1))
         transition_row = merge(1, equation, model%common_transition)
         do regime = 2, model%regimes
            value(equation) = value(equation) + starvars_logistic(transition, &
               model%gamma(transition_row, regime - 1), &
               model%location(transition_row, regime - 1))* &
               dot_product(base, model%coefficients(:, equation, regime))
         end do
      end do
   end function forecast_mean

   pure function starvars_joint_linearity_test(series, transitions, order, &
      exogenous, constant, significance) result(out)
      !! Translate VLSTARjoint's third-order logistic Taylor linearity test.
      real(dp), intent(in) :: series(:, :) !! Observation-by-variable response matrix.
      real(dp), intent(in) :: transitions(:, :) !! Candidate transition variables.
      integer, intent(in), optional :: order !! Positive VAR lag order.
      real(dp), intent(in), optional :: exogenous(:, :) !! Contemporaneous exogenous regressors.
      logical, intent(in), optional :: constant !! Include an intercept.
      real(dp), intent(in), optional :: significance !! Significance level; defaults to 0.05.
      type(starvars_joint_test_t) :: out
      real(dp), allocatable :: response(:, :), design(:, :), aligned(:)
      real(dp), allocatable :: augmented(:, :), inverse(:, :), coefficient(:, :)
      real(dp), allocatable :: residual0(:, :), residual1(:, :), rss0(:, :)
      real(dp), allocatable :: rss1(:, :), rss0_inverse(:, :), candidate(:)
      real(dp), allocatable :: covariance_ratio(:, :)
      real(dp) :: selected_significance, trace_value
      integer :: selected_order, variables, observations, regressors
      integer :: candidate_index, power, first, status, equation
      logical :: include_constant

      selected_order = 1
      if (present(order)) selected_order = order
      include_constant = .true.
      if (present(constant)) include_constant = constant
      selected_significance = 0.05_dp
      if (present(significance)) selected_significance = significance
      if (size(transitions, 1) /= size(series, 1) .or. &
         size(transitions, 2) < 1 .or. selected_significance <= 0.0_dp .or. &
         selected_significance >= 1.0_dp) then
         out%info = 1
         return
      end if
      call prepare_vlstar_data(series, transitions(:, 1), selected_order, &
         exogenous, include_constant, response, design, aligned)
      observations = size(response, 1)
      variables = size(response, 2)
      regressors = size(design, 2)
      allocate(inverse(regressors, regressors))
      call symmetric_pseudoinverse(matmul(transpose(design), design), &
         inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      coefficient = matmul(inverse, matmul(transpose(design), response))
      residual0 = response - matmul(design, coefficient)
      rss0 = matmul(transpose(residual0), residual0)
      allocate(rss0_inverse(variables, variables))
      call symmetric_pseudoinverse(rss0, rss0_inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      allocate(out%statistic(size(transitions, 2)), &
         out%p_value(size(transitions, 2)))
      allocate(augmented(observations, 3*regressors))
      do candidate_index = 1, size(transitions, 2)
         candidate = transitions(selected_order + 1:, candidate_index)
         do power = 1, 3
            first = (power - 1)*regressors + 1
            augmented(:, first:first + regressors - 1) = design* &
               spread(candidate**power, 2, regressors)
         end do
         deallocate(inverse)
         allocate(inverse(3*regressors, 3*regressors))
         call symmetric_pseudoinverse(matmul(transpose(augmented), augmented), &
            inverse, status)
         if (status /= 0) then
            out%statistic(candidate_index) = quiet_nan()
            out%p_value(candidate_index) = quiet_nan()
            cycle
         end if
         coefficient = matmul(inverse, matmul(transpose(augmented), residual0))
         residual1 = residual0 - matmul(augmented, coefficient)
         rss1 = matmul(transpose(residual1), residual1)
         covariance_ratio = matmul(rss0_inverse, rss1)
         trace_value = 0.0_dp
         do equation = 1, variables
            trace_value = trace_value + covariance_ratio(equation, equation)
         end do
         out%statistic(candidate_index) = max(0.0_dp, &
            real(observations, dp)*(real(variables, dp) - trace_value))
      end do
      out%degrees_of_freedom = 3*variables + regressors
      do candidate_index = 1, size(out%p_value)
         out%p_value(candidate_index) = regularized_gamma_q(0.5_dp* &
            real(out%degrees_of_freedom, dp), &
            0.5_dp*out%statistic(candidate_index))
      end do
      out%critical_value = chi_square_quantile_approx(1.0_dp - &
         0.5_dp*selected_significance, out%degrees_of_freedom)
      out%selected = minloc(out%p_value, dim=1)
   end function starvars_joint_linearity_test

   pure elemental real(dp) function chi_square_quantile_approx(probability, &
      degrees_of_freedom) result(value)
      !! Approximate a chi-square quantile with Wilson-Hilferty transformation.
      real(dp), intent(in) :: probability !! Cumulative probability.
      integer, intent(in) :: degrees_of_freedom !! Positive degrees of freedom.
      real(dp) :: degrees

      degrees = real(max(1, degrees_of_freedom), dp)
      value = degrees*(1.0_dp - 2.0_dp/(9.0_dp*degrees) + &
         normal_quantile(probability)*sqrt(2.0_dp/(9.0_dp*degrees)))**3
   end function chi_square_quantile_approx

   pure function starvars_long_run_variance(values) result(out)
      !! Translate lrvarbart's data-selected Bartlett long-run variance.
      real(dp), intent(in) :: values(:) !! Scalar observations.
      type(starvars_long_run_variance_t) :: out
      real(dp), allocatable :: autocorrelation(:), centered(:)
      real(dp) :: threshold, covariance_lag, weight
      integer :: observations, run_length, lag_max, lag, lookahead, stop_lag

      observations = size(values)
      if (observations < 3 .or. .not. all(ieee_is_finite(values))) then
         out%info = 1
         return
      end if
      centered = values - sum(values)/real(observations, dp)
      run_length = max(1, int(real(observations, dp)**(1.0_dp/3.0_dp)))
      lag_max = min(observations - 1, max(run_length, &
         int(2.0_dp*real(observations, dp)**(2.0_dp/3.0_dp))))
      autocorrelation = acf_values(values, lag_max)
      threshold = 1.4_dp*sqrt(log10(real(observations, dp))/ &
         real(observations, dp))
      stop_lag = max(1, run_length)
      do lag = 1, max(1, lag_max - run_length + 1)
         lookahead = min(lag_max, lag + run_length - 1)
         if (maxval(abs(autocorrelation(lag + 1:lookahead + 1))) < &
            threshold) then
            stop_lag = lag
            exit
         end if
      end do
      out%bandwidth = stop_lag
      out%variance = sum(centered**2)/real(observations, dp)
      do lag = 1, min(lag_max, 2*out%bandwidth)
         weight = real(2*out%bandwidth - lag + 1, dp)/ &
            real(2*out%bandwidth, dp)
         covariance_lag = sum(centered(:observations - lag)* &
            centered(lag + 1:))/real(observations, dp)
         out%variance = out%variance + 2.0_dp*weight*covariance_lag
      end do
      out%variance = max(out%variance, tiny(1.0_dp))
   end function starvars_long_run_variance

   pure function starvars_multiple_cumsum(data, confidence, max_breaks) &
      result(out)
      !! Translate multiCUMSUM's recursive common covariance-break tests.
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable data matrix.
      real(dp), intent(in), optional :: confidence !! Confidence level: 0.90, 0.95, or 0.99.
      integer, intent(in), optional :: max_breaks !! Maximum reported breaks from one through seven.
      type(starvars_cumsum_t) :: out
      real(dp), allocatable :: segment_lambda(:), segment_omega(:)
      integer, allocatable :: segment_start(:), segment_end(:), segment_break(:)
      real(dp) :: selected_confidence, best_value
      integer :: limit, segment_count, selected_segment, iteration, variables
      integer :: covariance_dimension, confidence_row, dimension_column

      selected_confidence = 0.95_dp
      if (present(confidence)) selected_confidence = confidence
      limit = 7
      if (present(max_breaks)) limit = max_breaks
      variables = size(data, 2)
      covariance_dimension = variables*(variables + 1)/2
      if (size(data, 1) < 8 .or. variables < 1 .or. limit < 1 .or. &
         limit > 7 .or. .not. all(ieee_is_finite(data))) then
         out%info = 1
         return
      end if
      confidence_row = nearest_confidence_row(selected_confidence)
      if (confidence_row == 0) then
         out%info = 2
         return
      end if
      dimension_column = cumsum_dimension_column(covariance_dimension)
      out%lambda_critical = cumsum_lambda_critical(confidence_row, &
         dimension_column)
      out%omega_critical = cumsum_omega_critical(confidence_row, &
         dimension_column)
      allocate(out%lambda(limit), out%omega(limit), &
         out%break_location(limit))
      allocate(segment_start(limit + 1), segment_end(limit + 1), &
         segment_break(limit + 1), segment_lambda(limit + 1), &
         segment_omega(limit + 1))
      out%lambda = quiet_nan()
      out%omega = quiet_nan()
      out%break_location = 0
      segment_count = 1
      segment_start(1) = 1
      segment_end(1) = size(data, 1)
      call cumsum_segment(data, segment_start(1), segment_end(1), &
         segment_lambda(1), segment_omega(1), segment_break(1))
      do iteration = 1, limit
         selected_segment = 0
         best_value = -huge(1.0_dp)
         do variables = 1, segment_count
            if (segment_lambda(variables) > best_value .and. &
               segment_break(variables) > segment_start(variables) .and. &
               segment_break(variables) < segment_end(variables)) then
               selected_segment = variables
               best_value = segment_lambda(variables)
            end if
         end do
         if (selected_segment == 0) exit
         out%lambda(iteration) = segment_lambda(selected_segment)
         out%omega(iteration) = segment_omega(selected_segment)
         out%break_location(iteration) = segment_break(selected_segment)
         segment_count = segment_count + 1
         segment_start(segment_count) = segment_break(selected_segment) + 1
         segment_end(segment_count) = segment_end(selected_segment)
         segment_end(selected_segment) = segment_break(selected_segment) - 1
         call cumsum_segment(data, segment_start(selected_segment), &
            segment_end(selected_segment), segment_lambda(selected_segment), &
            segment_omega(selected_segment), segment_break(selected_segment))
         call cumsum_segment(data, segment_start(segment_count), &
            segment_end(segment_count), segment_lambda(segment_count), &
            segment_omega(segment_count), segment_break(segment_count))
      end do
   end function starvars_multiple_cumsum

   pure subroutine cumsum_segment(data, first, last, lambda, omega, &
      break_location)
      !! Evaluate the covariance CUMSUM process on one inclusive segment.
      real(dp), intent(in) :: data(:, :) !! Complete data matrix.
      integer, intent(in) :: first !! First segment observation.
      integer, intent(in) :: last !! Last segment observation.
      real(dp), intent(out) :: lambda !! Maximum quadratic CUMSUM statistic.
      real(dp), intent(out) :: omega !! Mean quadratic CUMSUM statistic.
      integer, intent(out) :: break_location !! Maximizing absolute observation index.
      real(dp), allocatable :: products(:, :), total(:), cumulative(:)
      real(dp), allocatable :: inverse_variance(:), statistic(:)
      type(starvars_long_run_variance_t) :: long_run
      integer :: observations, variables, dimension, row, column

      observations = last - first + 1
      variables = size(data, 2)
      dimension = variables*(variables + 1)/2
      if (observations < 4) then
         lambda = -huge(1.0_dp)
         omega = quiet_nan()
         break_location = first
         return
      end if
      allocate(products(observations, dimension))
      do row = 1, observations
         products(row, :) = lower_outer_product(data(first + row - 1, :))
      end do
      total = sum(products, dim=1)
      allocate(inverse_variance(dimension))
      do column = 1, dimension
         long_run = starvars_long_run_variance([0.0_dp, products(:, column)])
         inverse_variance(column) = 1.0_dp/max(long_run%variance, tiny(1.0_dp))
      end do
      allocate(statistic(observations), cumulative(dimension))
      cumulative = 0.0_dp
      do row = 1, observations
         cumulative = cumulative + products(row, :)
         statistic(row) = sum(((cumulative - real(row, dp)/ &
            real(observations, dp)*total)**2)*inverse_variance)/ &
            real(observations, dp)
      end do
      lambda = maxval(statistic)
      omega = sum(statistic)/real(observations, dp)
      break_location = first - 1 + maxloc(statistic, dim=1)
   end subroutine cumsum_segment

   pure function lower_outer_product(values) result(packed)
      !! Pack the lower triangle of one vector outer product by columns.
      real(dp), intent(in) :: values(:) !! Input vector.
      real(dp) :: packed(size(values)*(size(values) + 1)/2)
      integer :: row, column, index

      index = 0
      do column = 1, size(values)
         do row = column, size(values)
            index = index + 1
            packed(index) = values(row)*values(column)
         end do
      end do
   end function lower_outer_product

   pure integer function nearest_confidence_row(confidence) result(row)
      !! Map a supported confidence level to the STARVARS table row.
      real(dp), intent(in) :: confidence !! Confidence level.

      if (abs(confidence - 0.90_dp) < 1.0e-8_dp) then
         row = 1
      else if (abs(confidence - 0.95_dp) < 1.0e-8_dp) then
         row = 2
      else if (abs(confidence - 0.99_dp) < 1.0e-8_dp) then
         row = 3
      else
         row = 0
      end if
   end function nearest_confidence_row

   pure integer function cumsum_dimension_column(dimension) result(column)
      !! Map a covariance-vector dimension to the STARVARS table column.
      integer, intent(in) :: dimension !! Packed covariance dimension.

      if (dimension <= 10) then
         column = 1
      else if (dimension <= 15) then
         column = 2
      else if (dimension <= 20) then
         column = 3
      else if (dimension <= 50) then
         column = 4
      else if (dimension <= 100) then
         column = 5
      else if (dimension <= 200) then
         column = 6
      else if (dimension <= 500) then
         column = 7
      else
         column = 8
      end if
   end function cumsum_dimension_column

   pure real(dp) function cumsum_lambda_critical(row, column) result(value)
      !! Return the STARVARS maximum-CUMSUM critical value.
      integer, intent(in) :: row !! Confidence table row.
      integer, intent(in) :: column !! Dimension table column.
      real(dp), parameter :: table(3, 8) = reshape([ &
         2.64_dp, 3.17_dp, 4.28_dp, 2.53_dp, 3.02_dp, 4.04_dp, &
         2.46_dp, 2.92_dp, 3.89_dp, 2.27_dp, 2.69_dp, 3.53_dp, &
         2.16_dp, 2.55_dp, 3.33_dp, 2.06_dp, 2.44_dp, 3.18_dp, &
         1.96_dp, 2.33_dp, 3.04_dp, 1.28_dp, 1.64_dp, 2.33_dp], [3, 8])

      value = table(row, column)
   end function cumsum_lambda_critical

   pure real(dp) function cumsum_omega_critical(row, column) result(value)
      !! Return the STARVARS integrated-CUMSUM critical value.
      integer, intent(in) :: row !! Confidence table row.
      integer, intent(in) :: column !! Dimension table column.
      real(dp), parameter :: table(3, 8) = reshape([ &
         1.33_dp, 1.84_dp, 2.90_dp, 1.33_dp, 1.81_dp, 2.80_dp, &
         1.32_dp, 1.79_dp, 2.74_dp, 1.31_dp, 1.74_dp, 2.59_dp, &
         1.31_dp, 1.71_dp, 2.51_dp, 1.30_dp, 1.69_dp, 2.46_dp, &
         1.29_dp, 1.68_dp, 2.41_dp, 1.28_dp, 1.64_dp, 2.33_dp], [3, 8])

      value = table(row, column)
   end function cumsum_omega_critical

   pure function starvars_realized_covariance_groups(data, groups, make_returns, &
      return_scale, compute_cholesky) result(out)
      !! Translate rcov aggregation using caller-supplied integer periods.
      real(dp), intent(in) :: data(:, :) !! Prices or returns by observation and asset.
      integer, intent(in) :: groups(:) !! Period identifier for each observation.
      logical, intent(in), optional :: make_returns !! Convert positive prices to log returns.
      real(dp), intent(in), optional :: return_scale !! Multiplier applied to log returns.
      logical, intent(in), optional :: compute_cholesky !! Return packed lower Cholesky factors.
      type(starvars_realized_covariance_t) :: out
      real(dp), allocatable :: covariance_matrix(:, :), lower(:, :)
      real(dp) :: scale
      integer :: observations, variables, dimension, group_count
      integer :: observation, group_index, row, column, index, status
      logical :: convert, factorize, found

      observations = size(data, 1)
      variables = size(data, 2)
      convert = .true.
      if (present(make_returns)) convert = make_returns
      scale = 100.0_dp
      if (present(return_scale)) scale = return_scale
      factorize = .false.
      if (present(compute_cholesky)) factorize = compute_cholesky
      if (observations < 1 .or. variables < 1 .or. &
         size(groups) /= observations .or. scale <= 0.0_dp .or. &
         .not. all(ieee_is_finite(data)) .or. &
         (convert .and. any(data <= 0.0_dp))) then
         out%info = 1
         return
      end if
      allocate(out%returns(observations, variables))
      if (convert) then
         out%returns(1, :) = 0.0_dp
         do observation = 2, observations
            out%returns(observation, :) = scale*log(data(observation, :)/ &
               data(observation - 1, :))
         end do
      else
         out%returns = data
      end if
      allocate(out%groups(observations))
      group_count = 0
      do observation = 1, observations
         found = .false.
         if (group_count > 0) &
            found = any(out%groups(:group_count) == groups(observation))
         if (.not. found) then
            group_count = group_count + 1
            out%groups(group_count) = groups(observation)
         end if
      end do
      out%groups = out%groups(:group_count)
      dimension = variables*(variables + 1)/2
      allocate(out%covariance(group_count, dimension), source=0.0_dp)
      do observation = 1, observations
         group_index = find_group(out%groups, groups(observation))
         out%covariance(group_index, :) = out%covariance(group_index, :) + &
            lower_outer_product(out%returns(observation, :))
      end do
      if (factorize) then
         allocate(out%cholesky(group_count, dimension), source=0.0_dp)
         allocate(covariance_matrix(variables, variables), &
            lower(variables, variables))
         do group_index = 1, group_count
            covariance_matrix = unpack_lower_triangle( &
               out%covariance(group_index, :), variables)
            call cholesky_lower(covariance_matrix, lower, status)
            if (status /= 0) then
               out%info = 2
               return
            end if
            index = 0
            do column = 1, variables
               do row = column, variables
                  index = index + 1
                  out%cholesky(group_index, index) = lower(row, column)
               end do
            end do
         end do
      else
         allocate(out%cholesky(0, 0))
      end if
   end function starvars_realized_covariance_groups

   pure function starvars_realized_covariance_dates(data, dates, frequency, &
      make_returns, return_scale, compute_cholesky) result(out)
      !! Aggregate rcov values by Gregorian day, month, quarter, or year.
      real(dp), intent(in) :: data(:, :) !! Prices or returns by observation and asset.
      type(date_t), intent(in) :: dates(:) !! Gregorian date for each observation.
      integer, intent(in), optional :: frequency !! Calendar aggregation frequency code.
      logical, intent(in), optional :: make_returns !! Convert positive prices to log returns.
      real(dp), intent(in), optional :: return_scale !! Multiplier applied to log returns.
      logical, intent(in), optional :: compute_cholesky !! Return packed lower Cholesky factors.
      type(starvars_realized_covariance_t) :: out
      integer, allocatable :: groups(:)
      integer :: selected_frequency, observation

      selected_frequency = starvars_frequency_daily
      if (present(frequency)) selected_frequency = frequency
      if (size(dates) /= size(data, 1) .or. .not. all(date_valid(dates)) .or. &
         selected_frequency < starvars_frequency_daily .or. &
         selected_frequency > starvars_frequency_yearly) then
         out%info = 1
         return
      end if
      allocate(groups(size(dates)))
      do observation = 1, size(dates)
         select case (selected_frequency)
         case (starvars_frequency_daily)
            groups(observation) = date_day_number(dates(observation))
         case (starvars_frequency_monthly)
            groups(observation) = 12*dates(observation)%year + &
               dates(observation)%month
         case (starvars_frequency_quarterly)
            groups(observation) = 4*dates(observation)%year + &
               (dates(observation)%month - 1)/3
         case (starvars_frequency_yearly)
            groups(observation) = dates(observation)%year
         end select
      end do
      out = starvars_realized_covariance_groups(data, groups, make_returns, &
         return_scale, compute_cholesky)
   end function starvars_realized_covariance_dates

   pure integer function find_group(groups, target) result(location)
      !! Locate a period identifier in the compact group vector.
      integer, intent(in) :: groups(:) !! Unique group identifiers.
      integer, intent(in) :: target !! Sought identifier.
      integer :: index

      location = 0
      do index = 1, size(groups)
         if (groups(index) == target) then
            location = index
            return
         end if
      end do
   end function find_group

   pure function unpack_lower_triangle(packed, order) result(matrix)
      !! Expand a packed symmetric lower triangle into a full matrix.
      real(dp), intent(in) :: packed(:) !! Packed lower triangle.
      integer, intent(in) :: order !! Matrix order.
      real(dp) :: matrix(order, order)
      integer :: row, column, index

      matrix = 0.0_dp
      index = 0
      do column = 1, order
         do row = column, order
            index = index + 1
            matrix(row, column) = packed(index)
            matrix(column, row) = packed(index)
         end do
      end do
   end function unpack_lower_triangle

   subroutine display_starvars_start(value, unit)
      !! Display selected VLSTAR transition starting values.
      type(starvars_start_t), intent(in) :: value !! Selected starting values.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: selected_unit, transition

      selected_unit = output_unit
      if (present(unit)) selected_unit = unit
      write(selected_unit, '(a)') 'VLSTAR starting values'
      if (allocated(value%gamma)) then
         do transition = 1, size(value%gamma, 2)
            write(selected_unit, '(a,i0)') '  transition: ', transition
            write(selected_unit, '(a,*(1x,es14.6))') '    gamma:', &
               value%gamma(:, transition)
            write(selected_unit, '(a,*(1x,es14.6))') '    location:', &
               value%location(:, transition)
         end do
      end if
      write(selected_unit, '(a,es14.6)') '  objective: ', value%objective
      write(selected_unit, '(a,i0)') '  status: ', value%info
   end subroutine display_starvars_start

   subroutine display_starvars_fit(value, unit)
      !! Display the principal VLSTAR fit summary.
      type(starvars_fit_t), intent(in) :: value !! Fitted VLSTAR model.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: selected_unit

      selected_unit = output_unit
      if (present(unit)) selected_unit = unit
      write(selected_unit, '(a)') 'Vector logistic STAR fit'
      write(selected_unit, '(a,i0)') '  order: ', value%order
      write(selected_unit, '(a,i0)') '  regimes: ', value%regimes
      write(selected_unit, '(a,l1)') '  common transition: ', &
         value%common_transition
      write(selected_unit, '(a,es14.6)') '  residual sum of squares: ', &
         value%rss
      write(selected_unit, '(a,es14.6)') '  log likelihood: ', &
         value%log_likelihood
      write(selected_unit, '(a,l1)') '  converged: ', value%converged
      write(selected_unit, '(a,i0)') '  status: ', value%info
   end subroutine display_starvars_fit

   subroutine display_starvars_forecast(value, unit, print_paths)
      !! Display a VLSTAR forecast and optionally its retained simulation paths.
      type(starvars_forecast_t), intent(in) :: value !! Forecast result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_paths !! Print retained simulation paths.
      integer :: selected_unit, row
      logical :: show_paths

      selected_unit = output_unit
      if (present(unit)) selected_unit = unit
      show_paths = .false.
      if (present(print_paths)) show_paths = print_paths
      write(selected_unit, '(a)') 'VLSTAR forecast'
      write(selected_unit, '(a,i0)') '  method: ', value%method
      write(selected_unit, '(a,i0)') '  simulations: ', value%simulations
      if (allocated(value%point)) then
         write(selected_unit, '(a)') '  point forecasts:'
         do row = 1, size(value%point, 1)
            write(selected_unit, '(*(1x,es14.6))') value%point(row, :)
         end do
      end if
      if (show_paths .and. allocated(value%paths)) &
         write(selected_unit, '(*(1x,es14.6))') value%paths
   end subroutine display_starvars_forecast

   subroutine display_starvars_joint_test(value, unit)
      !! Display joint-linearity statistics and p-values.
      type(starvars_joint_test_t), intent(in) :: value !! Joint-linearity test result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: selected_unit, candidate

      selected_unit = output_unit
      if (present(unit)) selected_unit = unit
      write(selected_unit, '(a)') 'VLSTAR joint linearity test'
      write(selected_unit, '(a,i0)') '  degrees of freedom: ', &
         value%degrees_of_freedom
      if (allocated(value%statistic)) then
         do candidate = 1, size(value%statistic)
            write(selected_unit, '(a,i0,a,es14.6,a,es14.6)') '  candidate ', &
               candidate, ': statistic=', value%statistic(candidate), &
               ', p=', value%p_value(candidate)
         end do
      end if
   end subroutine display_starvars_joint_test

   subroutine display_starvars_long_run_variance(value, unit)
      !! Display a Bartlett long-run variance estimate and bandwidth.
      type(starvars_long_run_variance_t), intent(in) :: value !! Long-run variance result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: selected_unit

      selected_unit = output_unit
      if (present(unit)) selected_unit = unit
      write(selected_unit, '(a)') 'Bartlett long-run variance'
      write(selected_unit, '(a,es14.6)') '  variance: ', value%variance
      write(selected_unit, '(a,i0)') '  bandwidth: ', value%bandwidth
      write(selected_unit, '(a,i0)') '  status: ', value%info
   end subroutine display_starvars_long_run_variance

   subroutine display_starvars_cumsum(value, unit)
      !! Display covariance-CUMSUM break locations and statistics.
      type(starvars_cumsum_t), intent(in) :: value !! Multiple-CUMSUM result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: selected_unit, break_index

      selected_unit = output_unit
      if (present(unit)) selected_unit = unit
      write(selected_unit, '(a)') 'Multiple covariance CUMSUM'
      if (allocated(value%break_location)) then
         do break_index = 1, size(value%break_location)
            if (value%break_location(break_index) <= 0) cycle
            write(selected_unit, '(a,i0,a,i0,a,es14.6,a,es14.6)') &
               '  break ', break_index, ': index=', &
               value%break_location(break_index), ', lambda=', &
               value%lambda(break_index), ', omega=', value%omega(break_index)
         end do
      end if
   end subroutine display_starvars_cumsum

   subroutine display_starvars_realized_covariance(value, unit, print_values)
      !! Display realized-covariance dimensions and optionally packed values.
      type(starvars_realized_covariance_t), intent(in) :: value !! Realized covariance result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_values !! Print packed covariance values.
      integer :: selected_unit, group
      logical :: show_values

      selected_unit = output_unit
      if (present(unit)) selected_unit = unit
      show_values = .false.
      if (present(print_values)) show_values = print_values
      write(selected_unit, '(a)') 'Realized covariance aggregation'
      if (allocated(value%covariance)) then
         write(selected_unit, '(a,i0)') '  periods: ', size(value%covariance, 1)
         write(selected_unit, '(a,i0)') '  packed dimension: ', &
            size(value%covariance, 2)
         if (show_values) then
            do group = 1, size(value%covariance, 1)
               write(selected_unit, '(a,i0,a,*(1x,es14.6))') '  period ', &
                  value%groups(group), ':', value%covariance(group, :)
            end do
         end if
      end if
      write(selected_unit, '(a,i0)') '  status: ', value%info
   end subroutine display_starvars_realized_covariance

end module starvars_mod
