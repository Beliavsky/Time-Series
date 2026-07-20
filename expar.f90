! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R EXPAR package.
module expar_mod
   !! Amplitude-dependent exponential autoregression translated from EXPAR.
   use kind_mod, only: dp
   use stats_mod, only: ols_fit
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd
   use utils_mod, only: lowercase
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan
   implicit none
   private

   type, public :: expar_fit_t
      !! Fitted exponential autoregression and conditional residual diagnostics.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: phi(:)
      real(dp), allocatable :: pi(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: gamma = 0.0_dp
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: order = 0
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: converged = .false.
   end type expar_fit_t

   type, public :: expar_selection_t
      !! EXPAR order-selection scores and selected fitted model.
      type(expar_fit_t) :: model
      real(dp), allocatable :: aic(:)
      real(dp), allocatable :: aicc(:)
      real(dp), allocatable :: bic(:)
      character(len=4) :: criterion = 'AIC'
      integer :: selected_order = 0
      integer :: info = 0
   end type expar_selection_t

   type, public :: expar_forecast_t
      !! Recursive EXPAR point forecasts from the end of the fitted sample.
      real(dp), allocatable :: mean(:)
      integer :: info = 0
   end type expar_forecast_t

   public :: expar_evaluate, expar_initial_parameters, expar_fit
   public :: expar_select, expar_forecast, expar_transition

contains

   pure function expar_evaluate(series, parameters) result(out)
      !! Evaluate an EXPAR parameter vector and its information criteria.
      real(dp), intent(in) :: series(:) !! Observed time series.
      real(dp), intent(in) :: parameters(:) !! Phi, pi, and gamma parameters.
      type(expar_fit_t) :: out
      real(dp) :: transition, prediction, nan_value
      integer :: observations, parameter_count, order, time, lag

      observations = size(series)
      parameter_count = size(parameters)
      if (parameter_count < 3 .or. mod(parameter_count - 1, 2) /= 0 .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      order = (parameter_count - 1)/2
      if (observations <= order) then
         out%info = 1
         return
      end if
      out%data = series
      out%phi = parameters(:order)
      out%pi = parameters(order + 1:2*order)
      out%gamma = parameters(parameter_count)
      out%order = order
      allocate(out%fitted(observations), out%residuals(observations))
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      out%fitted = nan_value
      out%residuals = 0.0_dp
      do time = order + 1, observations
         transition = expar_transition(out%gamma, series(time - 1))
         if (.not. ieee_is_finite(transition)) then
            out%info = 2
            return
         end if
         prediction = 0.0_dp
         do lag = 1, order
            prediction = prediction + (out%phi(lag) + &
               out%pi(lag)*transition)*series(time - lag)
         end do
         out%fitted(time) = prediction
         out%residuals(time) = series(time) - prediction
      end do
      out%rss = sum(out%residuals**2)
      if (.not. ieee_is_finite(out%rss) .or. out%rss <= 0.0_dp) then
         out%info = 2
         return
      end if
      out%aic = 2.0_dp*real(parameter_count, dp) + &
         real(observations, dp)*log(out%rss/real(observations, dp))
      if (observations > parameter_count + 1) then
         out%aicc = out%aic + &
            2.0_dp*real(parameter_count*(parameter_count + 1), dp)/ &
            real(observations - parameter_count - 1, dp)
      end if
      out%bic = real(parameter_count, dp)*log(real(observations, dp)) + &
         real(observations, dp)*log(out%rss/real(observations, dp))
   end function expar_evaluate

   pure function expar_initial_parameters(series, order) result(parameters)
      !! Construct EXPAR starting values from an ordinary AR regression.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: order !! Autoregressive order.
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: design(:, :), response(:), coefficient(:)
      real(dp), allocatable :: standard_error(:), residual(:)
      real(dp) :: rss
      integer :: rows, row, time, lag, status

      if (order < 1 .or. size(series) <= 2*order + 2 .or. &
         .not. all(ieee_is_finite(series))) then
         allocate(parameters(0))
         return
      end if
      rows = size(series) - order
      allocate(design(rows, order), response(rows))
      do row = 1, rows
         time = order + row
         response(row) = series(time)
         do lag = 1, order
            design(row, lag) = series(time - lag)
         end do
      end do
      call ols_fit(design, response, coefficient, standard_error, residual, &
         rss, status)
      allocate(parameters(2*order + 1))
      parameters = 0.5_dp
      if (status == 0 .and. size(coefficient) == order) then
         parameters(:order) = coefficient
      else
         parameters(:order) = 0.0_dp
      end if
      parameters(2*order + 1) = 0.5_dp
   end function expar_initial_parameters

   pure function expar_fit(series, order, initial, max_iterations, &
      tolerance) result(out)
      !! Estimate an EXPAR model by finite-difference BFGS minimization of RSS.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: order !! Autoregressive order.
      real(dp), intent(in), optional :: initial(:) !! Initial phi, pi, and gamma values.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(expar_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: starting(:)
      real(dp) :: selected_tolerance
      integer :: limit

      limit = 300
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (order < 1 .or. size(series) <= 2*order + 2 .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (present(initial)) then
         if (size(initial) /= 2*order + 1 .or. &
            .not. all(ieee_is_finite(initial))) then
            out%info = 2
            return
         end if
         starting = initial
      else
         starting = expar_initial_parameters(series, order)
      end if
      if (size(starting) /= 2*order + 1) then
         out%info = 2
         return
      end if
      optimization = bfgs_minimize_fd(objective, starting, limit, &
         selected_tolerance)
      if (.not. allocated(optimization%parameters) .or. &
         .not. ieee_is_finite(optimization%objective)) then
         out%info = 3
         return
      end if
      out = expar_evaluate(series, optimization%parameters)
      out%iterations = optimization%iterations
      out%optimizer_info = optimization%info
      out%converged = optimization%converged

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Evaluate residual sum of squares for the optimizer.
         real(dp), intent(in) :: parameters(:) !! Candidate EXPAR parameters.
         type(expar_fit_t) :: evaluated

         evaluated = expar_evaluate(series, parameters)
         value = huge(1.0_dp)
         if (evaluated%info == 0) value = evaluated%rss
      end function objective

   end function expar_fit

   pure function expar_select(series, maximum_order, criterion, &
      max_iterations, tolerance) result(out)
      !! Select an EXPAR order by AIC, corrected AIC, or BIC.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: maximum_order !! Largest autoregressive order considered.
      character(len=*), intent(in), optional :: criterion !! `AIC`, `AICc`, or `BIC`.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations per order.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(expar_selection_t) :: out
      type(expar_fit_t) :: candidate
      character(len=:), allocatable :: selected_criterion
      real(dp) :: score, best_score
      integer :: order

      selected_criterion = 'aic'
      if (present(criterion)) selected_criterion = lowercase(trim(criterion))
      if (maximum_order < 1 .or. &
         (selected_criterion /= 'aic' .and. selected_criterion /= 'aicc' .and. &
         selected_criterion /= 'aic_c' .and. selected_criterion /= 'bic')) then
         out%info = 1
         return
      end if
      allocate(out%aic(maximum_order), out%aicc(maximum_order), &
         out%bic(maximum_order))
      out%aic = huge(1.0_dp)
      out%aicc = huge(1.0_dp)
      out%bic = huge(1.0_dp)
      best_score = huge(1.0_dp)
      do order = 1, maximum_order
         candidate = expar_fit(series, order, max_iterations=max_iterations, &
            tolerance=tolerance)
         if (candidate%info /= 0) cycle
         out%aic(order) = candidate%aic
         out%aicc(order) = candidate%aicc
         out%bic(order) = candidate%bic
         select case (selected_criterion)
         case ('aic')
            score = candidate%aic
         case ('aicc', 'aic_c')
            score = candidate%aicc
         case default
            score = candidate%bic
         end select
         if (score < best_score) then
            best_score = score
            out%selected_order = order
            out%model = candidate
         end if
      end do
      if (out%selected_order == 0) then
         out%info = 2
         return
      end if
      select case (selected_criterion)
      case ('aic')
         out%criterion = 'AIC'
      case ('aicc', 'aic_c')
         out%criterion = 'AICc'
      case default
         out%criterion = 'BIC'
      end select
   end function expar_select

   pure function expar_forecast(model, horizon) result(out)
      !! Compute recursive EXPAR forecasts from the observed sample endpoint.
      type(expar_fit_t), intent(in) :: model !! Fitted EXPAR model.
      integer, intent(in) :: horizon !! Number of future observations.
      type(expar_forecast_t) :: out
      real(dp), allocatable :: work(:)
      real(dp) :: transition, prediction
      integer :: observations, step, time, lag

      if (model%info /= 0 .or. horizon < 1 .or. model%order < 1 .or. &
         .not. allocated(model%data) .or. .not. allocated(model%phi) .or. &
         .not. allocated(model%pi)) then
         out%info = 1
         return
      end if
      observations = size(model%data)
      allocate(work(observations + horizon), out%mean(horizon))
      work(:observations) = model%data
      do step = 1, horizon
         time = observations + step
         transition = expar_transition(model%gamma, work(time - 1))
         if (.not. ieee_is_finite(transition)) then
            out%info = 2
            return
         end if
         prediction = 0.0_dp
         do lag = 1, model%order
            prediction = prediction + (model%phi(lag) + &
               model%pi(lag)*transition)*work(time - lag)
         end do
         work(time) = prediction
         out%mean(step) = prediction
      end do
   end function expar_forecast

   pure elemental real(dp) function expar_transition(gamma, amplitude) &
      result(weight)
      !! Evaluate the exponential amplitude transition weight.
      real(dp), intent(in) :: gamma !! Exponential transition scale.
      real(dp), intent(in) :: amplitude !! Previous series amplitude.

      weight = exp(-gamma*amplitude**2)
   end function expar_transition

end module expar_mod
