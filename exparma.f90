! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R EXPARMA package.
module exparma_mod
   !! Exponential autoregressive moving-average algorithms from EXPARMA.
   use kind_mod, only: dp
   use expar_mod, only: expar_transition
   use itsmr_mod, only: itsmr_arma_model_t, hannan_rissanen_fit
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan
   implicit none
   private

   type, public :: exparma_fit_t
      !! Fitted amplitude-dependent EXPARMA model and residual diagnostics.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: phi(:)
      real(dp), allocatable :: pi(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: gamma = 0.0_dp
      real(dp) :: omega = 0.0_dp
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: ar_order = 0
      integer :: ma_order = 0
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: converged = .false.
   end type exparma_fit_t

   type, public :: exparma_selection_t
      !! EXPARMA AIC table and selected autoregressive and moving-average orders.
      type(exparma_fit_t) :: model
      real(dp), allocatable :: aic(:, :)
      integer :: selected_ar_order = 0
      integer :: selected_ma_order = 0
      integer :: info = 0
   end type exparma_selection_t

   public :: exparma_evaluate, exparma_initial_parameters
   public :: exparma_fit, exparma_select

contains

   pure function exparma_evaluate(series, ar_order, ma_order, parameters) &
      result(out)
      !! Evaluate fixed EXPARMA parameters, residuals, RSS, and package AIC.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      real(dp), intent(in) :: parameters(:) !! Phi, pi, theta, delta, gamma, and omega.
      type(exparma_fit_t) :: out
      real(dp) :: ar_weight, ma_weight, prediction, nan_value
      integer :: observations, parameter_count, start, time, lag, offset

      observations = size(series)
      parameter_count = 2*(ar_order + ma_order + 1)
      start = max(ar_order, ma_order) + 1
      if (ar_order < 1 .or. ma_order < 1 .or. &
         size(parameters) /= parameter_count .or. observations < start .or. &
         observations <= parameter_count .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      out%data = series
      out%phi = parameters(:ar_order)
      out%pi = parameters(ar_order + 1:2*ar_order)
      offset = 2*ar_order
      out%theta = parameters(offset + 1:offset + ma_order)
      offset = offset + ma_order
      out%delta = parameters(offset + 1:offset + ma_order)
      out%gamma = parameters(parameter_count - 1)
      out%omega = parameters(parameter_count)
      out%ar_order = ar_order
      out%ma_order = ma_order
      allocate(out%fitted(observations), out%residuals(observations))
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      out%fitted = nan_value
      out%residuals = 0.0_dp
      do time = start, observations
         ar_weight = expar_transition(out%gamma, series(time - 1))
         ma_weight = expar_transition(out%omega, out%residuals(time - 1))
         if (.not. ieee_is_finite(ar_weight) .or. &
            .not. ieee_is_finite(ma_weight)) then
            out%info = 2
            return
         end if
         prediction = 0.0_dp
         do lag = 1, ar_order
            prediction = prediction + (out%phi(lag) + &
               out%pi(lag)*ar_weight)*series(time - lag)
         end do
         do lag = 1, ma_order
            prediction = prediction + (out%theta(lag) + &
               out%delta(lag)*ma_weight)*out%residuals(time - lag)
         end do
         out%fitted(time) = prediction
         out%residuals(time) = series(time) - prediction
      end do
      out%rss = sum(out%residuals**2)
      if (.not. ieee_is_finite(out%rss) .or. out%rss <= 0.0_dp) then
         out%info = 2
         return
      end if
      out%aic = real(observations, dp)*log(out%rss/ &
         real(observations - parameter_count, dp)) + &
         2.0_dp*real(parameter_count, dp)
   end function exparma_evaluate

   pure function exparma_initial_parameters(series, ar_order, ma_order) &
      result(parameters)
      !! Construct EXPARMA starting values from a Hannan-Rissanen ARMA fit.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      real(dp), allocatable :: parameters(:)
      type(itsmr_arma_model_t) :: arma
      integer :: parameter_count, offset

      if (ar_order < 1 .or. ma_order < 1 .or. &
         size(series) <= 2*(ar_order + ma_order + 1) .or. &
         .not. all(ieee_is_finite(series))) then
         allocate(parameters(0))
         return
      end if
      parameter_count = 2*(ar_order + ma_order + 1)
      allocate(parameters(parameter_count))
      parameters = 0.5_dp
      arma = hannan_rissanen_fit(series, ar_order, ma_order)
      if (arma%info == 0) then
         parameters(:ar_order) = arma%ar
         offset = 2*ar_order
         parameters(offset + 1:offset + ma_order) = arma%ma
      else
         parameters(:ar_order) = 0.0_dp
         offset = 2*ar_order
         parameters(offset + 1:offset + ma_order) = 0.0_dp
      end if
   end function exparma_initial_parameters

   pure function exparma_fit(series, ar_order, ma_order, initial, &
      max_iterations, tolerance) result(out)
      !! Estimate an EXPARMA model by finite-difference BFGS RSS minimization.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      real(dp), intent(in), optional :: initial(:) !! Initial model parameters.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(exparma_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: starting(:)
      real(dp) :: selected_tolerance
      integer :: limit, parameter_count

      parameter_count = 2*(ar_order + ma_order + 1)
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (ar_order < 1 .or. ma_order < 1 .or. &
         size(series) <= parameter_count .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (present(initial)) then
         if (size(initial) /= parameter_count .or. &
            .not. all(ieee_is_finite(initial))) then
            out%info = 2
            return
         end if
         starting = initial
      else
         starting = exparma_initial_parameters(series, ar_order, ma_order)
      end if
      if (size(starting) /= parameter_count) then
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
      out = exparma_evaluate(series, ar_order, ma_order, &
         optimization%parameters)
      out%iterations = optimization%iterations
      out%optimizer_info = optimization%info
      out%converged = optimization%converged

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Evaluate EXPARMA RSS for the optimizer.
         real(dp), intent(in) :: parameters(:) !! Candidate EXPARMA parameters.
         type(exparma_fit_t) :: evaluated

         evaluated = exparma_evaluate(series, ar_order, ma_order, parameters)
         value = huge(1.0_dp)
         if (evaluated%info == 0) value = evaluated%rss
      end function objective

   end function exparma_fit

   pure function exparma_select(series, maximum_ar_order, maximum_ma_order, &
      max_iterations, tolerance) result(out)
      !! Select EXPARMA autoregressive and moving-average orders by minimum AIC.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: maximum_ar_order !! Largest AR order considered.
      integer, intent(in) :: maximum_ma_order !! Largest MA order considered.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations per model.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(exparma_selection_t) :: out
      type(exparma_fit_t) :: candidate
      real(dp) :: best_aic
      integer :: ar_order, ma_order

      if (maximum_ar_order < 1 .or. maximum_ma_order < 1) then
         out%info = 1
         return
      end if
      allocate(out%aic(maximum_ar_order, maximum_ma_order))
      out%aic = huge(1.0_dp)
      best_aic = huge(1.0_dp)
      do ar_order = 1, maximum_ar_order
         do ma_order = 1, maximum_ma_order
            candidate = exparma_fit(series, ar_order, ma_order, &
               max_iterations=max_iterations, tolerance=tolerance)
            if (candidate%info /= 0) cycle
            out%aic(ar_order, ma_order) = candidate%aic
            if (candidate%aic < best_aic) then
               best_aic = candidate%aic
               out%selected_ar_order = ar_order
               out%selected_ma_order = ma_order
               out%model = candidate
            end if
         end do
      end do
      if (out%selected_ar_order == 0) out%info = 2
   end function exparma_select

end module exparma_mod
