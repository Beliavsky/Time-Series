! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Markov-switching algorithms translated from the R MSwM package.
module mswm_mod
   !! Gaussian Markov-switching regression fitted by Hamilton filtering and EM.
   use kind_mod, only: dp
   use stats_mod, only: ols_fit, normal_quantile
   use utils_mod, only: lowercase
   use optimization_mod, only: finite_difference_hessian
   use linalg_mod, only: invert_matrix
   use random_mod, only: random_standard_normal, random_uniform
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   type, public :: mswm_filter_t
      !! Conditional moments, state probabilities, and likelihood contributions.
      real(dp), allocatable :: conditional_mean(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: likelihood(:, :)
      real(dp), allocatable :: predicted_probability(:, :)
      real(dp), allocatable :: filtered_probability(:, :)
      real(dp), allocatable :: marginal_log_likelihood(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: info = 0
   end type mswm_filter_t

   type, public :: mswm_fit_t
      !! Gaussian Markov-switching regression estimates and state inference.
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: standard_deviation(:)
      real(dp), allocatable :: transition(:, :)
      real(dp), allocatable :: initial_probability(:)
      real(dp), allocatable :: conditional_mean(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: filtered_probability(:, :)
      real(dp), allocatable :: smoothed_probability(:, :)
      real(dp), allocatable :: expected_transition(:, :, :)
      logical, allocatable :: switching(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      logical :: variance_switching = .true.
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: lag_order = 0
      real(dp), allocatable :: terminal_response(:)
      integer :: info = 0
   end type mswm_fit_t

   type, public :: mswm_glm_fit_t
      !! Generalized Markov-switching regression estimates and state inference.
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: transition(:, :)
      real(dp), allocatable :: initial_probability(:)
      real(dp), allocatable :: conditional_mean(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: filtered_probability(:, :)
      real(dp), allocatable :: smoothed_probability(:, :)
      real(dp), allocatable :: expected_transition(:, :, :)
      logical, allocatable :: switching(:)
      character(len=:), allocatable :: family
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: lag_order = 0
      real(dp), allocatable :: terminal_response(:)
      integer :: info = 0
   end type mswm_glm_fit_t

   type, public :: mswm_ar_data_t
      !! Response and predictor matrix aligned after autoregressive lagging.
      real(dp), allocatable :: response(:)
      real(dp), allocatable :: design(:, :)
      integer :: lag_order = 0
      logical :: includes_intercept = .true.
      integer :: info = 0
   end type mswm_ar_data_t

   type, public :: mswm_inference_t
      !! Hessian covariance and natural-scale standard errors for a Gaussian fit.
      real(dp), allocatable :: parameter(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: coefficient_standard_error(:, :)
      real(dp), allocatable :: deviation_standard_error(:)
      real(dp), allocatable :: transition_standard_error(:, :)
      integer :: info = 0
   end type mswm_inference_t

   type, public :: mswm_intervals_t
      !! Approximate normal confidence intervals for Gaussian switching parameters.
      real(dp), allocatable :: coefficient_lower(:, :)
      real(dp), allocatable :: coefficient_upper(:, :)
      real(dp), allocatable :: deviation_lower(:)
      real(dp), allocatable :: deviation_upper(:)
      real(dp), allocatable :: transition_lower(:, :)
      real(dp), allocatable :: transition_upper(:, :)
      real(dp) :: level = 0.95_dp
      integer :: info = 0
   end type mswm_intervals_t

   type, public :: mswm_multistart_t
      !! All supplied or random-start EM fits and the best likelihood solution.
      type(mswm_fit_t), allocatable :: fits(:)
      integer :: best = 0
      integer :: successful = 0
      integer :: info = 0
   end type mswm_multistart_t

   type, public :: mswm_glm_multistart_t
      !! All supplied or random-start switching GLM fits and the best solution.
      type(mswm_glm_fit_t), allocatable :: fits(:)
      integer :: best = 0
      integer :: successful = 0
      integer :: info = 0
   end type mswm_glm_multistart_t

   interface mswm_residuals
      module procedure mswm_gaussian_residuals
      module procedure mswm_glm_residuals
   end interface mswm_residuals

   interface mswm_states
      module procedure mswm_gaussian_states
      module procedure mswm_glm_states
   end interface mswm_states

   interface display
      module procedure display_mswm_filter
      module procedure display_mswm_fit
      module procedure display_mswm_glm_fit
      module procedure display_mswm_ar_data
      module procedure display_mswm_inference
      module procedure display_mswm_intervals
      module procedure display_mswm_multistart
      module procedure display_mswm_glm_multistart
   end interface display

   public :: mswm_gaussian_filter, mswm_gaussian_fit
   public :: mswm_glm_filter, mswm_glm_fit
   public :: mswm_ar_data, mswm_gaussian_ar_fit, mswm_glm_ar_fit
   public :: mswm_residuals, mswm_states
   public :: mswm_gaussian_inference, mswm_gaussian_intervals
   public :: mswm_glm_inference, mswm_glm_intervals
   public :: mswm_gaussian_multistart, mswm_gaussian_random_starts
   public :: mswm_glm_multistart, mswm_glm_random_starts
   public :: display, display_mswm_filter, display_mswm_fit
   public :: display_mswm_glm_fit, display_mswm_ar_data
   public :: display_mswm_inference, display_mswm_intervals
   public :: display_mswm_multistart, display_mswm_glm_multistart

contains

   subroutine display_mswm_filter(model, unit, print_obs)
      !! Display a switching-filter summary and optionally observation-level results.
      type(mswm_filter_t), intent(in) :: model !! Switching-filter result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print observation-level arrays.
      integer :: destination, observation
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Markov-switching filter'
      write(destination, '(a, i0)') 'Status: ', model%info
      if (model%info /= 0) return
      if (.not. allocated(model%conditional_mean)) then
         write(destination, '(a)') 'Filter results are not allocated.'
         return
      end if
      write(destination, '(a, i0)') 'Observations: ', size(model%conditional_mean, 1)
      write(destination, '(a, i0)') 'Regimes: ', size(model%conditional_mean, 2)
      write(destination, '(a, es14.6)') 'Log likelihood: ', model%log_likelihood
      if (show_observations) then
         write(destination, '(a)') 'Conditional means and filtered probabilities:'
         do observation = 1, size(model%conditional_mean, 1)
            write(destination, '(i8, 1x, *(es14.6, 1x))') observation, &
               model%conditional_mean(observation, :), &
               model%filtered_probability(observation, :)
         end do
      end if
   end subroutine display_mswm_filter

   subroutine display_mswm_fit(model, unit, print_obs)
      !! Display a Gaussian Markov-switching fit and optional observation results.
      type(mswm_fit_t), intent(in) :: model !! Gaussian switching fit to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print observation-level arrays.
      integer :: destination, regime, observation
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Gaussian Markov-switching regression'
      write(destination, '(a, i0)') 'Status: ', model%info
      if (model%info /= 0) return
      if (.not. allocated(model%coefficients)) then
         write(destination, '(a)') 'Fit results are not allocated.'
         return
      end if
      write(destination, '(a, i0)') 'Regimes: ', size(model%coefficients, 1)
      write(destination, '(a, i0)') 'Predictors: ', size(model%coefficients, 2)
      write(destination, '(a, i0)') 'Lag order: ', model%lag_order
      write(destination, '(a, l1)') 'Converged: ', model%converged
      write(destination, '(a, i0)') 'Iterations: ', model%iterations
      write(destination, '(a, es14.6)') 'Log likelihood: ', model%log_likelihood
      write(destination, '(a, es14.6)') 'AIC: ', model%aic
      write(destination, '(a, es14.6)') 'BIC: ', model%bic
      write(destination, '(a)') 'Coefficients by regime:'
      do regime = 1, size(model%coefficients, 1)
         write(destination, '(*(es14.6, 1x))') model%coefficients(regime, :)
      end do
      if (allocated(model%standard_deviation)) then
         write(destination, '(a)') 'Standard deviations:'
         write(destination, '(*(es14.6, 1x))') model%standard_deviation
      end if
      if (allocated(model%transition)) then
         write(destination, '(a)') 'Transition matrix:'
         do regime = 1, size(model%transition, 1)
            write(destination, '(*(es14.6, 1x))') model%transition(regime, :)
         end do
      end if
      if (show_observations .and. allocated(model%conditional_mean) .and. &
         allocated(model%smoothed_probability)) then
         write(destination, '(a)') 'Conditional means and smoothed probabilities:'
         do observation = 1, size(model%conditional_mean, 1)
            write(destination, '(i8, 1x, *(es14.6, 1x))') observation, &
               model%conditional_mean(observation, :), &
               model%smoothed_probability(observation, :)
         end do
      end if
   end subroutine display_mswm_fit

   subroutine display_mswm_glm_fit(model, unit, print_obs)
      !! Display a generalized Markov-switching fit and optional observation results.
      type(mswm_glm_fit_t), intent(in) :: model !! Generalized switching fit to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print observation-level arrays.
      integer :: destination, regime, observation
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Generalized Markov-switching regression'
      write(destination, '(a, i0)') 'Status: ', model%info
      if (model%info /= 0) return
      if (.not. allocated(model%coefficients)) then
         write(destination, '(a)') 'Fit results are not allocated.'
         return
      end if
      if (allocated(model%family)) write(destination, '(a, a)') &
         'Family: ', model%family
      write(destination, '(a, i0)') 'Regimes: ', size(model%coefficients, 1)
      write(destination, '(a, i0)') 'Predictors: ', size(model%coefficients, 2)
      write(destination, '(a, i0)') 'Lag order: ', model%lag_order
      write(destination, '(a, l1)') 'Converged: ', model%converged
      write(destination, '(a, i0)') 'Iterations: ', model%iterations
      write(destination, '(a, es14.6)') 'Log likelihood: ', model%log_likelihood
      write(destination, '(a, es14.6)') 'AIC: ', model%aic
      write(destination, '(a, es14.6)') 'BIC: ', model%bic
      write(destination, '(a)') 'Coefficients by regime:'
      do regime = 1, size(model%coefficients, 1)
         write(destination, '(*(es14.6, 1x))') model%coefficients(regime, :)
      end do
      if (allocated(model%transition)) then
         write(destination, '(a)') 'Transition matrix:'
         do regime = 1, size(model%transition, 1)
            write(destination, '(*(es14.6, 1x))') model%transition(regime, :)
         end do
      end if
      if (show_observations .and. allocated(model%conditional_mean) .and. &
         allocated(model%smoothed_probability)) then
         write(destination, '(a)') 'Conditional means and smoothed probabilities:'
         do observation = 1, size(model%conditional_mean, 1)
            write(destination, '(i8, 1x, *(es14.6, 1x))') observation, &
               model%conditional_mean(observation, :), &
               model%smoothed_probability(observation, :)
         end do
      end if
   end subroutine display_mswm_glm_fit

   subroutine display_mswm_ar_data(data, unit, print_obs)
      !! Display lagged-data dimensions and optionally the aligned observations.
      type(mswm_ar_data_t), intent(in) :: data !! Lagged data to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print aligned data rows.
      integer :: destination, observation
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Markov-switching autoregressive data'
      write(destination, '(a, i0)') 'Status: ', data%info
      if (data%info /= 0) return
      if (.not. allocated(data%response) .or. .not. allocated(data%design)) then
         write(destination, '(a)') 'Lagged data are not allocated.'
         return
      end if
      write(destination, '(a, i0)') 'Observations: ', size(data%response)
      write(destination, '(a, i0)') 'Predictors: ', size(data%design, 2)
      write(destination, '(a, i0)') 'Lag order: ', data%lag_order
      write(destination, '(a, l1)') 'Includes intercept: ', data%includes_intercept
      if (show_observations) then
         write(destination, '(a)') 'Response and design rows:'
         do observation = 1, size(data%response)
            write(destination, '(i8, 1x, *(es14.6, 1x))') observation, &
               data%response(observation), data%design(observation, :)
         end do
      end if
   end subroutine display_mswm_ar_data

   subroutine display_mswm_inference(inference, unit)
      !! Display standard errors from switching-model Hessian inference.
      type(mswm_inference_t), intent(in) :: inference !! Inference result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination, row

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Markov-switching parameter inference'
      write(destination, '(a, i0)') 'Status: ', inference%info
      if (inference%info /= 0) return
      if (allocated(inference%parameter)) write(destination, '(a, i0)') &
         'Unconstrained parameters: ', size(inference%parameter)
      if (allocated(inference%coefficient_standard_error)) then
         write(destination, '(a)') 'Coefficient standard errors by regime:'
         do row = 1, size(inference%coefficient_standard_error, 1)
            write(destination, '(*(es14.6, 1x))') &
               inference%coefficient_standard_error(row, :)
         end do
      end if
      if (allocated(inference%deviation_standard_error)) then
         if (size(inference%deviation_standard_error) > 0) then
            write(destination, '(a)') 'Standard-deviation standard errors:'
            write(destination, '(*(es14.6, 1x))') &
               inference%deviation_standard_error
         end if
      end if
      if (allocated(inference%transition_standard_error)) then
         write(destination, '(a)') 'Transition standard errors:'
         do row = 1, size(inference%transition_standard_error, 1)
            write(destination, '(*(es14.6, 1x))') &
               inference%transition_standard_error(row, :)
         end do
      end if
   end subroutine display_mswm_inference

   subroutine display_mswm_intervals(intervals, unit)
      !! Display approximate switching-model confidence intervals.
      type(mswm_intervals_t), intent(in) :: intervals !! Interval result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination, row

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Markov-switching confidence intervals'
      write(destination, '(a, i0)') 'Status: ', intervals%info
      if (intervals%info /= 0) return
      write(destination, '(a, f8.4)') 'Level: ', intervals%level
      if (allocated(intervals%coefficient_lower)) then
         write(destination, '(a)') 'Coefficient lower and upper limits by regime:'
         do row = 1, size(intervals%coefficient_lower, 1)
            write(destination, '(*(es14.6, 1x))') &
               intervals%coefficient_lower(row, :), &
               intervals%coefficient_upper(row, :)
         end do
      end if
      if (allocated(intervals%deviation_lower)) then
         if (size(intervals%deviation_lower) > 0) then
            write(destination, '(a)') 'Standard-deviation lower and upper limits:'
            write(destination, '(*(es14.6, 1x))') intervals%deviation_lower
            write(destination, '(*(es14.6, 1x))') intervals%deviation_upper
         end if
      end if
      if (allocated(intervals%transition_lower)) then
         write(destination, '(a)') 'Transition lower and upper limits by row:'
         do row = 1, size(intervals%transition_lower, 1)
            write(destination, '(*(es14.6, 1x))') &
               intervals%transition_lower(row, :), &
               intervals%transition_upper(row, :)
         end do
      end if
   end subroutine display_mswm_intervals

   subroutine display_mswm_multistart(result, unit)
      !! Display Gaussian multistart likelihoods and the selected fit.
      type(mswm_multistart_t), intent(in) :: result !! Multistart result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination, start

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Gaussian Markov-switching multistart fit'
      write(destination, '(a, i0)') 'Status: ', result%info
      write(destination, '(a, i0)') 'Successful starts: ', result%successful
      write(destination, '(a, i0)') 'Selected start: ', result%best
      if (.not. allocated(result%fits)) return
      do start = 1, size(result%fits)
         write(destination, '(a, i0, a, i0, a, es14.6)') 'Start ', start, &
            ': status=', result%fits(start)%info, ', log likelihood=', &
            result%fits(start)%log_likelihood
      end do
   end subroutine display_mswm_multistart

   subroutine display_mswm_glm_multistart(result, unit)
      !! Display generalized multistart likelihoods and the selected fit.
      type(mswm_glm_multistart_t), intent(in) :: result !! Multistart result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination, start

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Generalized Markov-switching multistart fit'
      write(destination, '(a, i0)') 'Status: ', result%info
      write(destination, '(a, i0)') 'Successful starts: ', result%successful
      write(destination, '(a, i0)') 'Selected start: ', result%best
      if (.not. allocated(result%fits)) return
      do start = 1, size(result%fits)
         write(destination, '(a, i0, a, i0, a, es14.6)') 'Start ', start, &
            ': status=', result%fits(start)%info, ', log likelihood=', &
            result%fits(start)%log_likelihood
      end do
   end subroutine display_mswm_glm_multistart

   pure function mswm_ar_data(response, lag_order, exogenous, include_intercept) &
      result(out)
      !! Align a response, optional predictors, and autoregressive response lags.
      real(dp), intent(in) :: response(:) !! Complete response series.
      integer, intent(in) :: lag_order !! Number of response lags.
      real(dp), intent(in), optional :: exogenous(:, :) !! Complete exogenous predictors.
      logical, intent(in), optional :: include_intercept !! Whether to prepend an intercept.
      type(mswm_ar_data_t) :: out
      logical :: selected_intercept
      integer :: observations, retained, exogenous_count, columns
      integer :: row, column, lag

      observations = size(response)
      selected_intercept = .true.
      if (present(include_intercept)) selected_intercept = include_intercept
      if (lag_order < 0 .or. observations <= lag_order .or. &
         .not. all(ieee_is_finite(response))) then
         out%info = 1
         return
      end if
      exogenous_count = 0
      if (present(exogenous)) then
         if (size(exogenous, 1) /= observations .or. &
            .not. all(ieee_is_finite(exogenous))) then
            out%info = 1
            return
         end if
         exogenous_count = size(exogenous, 2)
      end if
      columns = merge(1, 0, selected_intercept) + exogenous_count + lag_order
      if (columns < 1) then
         out%info = 1
         return
      end if
      retained = observations - lag_order
      allocate(out%response(retained), out%design(retained, columns))
      out%response = response(lag_order + 1:observations)
      column = 0
      if (selected_intercept) then
         column = 1
         out%design(:, column) = 1.0_dp
      end if
      if (present(exogenous)) then
         if (exogenous_count > 0) then
            out%design(:, column + 1:column + exogenous_count) = &
               exogenous(lag_order + 1:observations, :)
            column = column + exogenous_count
         end if
      end if
      do lag = 1, lag_order
         column = column + 1
         do row = 1, retained
            out%design(row, column) = response(lag_order + row - lag)
         end do
      end do
      out%lag_order = lag_order
      out%includes_intercept = selected_intercept
   end function mswm_ar_data

   pure function mswm_gaussian_ar_fit(response, lag_order, regimes, switching, &
      variance_switching, exogenous, include_intercept, initial_coefficients, &
      initial_standard_deviation, initial_transition, max_iterations, tolerance) &
      result(out)
      !! Fit Gaussian Markov-switching regression with automatic response lags.
      real(dp), intent(in) :: response(:) !! Complete response series.
      integer, intent(in) :: lag_order !! Number of response lags.
      integer, intent(in) :: regimes !! Number of latent regimes.
      logical, intent(in) :: switching(:) !! Flags for intercept, exogenous, then lag coefficients.
      logical, intent(in) :: variance_switching !! Whether innovation variance switches.
      real(dp), intent(in), optional :: exogenous(:, :) !! Complete exogenous predictors.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      real(dp), intent(in), optional :: initial_coefficients(:, :) !! Initial coefficients.
      real(dp), intent(in), optional :: initial_standard_deviation(:) !! Initial scales.
      real(dp), intent(in), optional :: initial_transition(:, :) !! Initial transition matrix.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative convergence tolerance.
      type(mswm_fit_t) :: out
      type(mswm_ar_data_t) :: data

      data = mswm_ar_data(response, lag_order, exogenous, include_intercept)
      if (data%info /= 0 .or. size(switching) /= size(data%design, 2)) then
         out%info = 1
         return
      end if
      out = mswm_gaussian_fit(data%response, data%design, regimes, switching, &
         variance_switching, initial_coefficients, initial_standard_deviation, &
         initial_transition, max_iterations, tolerance)
      if (out%info /= 0) return
      out%lag_order = lag_order
      if (lag_order > 0) then
         out%terminal_response = response(size(response) - lag_order + 1:)
      else
         allocate(out%terminal_response(0))
      end if
   end function mswm_gaussian_ar_fit

   pure function mswm_glm_ar_fit(response, lag_order, regimes, switching, family, &
      exogenous, include_intercept, initial_coefficients, initial_transition, &
      max_iterations, tolerance, irls_iterations) result(out)
      !! Fit a switching GLM with automatically constructed response lags.
      real(dp), intent(in) :: response(:) !! Complete response series.
      integer, intent(in) :: lag_order !! Number of response lags.
      integer, intent(in) :: regimes !! Number of latent regimes.
      logical, intent(in) :: switching(:) !! Flags for intercept, exogenous, then lag coefficients.
      character(len=*), intent(in) :: family !! Poisson, binomial, or Gamma family and link.
      real(dp), intent(in), optional :: exogenous(:, :) !! Complete exogenous predictors.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      real(dp), intent(in), optional :: initial_coefficients(:, :) !! Initial coefficients.
      real(dp), intent(in), optional :: initial_transition(:, :) !! Initial transition matrix.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative convergence tolerance.
      integer, intent(in), optional :: irls_iterations !! Maximum IRLS M-step iterations.
      type(mswm_glm_fit_t) :: out
      type(mswm_ar_data_t) :: data

      data = mswm_ar_data(response, lag_order, exogenous, include_intercept)
      if (data%info /= 0 .or. size(switching) /= size(data%design, 2)) then
         out%info = 1
         return
      end if
      out = mswm_glm_fit(data%response, data%design, regimes, switching, family, &
         initial_coefficients, initial_transition, max_iterations, tolerance, &
         irls_iterations)
      if (out%info /= 0) return
      out%lag_order = lag_order
      if (lag_order > 0) then
         out%terminal_response = response(size(response) - lag_order + 1:)
      else
         allocate(out%terminal_response(0))
      end if
   end function mswm_glm_ar_fit

   pure function mswm_gaussian_filter(response, design, coefficients, &
      standard_deviation, transition, initial_probability) result(out)
      !! Apply the scaled Hamilton filter to a Gaussian switching regression.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:, :) !! Regime-by-predictor coefficients.
      real(dp), intent(in) :: standard_deviation(:) !! Regime error standard deviations.
      real(dp), intent(in) :: transition(:, :) !! Current-state by next-state probabilities.
      real(dp), intent(in) :: initial_probability(:) !! Pre-sample state probabilities.
      type(mswm_filter_t) :: out
      real(dp), allocatable :: log_component(:), scaled(:)
      real(dp) :: maximum_log, total, pi_value
      integer :: observations, predictors, regimes, time, regime

      observations = size(response)
      predictors = size(design, 2)
      regimes = size(coefficients, 1)
      if (observations < 1 .or. predictors < 1 .or. &
         size(design, 1) /= observations .or. regimes < 1 .or. &
         size(coefficients, 2) /= predictors .or. &
         size(standard_deviation) /= regimes .or. &
         any(shape(transition) /= [regimes, regimes]) .or. &
         size(initial_probability) /= regimes .or. &
         any(standard_deviation <= 0.0_dp) .or. &
         any(transition < 0.0_dp) .or. any(initial_probability < 0.0_dp) .or. &
         maxval(abs(sum(transition, dim=2) - 1.0_dp)) > 1.0e-8_dp .or. &
         abs(sum(initial_probability) - 1.0_dp) > 1.0e-8_dp .or. &
         .not. all(ieee_is_finite(response)) .or. &
         .not. all(ieee_is_finite(design)) .or. &
         .not. all(ieee_is_finite(coefficients)) .or. &
         .not. all(ieee_is_finite(standard_deviation)) .or. &
         .not. all(ieee_is_finite(transition))) then
         out%info = 1
         return
      end if

      allocate(out%conditional_mean(observations, regimes))
      allocate(out%residuals(observations, regimes))
      allocate(out%likelihood(observations, regimes))
      allocate(out%predicted_probability(observations, regimes))
      allocate(out%filtered_probability(observations, regimes))
      allocate(out%marginal_log_likelihood(observations))
      allocate(log_component(regimes), scaled(regimes))
      out%conditional_mean = matmul(design, transpose(coefficients))
      out%residuals = spread(response, 2, regimes) - out%conditional_mean
      pi_value = acos(-1.0_dp)
      do regime = 1, regimes
         out%likelihood(:, regime) = exp( &
            -0.5_dp*(out%residuals(:, regime)/standard_deviation(regime))**2)/ &
            (sqrt(2.0_dp*pi_value)*standard_deviation(regime))
      end do

      out%predicted_probability(1, :) = &
         matmul(initial_probability, transition)
      do time = 1, observations
         if (time > 1) out%predicted_probability(time, :) = &
            matmul(out%filtered_probability(time - 1, :), transition)
         do regime = 1, regimes
            log_component(regime) = log(max( &
               out%predicted_probability(time, regime), tiny(1.0_dp))) - &
               log(standard_deviation(regime)) - 0.5_dp*log(2.0_dp*pi_value) - &
               0.5_dp*(out%residuals(time, regime)/standard_deviation(regime))**2
         end do
         maximum_log = maxval(log_component)
         scaled = exp(log_component - maximum_log)
         total = sum(scaled)
         if (total <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
         out%filtered_probability(time, :) = scaled/total
         out%marginal_log_likelihood(time) = maximum_log + log(total)
      end do
      out%log_likelihood = sum(out%marginal_log_likelihood)
   end function mswm_gaussian_filter

   pure function mswm_glm_filter(response, design, coefficients, transition, &
      initial_probability, family) result(out)
      !! Apply the scaled Hamilton filter to a generalized switching regression.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:, :) !! Regime-by-predictor coefficients.
      real(dp), intent(in) :: transition(:, :) !! Current-state by next-state probabilities.
      real(dp), intent(in) :: initial_probability(:) !! Pre-sample state probabilities.
      character(len=*), intent(in) :: family !! Poisson, binomial, or Gamma family and link.
      type(mswm_filter_t) :: out
      real(dp), allocatable :: log_component(:), scaled(:)
      real(dp) :: maximum_log, total
      integer :: observations, predictors, regimes, time, regime
      character(len=:), allocatable :: selected_family

      observations = size(response)
      predictors = size(design, 2)
      regimes = size(coefficients, 1)
      selected_family = trim(adjustl(lowercase(family)))
      if (observations < 1 .or. predictors < 1 .or. &
         size(design, 1) /= observations .or. regimes < 1 .or. &
         size(coefficients, 2) /= predictors .or. &
         any(shape(transition) /= [regimes, regimes]) .or. &
         size(initial_probability) /= regimes .or. &
         any(transition < 0.0_dp) .or. any(initial_probability < 0.0_dp) .or. &
         maxval(abs(sum(transition, dim=2) - 1.0_dp)) > 1.0e-8_dp .or. &
         abs(sum(initial_probability) - 1.0_dp) > 1.0e-8_dp .or. &
         .not. valid_glm_family(selected_family) .or. &
         .not. valid_glm_response(response, selected_family) .or. &
         .not. all(ieee_is_finite(design)) .or. &
         .not. all(ieee_is_finite(coefficients))) then
         out%info = 1
         return
      end if

      allocate(out%conditional_mean(observations, regimes))
      allocate(out%residuals(observations, regimes))
      allocate(out%likelihood(observations, regimes))
      allocate(out%predicted_probability(observations, regimes))
      allocate(out%filtered_probability(observations, regimes))
      allocate(out%marginal_log_likelihood(observations))
      allocate(log_component(regimes), scaled(regimes))
      do regime = 1, regimes
         do time = 1, observations
            out%conditional_mean(time, regime) = glm_mean( &
               dot_product(design(time, :), coefficients(regime, :)), &
               selected_family)
            out%residuals(time, regime) = response(time) - &
               out%conditional_mean(time, regime)
            log_component(regime) = glm_log_density(response(time), &
               out%conditional_mean(time, regime), selected_family)
            out%likelihood(time, regime) = exp(max(log_component(regime), &
               log(tiny(1.0_dp))))
         end do
      end do

      out%predicted_probability(1, :) = &
         matmul(initial_probability, transition)
      do time = 1, observations
         if (time > 1) out%predicted_probability(time, :) = &
            matmul(out%filtered_probability(time - 1, :), transition)
         do regime = 1, regimes
            log_component(regime) = log(max( &
               out%predicted_probability(time, regime), tiny(1.0_dp))) + &
               glm_log_density(response(time), &
               out%conditional_mean(time, regime), selected_family)
         end do
         maximum_log = maxval(log_component)
         scaled = exp(log_component - maximum_log)
         total = sum(scaled)
         if (total <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
         out%filtered_probability(time, :) = scaled/total
         out%marginal_log_likelihood(time) = maximum_log + log(total)
      end do
      out%log_likelihood = sum(out%marginal_log_likelihood)
   end function mswm_glm_filter

   pure function mswm_glm_fit(response, design, regimes, switching, family, &
      initial_coefficients, initial_transition, max_iterations, tolerance, &
      irls_iterations) result(out)
      !! Fit a generalized Markov-switching regression by EM and weighted IRLS.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      integer, intent(in) :: regimes !! Number of latent regimes.
      logical, intent(in) :: switching(:) !! Predictor coefficient switching flags.
      character(len=*), intent(in) :: family !! Poisson, binomial, or Gamma family and link.
      real(dp), intent(in), optional :: initial_coefficients(:, :) !! Initial coefficients.
      real(dp), intent(in), optional :: initial_transition(:, :) !! Initial transition matrix.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative convergence tolerance.
      integer, intent(in), optional :: irls_iterations !! Maximum IRLS steps per M-step.
      type(mswm_glm_fit_t) :: out
      type(mswm_filter_t) :: filtered
      real(dp), allocatable :: coefficients(:, :), transition(:, :)
      real(dp), allocatable :: initial_probability(:), smoothed(:, :)
      real(dp), allocatable :: expected_transition(:, :, :)
      real(dp), allocatable :: old_coefficients(:, :), probability(:, :)
      logical, allocatable :: common_switching(:)
      real(dp) :: selected_tolerance, old_log_likelihood
      real(dp) :: coefficient_change, likelihood_change, center
      integer :: selected_iterations, selected_irls, predictors, observations
      integer :: regime, iteration, status, parameter_count
      character(len=:), allocatable :: selected_family

      observations = size(response)
      predictors = size(design, 2)
      selected_family = trim(adjustl(lowercase(family)))
      selected_iterations = 100
      if (present(max_iterations)) selected_iterations = max_iterations
      selected_irls = 25
      if (present(irls_iterations)) selected_irls = irls_iterations
      selected_tolerance = 1.0e-8_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (observations < 2 .or. predictors < 1 .or. &
         size(design, 1) /= observations .or. regimes < 1 .or. &
         size(switching) /= predictors .or. selected_iterations < 1 .or. &
         selected_irls < 1 .or. selected_tolerance <= 0.0_dp .or. &
         .not. valid_glm_family(selected_family) .or. &
         .not. valid_glm_response(response, selected_family) .or. &
         .not. all(ieee_is_finite(design))) then
         out%info = 1
         return
      end if

      allocate(coefficients(regimes, predictors))
      allocate(transition(regimes, regimes), initial_probability(regimes))
      allocate(probability(observations, 1))
      allocate(common_switching(predictors))
      probability = 1.0_dp
      common_switching = .false.
      coefficients = 0.0_dp
      call maximize_glm(response, design, probability, common_switching, &
         selected_family, selected_irls, coefficients(:1, :), status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      coefficients = spread(coefficients(1, :), 1, regimes)
      do regime = 1, regimes
         center = real(regime, dp) - 0.5_dp*real(regimes + 1, dp)
         if (any(switching)) then
            coefficients(regime, minloc(merge(0, 1, switching), dim=1)) = &
               coefficients(regime, minloc(merge(0, 1, switching), dim=1)) + &
               0.4_dp*center
         end if
      end do
      if (regimes == 1) then
         transition = 1.0_dp
      else
         transition = 0.10_dp/real(regimes - 1, dp)
         do regime = 1, regimes
            transition(regime, regime) = 0.90_dp
         end do
      end if
      initial_probability = 1.0_dp/real(regimes, dp)
      if (present(initial_coefficients)) then
         if (any(shape(initial_coefficients) /= [regimes, predictors])) then
            out%info = 3
            return
         end if
         coefficients = initial_coefficients
      end if
      if (present(initial_transition)) then
         if (any(shape(initial_transition) /= [regimes, regimes]) .or. &
            any(initial_transition < 0.0_dp) .or. &
            maxval(abs(sum(initial_transition, dim=2) - 1.0_dp)) > 1.0e-8_dp) then
            out%info = 3
            return
         end if
         transition = initial_transition
      end if

      old_log_likelihood = -huge(1.0_dp)
      do iteration = 1, selected_iterations
         filtered = mswm_glm_filter(response, design, coefficients, transition, &
            initial_probability, selected_family)
         if (filtered%info /= 0) then
            out%info = 10 + filtered%info
            return
         end if
         call smooth_probabilities(filtered%filtered_probability, &
            filtered%predicted_probability, transition, smoothed, &
            expected_transition)
         old_coefficients = coefficients
         call maximize_glm(response, design, smoothed, switching, &
            selected_family, selected_irls, coefficients, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         call maximize_transition(smoothed, expected_transition, transition)
         initial_probability = smoothed(1, :)
         likelihood_change = abs(filtered%log_likelihood - old_log_likelihood)/ &
            (0.1_dp + abs(filtered%log_likelihood))
         coefficient_change = maxval(abs(coefficients - old_coefficients))/ &
            (0.1_dp + maxval(abs(coefficients)))
         out%iterations = iteration
         if (iteration > 1 .and. likelihood_change < selected_tolerance .and. &
            coefficient_change < selected_tolerance) then
            out%converged = .true.
            exit
         end if
         old_log_likelihood = filtered%log_likelihood
      end do

      filtered = mswm_glm_filter(response, design, coefficients, transition, &
         initial_probability, selected_family)
      if (filtered%info /= 0) then
         out%info = 30 + filtered%info
         return
      end if
      call smooth_probabilities(filtered%filtered_probability, &
         filtered%predicted_probability, transition, smoothed, &
         expected_transition)
      out%coefficients = coefficients
      out%transition = transition
      out%initial_probability = initial_probability
      out%conditional_mean = filtered%conditional_mean
      out%residuals = filtered%residuals
      out%filtered_probability = filtered%filtered_probability
      out%smoothed_probability = smoothed
      out%expected_transition = expected_transition
      out%switching = switching
      out%family = selected_family
      out%log_likelihood = filtered%log_likelihood
      parameter_count = count(.not. switching) + regimes*count(switching) + &
         regimes*(regimes - 1)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameter_count, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(observations, dp))* &
         real(parameter_count, dp)
   end function mswm_glm_fit

   pure function mswm_gaussian_fit(response, design, regimes, switching, &
      variance_switching, initial_coefficients, initial_standard_deviation, &
      initial_transition, max_iterations, tolerance) result(out)
      !! Fit a Gaussian Markov-switching regression by the MSwM EM algorithm.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      integer, intent(in) :: regimes !! Number of latent regimes.
      logical, intent(in) :: switching(:) !! Predictor coefficient switching flags.
      logical, intent(in), optional :: variance_switching !! Whether variance switches.
      real(dp), intent(in), optional :: initial_coefficients(:, :) !! Initial coefficients.
      real(dp), intent(in), optional :: initial_standard_deviation(:) !! Initial error scales.
      real(dp), intent(in), optional :: initial_transition(:, :) !! Initial transition matrix.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative convergence tolerance.
      type(mswm_fit_t) :: out
      type(mswm_filter_t) :: filtered
      real(dp), allocatable :: coefficients(:, :), standard_deviation(:)
      real(dp), allocatable :: transition(:, :), initial_probability(:)
      real(dp), allocatable :: smoothed(:, :), expected_transition(:, :, :)
      real(dp), allocatable :: old_coefficients(:, :), common_beta(:)
      real(dp), allocatable :: standard_error(:), common_residual(:)
      real(dp) :: common_rss, selected_tolerance, old_log_likelihood
      real(dp) :: coefficient_change, likelihood_change, scale, center
      integer :: selected_iterations, predictors, observations, regime, status
      integer :: iteration, parameter_count
      logical :: selected_variance_switching

      observations = size(response)
      predictors = size(design, 2)
      selected_variance_switching = .true.
      if (present(variance_switching)) &
         selected_variance_switching = variance_switching
      selected_iterations = 100
      if (present(max_iterations)) selected_iterations = max_iterations
      selected_tolerance = 1.0e-8_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (observations < 2 .or. predictors < 1 .or. &
         size(design, 1) /= observations .or. regimes < 1 .or. &
         size(switching) /= predictors .or. selected_iterations < 1 .or. &
         selected_tolerance <= 0.0_dp .or. &
         .not. all(ieee_is_finite(response)) .or. &
         .not. all(ieee_is_finite(design))) then
         out%info = 1
         return
      end if

      allocate(coefficients(regimes, predictors))
      allocate(standard_deviation(regimes), transition(regimes, regimes))
      allocate(initial_probability(regimes))
      call ols_fit(design, response, common_beta, standard_error, &
         common_residual, common_rss, status)
      if (status /= 0 .or. common_rss <= 0.0_dp) then
         out%info = 2
         return
      end if
      coefficients = spread(common_beta, 1, regimes)
      scale = sqrt(common_rss/real(max(1, observations - predictors), dp))
      do regime = 1, regimes
         center = real(regime, dp) - 0.5_dp*real(regimes + 1, dp)
         if (any(switching)) then
            coefficients(regime, minloc(merge(0, 1, switching), dim=1)) = &
               coefficients(regime, minloc(merge(0, 1, switching), dim=1)) + &
               0.5_dp*scale*center
         end if
         standard_deviation(regime) = scale*(0.8_dp + &
            0.4_dp*real(regime - 1, dp)/real(max(1, regimes - 1), dp))
      end do
      if (.not. selected_variance_switching) standard_deviation = scale
      if (regimes == 1) then
         transition = 1.0_dp
      else
         transition = (1.0_dp - 0.90_dp)/real(regimes - 1, dp)
         do regime = 1, regimes
            transition(regime, regime) = 0.90_dp
         end do
      end if
      initial_probability = 1.0_dp/real(regimes, dp)

      if (present(initial_coefficients)) then
         if (any(shape(initial_coefficients) /= [regimes, predictors])) then
            out%info = 3
            return
         end if
         coefficients = initial_coefficients
      end if
      if (present(initial_standard_deviation)) then
         if (size(initial_standard_deviation) /= regimes .or. &
            any(initial_standard_deviation <= 0.0_dp)) then
            out%info = 3
            return
         end if
         standard_deviation = initial_standard_deviation
      end if
      if (present(initial_transition)) then
         if (any(shape(initial_transition) /= [regimes, regimes]) .or. &
            any(initial_transition < 0.0_dp) .or. &
            maxval(abs(sum(initial_transition, dim=2) - 1.0_dp)) > 1.0e-8_dp) then
            out%info = 3
            return
         end if
         transition = initial_transition
      end if

      old_log_likelihood = -huge(1.0_dp)
      do iteration = 1, selected_iterations
         filtered = mswm_gaussian_filter(response, design, coefficients, &
            standard_deviation, transition, initial_probability)
         if (filtered%info /= 0) then
            out%info = 10 + filtered%info
            return
         end if
         call smooth_probabilities(filtered%filtered_probability, &
            filtered%predicted_probability, transition, smoothed, &
            expected_transition)
         old_coefficients = coefficients
         call maximize_gaussian(response, design, smoothed, switching, &
            selected_variance_switching, coefficients, standard_deviation, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         call maximize_transition(smoothed, expected_transition, transition)
         initial_probability = smoothed(1, :)

         likelihood_change = abs(filtered%log_likelihood - old_log_likelihood)/ &
            (0.1_dp + abs(filtered%log_likelihood))
         coefficient_change = maxval(abs(coefficients - old_coefficients))/ &
            (0.1_dp + maxval(abs(coefficients)))
         out%iterations = iteration
         if (iteration > 1 .and. likelihood_change < selected_tolerance .and. &
            coefficient_change < selected_tolerance) then
            out%converged = .true.
            exit
         end if
         old_log_likelihood = filtered%log_likelihood
      end do

      filtered = mswm_gaussian_filter(response, design, coefficients, &
         standard_deviation, transition, initial_probability)
      if (filtered%info /= 0) then
         out%info = 30 + filtered%info
         return
      end if
      call smooth_probabilities(filtered%filtered_probability, &
         filtered%predicted_probability, transition, smoothed, &
         expected_transition)
      out%coefficients = coefficients
      out%standard_deviation = standard_deviation
      out%transition = transition
      out%initial_probability = initial_probability
      out%conditional_mean = filtered%conditional_mean
      out%residuals = filtered%residuals
      out%filtered_probability = filtered%filtered_probability
      out%smoothed_probability = smoothed
      out%expected_transition = expected_transition
      out%switching = switching
      out%variance_switching = selected_variance_switching
      out%log_likelihood = filtered%log_likelihood
      parameter_count = count(.not. switching) + regimes*count(switching) + &
         merge(regimes, 1, selected_variance_switching) + regimes*(regimes - 1)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameter_count, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(observations, dp))* &
         real(parameter_count, dp)
   end function mswm_gaussian_fit

   pure function mswm_gaussian_inference(model, response, design, difference_step) &
      result(out)
      !! Compute finite-difference observed-information inference for a Gaussian fit.
      type(mswm_fit_t), intent(in) :: model !! Fitted Gaussian switching model.
      real(dp), intent(in) :: response(:) !! Regression response used for fitting.
      real(dp), intent(in) :: design(:, :) !! Regression design used for fitting.
      real(dp), intent(in), optional :: difference_step !! Hessian difference step.
      type(mswm_inference_t) :: out
      real(dp), allocatable :: inverse(:, :), gradient(:), block(:, :)
      real(dp) :: step, variance_value
      integer :: regimes, predictors, scales, transitions, position
      integer :: regime, predictor, target, source, status

      if (model%info /= 0 .or. .not. allocated(model%coefficients) .or. &
         size(response) < 1 .or. size(design, 1) /= size(response) .or. &
         size(design, 2) /= size(model%coefficients, 2)) then
         out%info = 1
         return
      end if
      step = 1.0e-4_dp
      if (present(difference_step)) step = difference_step
      if (step <= 0.0_dp) then
         out%info = 1
         return
      end if
      regimes = size(model%coefficients, 1)
      predictors = size(model%coefficients, 2)
      scales = merge(regimes, 1, model%variance_switching)
      transitions = regimes*(regimes - 1)
      call pack_gaussian_parameters(model, out%parameter)
      out%hessian = finite_difference_hessian(objective, out%parameter, step)
      call invert_matrix(out%hessian, inverse, status)
      if (status /= 0 .or. .not. all(ieee_is_finite(inverse))) then
         out%info = 10 + status
         return
      end if
      out%covariance = 0.5_dp*(inverse + transpose(inverse))
      allocate(out%deviation_standard_error(regimes))
      if (model%variance_switching) then
         do regime = 1, regimes
            out%deviation_standard_error(regime) = model%standard_deviation(regime)* &
               sqrt(max(0.0_dp, out%covariance(regime, regime)))
         end do
      else
         out%deviation_standard_error = model%standard_deviation(1)* &
            sqrt(max(0.0_dp, out%covariance(1, 1)))
      end if

      allocate(out%transition_standard_error(regimes, regimes))
      allocate(gradient(size(out%parameter)))
      out%transition_standard_error = 0.0_dp
      position = scales
      do regime = 1, regimes
         allocate(block(regimes, max(0, regimes - 1)))
         block = 0.0_dp
         do target = 1, regimes
            do source = 1, regimes - 1
               if (target == source) then
                  block(target, source) = model%transition(regime, target)* &
                     (1.0_dp - model%transition(regime, source))
               else
                  block(target, source) = -model%transition(regime, target)* &
                     model%transition(regime, source)
               end if
            end do
         end do
         do target = 1, regimes
            gradient = 0.0_dp
            if (regimes > 1) gradient(position + 1:position + regimes - 1) = &
               block(target, :)
            variance_value = dot_product(gradient, &
               matmul(out%covariance, gradient))
            out%transition_standard_error(regime, target) = &
               sqrt(max(0.0_dp, variance_value))
         end do
         position = position + regimes - 1
         deallocate(block)
      end do

      allocate(out%coefficient_standard_error(regimes, predictors))
      position = scales + transitions
      do predictor = 1, predictors
         if (model%switching(predictor)) then
            do regime = 1, regimes
               position = position + 1
               out%coefficient_standard_error(regime, predictor) = &
                  sqrt(max(0.0_dp, out%covariance(position, position)))
            end do
         else
            position = position + 1
            out%coefficient_standard_error(:, predictor) = &
               sqrt(max(0.0_dp, out%covariance(position, position)))
         end if
      end do

   contains

      pure function objective(parameters) result(value)
         !! Evaluate the negative Gaussian switching log likelihood.
         real(dp), intent(in) :: parameters(:) !! Unconstrained model parameters.
         real(dp) :: value
         real(dp), allocatable :: coefficients(:, :), standard_deviation(:)
         real(dp), allocatable :: transition(:, :)
         type(mswm_filter_t) :: filtered

         call unpack_gaussian_parameters(parameters, regimes, model%switching, &
            model%variance_switching, coefficients, standard_deviation, transition)
         filtered = mswm_gaussian_filter(response, design, coefficients, &
            standard_deviation, transition, model%initial_probability)
         if (filtered%info /= 0) then
            value = huge(1.0_dp)
         else
            value = -filtered%log_likelihood
         end if
      end function objective

   end function mswm_gaussian_inference

   pure function mswm_gaussian_intervals(model, inference, level) result(out)
      !! Form approximate normal confidence intervals on the natural parameter scale.
      type(mswm_fit_t), intent(in) :: model !! Fitted Gaussian switching model.
      type(mswm_inference_t), intent(in) :: inference !! Hessian inference result.
      real(dp), intent(in), optional :: level !! Confidence level.
      type(mswm_intervals_t) :: out
      real(dp) :: selected_level, critical

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      if (model%info /= 0 .or. inference%info /= 0 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      critical = normal_quantile(0.5_dp + 0.5_dp*selected_level)
      out%coefficient_lower = model%coefficients - &
         critical*inference%coefficient_standard_error
      out%coefficient_upper = model%coefficients + &
         critical*inference%coefficient_standard_error
      out%deviation_lower = max(0.0_dp, model%standard_deviation - &
         critical*inference%deviation_standard_error)
      out%deviation_upper = model%standard_deviation + &
         critical*inference%deviation_standard_error
      out%transition_lower = max(0.0_dp, model%transition - &
         critical*inference%transition_standard_error)
      out%transition_upper = min(1.0_dp, model%transition + &
         critical*inference%transition_standard_error)
      out%level = selected_level
   end function mswm_gaussian_intervals

   pure function mswm_glm_inference(model, response, design, difference_step) &
      result(out)
      !! Compute finite-difference observed-information inference for a GLM fit.
      type(mswm_glm_fit_t), intent(in) :: model !! Fitted switching GLM.
      real(dp), intent(in) :: response(:) !! Regression response used for fitting.
      real(dp), intent(in) :: design(:, :) !! Regression design used for fitting.
      real(dp), intent(in), optional :: difference_step !! Hessian difference step.
      type(mswm_inference_t) :: out
      real(dp), allocatable :: inverse(:, :), gradient(:), block(:, :)
      real(dp) :: step, variance_value
      integer :: regimes, predictors, transitions, position
      integer :: regime, predictor, target, source, status

      if (model%info /= 0 .or. .not. allocated(model%coefficients) .or. &
         .not. allocated(model%family) .or. size(response) < 1 .or. &
         size(design, 1) /= size(response) .or. &
         size(design, 2) /= size(model%coefficients, 2)) then
         out%info = 1
         return
      end if
      step = 1.0e-4_dp
      if (present(difference_step)) step = difference_step
      if (step <= 0.0_dp) then
         out%info = 1
         return
      end if
      regimes = size(model%coefficients, 1)
      predictors = size(model%coefficients, 2)
      transitions = regimes*(regimes - 1)
      call pack_glm_parameters(model, out%parameter)
      out%hessian = finite_difference_hessian(objective, out%parameter, step)
      call invert_matrix(out%hessian, inverse, status)
      if (status /= 0 .or. .not. all(ieee_is_finite(inverse))) then
         out%info = 10 + status
         return
      end if
      out%covariance = 0.5_dp*(inverse + transpose(inverse))
      allocate(out%deviation_standard_error(0))
      allocate(out%transition_standard_error(regimes, regimes))
      allocate(gradient(size(out%parameter)))
      out%transition_standard_error = 0.0_dp
      position = 0
      do regime = 1, regimes
         allocate(block(regimes, max(0, regimes - 1)))
         block = 0.0_dp
         do target = 1, regimes
            do source = 1, regimes - 1
               if (target == source) then
                  block(target, source) = model%transition(regime, target)* &
                     (1.0_dp - model%transition(regime, source))
               else
                  block(target, source) = -model%transition(regime, target)* &
                     model%transition(regime, source)
               end if
            end do
         end do
         do target = 1, regimes
            gradient = 0.0_dp
            if (regimes > 1) gradient(position + 1:position + regimes - 1) = &
               block(target, :)
            variance_value = dot_product(gradient, &
               matmul(out%covariance, gradient))
            out%transition_standard_error(regime, target) = &
               sqrt(max(0.0_dp, variance_value))
         end do
         position = position + regimes - 1
         deallocate(block)
      end do

      allocate(out%coefficient_standard_error(regimes, predictors))
      position = transitions
      do predictor = 1, predictors
         if (model%switching(predictor)) then
            do regime = 1, regimes
               position = position + 1
               out%coefficient_standard_error(regime, predictor) = &
                  sqrt(max(0.0_dp, out%covariance(position, position)))
            end do
         else
            position = position + 1
            out%coefficient_standard_error(:, predictor) = &
               sqrt(max(0.0_dp, out%covariance(position, position)))
         end if
      end do

   contains

      pure function objective(parameters) result(value)
         !! Evaluate the negative generalized switching log likelihood.
         real(dp), intent(in) :: parameters(:) !! Unconstrained model parameters.
         real(dp) :: value
         real(dp), allocatable :: coefficients(:, :), transition(:, :)
         type(mswm_filter_t) :: filtered

         call unpack_glm_parameters(parameters, regimes, model%switching, &
            coefficients, transition)
         filtered = mswm_glm_filter(response, design, coefficients, transition, &
            model%initial_probability, model%family)
         if (filtered%info /= 0) then
            value = huge(1.0_dp)
         else
            value = -filtered%log_likelihood
         end if
      end function objective

   end function mswm_glm_inference

   pure function mswm_glm_intervals(model, inference, level) result(out)
      !! Form approximate normal intervals for switching GLM parameters.
      type(mswm_glm_fit_t), intent(in) :: model !! Fitted switching GLM.
      type(mswm_inference_t), intent(in) :: inference !! Hessian inference result.
      real(dp), intent(in), optional :: level !! Confidence level.
      type(mswm_intervals_t) :: out
      real(dp) :: selected_level, critical

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      if (model%info /= 0 .or. inference%info /= 0 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      critical = normal_quantile(0.5_dp + 0.5_dp*selected_level)
      out%coefficient_lower = model%coefficients - &
         critical*inference%coefficient_standard_error
      out%coefficient_upper = model%coefficients + &
         critical*inference%coefficient_standard_error
      allocate(out%deviation_lower(0), out%deviation_upper(0))
      out%transition_lower = max(0.0_dp, model%transition - &
         critical*inference%transition_standard_error)
      out%transition_upper = min(1.0_dp, model%transition + &
         critical*inference%transition_standard_error)
      out%level = selected_level
   end function mswm_glm_intervals

   pure function mswm_gaussian_multistart(response, design, switching, &
      variance_switching, initial_coefficients, initial_standard_deviation, &
      initial_transition, max_iterations, tolerance) result(out)
      !! Fit all supplied Gaussian EM starts and select the largest likelihood.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      logical, intent(in) :: switching(:) !! Predictor coefficient switching flags.
      logical, intent(in) :: variance_switching !! Whether variance switches.
      real(dp), intent(in) :: initial_coefficients(:, :, :) !! Regime-predictor-start values.
      real(dp), intent(in) :: initial_standard_deviation(:, :) !! Regime-start scales.
      real(dp), intent(in) :: initial_transition(:, :, :) !! Regime transition starts.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations per start.
      real(dp), intent(in), optional :: tolerance !! Relative convergence tolerance.
      type(mswm_multistart_t) :: out
      integer :: regimes, starts, start

      regimes = size(initial_coefficients, 1)
      starts = size(initial_coefficients, 3)
      if (regimes < 1 .or. starts < 1 .or. &
         size(initial_coefficients, 2) /= size(design, 2) .or. &
         any(shape(initial_standard_deviation) /= [regimes, starts]) .or. &
         any(shape(initial_transition) /= [regimes, regimes, starts])) then
         out%info = 1
         return
      end if
      allocate(out%fits(starts))
      do start = 1, starts
         out%fits(start) = mswm_gaussian_fit(response, design, regimes, switching, &
            variance_switching, initial_coefficients(:, :, start), &
            initial_standard_deviation(:, start), initial_transition(:, :, start), &
            max_iterations, tolerance)
         if (out%fits(start)%info /= 0) cycle
         out%successful = out%successful + 1
         if (out%best == 0) then
            out%best = start
         else if (out%fits(start)%log_likelihood > &
            out%fits(out%best)%log_likelihood) then
            out%best = start
         end if
      end do
      if (out%best == 0) out%info = 2
   end function mswm_gaussian_multistart

   function mswm_gaussian_random_starts(response, design, regimes, switching, &
      variance_switching, starts, max_iterations, tolerance) result(out)
      !! Generate reproducible shared-RNG Gaussian EM starts and select the best fit.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      integer, intent(in) :: regimes !! Number of latent regimes.
      logical, intent(in) :: switching(:) !! Predictor coefficient switching flags.
      logical, intent(in) :: variance_switching !! Whether variance switches.
      integer, intent(in) :: starts !! Number of random initial fits.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations per start.
      real(dp), intent(in), optional :: tolerance !! Relative convergence tolerance.
      type(mswm_multistart_t) :: out
      real(dp), allocatable :: coefficients(:, :, :), deviation(:, :)
      real(dp), allocatable :: transition(:, :, :), beta(:), standard_error(:)
      real(dp), allocatable :: residual(:)
      real(dp) :: rss, scale, total
      integer :: start, regime, predictor, status

      if (regimes < 1 .or. starts < 1 .or. size(switching) /= size(design, 2)) then
         out%info = 1
         return
      end if
      call ols_fit(design, response, beta, standard_error, residual, rss, status)
      if (status /= 0 .or. rss <= 0.0_dp) then
         out%info = 2
         return
      end if
      scale = sqrt(rss/real(max(1, size(response) - size(design, 2)), dp))
      allocate(coefficients(regimes, size(design, 2), starts))
      allocate(deviation(regimes, starts), transition(regimes, regimes, starts))
      do start = 1, starts
         do regime = 1, regimes
            do predictor = 1, size(design, 2)
               coefficients(regime, predictor, start) = beta(predictor)
               if (switching(predictor)) coefficients(regime, predictor, start) = &
                  beta(predictor) + 0.5_dp*scale*random_standard_normal()
            end do
            deviation(regime, start) = scale*exp(0.25_dp*random_standard_normal())
            do predictor = 1, regimes
               transition(regime, predictor, start) = 0.1_dp + random_uniform()
            end do
            transition(regime, regime, start) = &
               transition(regime, regime, start) + real(regimes, dp)
            total = sum(transition(regime, :, start))
            transition(regime, :, start) = transition(regime, :, start)/total
         end do
         if (.not. variance_switching) deviation(:, start) = &
            sum(deviation(:, start))/real(regimes, dp)
      end do
      out = mswm_gaussian_multistart(response, design, switching, &
         variance_switching, coefficients, deviation, transition, &
         max_iterations, tolerance)
   end function mswm_gaussian_random_starts

   pure function mswm_glm_multistart(response, design, switching, family, &
      initial_coefficients, initial_transition, max_iterations, tolerance, &
      irls_iterations) result(out)
      !! Fit all supplied switching GLM starts and select the largest likelihood.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      logical, intent(in) :: switching(:) !! Predictor coefficient switching flags.
      character(len=*), intent(in) :: family !! Poisson, binomial, or Gamma family and link.
      real(dp), intent(in) :: initial_coefficients(:, :, :) !! Regime-predictor-start values.
      real(dp), intent(in) :: initial_transition(:, :, :) !! Regime transition starts.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations per start.
      real(dp), intent(in), optional :: tolerance !! Relative convergence tolerance.
      integer, intent(in), optional :: irls_iterations !! Maximum IRLS steps per M-step.
      type(mswm_glm_multistart_t) :: out
      integer :: regimes, starts, start

      regimes = size(initial_coefficients, 1)
      starts = size(initial_coefficients, 3)
      if (regimes < 1 .or. starts < 1 .or. &
         size(initial_coefficients, 2) /= size(design, 2) .or. &
         size(switching) /= size(design, 2) .or. &
         any(shape(initial_transition) /= [regimes, regimes, starts])) then
         out%info = 1
         return
      end if
      allocate(out%fits(starts))
      do start = 1, starts
         out%fits(start) = mswm_glm_fit(response, design, regimes, switching, &
            family, initial_coefficients(:, :, start), &
            initial_transition(:, :, start), max_iterations, tolerance, &
            irls_iterations)
         if (out%fits(start)%info /= 0) cycle
         out%successful = out%successful + 1
         if (out%best == 0) then
            out%best = start
         else if (out%fits(start)%log_likelihood > &
            out%fits(out%best)%log_likelihood) then
            out%best = start
         end if
      end do
      if (out%best == 0) out%info = 2
   end function mswm_glm_multistart

   function mswm_glm_random_starts(response, design, regimes, switching, family, &
      starts, max_iterations, tolerance, irls_iterations) result(out)
      !! Generate shared-RNG switching GLM starts and select the best fit.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      integer, intent(in) :: regimes !! Number of latent regimes.
      logical, intent(in) :: switching(:) !! Predictor coefficient switching flags.
      character(len=*), intent(in) :: family !! Poisson, binomial, or Gamma family and link.
      integer, intent(in) :: starts !! Number of random initial fits.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations per start.
      real(dp), intent(in), optional :: tolerance !! Relative convergence tolerance.
      integer, intent(in), optional :: irls_iterations !! Maximum IRLS steps per M-step.
      type(mswm_glm_multistart_t) :: out
      type(mswm_glm_fit_t) :: baseline
      real(dp), allocatable :: coefficients(:, :, :), transition(:, :, :)
      logical, allocatable :: baseline_switching(:)
      real(dp) :: total
      integer :: start, regime, predictor

      if (regimes < 1 .or. starts < 1 .or. size(switching) /= size(design, 2)) then
         out%info = 1
         return
      end if
      allocate(baseline_switching(size(design, 2)))
      baseline_switching = .false.
      baseline = mswm_glm_fit(response, design, 1, baseline_switching, family, &
         max_iterations=100, tolerance=1.0e-8_dp, irls_iterations=irls_iterations)
      if (baseline%info /= 0) then
         out%info = 2
         return
      end if
      allocate(coefficients(regimes, size(design, 2), starts))
      allocate(transition(regimes, regimes, starts))
      do start = 1, starts
         do regime = 1, regimes
            do predictor = 1, size(design, 2)
               coefficients(regime, predictor, start) = &
                  baseline%coefficients(1, predictor)
               if (switching(predictor)) then
                  coefficients(regime, predictor, start) = &
                     coefficients(regime, predictor, start) + &
                     0.5_dp*random_standard_normal()
               end if
            end do
            do predictor = 1, regimes
               transition(regime, predictor, start) = 0.1_dp + random_uniform()
            end do
            transition(regime, regime, start) = &
               transition(regime, regime, start) + real(regimes, dp)
            total = sum(transition(regime, :, start))
            transition(regime, :, start) = transition(regime, :, start)/total
         end do
      end do
      out = mswm_glm_multistart(response, design, switching, family, coefficients, &
         transition, max_iterations, tolerance, irls_iterations)
   end function mswm_glm_random_starts

   pure function mswm_gaussian_residuals(model, regime) result(residual)
      !! Return probability-weighted or regime-specific fitted residuals.
      type(mswm_fit_t), intent(in) :: model !! Fitted Markov-switching model.
      integer, intent(in), optional :: regime !! Optional selected regime.
      real(dp), allocatable :: residual(:)

      if (.not. allocated(model%residuals)) then
         allocate(residual(0))
         return
      end if
      if (present(regime)) then
         if (regime < 1 .or. regime > size(model%residuals, 2)) then
            allocate(residual(0))
            return
         end if
         residual = model%residuals(:, regime)
      else
         residual = sum(model%smoothed_probability*model%residuals, dim=2)
      end if
   end function mswm_gaussian_residuals

   pure function mswm_glm_residuals(model, regime) result(residual)
      !! Return weighted or regime-specific standardized GLM residuals.
      type(mswm_glm_fit_t), intent(in) :: model !! Fitted switching GLM.
      integer, intent(in), optional :: regime !! Optional selected regime.
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: standardized(:, :)
      integer :: selected_regime

      if (.not. allocated(model%residuals) .or. &
         .not. allocated(model%conditional_mean) .or. &
         .not. allocated(model%smoothed_probability) .or. &
         .not. allocated(model%family)) then
         allocate(residual(0))
         return
      end if
      if (any(shape(model%conditional_mean) /= shape(model%residuals)) .or. &
         any(shape(model%smoothed_probability) /= shape(model%residuals)) .or. &
         .not. valid_glm_family(model%family)) then
         allocate(residual(0))
         return
      end if
      allocate(standardized, mold=model%residuals)
      do selected_regime = 1, size(model%residuals, 2)
         standardized(:, selected_regime) = model%residuals(:, selected_regime)/ &
            sqrt(max(glm_variance(model%conditional_mean(:, selected_regime), &
            model%family), tiny(1.0_dp)))
      end do
      if (present(regime)) then
         if (regime < 1 .or. regime > size(standardized, 2)) then
            allocate(residual(0))
            return
         end if
         residual = standardized(:, regime)
      else
         residual = sum(model%smoothed_probability*standardized, dim=2)
      end if
   end function mswm_glm_residuals

   pure function mswm_gaussian_states(model) result(state)
      !! Return the maximum-smoothed-probability regime at each observation.
      type(mswm_fit_t), intent(in) :: model !! Fitted Markov-switching model.
      integer, allocatable :: state(:)
      integer :: time

      if (.not. allocated(model%smoothed_probability)) then
         allocate(state(0))
         return
      end if
      allocate(state(size(model%smoothed_probability, 1)))
      do time = 1, size(state)
         state(time) = maxloc(model%smoothed_probability(time, :), dim=1)
      end do
   end function mswm_gaussian_states

   pure function mswm_glm_states(model) result(state)
      !! Return the maximum-smoothed-probability GLM regime at each observation.
      type(mswm_glm_fit_t), intent(in) :: model !! Fitted switching GLM.
      integer, allocatable :: state(:)
      integer :: time

      if (.not. allocated(model%smoothed_probability)) then
         allocate(state(0))
         return
      end if
      allocate(state(size(model%smoothed_probability, 1)))
      do time = 1, size(state)
         state(time) = maxloc(model%smoothed_probability(time, :), dim=1)
      end do
   end function mswm_glm_states

   pure subroutine smooth_probabilities(filtered, predicted, transition, smoothed, &
      expected_transition)
      !! Apply the backward Kim smoother and compute transition probabilities.
      real(dp), intent(in) :: filtered(:, :) !! Hamilton filtered probabilities.
      real(dp), intent(in) :: predicted(:, :) !! One-step predicted probabilities.
      real(dp), intent(in) :: transition(:, :) !! Current-by-next transition matrix.
      real(dp), allocatable, intent(out) :: smoothed(:, :) !! Smoothed probabilities.
      real(dp), allocatable, intent(out) :: expected_transition(:, :, :) !! Pair probabilities.
      real(dp) :: total
      integer :: observations, regimes, time, current, next

      observations = size(filtered, 1)
      regimes = size(filtered, 2)
      allocate(smoothed(observations, regimes))
      allocate(expected_transition(regimes, regimes, max(0, observations - 1)))
      smoothed = 0.0_dp
      expected_transition = 0.0_dp
      smoothed(observations, :) = filtered(observations, :)
      do time = observations - 1, 1, -1
         do current = 1, regimes
            do next = 1, regimes
               expected_transition(current, next, time) = &
                  filtered(time, current)*transition(current, next)* &
                  smoothed(time + 1, next)/ &
                  max(predicted(time + 1, next), tiny(1.0_dp))
            end do
            smoothed(time, current) = &
               sum(expected_transition(current, :, time))
         end do
         total = sum(smoothed(time, :))
         if (total > tiny(1.0_dp)) then
            smoothed(time, :) = smoothed(time, :)/total
            expected_transition(:, :, time) = &
               expected_transition(:, :, time)/total
         end if
      end do
   end subroutine smooth_probabilities

   pure subroutine maximize_gaussian(response, design, probability, switching, &
      variance_switching, coefficients, standard_deviation, info)
      !! Maximize the Gaussian expected complete-data likelihood.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      real(dp), intent(in) :: probability(:, :) !! Smoothed regime probabilities.
      logical, intent(in) :: switching(:) !! Predictor switching flags.
      logical, intent(in) :: variance_switching !! Whether variance switches.
      real(dp), intent(inout) :: coefficients(:, :) !! Updated regime coefficients.
      real(dp), intent(inout) :: standard_deviation(:) !! Updated error scales.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: expanded_design(:, :), expanded_response(:)
      real(dp), allocatable :: beta(:), standard_error(:), residual(:)
      real(dp), allocatable :: model_residual(:, :)
      real(dp) :: weight, rss, total_weight, pooled
      integer :: observations, predictors, regimes, columns
      integer :: time, regime, predictor, row, column, status

      observations = size(response)
      predictors = size(design, 2)
      regimes = size(probability, 2)
      columns = count(.not. switching) + regimes*count(switching)
      allocate(expanded_design(observations*regimes, columns))
      allocate(expanded_response(observations*regimes))
      expanded_design = 0.0_dp
      do time = 1, observations
         do regime = 1, regimes
            row = (time - 1)*regimes + regime
            weight = sqrt(max(probability(time, regime), 0.0_dp))
            expanded_response(row) = weight*response(time)
            column = 0
            do predictor = 1, predictors
               if (switching(predictor)) then
                  column = column + regimes
                  expanded_design(row, column - regimes + regime) = &
                     weight*design(time, predictor)
               else
                  column = column + 1
                  expanded_design(row, column) = weight*design(time, predictor)
               end if
            end do
         end do
      end do
      call ols_fit(expanded_design, expanded_response, beta, standard_error, &
         residual, rss, status)
      if (status /= 0) then
         info = status
         return
      end if
      column = 0
      do predictor = 1, predictors
         if (switching(predictor)) then
            coefficients(:, predictor) = beta(column + 1:column + regimes)
            column = column + regimes
         else
            column = column + 1
            coefficients(:, predictor) = beta(column)
         end if
      end do
      allocate(model_residual(observations, regimes))
      model_residual = spread(response, 2, regimes) - &
         matmul(design, transpose(coefficients))
      if (variance_switching) then
         do regime = 1, regimes
            total_weight = sum(probability(:, regime))
            standard_deviation(regime) = sqrt(sum(probability(:, regime)* &
               model_residual(:, regime)**2)/max(total_weight, tiny(1.0_dp)))
         end do
      else
         pooled = sum(probability*model_residual**2)/real(observations, dp)
         standard_deviation = sqrt(max(pooled, tiny(1.0_dp)))
      end if
      standard_deviation = max(standard_deviation, sqrt(epsilon(1.0_dp)))
      info = 0
   end subroutine maximize_gaussian

   pure subroutine maximize_transition(smoothed, expected_transition, transition)
      !! Update each transition row from smoothed expected transition counts.
      real(dp), intent(in) :: smoothed(:, :) !! Smoothed state probabilities.
      real(dp), intent(in) :: expected_transition(:, :, :) !! Pair probabilities.
      real(dp), intent(inout) :: transition(:, :) !! Updated transition matrix.
      real(dp) :: denominator
      integer :: current

      do current = 1, size(transition, 1)
         denominator = sum(smoothed(:size(smoothed, 1) - 1, current))
         if (denominator > tiny(1.0_dp)) then
            transition(current, :) = &
               sum(expected_transition(current, :, :), dim=2)/denominator
            transition(current, :) = max(transition(current, :), 0.0_dp)
            transition(current, :) = transition(current, :)/ &
               max(sum(transition(current, :)), tiny(1.0_dp))
         end if
      end do
   end subroutine maximize_transition

   pure subroutine maximize_glm(response, design, probability, switching, family, &
      maximum_iterations, coefficients, info)
      !! Maximize a generalized expected complete-data likelihood by weighted IRLS.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      real(dp), intent(in) :: probability(:, :) !! Smoothed regime probabilities.
      logical, intent(in) :: switching(:) !! Predictor switching flags.
      character(len=*), intent(in) :: family !! Selected observation family and link.
      integer, intent(in) :: maximum_iterations !! Maximum IRLS iterations.
      real(dp), intent(inout) :: coefficients(:, :) !! Updated regime coefficients.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: expanded_design(:, :), weighted_design(:, :)
      real(dp), allocatable :: weighted_response(:), beta(:), old_beta(:)
      real(dp), allocatable :: standard_error(:), residual(:)
      real(dp) :: eta, mean_value, working_response, irls_weight, total_weight
      real(dp) :: rss, change
      integer :: observations, predictors, regimes, columns
      integer :: time, regime, predictor, row, column, iteration, status

      observations = size(response)
      predictors = size(design, 2)
      regimes = size(probability, 2)
      if (size(probability, 1) /= observations .or. &
         any(shape(coefficients) /= [regimes, predictors]) .or. &
         size(switching) /= predictors) then
         info = 1
         return
      end if
      columns = count(.not. switching) + regimes*count(switching)
      allocate(expanded_design(observations*regimes, columns))
      allocate(weighted_design(observations*regimes, columns))
      allocate(weighted_response(observations*regimes), beta(columns))
      expanded_design = 0.0_dp
      column = 0
      do predictor = 1, predictors
         if (switching(predictor)) then
            beta(column + 1:column + regimes) = coefficients(:, predictor)
            column = column + regimes
         else
            column = column + 1
            beta(column) = sum(coefficients(:, predictor))/real(regimes, dp)
         end if
      end do
      do time = 1, observations
         do regime = 1, regimes
            row = (time - 1)*regimes + regime
            column = 0
            do predictor = 1, predictors
               if (switching(predictor)) then
                  column = column + regimes
                  expanded_design(row, column - regimes + regime) = &
                     design(time, predictor)
               else
                  column = column + 1
                  expanded_design(row, column) = design(time, predictor)
               end if
            end do
         end do
      end do

      do iteration = 1, maximum_iterations
         old_beta = beta
         do time = 1, observations
            do regime = 1, regimes
               row = (time - 1)*regimes + regime
               eta = dot_product(expanded_design(row, :), beta)
               mean_value = glm_mean(eta, family)
               call glm_working_values(response(time), eta, mean_value, family, &
                  working_response, irls_weight)
               total_weight = sqrt(max(probability(time, regime)*irls_weight, &
                  tiny(1.0_dp)))
               weighted_design(row, :) = total_weight*expanded_design(row, :)
               weighted_response(row) = total_weight*working_response
            end do
         end do
         call ols_fit(weighted_design, weighted_response, beta, standard_error, &
            residual, rss, status)
         if (status /= 0 .or. .not. all(ieee_is_finite(beta))) then
            info = 10 + status
            return
         end if
         change = maxval(abs(beta - old_beta))/(0.1_dp + maxval(abs(beta)))
         if (change < 1.0e-10_dp) exit
      end do
      column = 0
      do predictor = 1, predictors
         if (switching(predictor)) then
            coefficients(:, predictor) = beta(column + 1:column + regimes)
            column = column + regimes
         else
            column = column + 1
            coefficients(:, predictor) = beta(column)
         end if
      end do
      info = 0
   end subroutine maximize_glm

   pure subroutine glm_working_values(response, eta, mean_value, family, &
      working_response, weight)
      !! Return one GLM IRLS working response and Fisher weight.
      real(dp), intent(in) :: response !! Observed response.
      real(dp), intent(in) :: eta !! Linear predictor.
      real(dp), intent(in) :: mean_value !! Conditional response mean.
      character(len=*), intent(in) :: family !! Selected family and link.
      real(dp), intent(out) :: working_response !! IRLS working response.
      real(dp), intent(out) :: weight !! IRLS Fisher weight.

      select case (family)
      case ('poisson')
         weight = max(mean_value, 1.0e-10_dp)
         working_response = eta + (response - mean_value)/weight
      case ('binomial')
         weight = max(mean_value*(1.0_dp - mean_value), 1.0e-10_dp)
         working_response = eta + (response - mean_value)/weight
      case ('gamma_log')
         weight = 1.0_dp
         working_response = eta + (response - mean_value)/mean_value
      case ('gamma_inverse')
         weight = max(mean_value**2, 1.0e-10_dp)
         working_response = eta - (response - mean_value)/weight
      case default
         weight = 1.0_dp
         working_response = response
      end select
   end subroutine glm_working_values

   pure real(dp) function glm_mean(eta, family) result(mean_value)
      !! Apply the selected inverse link to one linear predictor.
      real(dp), intent(in) :: eta !! Linear predictor.
      character(len=*), intent(in) :: family !! Selected family and link.

      select case (family)
      case ('poisson', 'gamma_log')
         mean_value = exp(max(-40.0_dp, min(40.0_dp, eta)))
      case ('binomial')
         if (eta >= 0.0_dp) then
            mean_value = 1.0_dp/(1.0_dp + exp(-min(40.0_dp, eta)))
         else
            mean_value = exp(max(-40.0_dp, eta))/ &
               (1.0_dp + exp(max(-40.0_dp, eta)))
         end if
         mean_value = max(1.0e-12_dp, min(1.0_dp - 1.0e-12_dp, mean_value))
      case ('gamma_inverse')
         mean_value = 1.0_dp/max(eta, 1.0e-10_dp)
      case default
         mean_value = 0.0_dp
      end select
   end function glm_mean

   pure elemental real(dp) function glm_variance(mean_value, family) &
      result(variance_value)
      !! Return the R family variance used to standardize one GLM residual.
      real(dp), intent(in) :: mean_value !! Conditional response mean.
      character(len=*), intent(in) :: family !! Selected family and link.

      select case (family)
      case ('poisson')
         variance_value = mean_value
      case ('binomial')
         variance_value = mean_value*(1.0_dp - mean_value)
      case ('gamma_log', 'gamma_inverse')
         variance_value = mean_value**2
      case default
         variance_value = 0.0_dp
      end select
   end function glm_variance

   pure real(dp) function glm_log_density(response, mean_value, family) &
      result(log_density)
      !! Evaluate one MSwM generalized observation log density.
      real(dp), intent(in) :: response !! Observed response.
      real(dp), intent(in) :: mean_value !! Conditional response mean.
      character(len=*), intent(in) :: family !! Selected family and link.

      select case (family)
      case ('poisson')
         log_density = response*log(mean_value) - mean_value - &
            log_gamma(response + 1.0_dp)
      case ('binomial')
         log_density = response*log(mean_value) + &
            (1.0_dp - response)*log(1.0_dp - mean_value)
      case ('gamma_log', 'gamma_inverse')
         log_density = (mean_value - 1.0_dp)*log(response) - response - &
            log_gamma(mean_value)
      case default
         log_density = -huge(1.0_dp)
      end select
   end function glm_log_density

   pure logical function valid_glm_family(family) result(valid)
      !! Return whether a generalized MSwM family and link is supported.
      character(len=*), intent(in) :: family !! Candidate family and link.

      valid = family == 'poisson' .or. family == 'binomial' .or. &
         family == 'gamma_log' .or. family == 'gamma_inverse'
   end function valid_glm_family

   pure logical function valid_glm_response(response, family) result(valid)
      !! Check the response support for a generalized observation family.
      real(dp), intent(in) :: response(:) !! Candidate response values.
      character(len=*), intent(in) :: family !! Selected family and link.

      valid = all(ieee_is_finite(response))
      if (.not. valid) return
      select case (family)
      case ('poisson')
         valid = all(response >= 0.0_dp .and. &
            abs(response - real(nint(response), dp)) < 1.0e-10_dp)
      case ('binomial')
         valid = all(response == 0.0_dp .or. response == 1.0_dp)
      case ('gamma_log', 'gamma_inverse')
         valid = all(response > 0.0_dp)
      case default
         valid = .false.
      end select
   end function valid_glm_response

   pure subroutine pack_gaussian_parameters(model, parameters)
      !! Pack a Gaussian fit into unconstrained Hessian coordinates.
      type(mswm_fit_t), intent(in) :: model !! Fitted Gaussian switching model.
      real(dp), allocatable, intent(out) :: parameters(:) !! Packed coordinates.
      integer :: regimes, predictors, scales, count_value, position
      integer :: regime, predictor, target

      regimes = size(model%coefficients, 1)
      predictors = size(model%coefficients, 2)
      scales = merge(regimes, 1, model%variance_switching)
      count_value = scales + regimes*(regimes - 1) + &
         count(.not. model%switching) + regimes*count(model%switching)
      allocate(parameters(count_value))
      position = 0
      if (model%variance_switching) then
         parameters(1:regimes) = log(model%standard_deviation)
         position = regimes
      else
         parameters(1) = log(model%standard_deviation(1))
         position = 1
      end if
      do regime = 1, regimes
         do target = 1, regimes - 1
            position = position + 1
            parameters(position) = log(max(model%transition(regime, target), &
               tiny(1.0_dp))/max(model%transition(regime, regimes), tiny(1.0_dp)))
         end do
      end do
      do predictor = 1, predictors
         if (model%switching(predictor)) then
            parameters(position + 1:position + regimes) = &
               model%coefficients(:, predictor)
            position = position + regimes
         else
            position = position + 1
            parameters(position) = model%coefficients(1, predictor)
         end if
      end do
   end subroutine pack_gaussian_parameters

   pure subroutine unpack_gaussian_parameters(parameters, regimes, switching, &
      variance_switching, coefficients, standard_deviation, transition)
      !! Unpack unconstrained Gaussian coordinates into natural parameters.
      real(dp), intent(in) :: parameters(:) !! Packed unconstrained coordinates.
      integer, intent(in) :: regimes !! Number of latent regimes.
      logical, intent(in) :: switching(:) !! Predictor switching flags.
      logical, intent(in) :: variance_switching !! Whether variance switches.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Natural coefficients.
      real(dp), allocatable, intent(out) :: standard_deviation(:) !! Positive scales.
      real(dp), allocatable, intent(out) :: transition(:, :) !! Row-stochastic transitions.
      real(dp), allocatable :: logits(:)
      real(dp) :: maximum_logit
      integer :: predictors, position, regime, predictor

      predictors = size(switching)
      allocate(coefficients(regimes, predictors), standard_deviation(regimes))
      allocate(transition(regimes, regimes), logits(regimes))
      position = 0
      if (variance_switching) then
         standard_deviation = exp(parameters(1:regimes))
         position = regimes
      else
         standard_deviation = exp(parameters(1))
         position = 1
      end if
      do regime = 1, regimes
         logits = 0.0_dp
         if (regimes > 1) logits(:regimes - 1) = &
            parameters(position + 1:position + regimes - 1)
         position = position + regimes - 1
         maximum_logit = maxval(logits)
         transition(regime, :) = exp(logits - maximum_logit)
         transition(regime, :) = transition(regime, :)/sum(transition(regime, :))
      end do
      do predictor = 1, predictors
         if (switching(predictor)) then
            coefficients(:, predictor) = parameters(position + 1:position + regimes)
            position = position + regimes
         else
            position = position + 1
            coefficients(:, predictor) = parameters(position)
         end if
      end do
   end subroutine unpack_gaussian_parameters

   pure subroutine pack_glm_parameters(model, parameters)
      !! Pack a switching GLM into unconstrained Hessian coordinates.
      type(mswm_glm_fit_t), intent(in) :: model !! Fitted switching GLM.
      real(dp), allocatable, intent(out) :: parameters(:) !! Packed coordinates.
      integer :: regimes, predictors, count_value, position
      integer :: regime, predictor, target

      regimes = size(model%coefficients, 1)
      predictors = size(model%coefficients, 2)
      count_value = regimes*(regimes - 1) + count(.not. model%switching) + &
         regimes*count(model%switching)
      allocate(parameters(count_value))
      position = 0
      do regime = 1, regimes
         do target = 1, regimes - 1
            position = position + 1
            parameters(position) = log(max(model%transition(regime, target), &
               tiny(1.0_dp))/max(model%transition(regime, regimes), tiny(1.0_dp)))
         end do
      end do
      do predictor = 1, predictors
         if (model%switching(predictor)) then
            parameters(position + 1:position + regimes) = &
               model%coefficients(:, predictor)
            position = position + regimes
         else
            position = position + 1
            parameters(position) = model%coefficients(1, predictor)
         end if
      end do
   end subroutine pack_glm_parameters

   pure subroutine unpack_glm_parameters(parameters, regimes, switching, &
      coefficients, transition)
      !! Unpack unconstrained switching GLM coordinates into natural parameters.
      real(dp), intent(in) :: parameters(:) !! Packed unconstrained coordinates.
      integer, intent(in) :: regimes !! Number of latent regimes.
      logical, intent(in) :: switching(:) !! Predictor switching flags.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Natural coefficients.
      real(dp), allocatable, intent(out) :: transition(:, :) !! Row-stochastic transitions.
      real(dp), allocatable :: logits(:)
      real(dp) :: maximum_logit
      integer :: predictors, position, regime, predictor

      predictors = size(switching)
      allocate(coefficients(regimes, predictors))
      allocate(transition(regimes, regimes), logits(regimes))
      position = 0
      do regime = 1, regimes
         logits = 0.0_dp
         if (regimes > 1) logits(:regimes - 1) = &
            parameters(position + 1:position + regimes - 1)
         position = position + regimes - 1
         maximum_logit = maxval(logits)
         transition(regime, :) = exp(logits - maximum_logit)
         transition(regime, :) = transition(regime, :)/sum(transition(regime, :))
      end do
      do predictor = 1, predictors
         if (switching(predictor)) then
            coefficients(:, predictor) = parameters(position + 1:position + regimes)
            position = position + regimes
         else
            position = position + 1
            coefficients(:, predictor) = parameters(position)
         end if
      end do
   end subroutine unpack_glm_parameters

end module mswm_mod
