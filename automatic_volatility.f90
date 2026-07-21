! SPDX-License-Identifier: MIT
! SPDX-FileComment: Automatic univariate conditional-volatility model selection.
module automatic_volatility_mod
   !! Fit, compare, and forecast a focused set of volatility models.
   use kind_mod, only: dp
   use rugarch_mod, only: rugarch_spec_t, rugarch_fit_t, &
      rugarch_forecast_t, rugarch_spec, rugarch_fit, rugarch_forecast, &
      rugarch_model_sgarch, rugarch_distribution_normal, &
      rugarch_distribution_student
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter :: model_constant = 1
   integer, parameter :: model_arch = 2
   integer, parameter :: model_garch_normal = 3
   integer, parameter :: model_garch_student = 4
   character(len=*), parameter :: selection_validation = "validation"
   character(len=*), parameter :: selection_aicc = "aicc"
   character(len=*), parameter :: selection_bic = "bic"

   type, public :: automatic_volatility_options_t
      !! Controls volatility-model validation, forecasting, and timing.
      integer :: horizon = 10
      integer :: validation_size = 0
      logical :: time_fits = .false.
      character(len=10) :: selection = selection_validation
   end type automatic_volatility_options_t

   type, public :: volatility_candidate_result_t
      !! Fit statistics and parameters for one volatility candidate.
      character(len=:), allocatable :: name
      integer :: model_code = 0
      real(dp) :: validation_qlike = huge(1.0_dp)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: mean = 0.0_dp
      real(dp) :: omega = 0.0_dp
      real(dp) :: shape = 0.0_dp
      real(dp) :: persistence = 0.0_dp
      real(dp) :: unconditional_variance = 0.0_dp
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: beta(:)
      integer :: parameter_count = 0
      real(dp) :: full_fit_seconds = 0.0_dp
      real(dp) :: validation_fit_seconds = 0.0_dp
      logical :: optimizer_converged = .false.
      logical :: selected = .false.
      integer :: info = 0
   end type volatility_candidate_result_t

   type, public :: automatic_volatility_result_t
      !! Ranked volatility candidates and selected conditional-sigma forecasts.
      type(volatility_candidate_result_t), allocatable :: candidates(:)
      real(dp), allocatable :: sigma_forecast(:)
      integer :: validation_size = 0
      character(len=10) :: selection_criterion = selection_validation
      logical :: fit_times_recorded = .false.
      integer :: selected = 0
      integer :: info = 0
   end type automatic_volatility_result_t

   type :: volatility_fit_t
      type(rugarch_fit_t) :: rugarch
      real(dp), allocatable :: mean_forecast(:)
      real(dp), allocatable :: variance_forecast(:)
      real(dp) :: mean = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: parameter_count = 0
      integer :: info = 0
   end type volatility_fit_t

   interface display
      module procedure display_automatic_volatility_result
   end interface display

   public :: automatic_volatility, display

contains

   function automatic_volatility(series, options) result(out)
      !! Fit volatility candidates and rank them by validation or information criterion.
      real(dp), intent(in) :: series(:) !! Finite transformed observations.
      type(automatic_volatility_options_t), intent(in), optional :: options !! Search options.
      type(automatic_volatility_result_t) :: out
      type(automatic_volatility_options_t) :: selected_options
      type(volatility_fit_t) :: fit, validation_fit
      type(volatility_fit_t), allocatable :: full_fits(:)
      real(dp), allocatable :: full_times(:), validation_times(:)
      character(len=10) :: criterion
      integer :: clock_end, clock_max, clock_rate, clock_start
      integer :: code, training_size, validation_size

      selected_options = automatic_volatility_options_t()
      if (present(options)) selected_options = options
      criterion = lower_ascii(trim(selected_options%selection))
      if (size(series) < 30 .or. selected_options%horizon < 1 .or. &
         any(.not. ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (criterion /= selection_validation .and. &
         criterion /= selection_aicc .and. criterion /= selection_bic) then
         out%info = 2
         return
      end if
      validation_size = 0
      if (criterion == selection_validation) then
         validation_size = selected_options%validation_size
         if (validation_size <= 0) then
            validation_size = max(10, min(50, size(series)/5))
         end if
         validation_size = min(validation_size, size(series) - 20)
      end if
      training_size = size(series) - validation_size
      out%validation_size = validation_size
      out%selection_criterion = criterion
      out%fit_times_recorded = selected_options%time_fits
      allocate(out%candidates(4), full_fits(4), full_times(4), &
         validation_times(4))
      full_times = 0.0_dp
      validation_times = 0.0_dp
      do code = 1, 4
         if (selected_options%time_fits) then
            call system_clock(clock_start, clock_rate, clock_max)
         end if
         fit = fit_volatility_candidate(series, selected_options%horizon, code)
         full_fits(code) = fit
         if (selected_options%time_fits) then
            call system_clock(clock_end)
            full_times(code) = clock_seconds(clock_start, clock_end, &
               clock_rate, clock_max)
         end if
         out%candidates(code) = summarize_candidate(fit, code, size(series))
         out%candidates(code)%full_fit_seconds = full_times(code)
         if (criterion == selection_validation) then
            if (selected_options%time_fits) then
               call system_clock(clock_start, clock_rate, clock_max)
            end if
            validation_fit = fit_volatility_candidate( &
               series(:training_size), validation_size, code)
            out%candidates(code)%validation_qlike = validation_qlike( &
               validation_fit, series(training_size + 1:))
            if (selected_options%time_fits) then
               call system_clock(clock_end)
               validation_times(code) = clock_seconds(clock_start, &
                  clock_end, clock_rate, clock_max)
               out%candidates(code)%validation_fit_seconds = &
                  validation_times(code)
            end if
         end if
      end do
      call sort_candidates(out%candidates, criterion)
      if (out%candidates(1)%info /= 0) then
         out%info = 3
         return
      end if
      out%selected = 1
      out%candidates(1)%selected = .true.
      code = out%candidates(1)%model_code
      out%sigma_forecast = sqrt(max(full_fits(code)%variance_forecast, &
         0.0_dp))
   end function automatic_volatility

   pure function fit_volatility_candidate(series, horizon, code) result(out)
      !! Fit one constant-variance, ARCH, or GARCH candidate.
      real(dp), intent(in) :: series(:) !! Transformed observations.
      integer, intent(in) :: horizon !! Forecast horizon.
      integer, intent(in) :: code !! Internal volatility model code.
      type(volatility_fit_t) :: out
      type(rugarch_spec_t) :: specification
      type(rugarch_forecast_t) :: forecast
      real(dp) :: pi, residual_sum

      if (code == model_constant) then
         out%mean = sum(series)/real(size(series), dp)
         residual_sum = sum((series - out%mean)**2)
         out%variance = max(residual_sum/real(size(series), dp), &
            tiny(1.0_dp))
         pi = acos(-1.0_dp)
         out%log_likelihood = -0.5_dp*real(size(series), dp)* &
            (log(2.0_dp*pi*out%variance) + 1.0_dp)
         out%parameter_count = 2
         out%aic = -2.0_dp*out%log_likelihood + 4.0_dp
         out%bic = -2.0_dp*out%log_likelihood + &
            2.0_dp*log(real(size(series), dp))
         allocate(out%mean_forecast(horizon), out%variance_forecast(horizon))
         out%mean_forecast = out%mean
         out%variance_forecast = out%variance
         return
      end if
      select case (code)
      case (model_arch)
         specification = rugarch_spec(variance_model=rugarch_model_sgarch, &
            distribution=rugarch_distribution_normal, arch_order=1, &
            garch_order=0, include_mean=.true.)
      case (model_garch_normal)
         specification = rugarch_spec(variance_model=rugarch_model_sgarch, &
            distribution=rugarch_distribution_normal, arch_order=1, &
            garch_order=1, include_mean=.true.)
      case (model_garch_student)
         specification = rugarch_spec(variance_model=rugarch_model_sgarch, &
            distribution=rugarch_distribution_student, arch_order=1, &
            garch_order=1, include_mean=.true.)
      end select
      out%rugarch = rugarch_fit(series, specification, max_iterations=200, &
         tolerance=1.0e-5_dp)
      if (out%rugarch%info /= 0) then
         out%info = 10 + out%rugarch%info
         return
      end if
      forecast = rugarch_forecast(out%rugarch, horizon)
      if (forecast%info /= 0) then
         out%info = 20 + forecast%info
         return
      end if
      out%mean_forecast = forecast%mean
      out%variance_forecast = forecast%variance
      out%log_likelihood = out%rugarch%log_likelihood
      out%aic = out%rugarch%aic
      out%bic = out%rugarch%bic
      out%parameter_count = size(out%rugarch%coefficients)
   end function fit_volatility_candidate

   pure function summarize_candidate(fit, code, observations) result(out)
      !! Convert one internal volatility fit to a public candidate summary.
      type(volatility_fit_t), intent(in) :: fit !! Internal fitted candidate.
      integer, intent(in) :: code !! Internal volatility model code.
      integer, intent(in) :: observations !! Number of fitted observations.
      type(volatility_candidate_result_t) :: out

      out%name = model_name(code)
      out%model_code = code
      out%info = fit%info
      if (fit%info /= 0) return
      out%log_likelihood = fit%log_likelihood
      out%aic = fit%aic
      out%bic = fit%bic
      out%aicc = fit%aic
      if (observations > fit%parameter_count + 1) then
         out%aicc = fit%aic + &
            2.0_dp*real(fit%parameter_count*(fit%parameter_count + 1), dp)/ &
            real(observations - fit%parameter_count - 1, dp)
      end if
      out%parameter_count = fit%parameter_count
      if (code == model_constant) then
         out%mean = fit%mean
         out%omega = fit%variance
         out%unconditional_variance = fit%variance
         out%optimizer_converged = .true.
         allocate(out%alpha(0), out%beta(0))
         return
      end if
      out%mean = fit%rugarch%parameters%mean
      out%omega = fit%rugarch%parameters%omega
      out%shape = fit%rugarch%parameters%shape
      out%persistence = fit%rugarch%persistence
      out%unconditional_variance = fit%rugarch%unconditional_variance
      out%alpha = fit%rugarch%parameters%alpha
      out%beta = fit%rugarch%parameters%beta
      out%optimizer_converged = fit%rugarch%converged
   end function summarize_candidate

   pure real(dp) function validation_qlike(fit, actual) result(score)
      !! Compute mean Gaussian QLIKE loss on a held-out tail.
      type(volatility_fit_t), intent(in) :: fit !! Training fit and forecasts.
      real(dp), intent(in) :: actual(:) !! Held-out observations.
      real(dp), allocatable :: error(:), variance(:)

      score = huge(1.0_dp)
      if (fit%info /= 0 .or. .not. allocated(fit%mean_forecast) .or. &
         .not. allocated(fit%variance_forecast)) return
      if (size(fit%mean_forecast) /= size(actual) .or. &
         size(fit%variance_forecast) /= size(actual)) return
      variance = max(fit%variance_forecast, tiny(1.0_dp))
      error = actual - fit%mean_forecast
      score = sum(log(variance) + error**2/variance)/real(size(actual), dp)
   end function validation_qlike

   pure subroutine sort_candidates(candidates, criterion)
      !! Sort volatility candidates by the requested selection statistic.
      type(volatility_candidate_result_t), intent(inout) :: candidates(:) !! Candidates to sort.
      character(len=*), intent(in) :: criterion !! Selection criterion.
      type(volatility_candidate_result_t) :: held
      integer :: i, j

      do i = 2, size(candidates)
         held = candidates(i)
         j = i - 1
         do while (j >= 1)
            if (candidate_score(candidates(j), criterion) <= &
               candidate_score(held, criterion)) exit
            candidates(j + 1) = candidates(j)
            j = j - 1
         end do
         candidates(j + 1) = held
      end do
   end subroutine sort_candidates

   pure real(dp) function candidate_score(candidate, criterion) result(score)
      !! Return the ranking value for one volatility candidate.
      type(volatility_candidate_result_t), intent(in) :: candidate !! Candidate to score.
      character(len=*), intent(in) :: criterion !! Selection criterion.

      if (candidate%info /= 0) then
         score = huge(1.0_dp)
      else if (trim(criterion) == selection_validation) then
         score = candidate%validation_qlike
      else if (trim(criterion) == selection_bic) then
         score = candidate%bic
      else
         score = candidate%aicc
      end if
   end function candidate_score

   pure function model_name(code) result(name)
      !! Return a readable volatility candidate name.
      integer, intent(in) :: code !! Internal volatility model code.
      character(len=:), allocatable :: name

      select case (code)
      case (model_constant)
         name = "Gaussian constant variance"
      case (model_arch)
         name = "Gaussian ARCH(1)"
      case (model_garch_normal)
         name = "Gaussian GARCH(1,1)"
      case (model_garch_student)
         name = "Student-t GARCH(1,1)"
      end select
   end function model_name

   subroutine display_automatic_volatility_result(result, print_parameters, &
      max_models, max_forecasts)
      !! Display ranked volatility candidates and conditional-sigma forecasts.
      type(automatic_volatility_result_t), intent(in) :: result !! Volatility result.
      logical, intent(in), optional :: print_parameters !! Print candidate parameters.
      integer, intent(in), optional :: max_models !! Maximum candidate summaries.
      integer, intent(in), optional :: max_forecasts !! Maximum sigma forecasts.
      integer :: display_count, forecast_count, i
      logical :: show_parameters

      write(*, '(a)') "Automatic volatility model comparison"
      write(*, '(a,i0)') "  status: ", result%info
      if (result%info /= 0) return
      write(*, '(a,a)') "  selection criterion: ", &
         trim(result%selection_criterion)
      if (result%selection_criterion == selection_validation) then
         write(*, '(a,i0)') "  validation observations: ", &
            result%validation_size
      end if
      show_parameters = .false.
      if (present(print_parameters)) show_parameters = print_parameters
      display_count = size(result%candidates)
      if (present(max_models)) then
         if (max_models > 0) display_count = min(max_models, display_count)
      end if
      do i = 1, display_count
         call display_volatility_candidate(result%candidates(i), &
            result%selection_criterion, result%fit_times_recorded, &
            show_parameters, i)
      end do
      write(*, '(a,a)') "Selected volatility model: ", &
         result%candidates(result%selected)%name
      forecast_count = size(result%sigma_forecast)
      if (present(max_forecasts)) then
         if (max_forecasts >= 0) then
            forecast_count = min(max_forecasts, forecast_count)
         end if
      end if
      if (forecast_count > 0) then
         write(*, '(a)') "Conditional standard-deviation forecasts"
         do i = 1, forecast_count
            write(*, '(2x,i5,2x,es16.8)') i, result%sigma_forecast(i)
         end do
      end if
   end subroutine display_automatic_volatility_result

   subroutine display_volatility_candidate(candidate, criterion, timed, &
      print_parameters, rank)
      !! Display one volatility candidate summary and optional parameters.
      type(volatility_candidate_result_t), intent(in) :: candidate !! Candidate to display.
      character(len=*), intent(in) :: criterion !! Selection criterion.
      logical, intent(in) :: timed !! Whether fit timings were requested.
      logical, intent(in) :: print_parameters !! Whether parameters are printed.
      integer, intent(in) :: rank !! Candidate rank.
      integer :: i

      write(*, '(2x,i3,2x,a)') rank, candidate%name
      if (candidate%info /= 0) then
         write(*, '(7x,a,i0)') "fit failed with status ", candidate%info
         write(*, '(a)') ""
         return
      end if
      if (trim(criterion) == selection_validation) then
         write(*, '(7x,a,es12.4)') "validation QLIKE ", &
            candidate%validation_qlike
      end if
      write(*, '(7x,a,es12.4,2x,a,es12.4,2x,a,es12.4,2x,a,es12.4)') &
         "logLik ", candidate%log_likelihood, "AIC ", candidate%aic, &
         "AICc ", candidate%aicc, "BIC ", candidate%bic
      write(*, '(7x,a,i0,2x,a,l1)') "parameter count ", &
         candidate%parameter_count, "optimizer converged ", &
         candidate%optimizer_converged
      if (timed) then
         write(*, '(7x,a,f10.4,a)') "full-data fit time ", &
            candidate%full_fit_seconds, " seconds"
         if (trim(criterion) == selection_validation) then
            write(*, '(7x,a,f10.4,a)') "validation fit time ", &
               candidate%validation_fit_seconds, " seconds"
         end if
      end if
      if (print_parameters) then
         write(*, '(7x,a,f16.8)') "mean ", candidate%mean
         write(*, '(7x,a,f16.8)') "omega/variance ", candidate%omega
         do i = 1, size(candidate%alpha)
            write(*, '(7x,a,i0,a,f16.8)') "alpha(", i, ") ", &
               candidate%alpha(i)
         end do
         do i = 1, size(candidate%beta)
            write(*, '(7x,a,i0,a,f16.8)') "beta(", i, ") ", &
               candidate%beta(i)
         end do
         if (candidate%model_code == model_garch_student) then
            write(*, '(7x,a,f16.8)') "shape ", candidate%shape
         end if
         write(*, '(7x,a,f16.8)') "persistence ", candidate%persistence
         write(*, '(7x,a,f16.8)') "unconditional variance ", &
            candidate%unconditional_variance
      end if
      write(*, '(a)') ""
   end subroutine display_volatility_candidate

   pure real(dp) function clock_seconds(start_count, end_count, count_rate, &
      maximum_count) result(seconds)
      !! Convert system-clock counts to elapsed seconds with wrap handling.
      integer, intent(in) :: start_count !! Starting clock count.
      integer, intent(in) :: end_count !! Ending clock count.
      integer, intent(in) :: count_rate !! Clock counts per second.
      integer, intent(in) :: maximum_count !! Maximum count before wrapping.
      integer :: elapsed

      if (count_rate <= 0) then
         seconds = 0.0_dp
      else
         elapsed = end_count - start_count
         if (elapsed < 0) elapsed = maximum_count - start_count + end_count + 1
         seconds = real(elapsed, dp)/real(count_rate, dp)
      end if
   end function clock_seconds

   pure function lower_ascii(value) result(lower)
      !! Convert ASCII letters to lowercase for option matching.
      character(len=*), intent(in) :: value !! Text to convert.
      character(len=len(value)) :: lower
      integer :: code, i

      lower = value
      do i = 1, len(value)
         code = iachar(lower(i:i))
         if (code >= iachar("A") .and. code <= iachar("Z")) then
            lower(i:i) = achar(code + iachar("a") - iachar("A"))
         end if
      end do
   end function lower_ascii

end module automatic_volatility_mod
