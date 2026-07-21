! SPDX-License-Identifier: MIT
! SPDX-FileComment: Command-line automatic modeling for a univariate series.
program ts_auto
   !! Read a numeric series, profile it, compare models, and forecast it.
   use kind_mod, only: dp
   use automatic_modeling_mod, only: automatic_model_options_t, &
      automatic_model_result_t, automatic_model, display
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: error_unit
   implicit none

   type :: series_read_t
      real(dp), allocatable :: values(:)
      integer :: skipped = 0
      logical :: truncated = .false.
      integer :: info = 0
   end type series_read_t

   type(automatic_model_options_t) :: options
   type(automatic_model_result_t) :: result
   type(series_read_t) :: input
   character(len=256) :: option_error
   character(len=1024) :: path
   integer :: clock_end, clock_max, clock_rate, clock_start
   integer :: display_lags, maximum_models, maximum_observations, status
   real(dp) :: elapsed_seconds
   logical :: print_all_ar, print_all_arma, print_parameters

   call system_clock(clock_start, clock_rate, clock_max)
   if (command_argument_count() < 1) then
      call print_usage()
      stop 1
   end if
   call get_command_argument(1, path, status=status)
   if (status /= 0 .or. len_trim(path) == 0) then
      call print_usage()
      stop 1
   end if
   call parse_options(options, maximum_observations, display_lags, &
      maximum_models, print_parameters, print_all_ar, print_all_arma, &
      status, option_error)
   if (status /= 0) then
      write(error_unit, '(a)') trim(option_error)
      call print_usage()
      stop 1
   end if
   input = read_series(trim(path), maximum_observations)
   if (input%info /= 0) then
      write(*, '(a,i0)') "Unable to read the input series; status ", input%info
      stop 1
   end if
   if (input%skipped > 0) then
      write(*, '(a,i0,a)') "Skipped ", input%skipped, " nonnumeric or nonfinite lines."
   end if
   write(*, '(a,i0,a)') "Read ", size(input%values), " valid observations."
   if (input%truncated) then
      write(*, '(a,i0,a)') "Input truncated at the requested limit of ", &
         maximum_observations, " observations."
   end if
   result = automatic_model(input%values, options)
   call display(result, display_lags, print_parameters, maximum_models, &
      print_all_ar, print_all_arma)
   call system_clock(clock_end)
   elapsed_seconds = elapsed_wall_time(clock_start, clock_end, &
      clock_rate, clock_max)
   write(*, '(a)') ""
   write(*, '(a,f12.3,a)') "Elapsed time: ", elapsed_seconds, " seconds"
   if (result%info /= 0) stop 1

contains

   pure function elapsed_wall_time(start_count, end_count, count_rate, &
      maximum_count) result(seconds)
      !! Convert system-clock counts to elapsed wall-clock seconds.
      integer, intent(in) :: start_count !! System-clock count at program start.
      integer, intent(in) :: end_count !! System-clock count after reporting results.
      integer, intent(in) :: count_rate !! System-clock counts per second.
      integer, intent(in) :: maximum_count !! Largest system-clock count before wrapping.
      real(dp) :: seconds
      integer :: elapsed_counts

      if (count_rate <= 0) then
         seconds = 0.0_dp
         return
      end if
      if (end_count >= start_count) then
         elapsed_counts = end_count - start_count
      else
         elapsed_counts = maximum_count - start_count + end_count + 1
      end if
      seconds = real(elapsed_counts, dp)/real(count_rate, dp)
   end function elapsed_wall_time

   subroutine parse_options(options, maximum_observations, display_lags, &
      maximum_models, print_parameters, print_all_ar, print_all_arma, &
      info, error_message)
      !! Parse command-line settings following the input path.
      type(automatic_model_options_t), intent(out) :: options !! Parsed modeling options.
      integer, intent(out) :: maximum_observations !! Maximum valid input values, or zero for all.
      integer, intent(out) :: display_lags !! Maximum positive correlation lag to print.
      integer, intent(out) :: maximum_models !! Maximum model summaries, or zero for all.
      logical, intent(out) :: print_parameters !! Print fitted model parameters when true.
      logical, intent(out) :: print_all_ar !! Print all tested autoregressive orders.
      logical, intent(out) :: print_all_arma !! Print all tested ARMA orders.
      integer, intent(out) :: info !! Zero on successful parsing.
      character(len=*), intent(out) :: error_message !! Explanation when parsing fails.
      character(len=128) :: argument, value
      integer :: i, read_status

      options = automatic_model_options_t()
      maximum_observations = 0
      display_lags = 5
      maximum_models = 0
      print_parameters = .false.
      print_all_ar = .false.
      print_all_arma = .false.
      info = 0
      error_message = ""
      i = 2
      do while (i <= command_argument_count())
         call get_command_argument(i, argument)
         if (trim(argument) == "--help" .or. trim(argument) == "-h") then
            call print_usage()
            stop
         end if
         if (trim(argument) == "--print-parameters" .or. &
            trim(argument) == "--print-param" .or. &
            trim(argument) == "--param") then
            print_parameters = .true.
            i = i + 1
            cycle
         end if
         if (trim(argument) == "--time-fits") then
            options%time_fits = .true.
            i = i + 1
            cycle
         end if
         if (trim(argument) == "--print-all-ar") then
            print_all_ar = .true.
            i = i + 1
            cycle
         end if
         if (trim(argument) == "--print-all-arma") then
            print_all_arma = .true.
            i = i + 1
            cycle
         end if
         if (.not. option_requires_value(trim(argument))) then
            info = 2
            error_message = "Unrecognized option: " // trim(argument)
            return
         end if
         if (i == command_argument_count()) then
            info = 1
            error_message = "Option " // trim(argument) // " requires a value."
            return
         end if
         call get_command_argument(i + 1, value)
         select case (trim(argument))
         case ("--frequency")
            read(value, *, iostat=read_status) options%frequency
         case ("--horizon")
            read(value, *, iostat=read_status) options%horizon
         case ("--max-lag")
            read(value, *, iostat=read_status) options%max_lag
         case ("--max-ar")
            read(value, *, iostat=read_status) options%max_ar_order
         case ("--validation")
            read(value, *, iostat=read_status) options%validation_size
         case ("--selection")
            options%selection = trim(value)
            read_status = 0
         case ("--observations", "--obs")
            read(value, *, iostat=read_status) maximum_observations
         case ("--display-lags")
            read(value, *, iostat=read_status) display_lags
         case ("--max-models")
            read(value, *, iostat=read_status) maximum_models
         case default
            info = 2
            error_message = "Unrecognized option: " // trim(argument)
            return
         end select
         if (read_status /= 0) then
            info = 3
            error_message = "Invalid value for " // trim(argument) // &
               ": " // trim(value)
            return
         end if
         if (maximum_observations < 0) then
            info = 4
            error_message = "The observation limit cannot be negative."
            return
         end if
         if (display_lags < 1) then
            info = 5
            error_message = "The number of displayed lags must be positive."
            return
         end if
         if (maximum_models < 0) then
            info = 7
            error_message = "The maximum number of displayed models cannot be negative."
            return
         end if
         i = i + 2
      end do
      if (options%max_lag < 1) then
         info = 6
         error_message = "The maximum correlation lag must be positive."
         return
      end if
      options%max_lag = max(options%max_lag, display_lags)
   end subroutine parse_options

   pure logical function option_requires_value(argument) result(required)
      !! Identify recognized command options that require a following value.
      character(len=*), intent(in) :: argument !! Command option to classify.

      select case (trim(argument))
      case ("--frequency", "--horizon", "--max-lag", "--max-ar", &
         "--validation", "--selection", "--observations", "--obs", &
         "--display-lags", "--max-models")
         required = .true.
      case default
         required = .false.
      end select
   end function option_requires_value

   function read_series(path, maximum_observations) result(out)
      !! Read the first numeric value from each nonempty text line.
      character(len=*), intent(in) :: path !! Input text or CSV path.
      integer, intent(in) :: maximum_observations !! Maximum valid values, or zero for all.
      type(series_read_t) :: out
      character(len=4096) :: line
      real(dp) :: value
      integer :: unit, io_status, count

      open(newunit=unit, file=path, status="old", action="read", iostat=io_status)
      if (io_status /= 0) then
         out%info = 1
         return
      end if
      count = 0
      do
         read(unit, '(a)', iostat=io_status) line
         if (io_status < 0) exit
         if (io_status > 0) then
            close(unit)
            out%info = 2
            return
         end if
         read(line, *, iostat=io_status) value
         if (io_status == 0 .and. ieee_is_finite(value)) then
            if (maximum_observations > 0 .and. &
               count >= maximum_observations) then
               out%truncated = .true.
               exit
            end if
            count = count + 1
         else if (len_trim(line) > 0) then
            out%skipped = out%skipped + 1
         end if
      end do
      if (count == 0) then
         close(unit)
         out%info = 3
         return
      end if
      rewind(unit)
      allocate(out%values(count))
      count = 0
      do
         read(unit, '(a)', iostat=io_status) line
         if (io_status < 0) exit
         if (io_status > 0) then
            close(unit)
            out%info = 4
            return
         end if
         read(line, *, iostat=io_status) value
         if (io_status == 0 .and. ieee_is_finite(value)) then
            count = count + 1
            out%values(count) = value
            if (count == size(out%values)) exit
         end if
      end do
      close(unit)
   end function read_series

   subroutine print_usage()
      !! Display command-line syntax and supported options.

      write(*, '(a)') "Usage: ts_auto FILE [options]"
      write(*, '(a)') ""
      write(*, '(a)') "FILE contains one numeric observation per line; a header is allowed."
      write(*, '(a)') "Options:"
      write(*, '(a)') "  --frequency N    observations per seasonal cycle (default 1)"
      write(*, '(a)') "  --horizon N      forecast horizon (default 10)"
      write(*, '(a)') "  --max-lag N      maximum ACF and PACF lag (default 24)"
      write(*, '(a)') "  --display-lags N positive correlation lags to print (default 5)"
      write(*, '(a)') "  --max-ar N       maximum autoregressive order (default 12)"
      write(*, '(a)') "  --validation N   held-out tail length (default automatic)"
      write(*, '(a)') "  --selection NAME validation, aicc, or bic (default validation)"
      write(*, '(a)') "  --print-parameters print fitted parameters for every candidate"
      write(*, '(a)') "  --print-param     alias for --print-parameters"
      write(*, '(a)') "  --param           alias for --print-parameters"
      write(*, '(a)') "  --time-fits       report elapsed time for each candidate fit"
      write(*, '(a)') "  --print-all-ar    report every tested autoregressive order"
      write(*, '(a)') "  --print-all-arma  report every tested ARMA order"
      write(*, '(a)') "  --max-models N    maximum model summaries; zero displays all"
      write(*, '(a)') "  --observations N use at most the first N valid observations"
      write(*, '(a)') "  --obs N          alias for --observations"
   end subroutine print_usage

end program ts_auto
