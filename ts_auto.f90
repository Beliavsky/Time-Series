! SPDX-License-Identifier: MIT
! SPDX-FileComment: Command-line automatic modeling for a univariate series.
program ts_auto
   !! Read a numeric series, profile it, compare models, and forecast it.
   use kind_mod, only: dp
   use automatic_modeling_mod, only: automatic_model_options_t, &
      automatic_model_result_t, automatic_model, display
   use stats_mod, only: correlation_matrix
   use random_mod, only: random_uniform, set_random_seed
   use resampling_mod, only: resample
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: error_unit
   implicit none

   type :: table_read_t
      !! Numeric columns and optional row labels read from a text table.
      real(dp), allocatable :: values(:, :)
      character(len=128), allocatable :: names(:)
      character(len=128), allocatable :: row_labels(:)
      logical :: truncated = .false.
      integer :: error_line = 0
      integer :: info = 0
      character(len=256) :: error_message = ""
   end type table_read_t

   type :: csv_row_t
      !! Fields parsed from one comma-separated input line.
      character(len=:), allocatable :: fields(:)
   end type csv_row_t

   type :: bootstrap_plan_t
      !! Shared row indices for IID bootstrap replicates at one stride.
      integer, allocatable :: indices(:, :)
   end type bootstrap_plan_t

   type :: numeric_result_t
      !! Numeric values returned by a checked transformation.
      real(dp), allocatable :: values(:)
      integer :: info = 0
   end type numeric_result_t

   type(automatic_model_options_t) :: options
   type(automatic_model_result_t) :: result
   type(table_read_t) :: input
   type(numeric_result_t) :: modeled, restored, sampled
   character(len=256) :: option_error
   character(len=16) :: transformation
   character(len=1024) :: path
   integer :: clock_end, clock_max, clock_rate, clock_start
   integer :: display_lags, maximum_forecasts, maximum_models
   integer, allocatable :: strides(:)
   type(bootstrap_plan_t), allocatable :: bootstrap_plans(:)
   integer :: column_index, maximum_observations, overall_status
   integer :: replicate, resample_count, seed, status, stride_index
   real(dp) :: elapsed_seconds
   logical :: print_all_ar, print_all_arma, print_correlations
   logical :: print_parameters

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
      print_correlations, maximum_forecasts, strides, transformation, &
      resample_count, seed, status, option_error)
   if (status /= 0) then
      write(error_unit, '(a)') trim(option_error)
      call print_usage()
      stop 1
   end if
   input = read_table(trim(path), maximum_observations)
   if (input%info /= 0) then
      write(error_unit, '(a,i0,a,i0,a)') "Unable to read the input table; status ", &
         input%info, " at line ", input%error_line, "."
      if (len_trim(input%error_message) > 0) then
         write(error_unit, '(a)') trim(input%error_message)
      end if
      stop 1
   end if
   write(*, '(a,i0,a,i0,a)') "Read ", size(input%values, 1), &
      " valid rows and ", size(input%values, 2), " numeric columns."
   if (input%truncated) then
      write(*, '(a,i0,a)') "Input truncated at the requested limit of ", &
         maximum_observations, " observations."
   end if
   allocate(bootstrap_plans(size(strides)))
   if (resample_count > 0) then
      call initialize_bootstrap_plans(size(input%values, 1), strides, &
         transformation, resample_count, seed, bootstrap_plans)
      write(*, '(a,i0)') "IID bootstrap resamples: ", resample_count
      write(*, '(a,i0)') "Random-number seed: ", seed
   end if
   overall_status = 0
   do stride_index = 1, size(strides)
      if (stride_index > 1) write(*, '(a)') ""
      write(*, '(a)') "============================================================"
      write(*, '(a,i0)') "Stride analysis: ", strides(stride_index)
      do replicate = 1, max(1, resample_count)
         if (resample_count > 0) then
            write(*, '(a)') ""
            write(*, '(a,i0)') "IID resample: ", replicate
         end if
         if (print_correlations) then
            if (resample_count > 0) then
               call display_transformed_correlation(input%values, &
                  input%names, strides(stride_index), transformation, &
                  bootstrap_plans(stride_index)%indices(:, replicate))
            else
               call display_transformed_correlation(input%values, &
                  input%names, strides(stride_index), transformation)
            end if
         end if
         do column_index = 1, size(input%values, 2)
            write(*, '(a)') ""
            write(*, '(a)') "############################################################"
            write(*, '(a,i0,a,a)') "Column analysis: ", column_index, &
               " (", trim(input%names(column_index)) // ")"
            sampled = subsample_series(input%values(:, column_index), &
               strides(stride_index))
            write(*, '(a,i0)') "Retained observations: ", size(sampled%values)
            modeled = transform_series(sampled%values, transformation)
            if (modeled%info /= 0) then
               call report_transformation_error(modeled%info)
               overall_status = 1
               cycle
            end if
            if (resample_count > 0) then
               modeled%values = resample(modeled%values, &
                  bootstrap_plans(stride_index)%indices(:, replicate))
            end if
            write(*, '(a,a)') "Transformation: ", &
               transformation_label(transformation)
            write(*, '(a,i0)') "Modeling observations: ", size(modeled%values)
            write(*, '(a)') "Forecast scale: original observations"
            if (transformation == "log" .or. &
               transformation == "log-diff") then
               write(*, '(a)') &
                  "Back-transformed log forecasts are median forecasts on the original scale."
            end if
            result = automatic_model(modeled%values, options)
            if (result%info == 0) then
               restored = restore_forecasts(result%forecast, sampled%values, &
                  transformation)
               if (restored%info /= 0) then
                  write(error_unit, '(a,i0,a)') "Stride ", &
                     strides(stride_index), &
                     ": unable to restore forecasts to the original data scale."
                  overall_status = 1
                  cycle
               end if
               result%forecast = restored%values
            else
               overall_status = 1
            end if
            call display(result, display_lags, print_parameters, &
               maximum_models, print_all_ar, print_all_arma, &
               maximum_forecasts)
         end do
      end do
   end do
   call system_clock(clock_end)
   elapsed_seconds = elapsed_wall_time(clock_start, clock_end, &
      clock_rate, clock_max)
   write(*, '(a)') ""
   write(*, '(a,f12.3,a)') "Elapsed time: ", elapsed_seconds, " seconds"
   if (overall_status /= 0) stop 1

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
      print_correlations, maximum_forecasts, strides, transformation, &
      resample_count, seed, info, error_message)
      !! Parse command-line settings following the input path.
      type(automatic_model_options_t), intent(out) :: options !! Parsed modeling options.
      integer, intent(out) :: maximum_observations !! Maximum valid input values, or zero for all.
      integer, intent(out) :: display_lags !! Maximum positive correlation lag to print.
      integer, intent(out) :: maximum_models !! Maximum model summaries, or zero for all.
      logical, intent(out) :: print_parameters !! Print fitted model parameters when true.
      logical, intent(out) :: print_all_ar !! Print all tested autoregressive orders.
      logical, intent(out) :: print_all_arma !! Print all tested ARMA orders.
      logical, intent(out) :: print_correlations !! Print transformed correlation matrices.
      integer, intent(out) :: maximum_forecasts !! Maximum forecasts printed.
      integer, allocatable, intent(out) :: strides(:) !! Positive sampling strides.
      character(len=*), intent(out) :: transformation !! Requested input transformation.
      integer, intent(out) :: resample_count !! Number of IID bootstrap replicates.
      integer, intent(out) :: seed !! Random-number seed for bootstrap sampling.
      integer, intent(out) :: info !! Zero on successful parsing.
      character(len=*), intent(out) :: error_message !! Explanation when parsing fails.
      character(len=128) :: argument, value
      integer, allocatable :: stride_buffer(:)
      integer :: i, read_status, stride_count, stride_value
      logical :: transformation_set

      options = automatic_model_options_t()
      maximum_observations = 0
      display_lags = 5
      maximum_models = 0
      maximum_forecasts = huge(0)
      print_parameters = .false.
      print_all_ar = .false.
      print_all_arma = .false.
      print_correlations = .false.
      transformation = "none"
      resample_count = 0
      seed = 12345
      transformation_set = .false.
      allocate(stride_buffer(max(1, command_argument_count())))
      stride_count = 0
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
         if (trim(argument) == "--corr") then
            print_correlations = .true.
            i = i + 1
            cycle
         end if
         if (trim(argument) == "--resample") then
            resample_count = 1
            if (i < command_argument_count()) then
               call get_command_argument(i + 1, value)
               if (index(trim(value), "--") /= 1) then
                  read(value, *, iostat=read_status) resample_count
                  if (read_status /= 0 .or. resample_count < 1) then
                     info = 13
                     error_message = "Invalid resample count: " // &
                        trim(value) // ". The count must be a positive integer."
                     return
                  end if
                  i = i + 1
               end if
            end if
            i = i + 1
            cycle
         end if
         if (trim(argument) == "--log" .or. trim(argument) == "--diff" .or. &
            trim(argument) == "--log-diff") then
            call set_transformation(transformation, transformation_set, &
               trim(argument(3:)), info, error_message)
            if (info /= 0) return
            i = i + 1
            cycle
         end if
         if (trim(argument) == "--stride") then
            i = i + 1
            if (i > command_argument_count()) then
               info = 11
               error_message = "Option --stride requires at least one value."
               return
            end if
            call get_command_argument(i, value)
            if (index(trim(value), "--") == 1) then
               info = 11
               error_message = "Option --stride requires at least one value."
               return
            end if
            do while (i <= command_argument_count())
               call get_command_argument(i, value)
               if (index(trim(value), "--") == 1) exit
               read(value, *, iostat=read_status) stride_value
               if (read_status /= 0 .or. stride_value < 1) then
                  info = 12
                  error_message = "Invalid stride: " // trim(value) // &
                     ". Strides must be positive integers."
                  return
               end if
               stride_count = stride_count + 1
               stride_buffer(stride_count) = stride_value
               i = i + 1
            end do
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
         case ("--transform")
            call set_transformation(transformation, transformation_set, &
               lower_ascii(trim(value)), info, error_message)
            if (info /= 0) return
            read_status = 0
         case ("--observations", "--obs")
            read(value, *, iostat=read_status) maximum_observations
         case ("--display-lags")
            read(value, *, iostat=read_status) display_lags
         case ("--max-models")
            read(value, *, iostat=read_status) maximum_models
         case ("--print-forecasts")
            read(value, *, iostat=read_status) maximum_forecasts
         case ("--seed")
            read(value, *, iostat=read_status) seed
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
         if (maximum_forecasts < 0) then
            info = 10
            error_message = &
               "The number of printed forecasts cannot be negative."
            return
         end if
         if (seed < 0) then
            info = 14
            error_message = "The random-number seed cannot be negative."
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
      if (stride_count == 0) then
         allocate(strides(1))
         strides = 1
      else
         allocate(strides(stride_count))
         strides = stride_buffer(:stride_count)
      end if
   end subroutine parse_options

   pure logical function option_requires_value(argument) result(required)
      !! Identify recognized command options that require a following value.
      character(len=*), intent(in) :: argument !! Command option to classify.

      select case (trim(argument))
      case ("--frequency", "--horizon", "--max-lag", "--max-ar", &
         "--validation", "--selection", "--observations", "--obs", &
         "--display-lags", "--max-models", "--print-forecasts", &
         "--transform", "--seed")
         required = .true.
      case default
         required = .false.
      end select
   end function option_requires_value

   pure subroutine set_transformation(transformation, transformation_set, &
      requested, info, error_message)
      !! Validate one transformation request and reject conflicting modes.
      character(len=*), intent(inout) :: transformation !! Selected transformation name.
      logical, intent(inout) :: transformation_set !! Whether a mode was already requested.
      character(len=*), intent(in) :: requested !! Newly requested transformation name.
      integer, intent(out) :: info !! Zero when the request is valid and compatible.
      character(len=*), intent(out) :: error_message !! Explanation of an invalid request.

      info = 0
      error_message = ""
      select case (trim(requested))
      case ("none", "log", "diff", "log-diff")
      case default
         info = 8
         error_message = "Invalid transformation: " // trim(requested) // &
            ". Expected none, log, diff, or log-diff."
         return
      end select
      if (transformation_set .and. &
         trim(transformation) /= trim(requested)) then
         info = 9
         error_message = "Conflicting transformations: " // &
            trim(transformation) // " and " // trim(requested) // "."
         return
      end if
      transformation = trim(requested)
      transformation_set = .true.
   end subroutine set_transformation

   subroutine initialize_bootstrap_plans(input_rows, strides, transformation, &
      resample_count, seed, plans)
      !! Generate shared IID bootstrap row indices for every stride.
      integer, intent(in) :: input_rows !! Number of original data rows.
      integer, intent(in) :: strides(:) !! Requested sampling strides.
      character(len=*), intent(in) :: transformation !! Transformation name.
      integer, intent(in) :: resample_count !! Number of bootstrap replicates.
      integer, intent(in) :: seed !! Random-number seed.
      type(bootstrap_plan_t), intent(inout) :: plans(:) !! Generated index plans.
      integer :: draw, replicate, retained, rows, stride_index

      call set_random_seed(seed)
      do stride_index = 1, size(strides)
         retained = (input_rows - 1)/strides(stride_index) + 1
         rows = retained
         if (trim(transformation) == "diff" .or. &
            trim(transformation) == "log-diff") rows = rows - 1
         allocate(plans(stride_index)%indices(max(0, rows), resample_count))
         do replicate = 1, resample_count
            do draw = 1, rows
               plans(stride_index)%indices(draw, replicate) = &
                  min(rows, 1 + int(random_uniform()*real(rows, dp)))
            end do
         end do
      end do
   end subroutine initialize_bootstrap_plans

   subroutine display_transformed_correlation(values, names, stride, &
      transformation, indices)
      !! Print one transformed-data correlation matrix with optional resampling.
      real(dp), intent(in) :: values(:, :) !! Original numeric data columns.
      character(len=*), intent(in) :: names(:) !! Numeric column names.
      integer, intent(in) :: stride !! Sampling stride.
      character(len=*), intent(in) :: transformation !! Transformation name.
      integer, intent(in), optional :: indices(:) !! Shared IID bootstrap row indices.
      type(numeric_result_t) :: sampled_column, transformed_column
      real(dp), allocatable :: correlation(:, :), transformed(:, :)
      character(len=12) :: label
      integer :: column, failed_column, i, rows

      write(*, '(a)') ""
      write(*, '(a)') "Transformed-data correlation matrix"
      failed_column = 0
      rows = 0
      do column = 1, size(values, 2)
         sampled_column = subsample_series(values(:, column), stride)
         transformed_column = transform_series(sampled_column%values, &
            transformation)
         if (.not. allocated(transformed_column%values)) then
            failed_column = column
            exit
         end if
         if (size(transformed_column%values) < 2 .or. &
            transformed_column%info == 1 .or. &
            transformed_column%info == 2) then
            failed_column = column
            exit
         end if
         if (column == 1) then
            rows = size(transformed_column%values)
            allocate(transformed(rows, size(values, 2)))
         else if (size(transformed_column%values) /= rows) then
            failed_column = column
            exit
         end if
         transformed(:, column) = transformed_column%values
      end do
      if (failed_column > 0) then
         write(*, '(2x,a,a)') "unavailable: transformation failed for ", &
            trim(names(failed_column))
         return
      end if
      if (present(indices)) transformed = resample(transformed, indices)
      write(*, '(2x,a,i0)') "observations ", rows
      correlation = correlation_matrix(transformed)
      write(*, '(18x)', advance='no')
      do column = 1, size(names)
         label = names(column)
         write(*, '(1x,a12)', advance='no') label
      end do
      write(*, '(a)') ""
      do i = 1, size(names)
         label = names(i)
         write(*, '(2x,a12,4x)', advance='no') label
         do column = 1, size(names)
            write(*, '(1x,f12.6)', advance='no') correlation(i, column)
         end do
         write(*, '(a)') ""
      end do
   end subroutine display_transformed_correlation

   pure function subsample_series(series, stride) result(out)
      !! Retain the first observation and every stride-th observation after it.
      real(dp), intent(in) :: series(:) !! Original observations.
      integer, intent(in) :: stride !! Positive spacing between retained observations.
      type(numeric_result_t) :: out
      integer :: i, retained

      if (stride < 1 .or. size(series) < 1) then
         out%info = 1
         allocate(out%values(0))
         return
      end if
      retained = (size(series) - 1)/stride + 1
      allocate(out%values(retained))
      do i = 1, retained
         out%values(i) = series(1 + (i - 1)*stride)
      end do
   end function subsample_series

   subroutine report_transformation_error(info)
      !! Report why one strided series cannot be transformed and modeled.
      integer, intent(in) :: info !! Transformation status code.

      select case (info)
      case (1)
         write(error_unit, '(a)') &
            "Logarithmic transformations require strictly positive observations."
      case (2)
         write(error_unit, '(a)') &
            "Differencing requires at least two observations."
      case default
         write(error_unit, '(a)') &
            "The transformation leaves fewer than eight observations for modeling."
      end select
   end subroutine report_transformation_error

   pure function transform_series(series, transformation) result(out)
      !! Transform observations before profiling and model fitting.
      real(dp), intent(in) :: series(:) !! Original-scale observations.
      character(len=*), intent(in) :: transformation !! Transformation name.
      type(numeric_result_t) :: out
      real(dp), allocatable :: logged(:)

      select case (trim(transformation))
      case ("none")
         out%values = series
      case ("log")
         if (any(series <= 0.0_dp)) then
            out%info = 1
            return
         end if
         out%values = log(series)
      case ("diff")
         if (size(series) < 2) then
            out%info = 2
            return
         end if
         out%values = series(2:) - series(:size(series) - 1)
      case ("log-diff")
         if (any(series <= 0.0_dp)) then
            out%info = 1
            return
         end if
         if (size(series) < 2) then
            out%info = 2
            return
         end if
         logged = log(series)
         out%values = logged(2:) - logged(:size(logged) - 1)
      end select
      if (size(out%values) < 8) out%info = 3
   end function transform_series

   pure function restore_forecasts(forecast, series, transformation) result(out)
      !! Restore transformed forecasts to the original observation scale.
      real(dp), intent(in) :: forecast(:) !! Forecasts on the modeling scale.
      real(dp), intent(in) :: series(:) !! Original observations used as anchors.
      character(len=*), intent(in) :: transformation !! Transformation name.
      type(numeric_result_t) :: out
      real(dp) :: accumulated
      integer :: i

      allocate(out%values(size(forecast)))
      select case (trim(transformation))
      case ("none")
         out%values = forecast
      case ("log")
         out%values = exp(forecast)
      case ("diff")
         accumulated = series(size(series))
         do i = 1, size(forecast)
            accumulated = accumulated + forecast(i)
            out%values(i) = accumulated
         end do
      case ("log-diff")
         accumulated = log(series(size(series)))
         do i = 1, size(forecast)
            accumulated = accumulated + forecast(i)
            out%values(i) = exp(accumulated)
         end do
      end select
      if (any(.not. ieee_is_finite(out%values))) out%info = 1
   end function restore_forecasts

   pure function transformation_label(transformation) result(label)
      !! Return a readable description of an input transformation.
      character(len=*), intent(in) :: transformation !! Transformation name.
      character(len=:), allocatable :: label

      select case (trim(transformation))
      case ("none")
         label = "none"
      case ("log")
         label = "natural logarithm"
      case ("diff")
         label = "first difference"
      case ("log-diff")
         label = "first difference of natural logarithms"
      end select
   end function transformation_label

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

   function read_table(path, maximum_observations) result(out)
      !! Read numeric columns after an optional header and leading index fields.
      character(len=*), intent(in) :: path !! Input text or CSV path.
      integer, intent(in) :: maximum_observations !! Maximum data rows, or zero for all.
      type(table_read_t) :: out
      type(csv_row_t) :: first_row, header_row, row, sample_row
      type(numeric_result_t) :: numeric
      character(len=4096) :: line
      character(len=128) :: generated_name
      integer :: columns, count, data_start, first_line, i, io_status
      integer :: line_number, sample_line, unit
      logical :: has_header

      open(newunit=unit, file=path, status="old", action="read", iostat=io_status)
      if (io_status /= 0) then
         out%info = 1
         out%error_message = "The input file could not be opened."
         return
      end if
      call read_next_nonempty_line(unit, line, first_line, io_status)
      if (io_status /= 0) then
         close(unit)
         out%info = 3
         out%error_message = "The input file contains no data."
         return
      end if
      first_row = split_csv_line(line)
      has_header = .not. row_has_numeric_data(first_row)
      if (has_header) then
         header_row = first_row
         call read_next_nonempty_line(unit, line, sample_line, io_status, &
            first_line)
         if (io_status /= 0) then
            close(unit)
            out%info = 3
            out%error_line = first_line
            out%error_message = "A header was found without any data rows."
            return
         end if
         sample_row = split_csv_line(line)
         if (size(header_row%fields) /= size(sample_row%fields)) then
            close(unit)
            out%info = 4
            out%error_line = sample_line
            out%error_message = &
               "The header and first data row have different field counts."
            return
         end if
      else
         sample_row = first_row
         sample_line = first_line
      end if
      data_start = first_numeric_field(sample_row)
      if (data_start == 0 .or. &
         .not. numeric_fields_are_valid(sample_row, data_start)) then
         close(unit)
         out%info = 5
         out%error_line = sample_line
         out%error_message = &
            "No contiguous numeric data columns were found after the index fields."
         return
      end if
      columns = size(sample_row%fields) - data_start + 1
      allocate(out%names(columns))
      out%names = ""
      do i = 1, columns
         if (has_header) then
            out%names(i) = trim(header_row%fields(data_start + i - 1))
         end if
         if (len_trim(out%names(i)) == 0) then
            write(generated_name, '(a,i0)') "Column ", i
            out%names(i) = generated_name
         end if
      end do

      rewind(unit)
      count = 0
      line_number = 0
      do
         read(unit, '(a)', iostat=io_status) line
         if (io_status < 0) exit
         line_number = line_number + 1
         if (io_status > 0) then
            close(unit)
            out%info = 2
            out%error_line = line_number
            out%error_message = "The input file could not be read."
            return
         end if
         if (len_trim(line) == 0) cycle
         if (has_header .and. line_number == first_line) cycle
         row = split_csv_line(line)
         if (size(row%fields) /= size(sample_row%fields)) then
            close(unit)
            out%info = 6
            out%error_line = line_number
            out%error_message = "The row has an inconsistent number of fields."
            return
         end if
         numeric = numeric_row_values(row, data_start)
         if (numeric%info /= 0) then
            close(unit)
            out%info = 7
            out%error_line = line_number
            out%error_message = &
               "A numeric data column is missing, nonnumeric, or nonfinite."
            return
         end if
         if (maximum_observations > 0 .and. &
            count >= maximum_observations) then
            out%truncated = .true.
            exit
         end if
         count = count + 1
      end do
      if (count == 0) then
         close(unit)
         out%info = 8
         out%error_message = "No valid numeric data rows were found."
         return
      end if

      allocate(out%values(count, columns), out%row_labels(count))
      out%row_labels = ""
      rewind(unit)
      count = 0
      line_number = 0
      do
         read(unit, '(a)', iostat=io_status) line
         if (io_status /= 0) exit
         line_number = line_number + 1
         if (len_trim(line) == 0) cycle
         if (has_header .and. line_number == first_line) cycle
         row = split_csv_line(line)
         numeric = numeric_row_values(row, data_start)
         count = count + 1
         out%values(count, :) = numeric%values
         if (data_start > 1) out%row_labels(count) = trim(row%fields(1))
         if (count == size(out%values, 1)) exit
      end do
      close(unit)
   end function read_table

   subroutine read_next_nonempty_line(unit, line, line_number, info, &
      starting_line)
      !! Read the next nonempty line and return its physical line number.
      integer, intent(in) :: unit !! Connected input unit.
      character(len=*), intent(out) :: line !! Next nonempty line.
      integer, intent(out) :: line_number !! Physical line number read.
      integer, intent(out) :: info !! Zero on success and nonzero at end or error.
      integer, intent(in), optional :: starting_line !! Previously consumed line count.

      line_number = 0
      if (present(starting_line)) line_number = starting_line
      do
         read(unit, '(a)', iostat=info) line
         if (info /= 0) return
         line_number = line_number + 1
         if (len_trim(line) > 0) return
      end do
   end subroutine read_next_nonempty_line

   pure function split_csv_line(line) result(out)
      !! Split a CSV record while respecting quoted fields and escaped quotes.
      character(len=*), intent(in) :: line !! One comma-separated record.
      type(csv_row_t) :: out
      integer :: field, field_count, i, length, position
      logical :: quoted

      length = len_trim(line)
      field_count = 1
      quoted = .false.
      i = 1
      do while (i <= length)
         if (line(i:i) == '"') then
            if (quoted .and. i < length .and. line(i + 1:i + 1) == '"') then
               i = i + 2
               cycle
            end if
            quoted = .not. quoted
         else if (line(i:i) == "," .and. .not. quoted) then
            field_count = field_count + 1
         end if
         i = i + 1
      end do
      allocate(character(len=max(1, length)) :: out%fields(field_count))
      out%fields = ""
      field = 1
      position = 0
      quoted = .false.
      i = 1
      do while (i <= length)
         if (line(i:i) == '"') then
            if (quoted .and. i < length .and. line(i + 1:i + 1) == '"') then
               position = position + 1
               out%fields(field)(position:position) = '"'
               i = i + 2
               cycle
            end if
            quoted = .not. quoted
         else if (line(i:i) == "," .and. .not. quoted) then
            out%fields(field) = adjustl(out%fields(field))
            field = field + 1
            position = 0
         else
            position = position + 1
            out%fields(field)(position:position) = line(i:i)
         end if
         i = i + 1
      end do
      out%fields(field) = adjustl(out%fields(field))
   end function split_csv_line

   logical function row_has_numeric_data(row) result(has_data)
      !! Identify a row with leading labels followed by numeric fields.
      type(csv_row_t), intent(in) :: row !! Parsed input row.
      integer :: first

      first = first_numeric_field(row)
      has_data = first > 0
      if (has_data) has_data = numeric_fields_are_valid(row, first)
   end function row_has_numeric_data

   integer function first_numeric_field(row) result(first)
      !! Locate the first finite numeric field in a parsed row.
      type(csv_row_t), intent(in) :: row !! Parsed input row.
      real(dp) :: value
      integer :: i, read_status

      first = 0
      do i = 1, size(row%fields)
         read(row%fields(i), *, iostat=read_status) value
         if (read_status == 0 .and. ieee_is_finite(value)) then
            first = i
            return
         end if
      end do
   end function first_numeric_field

   logical function numeric_fields_are_valid(row, first) result(valid)
      !! Check that all fields from a starting column are finite numbers.
      type(csv_row_t), intent(in) :: row !! Parsed input row.
      integer, intent(in) :: first !! First numeric data field.
      real(dp) :: value
      integer :: i, read_status

      valid = first >= 1 .and. first <= size(row%fields)
      if (.not. valid) return
      do i = first, size(row%fields)
         read(row%fields(i), *, iostat=read_status) value
         if (read_status /= 0 .or. .not. ieee_is_finite(value)) then
            valid = .false.
            return
         end if
      end do
   end function numeric_fields_are_valid

   function numeric_row_values(row, first) result(out)
      !! Convert the numeric suffix of a parsed row to real values.
      type(csv_row_t), intent(in) :: row !! Parsed input row.
      integer, intent(in) :: first !! First numeric data field.
      type(numeric_result_t) :: out
      integer :: i, read_status

      allocate(out%values(size(row%fields) - first + 1))
      do i = first, size(row%fields)
         read(row%fields(i), *, iostat=read_status) out%values(i - first + 1)
         if (read_status /= 0 .or. &
            .not. ieee_is_finite(out%values(i - first + 1))) then
            out%info = 1
            return
         end if
      end do
   end function numeric_row_values

   subroutine print_usage()
      !! Display command-line syntax and supported options.

      write(*, '(a)') "Usage: ts_auto FILE [options]"
      write(*, '(a)') ""
      write(*, '(a)') &
         "FILE contains numeric columns with an optional header and leading index column."
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
      write(*, '(a)') "  --corr            print transformed-data correlation matrices"
      write(*, '(a)') "  --resample [N]    analyze N IID row-bootstrap samples (default 1)"
      write(*, '(a)') "  --seed N          bootstrap random-number seed (default 12345)"
      write(*, '(a)') "  --transform NAME  none, log, diff, or log-diff (default none)"
      write(*, '(a)') "  --log             alias for --transform log"
      write(*, '(a)') "  --diff            alias for --transform diff"
      write(*, '(a)') "  --log-diff        alias for --transform log-diff"
      write(*, '(a)') "  --stride N...     analyze the series separately at each stride"
      write(*, '(a)') "  --max-models N    maximum model summaries; zero displays all"
      write(*, '(a)') "  --print-forecasts N maximum forecasts printed; zero prints none"
      write(*, '(a)') "  --observations N use at most the first N valid observations"
      write(*, '(a)') "  --obs N          alias for --observations"
   end subroutine print_usage

end program ts_auto
