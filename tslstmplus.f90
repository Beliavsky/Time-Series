! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R TSLSTMplus package.
module tslstmplus_mod
   !! LSTM time-series regression with lagged responses and exogenous inputs.
   use kind_mod, only: dp
   use recurrent_network_mod, only: recurrent_lstm_model_t, &
      recurrent_lstm_fit, recurrent_lstm_predict, &
      recurrent_activation_linear, recurrent_activation_sigmoid, &
      recurrent_activation_relu, recurrent_activation_tanh, &
      recurrent_optimizer_sgd, recurrent_optimizer_adam, &
      recurrent_optimizer_rmsprop, recurrent_loss_mse, recurrent_loss_mae
   use utils_mod, only: quiet_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: tslstmplus_scale_none = 0
   integer, parameter, public :: tslstmplus_scale_standard = 1
   integer, parameter, public :: tslstmplus_scale_minmax = 2
   integer, parameter, public :: tslstmplus_activation_linear = &
      recurrent_activation_linear
   integer, parameter, public :: tslstmplus_activation_sigmoid = &
      recurrent_activation_sigmoid
   integer, parameter, public :: tslstmplus_activation_relu = &
      recurrent_activation_relu
   integer, parameter, public :: tslstmplus_activation_tanh = &
      recurrent_activation_tanh
   integer, parameter, public :: tslstmplus_optimizer_sgd = &
      recurrent_optimizer_sgd
   integer, parameter, public :: tslstmplus_optimizer_adam = &
      recurrent_optimizer_adam
   integer, parameter, public :: tslstmplus_optimizer_rmsprop = &
      recurrent_optimizer_rmsprop
   integer, parameter, public :: tslstmplus_loss_mse = recurrent_loss_mse
   integer, parameter, public :: tslstmplus_loss_mae = recurrent_loss_mae

   type, public :: tslstmplus_scale_t
      !! Retained columnwise location and scale transformation.
      real(dp), allocatable :: location(:)
      real(dp), allocatable :: scale(:)
      integer :: method = tslstmplus_scale_none
   end type tslstmplus_scale_t

   type, public :: tslstmplus_prepared_t
      !! Prepared response targets and sample-by-time-by-feature inputs.
      real(dp), allocatable :: sequences(:, :, :)
      real(dp), allocatable :: target(:)
      integer, allocatable :: xreg_lags(:)
      integer :: response_lag = 0
      integer :: maximum_lag = 0
      integer :: features = 0
      integer :: info = 0
      logical :: lags_as_sequences = .true.
   end type tslstmplus_prepared_t

   type, public :: tslstmplus_model_t
      !! Fitted TSLSTMplus recurrent forecasting model.
      type(recurrent_lstm_model_t) :: network
      type(tslstmplus_scale_t) :: output_scale
      type(tslstmplus_scale_t) :: input_scale
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      integer, allocatable :: xreg_lags(:)
      integer :: response_lag = 0
      integer :: maximum_lag = 0
      integer :: regressors = 0
      integer :: batch_size = 1
      integer :: info = 0
      logical :: lags_as_sequences = .true.
      logical :: stateful = .false.
   end type tslstmplus_model_t

   interface display
      module procedure display_tslstmplus_model
   end interface display

   public :: tslstmplus_prepare_data, tslstmplus_prepare_sequences
   public :: tslstmplus_fit, tslstmplus_predict, tslstmplus_forecast
   public :: display

contains

   pure function tslstmplus_prepare_data(series, xreg, response_lag, &
      xreg_lags) result(data)
      !! Construct complete current-response and lagged-feature rows.
      real(dp), intent(in) :: series(:) !! Response observations.
      real(dp), intent(in), optional :: xreg(:, :) !! Exogenous observations.
      integer, intent(in), optional :: response_lag !! Response lag order.
      integer, intent(in), optional :: xreg_lags(:) !! Exogenous lag orders.
      real(dp), allocatable :: data(:, :)
      integer, allocatable :: selected_lags(:)
      integer :: p, variables, maximum_lag, columns, row, time, lag, variable
      integer :: column

      p = 0
      if (present(response_lag)) p = response_lag
      variables = 0
      if (present(xreg)) variables = size(xreg, 2)
      call resolve_xreg_lags(variables, xreg_lags, selected_lags)
      if (size(series) < 1 .or. p < 0 .or. any(selected_lags < 0) .or. &
         (present(xreg) .and. size(xreg, 1) /= size(series))) then
         allocate(data(0, 0))
         return
      end if
      maximum_lag = p
      if (variables > 0) maximum_lag = max(maximum_lag, maxval(selected_lags))
      if (size(series) <= maximum_lag) then
         allocate(data(0, 0))
         return
      end if
      columns = 1 + p
      if (variables > 0) columns = columns + sum(selected_lags + 1)
      allocate(data(size(series) - maximum_lag, columns))
      do row = 1, size(data, 1)
         time = maximum_lag + row
         data(row, 1) = series(time)
         column = 1
         do lag = 1, p
            column = column + 1
            data(row, column) = series(time - lag)
         end do
         do variable = 1, variables
            do lag = 0, selected_lags(variable)
               column = column + 1
               data(row, column) = xreg(time - lag, variable)
            end do
         end do
      end do
   end function tslstmplus_prepare_data

   pure function tslstmplus_prepare_sequences(series, xreg, response_lag, &
      xreg_lags, lags_as_sequences) result(prepared)
      !! Prepare LSTM tensors using lags as timesteps or flattened features.
      real(dp), intent(in) :: series(:) !! Response observations.
      real(dp), intent(in), optional :: xreg(:, :) !! Exogenous observations.
      integer, intent(in), optional :: response_lag !! Response lag order.
      integer, intent(in), optional :: xreg_lags(:) !! Exogenous lag orders.
      logical, intent(in), optional :: lags_as_sequences !! Treat lags as timesteps.
      type(tslstmplus_prepared_t) :: prepared
      integer, allocatable :: selected_lags(:)
      integer :: p, variables, maximum_lag, samples, timesteps, features
      integer :: sample, time, step, variable, lag, column
      logical :: use_sequences

      p = 0
      if (present(response_lag)) p = response_lag
      variables = 0
      if (present(xreg)) variables = size(xreg, 2)
      call resolve_xreg_lags(variables, xreg_lags, selected_lags)
      use_sequences = .true.
      if (present(lags_as_sequences)) use_sequences = lags_as_sequences
      prepared%response_lag = p
      prepared%xreg_lags = selected_lags
      prepared%lags_as_sequences = use_sequences
      if (size(series) < 1 .or. p < 0 .or. any(selected_lags < 0) .or. &
         (present(xreg) .and. size(xreg, 1) /= size(series)) .or. &
         .not. all(ieee_is_finite(series))) then
         prepared%info = 1
         return
      end if
      if (present(xreg)) then
         if (.not. all(ieee_is_finite(xreg))) then
            prepared%info = 1
            return
         end if
      end if
      maximum_lag = p
      if (variables > 0) maximum_lag = max(maximum_lag, maxval(selected_lags))
      if (p == 0 .and. variables == 0) then
         prepared%info = 2
         return
      end if
      samples = size(series) - maximum_lag
      if (samples < 2) then
         prepared%info = 3
         return
      end if
      prepared%maximum_lag = maximum_lag
      allocate(prepared%target(samples))
      prepared%target = series(maximum_lag + 1:)
      if (use_sequences) then
         timesteps = max(1, maximum_lag)
         features = variables + merge(1, 0, p > 0)
         prepared%features = features
         allocate(prepared%sequences(samples, timesteps, features))
         do sample = 1, samples
            time = maximum_lag + sample
            column = 0
            if (p > 0) then
               column = column + 1
               do step = 1, timesteps
                  prepared%sequences(sample, step, column) = &
                     series(time - timesteps + step - 1)
               end do
            end if
            do variable = 1, variables
               column = column + 1
               do step = 1, timesteps
                  prepared%sequences(sample, step, column) = &
                     xreg(time - timesteps + step, variable)
               end do
            end do
         end do
      else
         features = p
         if (variables > 0) features = features + sum(selected_lags + 1)
         prepared%features = features
         allocate(prepared%sequences(samples, 1, features))
         do sample = 1, samples
            time = maximum_lag + sample
            column = 0
            do lag = 1, p
               column = column + 1
               prepared%sequences(sample, 1, column) = series(time - lag)
            end do
            do variable = 1, variables
               do lag = 0, selected_lags(variable)
                  column = column + 1
                  prepared%sequences(sample, 1, column) = &
                     xreg(time - lag, variable)
               end do
            end do
         end do
      end if
   end function tslstmplus_prepare_sequences

   pure function tslstmplus_fit(series, lstm_units, xreg, response_lag, &
      xreg_lags, dense_units, dropout_rate, epochs, learning_rate, optimizer, &
      loss, scale_output, scale_input, batch_size, lstm_activation, &
      recurrent_activation, dense_activation, validation_split, patience, &
      minimum_delta, lags_as_sequences, stateful, seed) result(model)
      !! Fit a configurable stand-alone TSLSTMplus forecasting model.
      real(dp), intent(in) :: series(:) !! Response observations.
      integer, intent(in) :: lstm_units(:) !! Units in successive LSTM layers.
      real(dp), intent(in), optional :: xreg(:, :) !! Exogenous observations.
      integer, intent(in), optional :: response_lag !! Response lag order.
      integer, intent(in), optional :: xreg_lags(:) !! Exogenous lag orders.
      integer, intent(in), optional :: dense_units(:) !! Hidden dense-layer sizes.
      real(dp), intent(in), optional :: dropout_rate !! LSTM input dropout fraction.
      integer, intent(in), optional :: epochs !! Maximum training epochs.
      real(dp), intent(in), optional :: learning_rate !! Optimizer learning rate.
      integer, intent(in), optional :: optimizer !! SGD, Adam, or RMSprop code.
      integer, intent(in), optional :: loss !! MSE or MAE loss code.
      integer, intent(in), optional :: scale_output !! Response scaling code.
      integer, intent(in), optional :: scale_input !! Exogenous scaling code.
      integer, intent(in), optional :: batch_size !! Mini-batch size.
      integer, intent(in), optional :: lstm_activation !! Candidate and cell activation.
      integer, intent(in), optional :: recurrent_activation !! Gate activation.
      integer, intent(in), optional :: dense_activation !! Hidden dense activation.
      real(dp), intent(in), optional :: validation_split !! Terminal validation fraction.
      integer, intent(in), optional :: patience !! Early-stopping patience.
      real(dp), intent(in), optional :: minimum_delta !! Minimum validation improvement.
      logical, intent(in), optional :: lags_as_sequences !! Treat lags as timesteps.
      logical, intent(in), optional :: stateful !! Carry states between samples.
      integer, intent(in), optional :: seed !! Deterministic initialization seed.
      type(tslstmplus_model_t) :: model
      type(tslstmplus_prepared_t) :: prepared
      real(dp), allocatable :: scaled_series(:), scaled_xreg(:, :)
      real(dp), allocatable :: targets(:, :), prediction(:, :)
      integer, allocatable :: selected_dense(:)
      integer :: p, input_method, output_method
      logical :: sequence_mode

      p = 0
      if (present(response_lag)) p = response_lag
      input_method = tslstmplus_scale_none
      output_method = tslstmplus_scale_none
      if (present(scale_input)) input_method = scale_input
      if (present(scale_output)) output_method = scale_output
      sequence_mode = .true.
      if (present(lags_as_sequences)) sequence_mode = lags_as_sequences
      model%stateful = .false.
      if (present(stateful)) model%stateful = stateful
      if (model%stateful) sequence_mode = .true.
      model%lags_as_sequences = sequence_mode
      model%response_lag = p
      model%batch_size = 1
      if (present(batch_size)) model%batch_size = batch_size
      model%regressors = 0
      if (present(xreg)) model%regressors = size(xreg, 2)
      if (size(series) < 3 .or. .not. all(ieee_is_finite(series)) .or. &
         input_method < tslstmplus_scale_none .or. &
         input_method > tslstmplus_scale_minmax .or. &
         output_method < tslstmplus_scale_none .or. &
         output_method > tslstmplus_scale_minmax) then
         model%info = 1
         return
      end if
      call fit_scale_vector(series, output_method, scaled_series, &
         model%output_scale)
      if (present(xreg)) then
         if (size(xreg, 1) /= size(series) .or. &
            .not. all(ieee_is_finite(xreg))) then
            model%info = 1
            return
         end if
         call fit_scale_matrix(xreg, input_method, scaled_xreg, &
            model%input_scale)
         prepared = tslstmplus_prepare_sequences(scaled_series, scaled_xreg, &
            p, xreg_lags, sequence_mode)
      else
         allocate(scaled_xreg(size(series), 0))
         call empty_scale(input_method, model%input_scale)
         prepared = tslstmplus_prepare_sequences(scaled_series, &
            response_lag=p, xreg_lags=xreg_lags, &
            lags_as_sequences=sequence_mode)
      end if
      if (prepared%info /= 0) then
         model%info = 10 + prepared%info
         return
      end if
      model%xreg_lags = prepared%xreg_lags
      model%maximum_lag = prepared%maximum_lag
      allocate(targets(size(prepared%target), 1))
      targets(:, 1) = prepared%target
      allocate(selected_dense(0))
      if (present(dense_units)) selected_dense = dense_units
      model%network = recurrent_lstm_fit(prepared%sequences, targets, &
         lstm_units, selected_dense, epochs, model%batch_size, learning_rate, &
         optimizer, loss, dropout_rate, dense_activation, validation_split, &
         patience, minimum_delta, model%stateful, seed, lstm_activation, &
         recurrent_activation)
      if (model%network%info /= 0) then
         model%info = 20 + model%network%info
         return
      end if
      prediction = recurrent_lstm_predict(model%network, prepared%sequences)
      allocate(model%fitted(size(series)), model%residuals(size(series)))
      model%fitted = quiet_nan()
      model%residuals = quiet_nan()
      model%fitted(model%maximum_lag + 1:) = &
         undo_scale_vector(prediction(:, 1), model%output_scale)
      model%residuals(model%maximum_lag + 1:) = &
         series(model%maximum_lag + 1:) - &
         model%fitted(model%maximum_lag + 1:)
   end function tslstmplus_fit

   pure function tslstmplus_predict(model, series, xreg) result(prediction)
      !! Return rolling one-step predictions from a fitted LSTM model.
      type(tslstmplus_model_t), intent(in) :: model !! Fitted TSLSTMplus model.
      real(dp), intent(in) :: series(:) !! Response history.
      real(dp), intent(in), optional :: xreg(:, :) !! Exogenous history.
      real(dp), allocatable :: prediction(:)
      type(tslstmplus_prepared_t) :: prepared
      real(dp), allocatable :: scaled_series(:), scaled_xreg(:, :)
      real(dp), allocatable :: normalized(:, :)

      if (model%info /= 0 .or. .not. all(ieee_is_finite(series))) then
         allocate(prediction(0))
         return
      end if
      scaled_series = apply_scale_vector(series, model%output_scale)
      if (model%regressors > 0) then
         if (.not. present(xreg)) then
            allocate(prediction(0))
            return
         end if
         if (size(xreg, 1) /= size(series) .or. &
            size(xreg, 2) /= model%regressors .or. &
            .not. all(ieee_is_finite(xreg))) then
            allocate(prediction(0))
            return
         end if
         scaled_xreg = apply_scale_matrix(xreg, model%input_scale)
         prepared = tslstmplus_prepare_sequences(scaled_series, scaled_xreg, &
            model%response_lag, model%xreg_lags, model%lags_as_sequences)
      else
         prepared = tslstmplus_prepare_sequences(scaled_series, &
            response_lag=model%response_lag, xreg_lags=model%xreg_lags, &
            lags_as_sequences=model%lags_as_sequences)
      end if
      if (prepared%info /= 0) then
         allocate(prediction(0))
         return
      end if
      normalized = recurrent_lstm_predict(model%network, prepared%sequences)
      prediction = undo_scale_vector(normalized(:, 1), model%output_scale)
   end function tslstmplus_predict

   pure function tslstmplus_forecast(model, series, horizon, xreg, xreg_new) &
      result(forecast)
      !! Recursively forecast with optional known future exogenous variables.
      type(tslstmplus_model_t), intent(in) :: model !! Fitted TSLSTMplus model.
      real(dp), intent(in) :: series(:) !! Response history.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: xreg(:, :) !! Historical exogenous values.
      real(dp), intent(in), optional :: xreg_new(:, :) !! Future exogenous values.
      real(dp), allocatable :: forecast(:)
      type(tslstmplus_prepared_t) :: prepared
      real(dp), allocatable :: working_series(:), combined_xreg(:, :)
      real(dp), allocatable :: normalized(:, :), scaled_history(:)
      integer :: step, observations

      if (model%info /= 0 .or. horizon < 1 .or. &
         .not. all(ieee_is_finite(series))) then
         allocate(forecast(0))
         return
      end if
      observations = size(series)
      allocate(working_series(observations + horizon), source=0.0_dp)
      working_series(:observations) = apply_scale_vector(series, &
         model%output_scale)
      if (model%regressors > 0) then
         if (.not. present(xreg) .or. .not. present(xreg_new)) then
            allocate(forecast(0))
            return
         end if
         if (size(xreg, 1) /= observations .or. &
            size(xreg, 2) /= model%regressors .or. &
            size(xreg_new, 1) < horizon .or. &
            size(xreg_new, 2) /= model%regressors .or. &
            .not. all(ieee_is_finite(xreg)) .or. &
            .not. all(ieee_is_finite(xreg_new(:horizon, :)))) then
            allocate(forecast(0))
            return
         end if
         allocate(combined_xreg(observations + horizon, model%regressors))
         combined_xreg(:observations, :) = apply_scale_matrix(xreg, &
            model%input_scale)
         combined_xreg(observations + 1:, :) = apply_scale_matrix( &
            xreg_new(:horizon, :), model%input_scale)
      else
         allocate(combined_xreg(observations + horizon, 0))
      end if
      allocate(forecast(horizon))
      do step = 1, horizon
         if (model%regressors > 0) then
            prepared = tslstmplus_prepare_sequences( &
               working_series(:observations + step), &
               combined_xreg(:observations + step, :), model%response_lag, &
               model%xreg_lags, model%lags_as_sequences)
         else
            prepared = tslstmplus_prepare_sequences( &
               working_series(:observations + step), &
               response_lag=model%response_lag, xreg_lags=model%xreg_lags, &
               lags_as_sequences=model%lags_as_sequences)
         end if
         if (prepared%info /= 0) then
            deallocate(forecast)
            allocate(forecast(0))
            return
         end if
         normalized = recurrent_lstm_predict(model%network, prepared%sequences)
         working_series(observations + step) = &
            normalized(size(normalized, 1), 1)
         scaled_history = undo_scale_vector( &
            [working_series(observations + step)], model%output_scale)
         forecast(step) = scaled_history(1)
      end do
   end function tslstmplus_forecast

   subroutine display_tslstmplus_model(model, unit, print_obs)
      !! Display an LSTM fit and optionally its aligned observations.
      type(tslstmplus_model_t), intent(in) :: model !! Fitted TSLSTMplus model.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Print fitted values and residuals.
      integer :: destination, layer, observation
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'TSLSTMplus fit'
      write(destination, '(a,i0)') '  response lag: ', model%response_lag
      write(destination, '(a,i0)') '  regressors: ', model%regressors
      write(destination, '(a,l1)') '  lags as sequences: ', &
         model%lags_as_sequences
      write(destination, '(a,l1)') '  stateful: ', model%stateful
      write(destination, '(a,i0)') '  epochs trained: ', &
         model%network%epochs_trained
      do layer = 1, size(model%network%lstm_layers)
         write(destination, '(a,i0,a,i0)') '  LSTM layer ', layer, &
            ' units: ', model%network%lstm_layers(layer)%units
      end do
      if (show_observations) then
         write(destination, '(a)') '  fitted residual'
         do observation = 1, size(model%fitted)
            write(destination, '(2(es14.6,1x))') model%fitted(observation), &
               model%residuals(observation)
         end do
      end if
   end subroutine display_tslstmplus_model

   pure subroutine resolve_xreg_lags(variables, supplied, selected)
      !! Expand one supplied exogenous lag or retain per-variable lags.
      integer, intent(in) :: variables !! Number of exogenous variables.
      integer, intent(in), optional :: supplied(:) !! Supplied lag orders.
      integer, allocatable, intent(out) :: selected(:) !! Expanded lag orders.

      allocate(selected(variables), source=0)
      if (.not. present(supplied) .or. variables == 0) return
      if (size(supplied) == 1) then
         selected = supplied(1)
      else if (size(supplied) == variables) then
         selected = supplied
      else
         selected = -1
      end if
   end subroutine resolve_xreg_lags

   pure subroutine fit_scale_vector(values, method, transformed, scale)
      !! Fit and apply a response scaling transformation.
      real(dp), intent(in) :: values(:) !! Original values.
      integer, intent(in) :: method !! Scaling method code.
      real(dp), allocatable, intent(out) :: transformed(:) !! Scaled values.
      type(tslstmplus_scale_t), intent(out) :: scale !! Retained transformation.
      real(dp) :: center, spread_value

      scale%method = method
      allocate(scale%location(1), scale%scale(1))
      center = 0.0_dp
      spread_value = 1.0_dp
      if (method == tslstmplus_scale_standard) then
         center = sum(values)/real(size(values), dp)
         spread_value = sqrt(sum((values - center)**2)/ &
            real(max(1, size(values) - 1), dp))
      else if (method == tslstmplus_scale_minmax) then
         center = minval(values)
         spread_value = maxval(values) - center
      end if
      spread_value = max(spread_value, sqrt(epsilon(1.0_dp)))
      scale%location(1) = center
      scale%scale(1) = spread_value
      transformed = (values - center)/spread_value
   end subroutine fit_scale_vector

   pure subroutine fit_scale_matrix(values, method, transformed, scale)
      !! Fit and apply columnwise exogenous scaling transformations.
      real(dp), intent(in) :: values(:, :) !! Original exogenous values.
      integer, intent(in) :: method !! Scaling method code.
      real(dp), allocatable, intent(out) :: transformed(:, :) !! Scaled values.
      type(tslstmplus_scale_t), intent(out) :: scale !! Retained transformations.
      integer :: variable

      scale%method = method
      allocate(scale%location(size(values, 2)), scale%scale(size(values, 2)))
      allocate(transformed, source=values)
      scale%location = 0.0_dp
      scale%scale = 1.0_dp
      do variable = 1, size(values, 2)
         if (method == tslstmplus_scale_standard) then
            scale%location(variable) = sum(values(:, variable))/ &
               real(size(values, 1), dp)
            scale%scale(variable) = sqrt(sum((values(:, variable) - &
               scale%location(variable))**2)/ &
               real(max(1, size(values, 1) - 1), dp))
         else if (method == tslstmplus_scale_minmax) then
            scale%location(variable) = minval(values(:, variable))
            scale%scale(variable) = maxval(values(:, variable)) - &
               scale%location(variable)
         end if
         scale%scale(variable) = max(scale%scale(variable), &
            sqrt(epsilon(1.0_dp)))
         transformed(:, variable) = (values(:, variable) - &
            scale%location(variable))/scale%scale(variable)
      end do
   end subroutine fit_scale_matrix

   pure subroutine empty_scale(method, scale)
      !! Construct an empty retained scaling transformation.
      integer, intent(in) :: method !! Scaling method code.
      type(tslstmplus_scale_t), intent(out) :: scale !! Empty transformation.

      scale%method = method
      allocate(scale%location(0), scale%scale(0))
   end subroutine empty_scale

   pure function apply_scale_vector(values, scale) result(transformed)
      !! Apply a retained univariate scaling transformation.
      real(dp), intent(in) :: values(:) !! Original values.
      type(tslstmplus_scale_t), intent(in) :: scale !! Retained transformation.
      real(dp) :: transformed(size(values))

      transformed = (values - scale%location(1))/scale%scale(1)
   end function apply_scale_vector

   pure function undo_scale_vector(values, scale) result(transformed)
      !! Undo a retained univariate scaling transformation.
      real(dp), intent(in) :: values(:) !! Scaled values.
      type(tslstmplus_scale_t), intent(in) :: scale !! Retained transformation.
      real(dp) :: transformed(size(values))

      transformed = values*scale%scale(1) + scale%location(1)
   end function undo_scale_vector

   pure function apply_scale_matrix(values, scale) result(transformed)
      !! Apply retained columnwise exogenous scaling transformations.
      real(dp), intent(in) :: values(:, :) !! Original values.
      type(tslstmplus_scale_t), intent(in) :: scale !! Retained transformations.
      real(dp) :: transformed(size(values, 1), size(values, 2))
      integer :: variable

      do variable = 1, size(values, 2)
         transformed(:, variable) = (values(:, variable) - &
            scale%location(variable))/scale%scale(variable)
      end do
   end function apply_scale_matrix

end module tslstmplus_mod
