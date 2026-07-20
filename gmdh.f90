! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R GMDH package.
module gmdh_mod
   !! Group-method-of-data-handling networks for short-term forecasting.
   use kind_mod, only: dp
   use linalg_mod, only: symmetric_pseudoinverse
   use stats_mod, only: normal_quantile
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: gmdh_method_classic = 1
   integer, parameter, public :: gmdh_method_revised = 2
   integer, parameter, public :: gmdh_transfer_polynomial = 101
   integer, parameter, public :: gmdh_transfer_sigmoid = 102
   integer, parameter, public :: gmdh_transfer_rbf = 103
   integer, parameter, public :: gmdh_transfer_tangent = 104
   integer, parameter, public :: gmdh_node_quadratic = 1
   integer, parameter, public :: gmdh_node_feedback = 2

   type, public :: gmdh_ridge_fit_t
      !! Validation-selected ridge-regression result.
      real(dp), allocatable :: coefficients(:)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: validation_mse = huge(1.0_dp)
      integer :: info = 0
   end type gmdh_ridge_fit_t

   type, public :: gmdh_node_t
      !! One selected quadratic or revised feedback neuron.
      real(dp), allocatable :: coefficients(:)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: mse = huge(1.0_dp)
      integer :: node_type = gmdh_node_quadratic
      integer :: first_input = 0
      integer :: second_input = 0
      integer :: feedback_inputs = 0
      integer :: transfer = gmdh_transfer_polynomial
      integer :: info = 0
   end type gmdh_node_t

   type, public :: gmdh_layer_t
      !! Selected neurons in one self-organized GMDH layer.
      type(gmdh_node_t), allocatable :: nodes(:)
   end type gmdh_layer_t

   type, public :: gmdh_model_t
      !! Fitted classical or revised GMDH forecasting network.
      type(gmdh_layer_t), allocatable :: layers(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: shift = 0.0_dp
      real(dp) :: scale = 1.0_dp
      real(dp) :: validation_weight = 0.7_dp
      integer :: input_count = 0
      integer :: layer_count = 0
      integer :: method = gmdh_method_classic
      integer :: observations = 0
      integer :: info = 0
   end type gmdh_model_t

   type, public :: gmdh_forecast_t
      !! GMDH point forecasts and rolling-error prediction intervals.
      type(gmdh_model_t) :: model
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: level = 95.0_dp
      integer :: info = 0
   end type gmdh_forecast_t

   type :: node_fit_t
      !! Internal neuron fit and its training-sample output.
      type(gmdh_node_t) :: node
      real(dp), allocatable :: fitted(:)
   end type node_fit_t

   interface display
      module procedure display_gmdh_model
      module procedure display_gmdh_forecast
   end interface display

   public :: gmdh_transform, gmdh_inverse_transform
   public :: gmdh_ridge_fit, gmdh_fit, gmdh_predict, gmdh_forecast
   public :: gmdh_coefficients, gmdh_fitted_values, gmdh_residuals
   public :: display, display_gmdh_model, display_gmdh_forecast

contains

   pure elemental real(dp) function gmdh_transform(value, transfer) result(out)
      !! Transform a normalized response for one GMDH neuron family.
      real(dp), intent(in) :: value !! Normalized response value.
      integer, intent(in) :: transfer !! Transfer-function selector.

      select case (transfer)
      case (gmdh_transfer_polynomial)
         out = value
      case (gmdh_transfer_sigmoid)
         out = log(value/(1.0_dp - value))
      case (gmdh_transfer_rbf)
         out = sqrt(-log(value))
      case (gmdh_transfer_tangent)
         out = atan(value)*180.0_dp/acos(-1.0_dp)
      case default
         out = value
      end select
   end function gmdh_transform

   pure elemental real(dp) function gmdh_inverse_transform(value, transfer) &
      result(out)
      !! Map a neuron linear predictor back to normalized response space.
      real(dp), intent(in) :: value !! Transformed response value.
      integer, intent(in) :: transfer !! Transfer-function selector.

      select case (transfer)
      case (gmdh_transfer_polynomial)
         out = value
      case (gmdh_transfer_sigmoid)
         if (value >= 0.0_dp) then
            out = 1.0_dp/(1.0_dp + exp(-value))
         else
            out = exp(value)/(1.0_dp + exp(value))
         end if
      case (gmdh_transfer_rbf)
         out = exp(-value**2)
      case (gmdh_transfer_tangent)
         out = tan(value*acos(-1.0_dp)/180.0_dp)
      case default
         out = value
      end select
   end function gmdh_inverse_transform

   pure function gmdh_ridge_fit(predictors, response, lambdas, &
      validation_weight) result(out)
      !! Select a ridge penalty by holdout error and refit using all observations.
      real(dp), intent(in) :: predictors(:, :) !! Design matrix including its intercept column.
      real(dp), intent(in) :: response(:) !! Transformed training response.
      real(dp), intent(in) :: lambdas(:) !! Candidate nonnegative ridge penalties.
      real(dp), intent(in) :: validation_weight !! Fraction assigned to penalty training.
      type(gmdh_ridge_fit_t) :: out
      real(dp), allocatable :: gram(:, :), inverse(:, :), right(:)
      real(dp), allocatable :: candidate(:), prediction(:)
      real(dp) :: cost
      integer :: observations, parameters, training, candidate_index, status
      integer :: parameter, best_index

      observations = size(predictors, 1)
      parameters = size(predictors, 2)
      training = nint(real(observations, dp)*validation_weight)
      if (observations /= size(response) .or. observations < 2 .or. &
         parameters < 1 .or. size(lambdas) < 1 .or. training < 1 .or. &
         training >= observations .or. validation_weight <= 0.0_dp .or. &
         validation_weight >= 1.0_dp .or. any(lambdas < 0.0_dp) .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         .not. all(ieee_is_finite(response)) .or. &
         .not. all(ieee_is_finite(lambdas))) then
         out%info = 1
         return
      end if
      allocate(gram(parameters, parameters), inverse(parameters, parameters))
      allocate(right(parameters), candidate(parameters))
      allocate(prediction(observations - training))
      best_index = 0
      do candidate_index = 1, size(lambdas)
         gram = matmul(transpose(predictors(:training, :)), &
            predictors(:training, :))
         do parameter = 2, parameters
            gram(parameter, parameter) = gram(parameter, parameter) + &
               lambdas(candidate_index)
         end do
         call symmetric_pseudoinverse(gram, inverse, status)
         if (status /= 0) cycle
         right = matmul(transpose(predictors(:training, :)), &
            response(:training))
         candidate = matmul(inverse, right)
         prediction = matmul(predictors(training + 1:, :), candidate)
         cost = sum((prediction - response(training + 1:))**2)/ &
            real(observations - training, dp)
         if (cost < out%validation_mse) then
            out%validation_mse = cost
            out%lambda = lambdas(candidate_index)
            best_index = candidate_index
         end if
      end do
      if (best_index == 0) then
         out%info = 2
         return
      end if
      gram = matmul(transpose(predictors), predictors)
      do parameter = 2, parameters
         gram(parameter, parameter) = gram(parameter, parameter) + out%lambda
      end do
      call symmetric_pseudoinverse(gram, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      right = matmul(transpose(predictors), response)
      out%coefficients = matmul(inverse, right)
   end function gmdh_ridge_fit

   pure function gmdh_fit(series, input_count, layer_count, method, &
      transfers, validation_weight, lambdas) result(model)
      !! Fit the classical or revised GMDH algorithm used by GMDH::fcast.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in), optional :: input_count !! Number of autoregressive inputs.
      integer, intent(in), optional :: layer_count !! Number of self-organized layers.
      integer, intent(in), optional :: method !! Classical or revised GMDH selector.
      integer, intent(in), optional :: transfers(:) !! Permitted transfer-function selectors.
      real(dp), intent(in), optional :: validation_weight !! Ridge-training fraction.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge penalties.
      type(gmdh_model_t) :: model
      type(node_fit_t) :: candidate_fit
      type(gmdh_node_t), allocatable :: candidates(:)
      real(dp), allocatable :: normalized(:), response(:), current(:, :)
      real(dp), allocatable :: outputs(:, :), scores(:), selected_lambdas(:)
      integer, allocatable :: selected_transfers(:)
      logical, allocatable :: selected(:)
      integer :: inputs, layers, selected_method, observations, rows
      integer :: layer, first, second, candidate, candidate_count, keep
      integer :: selected_index, retained, position, lag
      real(dp) :: weight

      inputs = 4
      if (present(input_count)) inputs = input_count
      layers = 3
      if (present(layer_count)) layers = layer_count
      selected_method = gmdh_method_classic
      if (present(method)) selected_method = method
      weight = 0.7_dp
      if (present(validation_weight)) weight = validation_weight
      if (present(transfers)) then
         selected_transfers = transfers
      else
         selected_transfers = [gmdh_transfer_polynomial, &
            gmdh_transfer_sigmoid, gmdh_transfer_rbf, gmdh_transfer_tangent]
      end if
      if (present(lambdas)) then
         selected_lambdas = lambdas
      else
         selected_lambdas = [0.0_dp, 0.01_dp, 0.02_dp, 0.04_dp, 0.08_dp, &
            0.16_dp, 0.32_dp, 0.64_dp, 1.28_dp, 2.56_dp, 5.12_dp, 10.24_dp]
      end if
      observations = size(series)
      rows = observations - inputs
      if (inputs < 2 .or. layers < 1 .or. rows < 3 .or. &
         (selected_method /= gmdh_method_classic .and. &
         selected_method /= gmdh_method_revised) .or. weight <= 0.0_dp .or. &
         weight >= 1.0_dp .or. size(selected_transfers) < 1 .or. &
         any(selected_transfers < gmdh_transfer_polynomial) .or. &
         any(selected_transfers > gmdh_transfer_tangent) .or. &
         size(selected_lambdas) < 1 .or. any(selected_lambdas < 0.0_dp) .or. &
         .not. all(ieee_is_finite(series))) then
         model%info = 1
         return
      end if
      model%input_count = inputs
      model%layer_count = layers
      model%method = selected_method
      model%validation_weight = weight
      model%observations = observations
      if (minval(series) <= 0.0_dp) model%shift = abs(minval(series)) + 1.0_dp
      model%scale = maxval(series + model%shift) + 1.0_dp
      normalized = (series + model%shift)/model%scale
      if (any(normalized <= 0.0_dp) .or. any(normalized >= 1.0_dp)) then
         model%info = 1
         return
      end if
      allocate(response(rows), current(rows, inputs))
      response = normalized(inputs + 1:)
      do lag = 1, inputs
         current(:, lag) = normalized(inputs + 1 - lag:observations - lag)
      end do
      allocate(model%layers(layers))
      do layer = 1, layers
         candidate_count = inputs*(inputs - 1)/2
         if (selected_method == gmdh_method_revised) &
            candidate_count = candidate_count + inputs
         allocate(candidates(candidate_count), outputs(rows, candidate_count))
         allocate(scores(candidate_count), source=huge(1.0_dp))
         candidate = 0
         do first = 1, inputs - 1
            do second = first + 1, inputs
               candidate = candidate + 1
               candidate_fit = fit_node(current, response, &
                  gmdh_node_quadratic, first, second, 0, selected_transfers, &
                  selected_lambdas, weight)
               candidates(candidate) = candidate_fit%node
               if (allocated(candidate_fit%fitted)) then
                  outputs(:, candidate) = candidate_fit%fitted
                  scores(candidate) = candidate_fit%node%mse
               end if
            end do
         end do
         if (selected_method == gmdh_method_revised) then
            do first = 1, inputs
               candidate = candidate + 1
               candidate_fit = fit_node(current, response, &
                  gmdh_node_feedback, 0, 0, first, selected_transfers, &
                  selected_lambdas, weight)
               candidates(candidate) = candidate_fit%node
               if (allocated(candidate_fit%fitted)) then
                  outputs(:, candidate) = candidate_fit%fitted
                  scores(candidate) = candidate_fit%node%mse
               end if
            end do
         end if
         keep = merge(1, inputs, layer == layers)
         if (count(candidates%info == 0) < keep) then
            model%info = 2
            return
         end if
         allocate(selected(candidate_count), source=.false.)
         do retained = 1, keep
            selected_index = minloc(scores, dim=1, mask=.not. selected)
            selected(selected_index) = .true.
         end do
         allocate(model%layers(layer)%nodes(keep))
         if (layer < layers) then
            deallocate(current)
            allocate(current(rows, keep))
         end if
         position = 0
         do candidate = 1, candidate_count
            if (.not. selected(candidate)) cycle
            position = position + 1
            model%layers(layer)%nodes(position) = candidates(candidate)
            if (layer < layers) current(:, position) = outputs(:, candidate)
         end do
         if (layer == layers) then
            model%fitted = outputs(:, selected_index)*model%scale - model%shift
            model%residuals = series(inputs + 1:) - model%fitted
         end if
         deallocate(candidates, outputs, scores, selected)
      end do
   end function gmdh_fit

   pure function gmdh_predict(model, history, horizon) result(forecast)
      !! Recursively forecast from a fitted GMDH network and supplied history.
      type(gmdh_model_t), intent(in) :: model !! Fitted GMDH network.
      real(dp), intent(in) :: history(:) !! Observed history ending at the forecast origin.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), allocatable :: forecast(:)
      real(dp), allocatable :: work(:), layer_input(:), layer_output(:)
      integer :: step, lag, layer, inputs

      inputs = model%input_count
      if (model%info /= 0 .or. horizon < 1 .or. size(history) < inputs .or. &
         .not. allocated(model%layers) .or. &
         .not. all(ieee_is_finite(history))) then
         allocate(forecast(0))
         return
      end if
      allocate(forecast(horizon), work(inputs + horizon))
      work(:inputs) = (history(size(history) - inputs + 1:) + model%shift)/ &
         model%scale
      do step = 1, horizon
         allocate(layer_input(inputs))
         do lag = 1, inputs
            layer_input(lag) = work(inputs + step - lag)
         end do
         do layer = 1, model%layer_count
            layer_output = evaluate_layer(model%layers(layer), layer_input)
            call move_alloc(layer_output, layer_input)
         end do
         work(inputs + step) = layer_input(1)
         forecast(step) = layer_input(1)*model%scale - model%shift
         deallocate(layer_input)
      end do
   end function gmdh_predict

   pure function gmdh_forecast(series, horizon, input_count, layer_count, &
      method, level, transfers, validation_weight, lambdas) result(out)
      !! Fit GMDH and form rolling-origin normal prediction intervals.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in), optional :: horizon !! Number of forecasts, limited to five.
      integer, intent(in), optional :: input_count !! Number of autoregressive inputs.
      integer, intent(in), optional :: layer_count !! Number of self-organized layers.
      integer, intent(in), optional :: method !! Classical or revised GMDH selector.
      real(dp), intent(in), optional :: level !! Prediction-interval confidence percentage.
      integer, intent(in), optional :: transfers(:) !! Permitted transfer-function selectors.
      real(dp), intent(in), optional :: validation_weight !! Ridge-training fraction.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge penalties.
      type(gmdh_forecast_t) :: out
      real(dp), allocatable :: errors(:, :), prediction(:)
      real(dp) :: probability, critical, error_mean
      integer, parameter :: interval_horizon = 5
      integer :: steps, inputs, origins, row, origin, step

      steps = 5
      if (present(horizon)) steps = horizon
      inputs = 4
      if (present(input_count)) inputs = input_count
      if (present(level)) out%level = level
      if (steps < 1 .or. steps > 5 .or. out%level < 0.0_dp .or. &
         out%level > 100.0_dp) then
         out%info = 1
         return
      end if
      out%model = gmdh_fit(series, input_count, layer_count, method, &
         transfers, validation_weight, lambdas)
      if (out%model%info /= 0) then
         out%info = out%model%info
         return
      end if
      out%mean = gmdh_predict(out%model, series, steps)
      out%fitted = out%model%fitted
      out%residuals = out%model%residuals
      origins = size(series) - inputs - interval_horizon + 1
      if (origins < 2) then
         out%info = 3
         return
      end if
      allocate(errors(origins, interval_horizon))
      do row = 1, origins
         origin = inputs + row - 1
         prediction = gmdh_predict(out%model, series(:origin), interval_horizon)
         errors(row, :) = prediction - &
            series(origin + 1:origin + interval_horizon)
      end do
      allocate(out%standard_error(steps))
      do step = 1, steps
         error_mean = sum(errors(:, step))/real(origins, dp)
         out%standard_error(step) = sqrt(sum((errors(:, step) - &
            error_mean)**2)/real(origins - 1, dp))
      end do
      probability = 1.0_dp - (1.0_dp - out%level/100.0_dp)/2.0_dp
      critical = normal_quantile(probability)
      out%lower = out%mean - critical*out%standard_error
      out%upper = out%mean + critical*out%standard_error
   end function gmdh_forecast

   pure function gmdh_coefficients(model) result(values)
      !! Return all selected neuron coefficients in layer and node order.
      type(gmdh_model_t), intent(in) :: model !! Fitted GMDH network.
      real(dp), allocatable :: values(:)
      integer :: layer, node, count, position, width

      count = 0
      if (allocated(model%layers)) then
         do layer = 1, size(model%layers)
            if (.not. allocated(model%layers(layer)%nodes)) cycle
            do node = 1, size(model%layers(layer)%nodes)
               if (allocated(model%layers(layer)%nodes(node)%coefficients)) &
                  count = count + &
                     size(model%layers(layer)%nodes(node)%coefficients)
            end do
         end do
      end if
      allocate(values(count))
      if (.not. allocated(model%layers)) return
      position = 0
      do layer = 1, size(model%layers)
         if (.not. allocated(model%layers(layer)%nodes)) cycle
         do node = 1, size(model%layers(layer)%nodes)
            if (.not. allocated( &
               model%layers(layer)%nodes(node)%coefficients)) cycle
            width = size(model%layers(layer)%nodes(node)%coefficients)
            values(position + 1:position + width) = &
               model%layers(layer)%nodes(node)%coefficients
            position = position + width
         end do
      end do
   end function gmdh_coefficients

   pure function gmdh_fitted_values(model) result(values)
      !! Return in-sample fitted values from a GMDH model.
      type(gmdh_model_t), intent(in) :: model !! Fitted GMDH network.
      real(dp), allocatable :: values(:)

      if (allocated(model%fitted)) then
         values = model%fitted
      else
         allocate(values(0))
      end if
   end function gmdh_fitted_values

   pure function gmdh_residuals(model) result(values)
      !! Return in-sample residuals from a GMDH model.
      type(gmdh_model_t), intent(in) :: model !! Fitted GMDH network.
      real(dp), allocatable :: values(:)

      if (allocated(model%residuals)) then
         values = model%residuals
      else
         allocate(values(0))
      end if
   end function gmdh_residuals

   subroutine display_gmdh_model(model, unit, print_obs)
      !! Display a fitted GMDH network and optionally its fitted observations.
      type(gmdh_model_t), intent(in) :: model !! Fitted GMDH network to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted values and residuals.
      integer :: destination, layer, node, index
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Group method of data handling model'
      write(destination, '(a, i0)') 'Status: ', model%info
      write(destination, '(a, i0)') 'Method: ', model%method
      write(destination, '(a, i0)') 'Inputs: ', model%input_count
      write(destination, '(a, i0)') 'Layers: ', model%layer_count
      write(destination, '(a, i0)') 'Observations: ', model%observations
      if (allocated(model%layers)) then
         do layer = 1, size(model%layers)
            if (.not. allocated(model%layers(layer)%nodes)) then
               write(destination, '(a, i0, a)') 'Layer ', layer, &
                  ' is not initialized.'
               cycle
            end if
            write(destination, '(a, i0, a, i0)') 'Layer ', layer, &
               ' selected nodes: ', size(model%layers(layer)%nodes)
            do node = 1, size(model%layers(layer)%nodes)
               write(destination, '(a, i0, a, i0, a, es14.6)') '  Node ', node, &
                  ', transfer ', model%layers(layer)%nodes(node)%transfer, &
                  ', MSE ', model%layers(layer)%nodes(node)%mse
            end do
         end do
      end if
      if (show_observations .and. allocated(model%fitted)) then
         write(destination, '(a)') 'Index, fitted value, residual:'
         do index = 1, size(model%fitted)
            write(destination, '(i8, 2(1x, es14.6))') &
               model%input_count + index, model%fitted(index), &
               model%residuals(index)
         end do
      end if
   end subroutine display_gmdh_model

   subroutine display_gmdh_forecast(forecast, unit, print_obs)
      !! Display GMDH forecasts and optionally the fitted observations.
      type(gmdh_forecast_t), intent(in) :: forecast !! Forecast result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted values and residuals.
      integer :: destination, step

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'GMDH short-term forecast'
      write(destination, '(a, i0)') 'Status: ', forecast%info
      write(destination, '(a, f8.3)') 'Prediction interval level: ', forecast%level
      if (allocated(forecast%mean) .and. allocated(forecast%lower) .and. &
         allocated(forecast%upper) .and. &
         allocated(forecast%standard_error)) then
         write(destination, '(a)') 'Horizon, point, lower, upper, standard error:'
         do step = 1, size(forecast%mean)
            write(destination, '(i8, 4(1x, es14.6))') step, &
               forecast%mean(step), forecast%lower(step), forecast%upper(step), &
               forecast%standard_error(step)
         end do
      end if
      if (present(print_obs)) &
         call display_gmdh_model(forecast%model, destination, print_obs)
   end subroutine display_gmdh_forecast

   pure function fit_node(predictors, response, node_type, first_input, &
      second_input, feedback_inputs, transfers, lambdas, validation_weight) &
      result(out)
      !! Fit every permitted transfer function for one candidate neuron.
      real(dp), intent(in) :: predictors(:, :) !! Current layer inputs.
      real(dp), intent(in) :: response(:) !! Normalized training response.
      integer, intent(in) :: node_type !! Quadratic or feedback node selector.
      integer, intent(in) :: first_input !! First quadratic input index.
      integer, intent(in) :: second_input !! Second quadratic input index.
      integer, intent(in) :: feedback_inputs !! Number of feedback-node inputs.
      integer, intent(in) :: transfers(:) !! Permitted transfer functions.
      real(dp), intent(in) :: lambdas(:) !! Candidate ridge penalties.
      real(dp), intent(in) :: validation_weight !! Ridge-training fraction.
      type(node_fit_t) :: out
      type(gmdh_ridge_fit_t) :: ridge
      real(dp), allocatable :: design(:, :), transformed(:), linear(:)
      real(dp), allocatable :: fitted(:)
      real(dp) :: mse
      integer :: transfer_index

      if (node_type == gmdh_node_quadratic) then
         allocate(design(size(response), 6))
         design(:, 1) = 1.0_dp
         design(:, 2) = predictors(:, first_input)
         design(:, 3) = predictors(:, second_input)
         design(:, 4) = predictors(:, first_input)*predictors(:, second_input)
         design(:, 5) = predictors(:, first_input)**2
         design(:, 6) = predictors(:, second_input)**2
      else
         allocate(design(size(response), feedback_inputs + 1))
         design(:, 1) = 1.0_dp
         design(:, 2:) = predictors(:, :feedback_inputs)
      end if
      do transfer_index = 1, size(transfers)
         transformed = gmdh_transform(response, transfers(transfer_index))
         ridge = gmdh_ridge_fit(design, transformed, lambdas, &
            validation_weight)
         if (ridge%info /= 0) cycle
         linear = matmul(design, ridge%coefficients)
         fitted = gmdh_inverse_transform(linear, transfers(transfer_index))
         mse = sum((fitted - response)**2)/real(size(response), dp)
         if (mse < out%node%mse) then
            out%node%coefficients = ridge%coefficients
            out%node%lambda = ridge%lambda
            out%node%mse = mse
            out%node%node_type = node_type
            out%node%first_input = first_input
            out%node%second_input = second_input
            out%node%feedback_inputs = feedback_inputs
            out%node%transfer = transfers(transfer_index)
            out%fitted = fitted
         end if
      end do
      if (.not. allocated(out%node%coefficients)) out%node%info = 1
   end function fit_node

   pure function evaluate_layer(layer, inputs) result(outputs)
      !! Evaluate every selected node in one GMDH layer.
      type(gmdh_layer_t), intent(in) :: layer !! Selected GMDH layer.
      real(dp), intent(in) :: inputs(:) !! Layer input vector.
      real(dp), allocatable :: outputs(:)
      real(dp) :: design(6), linear
      integer :: node, count

      allocate(outputs(size(layer%nodes)))
      do node = 1, size(layer%nodes)
         design = 0.0_dp
         design(1) = 1.0_dp
         if (layer%nodes(node)%node_type == gmdh_node_quadratic) then
            design(2) = inputs(layer%nodes(node)%first_input)
            design(3) = inputs(layer%nodes(node)%second_input)
            design(4) = design(2)*design(3)
            design(5) = design(2)**2
            design(6) = design(3)**2
            count = 6
         else
            count = layer%nodes(node)%feedback_inputs + 1
            design(2:count) = inputs(:count - 1)
         end if
         linear = dot_product(design(:count), &
            layer%nodes(node)%coefficients)
         outputs(node) = gmdh_inverse_transform(linear, &
            layer%nodes(node)%transfer)
      end do
   end function evaluate_layer

end module gmdh_mod
