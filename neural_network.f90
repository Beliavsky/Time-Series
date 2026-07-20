! SPDX-License-Identifier: MIT
! SPDX-FileComment: Reusable feed-forward neural-network regression algorithms.
module neural_network_mod
   !! Feed-forward neural regression with analytic-gradient BFGS training.
   use kind_mod, only: dp
   use optimization_mod, only: optimization_result_t, bfgs_minimize
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   type, public :: neural_network_layer_t
      !! Bias vector and incoming weights for one dense network layer.
      real(dp), allocatable :: bias(:)
      real(dp), allocatable :: weights(:, :)
   end type neural_network_layer_t

   type :: matrix_t
      !! Internal allocatable matrix container used during backpropagation.
      real(dp), allocatable :: values(:, :)
   end type matrix_t

   type, public :: neural_network_t
      !! Trained sigmoid-hidden, linear-output feed-forward regression network.
      real(dp), allocatable :: input_mean(:)
      real(dp), allocatable :: input_scale(:)
      real(dp), allocatable :: hidden_bias(:)
      real(dp), allocatable :: hidden_weights(:, :)
      real(dp), allocatable :: output_bias(:)
      real(dp), allocatable :: output_weights(:, :)
      real(dp), allocatable :: direct_weights(:, :)
      type(neural_network_layer_t), allocatable :: layers(:)
      integer, allocatable :: hidden_counts(:)
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: decay = 0.0_dp
      integer :: input_count = 0
      integer :: hidden_count = 0
      integer :: layer_count = 0
      integer :: output_count = 0
      integer :: iterations = 0
      logical :: converged = .false.
      logical :: direct = .false.
      integer :: info = 0
   end type neural_network_t

   public :: neural_network_fit, neural_network_predict
   public :: neural_network_from_parameters
   public :: neural_network_parameter_count

   interface neural_network_fit
      module procedure neural_network_fit_single
      module procedure neural_network_fit_layers
   end interface neural_network_fit

   interface neural_network_parameter_count
      module procedure neural_network_parameter_count_single
      module procedure neural_network_parameter_count_layers
   end interface neural_network_parameter_count

contains

   pure function neural_network_from_parameters(inputs, hidden, outputs, &
      parameters, direct, input_mean, input_scale) result(network)
      !! Construct a one-hidden-layer network from packed nnet-order weights.
      integer, intent(in) :: inputs !! Number of input variables.
      integer, intent(in) :: hidden !! Number of hidden units.
      integer, intent(in) :: outputs !! Number of output variables.
      real(dp), intent(in) :: parameters(:) !! Packed network parameters.
      logical, intent(in), optional :: direct !! Include trailing direct weights.
      real(dp), intent(in), optional :: input_mean(:) !! Predictor centers.
      real(dp), intent(in), optional :: input_scale(:) !! Predictor scales.
      type(neural_network_t) :: network
      integer :: base_parameters, expected
      logical :: use_direct

      use_direct = .false.
      if (present(direct)) use_direct = direct
      base_parameters = neural_network_parameter_count(inputs, hidden, outputs)
      expected = base_parameters + merge(inputs*outputs, 0, use_direct)
      if (inputs < 1 .or. hidden < 1 .or. outputs < 1 .or. &
         size(parameters) /= expected .or. &
         .not. all(ieee_is_finite(parameters))) then
         network%info = 1
         return
      end if
      if (present(input_mean)) then
         if (size(input_mean) /= inputs .or. &
            .not. all(ieee_is_finite(input_mean))) then
            network%info = 2
            return
         end if
      end if
      if (present(input_scale)) then
         if (size(input_scale) /= inputs .or. any(input_scale <= 0.0_dp) .or. &
            .not. all(ieee_is_finite(input_scale))) then
            network%info = 2
            return
         end if
      end if
      network%input_count = inputs
      network%hidden_count = hidden
      network%hidden_counts = [hidden]
      network%layer_count = 1
      network%output_count = outputs
      network%direct = use_direct
      allocate(network%input_mean(inputs), source=0.0_dp)
      allocate(network%input_scale(inputs), source=1.0_dp)
      if (present(input_mean)) network%input_mean = input_mean
      if (present(input_scale)) network%input_scale = input_scale
      call unpack_parameters(parameters(:base_parameters), inputs, hidden, &
         outputs, network%hidden_bias, network%hidden_weights, &
         network%output_bias, network%output_weights)
      allocate(network%layers(2))
      network%layers(1)%bias = network%hidden_bias
      network%layers(1)%weights = network%hidden_weights
      network%layers(2)%bias = network%output_bias
      network%layers(2)%weights = network%output_weights
      if (use_direct) then
         network%direct_weights = reshape(parameters(base_parameters + 1:), &
            [outputs, inputs])
      else
         allocate(network%direct_weights(outputs, 0))
      end if
      network%converged = .true.
   end function neural_network_from_parameters

   pure function neural_network_fit_single(predictors, response, hidden_count, &
      max_iterations, tolerance, decay, initial_parameters, direct, &
      scale_predictors) result(network)
      !! Fit a one-hidden-layer neural regression using analytic-gradient BFGS.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-input predictor matrix.
      real(dp), intent(in) :: response(:, :) !! Observation-by-output response matrix.
      integer, intent(in) :: hidden_count !! Number of sigmoid hidden units.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight-decay coefficient.
      real(dp), intent(in), optional :: initial_parameters(:) !! Initial packed network weights.
      logical, intent(in), optional :: direct !! Include jointly optimized input-output weights.
      logical, intent(in), optional :: scale_predictors !! Standardize predictor columns internally.
      type(neural_network_t) :: network
      type(optimization_result_t) :: optimum
      real(dp), allocatable :: scaled(:, :), initial(:), fitted(:, :)
      real(dp) :: selected_tolerance, selected_decay, variance
      integer :: observations, inputs, outputs, parameters, base_parameters
      integer :: limit, input, index
      logical :: use_direct, use_predictor_scaling

      observations = size(predictors, 1)
      inputs = size(predictors, 2)
      outputs = size(response, 2)
      limit = 500
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_decay = 0.0_dp
      if (present(decay)) selected_decay = decay
      use_direct = .false.
      if (present(direct)) use_direct = direct
      use_predictor_scaling = .true.
      if (present(scale_predictors)) &
         use_predictor_scaling = scale_predictors
      if (observations /= size(response, 1) .or. observations < 2 .or. &
         inputs < 1 .or. outputs < 1 .or. hidden_count < 1 .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. selected_decay < 0.0_dp .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         .not. all(ieee_is_finite(response))) then
         network%info = 1
         return
      end if
      network%input_count = inputs
      network%hidden_count = hidden_count
      network%layer_count = 1
      network%hidden_counts = [hidden_count]
      network%output_count = outputs
      network%decay = selected_decay
      network%direct = use_direct
      allocate(network%input_mean(inputs), network%input_scale(inputs))
      allocate(scaled(observations, inputs))
      if (use_predictor_scaling) then
         do input = 1, inputs
            network%input_mean(input) = sum(predictors(:, input))/ &
               real(observations, dp)
            variance = sum((predictors(:, input) - &
               network%input_mean(input))**2)/ &
               real(max(1, observations - 1), dp)
            network%input_scale(input) = sqrt(max(variance, epsilon(1.0_dp)))
            scaled(:, input) = (predictors(:, input) - &
               network%input_mean(input))/network%input_scale(input)
         end do
      else
         network%input_mean = 0.0_dp
         network%input_scale = 1.0_dp
         scaled = predictors
      end if
      base_parameters = neural_network_parameter_count(inputs, hidden_count, &
         outputs)
      parameters = base_parameters
      if (use_direct) parameters = parameters + inputs*outputs
      allocate(initial(parameters))
      if (present(initial_parameters)) then
         if (size(initial_parameters) /= parameters .or. &
            .not. all(ieee_is_finite(initial_parameters))) then
            network%info = 2
            return
         end if
         initial = initial_parameters
      else
         do index = 1, parameters
            initial(index) = 0.15_dp*sin(1.6180339887498948_dp*real(index, dp))/ &
               sqrt(real(max(1, inputs), dp))
         end do
      end if
      optimum = bfgs_minimize(objective_function, gradient_function, initial, &
         limit, selected_tolerance)
      if (optimum%info == 1 .or. optimum%info == 2 .or. &
         .not. ieee_is_finite(optimum%objective)) then
         network%info = 3
         return
      end if
      call unpack_parameters(optimum%parameters(:base_parameters), inputs, &
         hidden_count, outputs, &
         network%hidden_bias, network%hidden_weights, network%output_bias, &
         network%output_weights)
      if (use_direct) then
         network%direct_weights = reshape( &
            optimum%parameters(base_parameters + 1:), [outputs, inputs])
      else
         allocate(network%direct_weights(outputs, 0))
      end if
      allocate(network%layers(2))
      network%layers(1)%bias = network%hidden_bias
      network%layers(1)%weights = network%hidden_weights
      network%layers(2)%bias = network%output_bias
      network%layers(2)%weights = network%output_weights
      network%objective = optimum%objective
      network%iterations = optimum%iterations
      network%converged = optimum%converged
      fitted = neural_network_predict(network, predictors)
      network%rss = sum((fitted - response)**2)

   contains

      pure function objective_function(parameters_in) result(value)
         !! Evaluate penalized squared error for the host training data.
         real(dp), intent(in) :: parameters_in(:) !! Packed network weights.
         real(dp) :: value
         real(dp), allocatable :: prediction(:, :)

         call evaluate_parameters(parameters_in(:base_parameters), scaled, &
            hidden_count, outputs, &
            prediction)
         if (use_direct) prediction = prediction + matmul(scaled, transpose( &
            reshape(parameters_in(base_parameters + 1:), [outputs, inputs])))
         value = 0.5_dp*sum((prediction - response)**2) + &
            0.5_dp*selected_decay*sum(parameters_in**2)
      end function objective_function

      pure function gradient_function(parameters_in) result(gradient)
         !! Evaluate the analytic gradient for the host training data.
         real(dp), intent(in) :: parameters_in(:) !! Packed network weights.
         real(dp) :: gradient(size(parameters_in))
         real(dp), allocatable :: hidden_bias(:), hidden_weights(:, :)
         real(dp), allocatable :: output_bias(:), output_weights(:, :)
         real(dp), allocatable :: hidden(:, :), prediction(:, :), error(:, :)
         real(dp), allocatable :: hidden_delta(:, :), gradient_hidden_bias(:)
         real(dp), allocatable :: gradient_hidden_weights(:, :)
         real(dp), allocatable :: gradient_output_bias(:)
         real(dp), allocatable :: gradient_output_weights(:, :)
         real(dp), allocatable :: base_gradient(:)

         call unpack_parameters(parameters_in(:base_parameters), inputs, &
            hidden_count, outputs, &
            hidden_bias, hidden_weights, output_bias, output_weights)
         hidden = sigmoid_matrix(matmul(scaled, transpose(hidden_weights)) + &
            spread(hidden_bias, 1, observations))
         prediction = matmul(hidden, transpose(output_weights)) + &
            spread(output_bias, 1, observations)
         if (use_direct) prediction = prediction + matmul(scaled, transpose( &
            reshape(parameters_in(base_parameters + 1:), [outputs, inputs])))
         error = prediction - response
         gradient_output_bias = sum(error, dim=1)
         gradient_output_weights = matmul(transpose(error), hidden)
         hidden_delta = matmul(error, output_weights)*hidden*(1.0_dp - hidden)
         gradient_hidden_bias = sum(hidden_delta, dim=1)
         gradient_hidden_weights = matmul(transpose(hidden_delta), scaled)
         allocate(base_gradient(base_parameters))
         call pack_parameters(gradient_hidden_bias, gradient_hidden_weights, &
            gradient_output_bias, gradient_output_weights, base_gradient)
         gradient(:base_parameters) = base_gradient
         if (use_direct) gradient(base_parameters + 1:) = reshape( &
            matmul(transpose(error), scaled), [inputs*outputs])
         gradient = gradient + selected_decay*parameters_in
      end function gradient_function

   end function neural_network_fit_single

   pure function neural_network_fit_layers(predictors, response, hidden_counts, &
      max_iterations, tolerance, decay, initial_parameters) result(network)
      !! Fit a multilayer neural regression using analytic-gradient BFGS.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-input predictor matrix.
      real(dp), intent(in) :: response(:, :) !! Observation-by-output response matrix.
      integer, intent(in) :: hidden_counts(:) !! Numbers of units in successive hidden layers.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight-decay coefficient.
      real(dp), intent(in), optional :: initial_parameters(:) !! Initial packed network weights.
      type(neural_network_t) :: network
      type(optimization_result_t) :: optimum
      real(dp), allocatable :: scaled(:, :), initial(:), fitted(:, :)
      real(dp) :: selected_tolerance, selected_decay, variance
      integer :: observations, inputs, outputs, parameters, limit
      integer :: input, index, layer_count

      observations = size(predictors, 1)
      inputs = size(predictors, 2)
      outputs = size(response, 2)
      layer_count = size(hidden_counts)
      limit = 500
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_decay = 0.0_dp
      if (present(decay)) selected_decay = decay
      if (observations /= size(response, 1) .or. observations < 2 .or. &
         inputs < 1 .or. outputs < 1 .or. layer_count < 1 .or. &
         any(hidden_counts < 1) .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. selected_decay < 0.0_dp .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         .not. all(ieee_is_finite(response))) then
         network%info = 1
         return
      end if
      network%input_count = inputs
      network%hidden_count = hidden_counts(1)
      network%hidden_counts = hidden_counts
      network%layer_count = layer_count
      network%output_count = outputs
      network%decay = selected_decay
      allocate(network%input_mean(inputs), network%input_scale(inputs))
      allocate(scaled(observations, inputs))
      do input = 1, inputs
         network%input_mean(input) = sum(predictors(:, input))/real(observations, dp)
         variance = sum((predictors(:, input) - network%input_mean(input))**2)/ &
            real(max(1, observations - 1), dp)
         network%input_scale(input) = sqrt(max(variance, epsilon(1.0_dp)))
         scaled(:, input) = (predictors(:, input) - network%input_mean(input))/ &
            network%input_scale(input)
      end do
      parameters = neural_network_parameter_count_layers(inputs, hidden_counts, &
         outputs)
      allocate(initial(parameters))
      if (present(initial_parameters)) then
         if (size(initial_parameters) /= parameters .or. &
            .not. all(ieee_is_finite(initial_parameters))) then
            network%info = 2
            return
         end if
         initial = initial_parameters
      else
         do index = 1, parameters
            initial(index) = 0.15_dp*sin(1.6180339887498948_dp*real(index, dp))/ &
               sqrt(real(max(1, inputs), dp))
         end do
      end if
      optimum = bfgs_minimize(objective_function, gradient_function, initial, &
         limit, selected_tolerance)
      if (optimum%info == 1 .or. optimum%info == 2 .or. &
         .not. ieee_is_finite(optimum%objective)) then
         network%info = 3
         return
      end if
      call unpack_multilayer_parameters(optimum%parameters, inputs, &
         hidden_counts, outputs, network%layers)
      network%hidden_bias = network%layers(1)%bias
      network%hidden_weights = network%layers(1)%weights
      network%output_bias = network%layers(layer_count + 1)%bias
      network%output_weights = network%layers(layer_count + 1)%weights
      network%objective = optimum%objective
      network%iterations = optimum%iterations
      network%converged = optimum%converged
      fitted = neural_network_predict(network, predictors)
      network%rss = sum((fitted - response)**2)

   contains

      pure function objective_function(parameters_in) result(value)
         !! Evaluate penalized squared error for the host training data.
         real(dp), intent(in) :: parameters_in(:) !! Packed weights for all layers.
         real(dp) :: value
         real(dp), allocatable :: prediction(:, :)

         call evaluate_multilayer_parameters(parameters_in, scaled, &
            hidden_counts, outputs, prediction)
         value = 0.5_dp*sum((prediction - response)**2) + &
            0.5_dp*selected_decay*sum(parameters_in**2)
      end function objective_function

      pure function gradient_function(parameters_in) result(gradient)
         !! Evaluate the analytic multilayer backpropagation gradient.
         real(dp), intent(in) :: parameters_in(:) !! Packed weights for all layers.
         real(dp) :: gradient(size(parameters_in))
         type(neural_network_layer_t), allocatable :: layers(:), derivatives(:)
         type(matrix_t), allocatable :: activations(:)
         real(dp), allocatable :: prediction(:, :), delta(:, :)
         integer :: layer

         call unpack_multilayer_parameters(parameters_in, inputs, hidden_counts, &
            outputs, layers)
         allocate(activations(layer_count + 1), derivatives(layer_count + 1))
         activations(1)%values = scaled
         do layer = 1, layer_count
            activations(layer + 1)%values = sigmoid_matrix(matmul( &
               activations(layer)%values, transpose(layers(layer)%weights)) + &
               spread(layers(layer)%bias, 1, observations))
         end do
         prediction = matmul(activations(layer_count + 1)%values, &
            transpose(layers(layer_count + 1)%weights)) + &
            spread(layers(layer_count + 1)%bias, 1, observations)
         delta = prediction - response
         derivatives(layer_count + 1)%bias = sum(delta, dim=1)
         derivatives(layer_count + 1)%weights = matmul(transpose(delta), &
            activations(layer_count + 1)%values)
         do layer = layer_count, 1, -1
            delta = matmul(delta, layers(layer + 1)%weights)* &
               activations(layer + 1)%values* &
               (1.0_dp - activations(layer + 1)%values)
            derivatives(layer)%bias = sum(delta, dim=1)
            derivatives(layer)%weights = matmul(transpose(delta), &
               activations(layer)%values)
         end do
         call pack_multilayer_parameters(derivatives, gradient)
         gradient = gradient + selected_decay*parameters_in
      end function gradient_function

   end function neural_network_fit_layers

   pure function neural_network_predict(network, predictors) result(prediction)
      !! Evaluate a trained network for new predictor rows.
      type(neural_network_t), intent(in) :: network !! Trained neural network.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-input predictor matrix.
      real(dp), allocatable :: prediction(:, :)
      real(dp), allocatable :: scaled(:, :), hidden(:, :)
      integer :: layer

      if (network%info /= 0 .or. size(predictors, 2) /= network%input_count) then
         allocate(prediction(0, 0))
         return
      end if
      scaled = (predictors - spread(network%input_mean, 1, size(predictors, 1)))/ &
         spread(network%input_scale, 1, size(predictors, 1))
      if (allocated(network%layers)) then
         hidden = scaled
         do layer = 1, network%layer_count
            hidden = sigmoid_matrix(matmul(hidden, &
               transpose(network%layers(layer)%weights)) + &
               spread(network%layers(layer)%bias, 1, size(predictors, 1)))
         end do
         prediction = matmul(hidden, &
            transpose(network%layers(network%layer_count + 1)%weights)) + &
            spread(network%layers(network%layer_count + 1)%bias, 1, &
            size(predictors, 1))
      else
         hidden = sigmoid_matrix(matmul(scaled, transpose(network%hidden_weights)) + &
            spread(network%hidden_bias, 1, size(predictors, 1)))
         prediction = matmul(hidden, transpose(network%output_weights)) + &
            spread(network%output_bias, 1, size(predictors, 1))
      end if
      if (network%direct .and. allocated(network%direct_weights)) &
         prediction = prediction + matmul(scaled, &
            transpose(network%direct_weights))
   end function neural_network_predict

   pure integer function neural_network_parameter_count_single(inputs, hidden, &
      outputs, direct) result(count)
      !! Return the number of biases and weights in a dense one-hidden-layer network.
      integer, intent(in) :: inputs !! Number of input variables.
      integer, intent(in) :: hidden !! Number of hidden units.
      integer, intent(in) :: outputs !! Number of output variables.
      logical, intent(in), optional :: direct !! Include input-output weights.

      count = hidden*(inputs + 1) + outputs*(hidden + 1)
      if (present(direct)) then
         if (direct) count = count + inputs*outputs
      end if
   end function neural_network_parameter_count_single

   pure integer function neural_network_parameter_count_layers(inputs, &
      hidden_counts, outputs) result(count)
      !! Return the parameter count for a dense multilayer network.
      integer, intent(in) :: inputs !! Number of input variables.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      integer, intent(in) :: outputs !! Number of output variables.
      integer :: layer, previous

      if (inputs < 1 .or. outputs < 1 .or. size(hidden_counts) < 1 .or. &
         any(hidden_counts < 1)) then
         count = 0
         return
      end if
      count = 0
      previous = inputs
      do layer = 1, size(hidden_counts)
         count = count + hidden_counts(layer)*(previous + 1)
         previous = hidden_counts(layer)
      end do
      count = count + outputs*(previous + 1)
   end function neural_network_parameter_count_layers

   pure subroutine evaluate_parameters(parameters, predictors, hidden_count, &
      output_count, prediction)
      !! Evaluate packed network parameters for standardized predictors.
      real(dp), intent(in) :: parameters(:) !! Packed network weights.
      real(dp), intent(in) :: predictors(:, :) !! Standardized predictor rows.
      integer, intent(in) :: hidden_count !! Number of hidden units.
      integer, intent(in) :: output_count !! Number of outputs.
      real(dp), allocatable, intent(out) :: prediction(:, :) !! Network output rows.
      real(dp), allocatable :: hidden_bias(:), hidden_weights(:, :)
      real(dp), allocatable :: output_bias(:), output_weights(:, :), hidden(:, :)

      call unpack_parameters(parameters, size(predictors, 2), hidden_count, &
         output_count, hidden_bias, hidden_weights, output_bias, output_weights)
      hidden = sigmoid_matrix(matmul(predictors, transpose(hidden_weights)) + &
         spread(hidden_bias, 1, size(predictors, 1)))
      prediction = matmul(hidden, transpose(output_weights)) + &
         spread(output_bias, 1, size(predictors, 1))
   end subroutine evaluate_parameters

   pure subroutine unpack_parameters(parameters, inputs, hidden, outputs, &
      hidden_bias, hidden_weights, output_bias, output_weights)
      !! Unpack a parameter vector into network biases and weight matrices.
      real(dp), intent(in) :: parameters(:) !! Packed network parameters.
      integer, intent(in) :: inputs !! Number of inputs.
      integer, intent(in) :: hidden !! Number of hidden units.
      integer, intent(in) :: outputs !! Number of outputs.
      real(dp), allocatable, intent(out) :: hidden_bias(:) !! Hidden-unit biases.
      real(dp), allocatable, intent(out) :: hidden_weights(:, :) !! Hidden input weights.
      real(dp), allocatable, intent(out) :: output_bias(:) !! Output biases.
      real(dp), allocatable, intent(out) :: output_weights(:, :) !! Hidden-to-output weights.
      integer :: first, last

      allocate(hidden_bias(hidden), hidden_weights(hidden, inputs))
      allocate(output_bias(outputs), output_weights(outputs, hidden))
      hidden_bias = parameters(:hidden)
      first = hidden + 1
      last = first + hidden*inputs - 1
      hidden_weights = reshape(parameters(first:last), [hidden, inputs])
      first = last + 1
      last = first + outputs - 1
      output_bias = parameters(first:last)
      first = last + 1
      output_weights = reshape(parameters(first:), [outputs, hidden])
   end subroutine unpack_parameters

   pure subroutine pack_parameters(hidden_bias, hidden_weights, output_bias, &
      output_weights, parameters)
      !! Pack network biases and weight matrices into one parameter vector.
      real(dp), intent(in) :: hidden_bias(:) !! Hidden-unit bias values.
      real(dp), intent(in) :: hidden_weights(:, :) !! Hidden input-weight matrix.
      real(dp), intent(in) :: output_bias(:) !! Output bias values.
      real(dp), intent(in) :: output_weights(:, :) !! Hidden-to-output weights.
      real(dp), intent(out) :: parameters(:) !! Packed network parameters.
      integer :: first, last

      parameters(:size(hidden_bias)) = hidden_bias
      first = size(hidden_bias) + 1
      last = first + size(hidden_weights) - 1
      parameters(first:last) = reshape(hidden_weights, [size(hidden_weights)])
      first = last + 1
      last = first + size(output_bias) - 1
      parameters(first:last) = output_bias
      first = last + 1
      parameters(first:) = reshape(output_weights, [size(output_weights)])
   end subroutine pack_parameters

   pure subroutine evaluate_multilayer_parameters(parameters, predictors, &
      hidden_counts, output_count, prediction)
      !! Evaluate packed multilayer parameters for standardized predictors.
      real(dp), intent(in) :: parameters(:) !! Packed weights for all layers.
      real(dp), intent(in) :: predictors(:, :) !! Standardized predictor rows.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      integer, intent(in) :: output_count !! Number of linear output units.
      real(dp), allocatable, intent(out) :: prediction(:, :) !! Network output rows.
      type(neural_network_layer_t), allocatable :: layers(:)
      real(dp), allocatable :: activation(:, :)
      integer :: layer

      call unpack_multilayer_parameters(parameters, size(predictors, 2), &
         hidden_counts, output_count, layers)
      activation = predictors
      do layer = 1, size(hidden_counts)
         activation = sigmoid_matrix(matmul(activation, &
            transpose(layers(layer)%weights)) + &
            spread(layers(layer)%bias, 1, size(predictors, 1)))
      end do
      prediction = matmul(activation, &
         transpose(layers(size(hidden_counts) + 1)%weights)) + &
         spread(layers(size(hidden_counts) + 1)%bias, 1, size(predictors, 1))
   end subroutine evaluate_multilayer_parameters

   pure subroutine unpack_multilayer_parameters(parameters, inputs, &
      hidden_counts, outputs, layers)
      !! Unpack biases and weights for every hidden and output layer.
      real(dp), intent(in) :: parameters(:) !! Packed network parameters.
      integer, intent(in) :: inputs !! Number of input variables.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      integer, intent(in) :: outputs !! Number of output variables.
      type(neural_network_layer_t), allocatable, intent(out) :: layers(:) !! Unpacked layers.
      integer :: layer, previous, units, first, last

      allocate(layers(size(hidden_counts) + 1))
      previous = inputs
      first = 1
      do layer = 1, size(layers)
         if (layer <= size(hidden_counts)) then
            units = hidden_counts(layer)
         else
            units = outputs
         end if
         last = first + units - 1
         layers(layer)%bias = parameters(first:last)
         first = last + 1
         last = first + units*previous - 1
         layers(layer)%weights = reshape(parameters(first:last), [units, previous])
         first = last + 1
         previous = units
      end do
   end subroutine unpack_multilayer_parameters

   pure subroutine pack_multilayer_parameters(layers, parameters)
      !! Pack all layer derivatives or parameters into one vector.
      type(neural_network_layer_t), intent(in) :: layers(:) !! Network layer values.
      real(dp), intent(out) :: parameters(:) !! Packed network values.
      integer :: layer, first, last

      first = 1
      do layer = 1, size(layers)
         last = first + size(layers(layer)%bias) - 1
         parameters(first:last) = layers(layer)%bias
         first = last + 1
         last = first + size(layers(layer)%weights) - 1
         parameters(first:last) = reshape(layers(layer)%weights, &
            [size(layers(layer)%weights)])
         first = last + 1
      end do
   end subroutine pack_multilayer_parameters

   pure elemental real(dp) function sigmoid(value) result(transformed)
      !! Evaluate a numerically stable logistic sigmoid.
      real(dp), intent(in) :: value !! Linear hidden-unit activation.

      if (value >= 0.0_dp) then
         transformed = 1.0_dp/(1.0_dp + exp(-min(value, 700.0_dp)))
      else
         transformed = exp(max(value, -700.0_dp))/ &
            (1.0_dp + exp(max(value, -700.0_dp)))
      end if
   end function sigmoid

   pure function sigmoid_matrix(values) result(transformed)
      !! Apply the logistic sigmoid elementally to a matrix.
      real(dp), intent(in) :: values(:, :) !! Linear hidden activations.
      real(dp) :: transformed(size(values, 1), size(values, 2))

      transformed = sigmoid(values)
   end function sigmoid_matrix

end module neural_network_mod
