! SPDX-License-Identifier: MIT
! SPDX-FileComment: Reusable recurrent neural-network algorithms.
module recurrent_network_mod
   !! Stand-alone stacked LSTM regression with analytic backpropagation.
   use kind_mod, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: recurrent_activation_linear = 0
   integer, parameter, public :: recurrent_activation_sigmoid = 1
   integer, parameter, public :: recurrent_activation_relu = 2
   integer, parameter, public :: recurrent_activation_tanh = 3
   integer, parameter, public :: recurrent_optimizer_sgd = 1
   integer, parameter, public :: recurrent_optimizer_adam = 2
   integer, parameter, public :: recurrent_optimizer_rmsprop = 3
   integer, parameter, public :: recurrent_loss_mse = 1
   integer, parameter, public :: recurrent_loss_mae = 2

   type, public :: recurrent_lstm_layer_t
      !! Parameters and optimizer state for one LSTM layer.
      real(dp), allocatable :: input_weights(:, :)
      real(dp), allocatable :: recurrent_weights(:, :)
      real(dp), allocatable :: bias(:)
      real(dp), allocatable :: input_first_moment(:, :)
      real(dp), allocatable :: recurrent_first_moment(:, :)
      real(dp), allocatable :: bias_first_moment(:)
      real(dp), allocatable :: input_second_moment(:, :)
      real(dp), allocatable :: recurrent_second_moment(:, :)
      real(dp), allocatable :: bias_second_moment(:)
      integer :: inputs = 0
      integer :: units = 0
   end type recurrent_lstm_layer_t

   type, public :: recurrent_dense_layer_t
      !! Parameters and optimizer state for one dense output layer.
      real(dp), allocatable :: weights(:, :)
      real(dp), allocatable :: bias(:)
      real(dp), allocatable :: weight_first_moment(:, :)
      real(dp), allocatable :: bias_first_moment(:)
      real(dp), allocatable :: weight_second_moment(:, :)
      real(dp), allocatable :: bias_second_moment(:)
      integer :: activation = recurrent_activation_linear
   end type recurrent_dense_layer_t

   type, public :: recurrent_lstm_model_t
      !! Fitted stacked LSTM and dense regression head.
      type(recurrent_lstm_layer_t), allocatable :: lstm_layers(:)
      type(recurrent_dense_layer_t), allocatable :: dense_layers(:)
      integer, allocatable :: units(:)
      integer, allocatable :: dense_units(:)
      real(dp), allocatable :: training_loss(:)
      real(dp), allocatable :: validation_loss(:)
      real(dp) :: learning_rate = 0.001_dp
      real(dp) :: dropout_rate = 0.0_dp
      real(dp) :: validation_split = 0.1_dp
      real(dp) :: minimum_delta = 0.0_dp
      integer :: optimizer = recurrent_optimizer_rmsprop
      integer :: loss = recurrent_loss_mse
      integer :: dense_activation = recurrent_activation_relu
      integer :: lstm_activation = recurrent_activation_tanh
      integer :: recurrent_activation = recurrent_activation_sigmoid
      integer :: epochs_requested = 0
      integer :: epochs_trained = 0
      integer :: batch_size = 1
      integer :: patience = 3
      integer :: input_count = 0
      integer :: output_count = 0
      logical :: stateful = .false.
      logical :: converged = .false.
      integer :: info = 0
   end type recurrent_lstm_model_t

   type :: vector_t
      !! Internal differently sized vector container.
      real(dp), allocatable :: value(:)
   end type vector_t

   type :: matrix_t
      !! Internal differently sized matrix container.
      real(dp), allocatable :: value(:, :)
   end type matrix_t

   type :: lstm_cache_t
      !! Forward-pass values retained for one LSTM layer.
      real(dp), allocatable :: input(:, :)
      real(dp), allocatable :: multiplier(:, :)
      real(dp), allocatable :: hidden(:, :)
      real(dp), allocatable :: cell(:, :)
      real(dp), allocatable :: input_gate(:, :)
      real(dp), allocatable :: forget_gate(:, :)
      real(dp), allocatable :: output_gate(:, :)
      real(dp), allocatable :: candidate(:, :)
   end type lstm_cache_t

   type :: dense_cache_t
      !! Forward-pass values retained for one dense layer.
      real(dp), allocatable :: input(:)
      real(dp), allocatable :: preactivation(:)
      real(dp), allocatable :: output(:)
   end type dense_cache_t

   type :: lstm_gradient_t
      !! Accumulated gradient for one LSTM layer.
      real(dp), allocatable :: input_weights(:, :)
      real(dp), allocatable :: recurrent_weights(:, :)
      real(dp), allocatable :: bias(:)
   end type lstm_gradient_t

   type :: dense_gradient_t
      !! Accumulated gradient for one dense layer.
      real(dp), allocatable :: weights(:, :)
      real(dp), allocatable :: bias(:)
   end type dense_gradient_t

   public :: recurrent_lstm_fit, recurrent_lstm_predict

contains

   pure function recurrent_lstm_fit(sequences, targets, units, dense_units, &
      epochs, batch_size, learning_rate, optimizer, loss, dropout_rate, &
      dense_activation, validation_split, patience, minimum_delta, stateful, &
      seed, lstm_activation, recurrent_activation) result(model)
      !! Fit stacked LSTM regression by mini-batch truncated BPTT.
      real(dp), intent(in) :: sequences(:, :, :) !! Sample-by-time-by-feature inputs.
      real(dp), intent(in) :: targets(:, :) !! Sample-by-output targets.
      integer, intent(in) :: units(:) !! Units in successive LSTM layers.
      integer, intent(in), optional :: dense_units(:) !! Units in hidden dense layers.
      integer, intent(in), optional :: epochs !! Maximum training epochs.
      integer, intent(in), optional :: batch_size !! Mini-batch size.
      real(dp), intent(in), optional :: learning_rate !! Optimizer learning rate.
      integer, intent(in), optional :: optimizer !! SGD, Adam, or RMSprop code.
      integer, intent(in), optional :: loss !! MSE or MAE loss code.
      real(dp), intent(in), optional :: dropout_rate !! LSTM input dropout fraction.
      integer, intent(in), optional :: dense_activation !! Hidden dense activation code.
      real(dp), intent(in), optional :: validation_split !! Terminal validation fraction.
      integer, intent(in), optional :: patience !! Early-stopping patience.
      real(dp), intent(in), optional :: minimum_delta !! Minimum validation improvement.
      logical, intent(in), optional :: stateful !! Carry states between sample sequences.
      integer, intent(in), optional :: seed !! Deterministic initialization seed.
      integer, intent(in), optional :: lstm_activation !! Candidate and cell activation code.
      integer, intent(in), optional :: recurrent_activation !! Gate activation code.
      type(recurrent_lstm_model_t) :: model
      integer, allocatable :: selected_dense_units(:)
      integer :: layer, inputs, outputs, selected_seed

      model%units = units
      allocate(selected_dense_units(0))
      if (present(dense_units)) selected_dense_units = dense_units
      model%dense_units = selected_dense_units
      model%epochs_requested = 10
      if (present(epochs)) model%epochs_requested = epochs
      model%batch_size = 1
      if (present(batch_size)) model%batch_size = batch_size
      model%learning_rate = 0.001_dp
      if (present(learning_rate)) model%learning_rate = learning_rate
      model%optimizer = recurrent_optimizer_rmsprop
      if (present(optimizer)) model%optimizer = optimizer
      model%loss = recurrent_loss_mse
      if (present(loss)) model%loss = loss
      model%dropout_rate = 0.0_dp
      if (present(dropout_rate)) model%dropout_rate = dropout_rate
      model%dense_activation = recurrent_activation_relu
      if (present(dense_activation)) model%dense_activation = dense_activation
      model%lstm_activation = recurrent_activation_tanh
      if (present(lstm_activation)) model%lstm_activation = lstm_activation
      model%recurrent_activation = recurrent_activation_sigmoid
      if (present(recurrent_activation)) &
         model%recurrent_activation = recurrent_activation
      model%validation_split = 0.1_dp
      if (present(validation_split)) model%validation_split = validation_split
      model%patience = 3
      if (present(patience)) model%patience = patience
      model%minimum_delta = 0.0_dp
      if (present(minimum_delta)) model%minimum_delta = minimum_delta
      model%stateful = .false.
      if (present(stateful)) model%stateful = stateful
      selected_seed = 1
      if (present(seed)) selected_seed = seed
      model%input_count = size(sequences, 3)
      model%output_count = size(targets, 2)
      if (size(sequences, 1) /= size(targets, 1) .or. &
         size(sequences, 1) < 2 .or. size(sequences, 2) < 1 .or. &
         size(sequences, 3) < 1 .or. size(targets, 2) < 1 .or. &
         size(units) < 1 .or. any(units < 1) .or. &
         any(selected_dense_units < 1) .or. model%epochs_requested < 1 .or. &
         model%batch_size < 1 .or. model%learning_rate <= 0.0_dp .or. &
         model%dropout_rate < 0.0_dp .or. model%dropout_rate >= 1.0_dp .or. &
         model%validation_split < 0.0_dp .or. &
         model%validation_split >= 1.0_dp .or. model%patience < 0 .or. &
         model%minimum_delta < 0.0_dp .or. &
         (model%optimizer < recurrent_optimizer_sgd .or. &
         model%optimizer > recurrent_optimizer_rmsprop) .or. &
         (model%loss /= recurrent_loss_mse .and. &
         model%loss /= recurrent_loss_mae) .or. &
         (model%dense_activation < recurrent_activation_linear .or. &
         model%dense_activation > recurrent_activation_tanh) .or. &
         (model%lstm_activation < recurrent_activation_linear .or. &
         model%lstm_activation > recurrent_activation_tanh) .or. &
         (model%recurrent_activation < recurrent_activation_linear .or. &
         model%recurrent_activation > recurrent_activation_tanh) .or. &
         .not. all(ieee_is_finite(sequences)) .or. &
         .not. all(ieee_is_finite(targets))) then
         model%info = 1
         return
      end if
      allocate(model%lstm_layers(size(units)))
      inputs = model%input_count
      do layer = 1, size(units)
         call initialize_lstm_layer(model%lstm_layers(layer), inputs, &
            units(layer), selected_seed + 1009*layer)
         inputs = units(layer)
      end do
      allocate(model%dense_layers(size(selected_dense_units) + 1))
      do layer = 1, size(model%dense_layers)
         if (layer <= size(selected_dense_units)) then
            outputs = selected_dense_units(layer)
            call initialize_dense_layer(model%dense_layers(layer), inputs, &
               outputs, model%dense_activation, selected_seed + 2017*layer)
         else
            outputs = model%output_count
            call initialize_dense_layer(model%dense_layers(layer), inputs, &
               outputs, recurrent_activation_linear, selected_seed + 2017*layer)
         end if
         inputs = outputs
      end do
      allocate(model%training_loss(model%epochs_requested))
      allocate(model%validation_loss(model%epochs_requested))
      model%training_loss = huge(1.0_dp)
      model%validation_loss = huge(1.0_dp)
      call train_model(model, sequences, targets, selected_seed)
   end function recurrent_lstm_fit

   pure function recurrent_lstm_predict(model, sequences) result(prediction)
      !! Predict all supplied sequences with fitted recurrent parameters.
      type(recurrent_lstm_model_t), intent(in) :: model !! Fitted recurrent model.
      real(dp), intent(in) :: sequences(:, :, :) !! Sample-by-time-by-feature inputs.
      real(dp), allocatable :: prediction(:, :)
      type(vector_t), allocatable :: hidden_state(:), cell_state(:)
      type(lstm_cache_t), allocatable :: lstm_cache(:)
      type(dense_cache_t), allocatable :: dense_cache(:)
      type(vector_t), allocatable :: final_hidden(:), final_cell(:)
      real(dp), allocatable :: sample_prediction(:)
      integer :: sample, slots, slot, first_state, last_state, layers

      if (model%info /= 0 .or. size(sequences, 3) /= model%input_count .or. &
         size(sequences, 2) < 1 .or. .not. all(ieee_is_finite(sequences))) then
         allocate(prediction(0, 0))
         return
      end if
      allocate(prediction(size(sequences, 1), model%output_count))
      layers = size(model%lstm_layers)
      slots = merge(model%batch_size, 1, model%stateful)
      call zero_states(model, hidden_state, cell_state, slots)
      do sample = 1, size(sequences, 1)
         if (.not. model%stateful) call clear_states(hidden_state, cell_state)
         slot = 1 + mod(sample - 1, slots)
         first_state = (slot - 1)*layers + 1
         last_state = slot*layers
         call forward_sample(model, sequences(sample, :, :), &
            hidden_state(first_state:last_state), &
            cell_state(first_state:last_state), .false., 0, sample, &
            lstm_cache, dense_cache, &
            sample_prediction, final_hidden, final_cell)
         prediction(sample, :) = sample_prediction
         if (model%stateful) then
            hidden_state(first_state:last_state) = final_hidden
            cell_state(first_state:last_state) = final_cell
         end if
      end do
   end function recurrent_lstm_predict

   pure subroutine initialize_lstm_layer(layer, inputs, units, seed)
      !! Initialize one LSTM layer and its optimizer state.
      type(recurrent_lstm_layer_t), intent(out) :: layer !! Initialized layer.
      integer, intent(in) :: inputs !! Input feature count.
      integer, intent(in) :: units !! Hidden-unit count.
      integer, intent(in) :: seed !! Deterministic seed.
      real(dp) :: input_scale, recurrent_scale
      integer :: row, column

      layer%inputs = inputs
      layer%units = units
      allocate(layer%input_weights(4*units, inputs))
      allocate(layer%recurrent_weights(4*units, units))
      allocate(layer%bias(4*units))
      input_scale = sqrt(2.0_dp/real(max(1, inputs + units), dp))
      recurrent_scale = 1.0_dp/sqrt(real(max(1, units), dp))
      do row = 1, 4*units
         do column = 1, inputs
            layer%input_weights(row, column) = input_scale* &
               deterministic_weight(seed, row, column)
         end do
         do column = 1, units
            layer%recurrent_weights(row, column) = recurrent_scale* &
               deterministic_weight(seed + 37, row, column)
         end do
      end do
      layer%bias = 0.0_dp
      layer%bias(units + 1:2*units) = 1.0_dp
      allocate(layer%input_first_moment(4*units, inputs), source=0.0_dp)
      allocate(layer%recurrent_first_moment(4*units, units), source=0.0_dp)
      allocate(layer%bias_first_moment(4*units), source=0.0_dp)
      allocate(layer%input_second_moment(4*units, inputs), source=0.0_dp)
      allocate(layer%recurrent_second_moment(4*units, units), source=0.0_dp)
      allocate(layer%bias_second_moment(4*units), source=0.0_dp)
   end subroutine initialize_lstm_layer

   pure subroutine initialize_dense_layer(layer, inputs, outputs, activation, seed)
      !! Initialize one dense layer and its optimizer state.
      type(recurrent_dense_layer_t), intent(out) :: layer !! Initialized layer.
      integer, intent(in) :: inputs !! Input width.
      integer, intent(in) :: outputs !! Output width.
      integer, intent(in) :: activation !! Activation code.
      integer, intent(in) :: seed !! Deterministic seed.
      real(dp) :: scale
      integer :: row, column

      allocate(layer%weights(outputs, inputs), layer%bias(outputs))
      scale = sqrt(2.0_dp/real(max(1, inputs + outputs), dp))
      do row = 1, outputs
         do column = 1, inputs
            layer%weights(row, column) = scale* &
               deterministic_weight(seed, row, column)
         end do
      end do
      layer%bias = 0.0_dp
      layer%activation = activation
      allocate(layer%weight_first_moment(outputs, inputs), source=0.0_dp)
      allocate(layer%bias_first_moment(outputs), source=0.0_dp)
      allocate(layer%weight_second_moment(outputs, inputs), source=0.0_dp)
      allocate(layer%bias_second_moment(outputs), source=0.0_dp)
   end subroutine initialize_dense_layer

   pure subroutine train_model(model, sequences, targets, seed)
      !! Run mini-batch BPTT with terminal validation and early stopping.
      type(recurrent_lstm_model_t), intent(inout) :: model !! Model updated in place.
      real(dp), intent(in) :: sequences(:, :, :) !! Training sequences.
      real(dp), intent(in) :: targets(:, :) !! Training targets.
      integer, intent(in) :: seed !! Deterministic dropout seed.
      type(lstm_gradient_t), allocatable :: lstm_gradient(:)
      type(dense_gradient_t), allocatable :: dense_gradient(:)
      type(recurrent_lstm_layer_t), allocatable :: best_lstm(:)
      type(recurrent_dense_layer_t), allocatable :: best_dense(:)
      type(vector_t), allocatable :: hidden_state(:), cell_state(:)
      type(vector_t), allocatable :: final_hidden(:), final_cell(:)
      type(lstm_cache_t), allocatable :: lstm_cache(:)
      type(dense_cache_t), allocatable :: dense_cache(:)
      real(dp), allocatable :: prediction(:)
      real(dp) :: best_loss, validation_loss
      integer :: training_samples, epoch, observation, row, batch_count
      integer :: update, stale_epochs, slots, slot, first_state, last_state
      integer :: layers

      training_samples = size(sequences, 1) - &
         int(floor(model%validation_split*real(size(sequences, 1), dp)))
      training_samples = max(1, min(size(sequences, 1), training_samples))
      call allocate_gradients(model, lstm_gradient, dense_gradient)
      layers = size(model%lstm_layers)
      slots = merge(model%batch_size, 1, model%stateful)
      call zero_states(model, hidden_state, cell_state, slots)
      best_lstm = model%lstm_layers
      best_dense = model%dense_layers
      best_loss = huge(1.0_dp)
      update = 0
      stale_epochs = 0
      do epoch = 1, model%epochs_requested
         call zero_gradients(lstm_gradient, dense_gradient)
         call clear_states(hidden_state, cell_state)
         batch_count = 0
         model%training_loss(epoch) = 0.0_dp
         do observation = 1, training_samples
            if (model%stateful) then
               row = observation
            else
               row = 1 + mod(observation + epoch - 2, training_samples)
               call clear_states(hidden_state, cell_state)
            end if
            slot = 1 + mod(observation - 1, slots)
            first_state = (slot - 1)*layers + 1
            last_state = slot*layers
            call forward_sample(model, sequences(row, :, :), &
               hidden_state(first_state:last_state), &
               cell_state(first_state:last_state), .true., seed + epoch, row, &
               lstm_cache, dense_cache, &
               prediction, final_hidden, final_cell)
            model%training_loss(epoch) = model%training_loss(epoch) + &
               sample_loss(prediction, targets(row, :), model%loss)
            call backward_sample(model, targets(row, :), prediction, &
               lstm_cache, dense_cache, lstm_gradient, dense_gradient)
            if (model%stateful) then
               hidden_state(first_state:last_state) = final_hidden
               cell_state(first_state:last_state) = final_cell
            end if
            batch_count = batch_count + 1
            if (batch_count < model%batch_size .and. &
               observation < training_samples) cycle
            update = update + 1
            call update_parameters(model, lstm_gradient, dense_gradient, &
               batch_count, update)
            call zero_gradients(lstm_gradient, dense_gradient)
            batch_count = 0
         end do
         model%training_loss(epoch) = model%training_loss(epoch)/ &
            real(training_samples, dp)
         if (training_samples < size(sequences, 1)) then
            validation_loss = evaluate_loss(model, &
               sequences(training_samples + 1:, :, :), &
               targets(training_samples + 1:, :))
         else
            validation_loss = model%training_loss(epoch)
         end if
         model%validation_loss(epoch) = validation_loss
         model%epochs_trained = epoch
         if (validation_loss < best_loss - model%minimum_delta) then
            best_loss = validation_loss
            best_lstm = model%lstm_layers
            best_dense = model%dense_layers
            stale_epochs = 0
         else
            stale_epochs = stale_epochs + 1
         end if
         if (model%patience > 0 .and. stale_epochs >= model%patience) exit
      end do
      model%lstm_layers = best_lstm
      model%dense_layers = best_dense
      model%converged = model%epochs_trained < model%epochs_requested .or. &
         model%epochs_trained == model%epochs_requested
   end subroutine train_model

   pure subroutine forward_sample(model, sequence, initial_hidden, initial_cell, &
      training, dropout_seed, sample, lstm_cache, dense_cache, prediction, &
      final_hidden, final_cell)
      !! Evaluate one sequence and retain values required by BPTT.
      type(recurrent_lstm_model_t), intent(in) :: model !! Recurrent model.
      real(dp), intent(in) :: sequence(:, :) !! Time-by-feature sequence.
      type(vector_t), intent(in) :: initial_hidden(:) !! Initial hidden states.
      type(vector_t), intent(in) :: initial_cell(:) !! Initial cell states.
      logical, intent(in) :: training !! Apply training dropout.
      integer, intent(in) :: dropout_seed !! Deterministic dropout seed.
      integer, intent(in) :: sample !! Sample index for dropout.
      type(lstm_cache_t), allocatable, intent(out) :: lstm_cache(:) !! LSTM caches.
      type(dense_cache_t), allocatable, intent(out) :: dense_cache(:) !! Dense caches.
      real(dp), allocatable, intent(out) :: prediction(:) !! Output prediction.
      type(vector_t), allocatable, intent(out) :: final_hidden(:) !! Final hidden states.
      type(vector_t), allocatable, intent(out) :: final_cell(:) !! Final cell states.
      real(dp), allocatable :: layer_input(:, :), gate(:), dense_value(:)
      integer :: layer, time, units

      allocate(lstm_cache(size(model%lstm_layers)))
      allocate(final_hidden(size(model%lstm_layers)))
      allocate(final_cell(size(model%lstm_layers)))
      layer_input = sequence
      do layer = 1, size(model%lstm_layers)
         units = model%lstm_layers(layer)%units
         allocate(lstm_cache(layer)%input(size(sequence, 1), size(layer_input, 2)))
         allocate(lstm_cache(layer)%multiplier(size(sequence, 1), &
            size(layer_input, 2)))
         lstm_cache(layer)%multiplier = 1.0_dp
         if (training .and. model%dropout_rate > 0.0_dp) &
            call dropout_multiplier(lstm_cache(layer)%multiplier, &
            model%dropout_rate, dropout_seed, sample, layer)
         lstm_cache(layer)%input = layer_input*lstm_cache(layer)%multiplier
         allocate(lstm_cache(layer)%hidden(0:size(sequence, 1), units))
         allocate(lstm_cache(layer)%cell(0:size(sequence, 1), units))
         allocate(lstm_cache(layer)%input_gate(size(sequence, 1), units))
         allocate(lstm_cache(layer)%forget_gate(size(sequence, 1), units))
         allocate(lstm_cache(layer)%output_gate(size(sequence, 1), units))
         allocate(lstm_cache(layer)%candidate(size(sequence, 1), units))
         lstm_cache(layer)%hidden(0, :) = initial_hidden(layer)%value
         lstm_cache(layer)%cell(0, :) = initial_cell(layer)%value
         do time = 1, size(sequence, 1)
            gate = matmul(model%lstm_layers(layer)%input_weights, &
               lstm_cache(layer)%input(time, :)) + &
               matmul(model%lstm_layers(layer)%recurrent_weights, &
               lstm_cache(layer)%hidden(time - 1, :)) + &
               model%lstm_layers(layer)%bias
            lstm_cache(layer)%input_gate(time, :) = &
               activation_value(gate(:units), model%recurrent_activation)
            lstm_cache(layer)%forget_gate(time, :) = &
               activation_value(gate(units + 1:2*units), &
               model%recurrent_activation)
            lstm_cache(layer)%candidate(time, :) = &
               activation_value(gate(2*units + 1:3*units), &
               model%lstm_activation)
            lstm_cache(layer)%output_gate(time, :) = &
               activation_value(gate(3*units + 1:4*units), &
               model%recurrent_activation)
            lstm_cache(layer)%cell(time, :) = &
               lstm_cache(layer)%forget_gate(time, :)* &
               lstm_cache(layer)%cell(time - 1, :) + &
               lstm_cache(layer)%input_gate(time, :)* &
               lstm_cache(layer)%candidate(time, :)
            lstm_cache(layer)%hidden(time, :) = &
               lstm_cache(layer)%output_gate(time, :)* &
               activation_value(lstm_cache(layer)%cell(time, :), &
               model%lstm_activation)
         end do
         final_hidden(layer)%value = &
            lstm_cache(layer)%hidden(size(sequence, 1), :)
         final_cell(layer)%value = &
            lstm_cache(layer)%cell(size(sequence, 1), :)
         layer_input = lstm_cache(layer)%hidden(1:, :)
      end do
      allocate(dense_cache(size(model%dense_layers)))
      dense_value = final_hidden(size(final_hidden))%value
      do layer = 1, size(model%dense_layers)
         dense_cache(layer)%input = dense_value
         dense_cache(layer)%preactivation = &
            matmul(model%dense_layers(layer)%weights, dense_value) + &
            model%dense_layers(layer)%bias
         dense_cache(layer)%output = activation_value( &
            dense_cache(layer)%preactivation, &
            model%dense_layers(layer)%activation)
         dense_value = dense_cache(layer)%output
      end do
      prediction = dense_value
   end subroutine forward_sample

   pure subroutine backward_sample(model, target, prediction, lstm_cache, &
      dense_cache, lstm_gradient, dense_gradient)
      !! Accumulate one sample's dense and recurrent gradients.
      type(recurrent_lstm_model_t), intent(in) :: model !! Recurrent model.
      real(dp), intent(in) :: target(:) !! Desired output.
      real(dp), intent(in) :: prediction(:) !! Predicted output.
      type(lstm_cache_t), intent(in) :: lstm_cache(:) !! Forward LSTM caches.
      type(dense_cache_t), intent(in) :: dense_cache(:) !! Forward dense caches.
      type(lstm_gradient_t), intent(inout) :: lstm_gradient(:) !! LSTM gradients.
      type(dense_gradient_t), intent(inout) :: dense_gradient(:) !! Dense gradients.
      type(matrix_t), allocatable :: upstream(:)
      real(dp), allocatable :: delta(:), pre_delta(:), input_delta(:)
      real(dp), allocatable :: hidden_delta(:), cell_delta(:), cell_total(:)
      real(dp), allocatable :: gate_delta(:), candidate_delta(:)
      real(dp), allocatable :: forget_delta(:), output_delta(:)
      real(dp), allocatable :: input_gate_delta(:), cell_tanh(:)
      integer :: layer, time, units

      if (model%loss == recurrent_loss_mae) then
         delta = merge(sign(1.0_dp, prediction - target), 0.0_dp, &
            prediction /= target)/real(size(target), dp)
      else
         delta = (prediction - target)/real(size(target), dp)
      end if
      do layer = size(model%dense_layers), 1, -1
         pre_delta = delta*activation_derivative( &
            dense_cache(layer)%preactivation, &
            model%dense_layers(layer)%activation)
         dense_gradient(layer)%weights = dense_gradient(layer)%weights + &
            outer_product(pre_delta, dense_cache(layer)%input)
         dense_gradient(layer)%bias = dense_gradient(layer)%bias + pre_delta
         delta = matmul(transpose(model%dense_layers(layer)%weights), pre_delta)
      end do
      allocate(upstream(size(model%lstm_layers)))
      do layer = 1, size(model%lstm_layers)
         allocate(upstream(layer)%value(size(lstm_cache(layer)%input, 1), &
            model%lstm_layers(layer)%units), source=0.0_dp)
      end do
      upstream(size(upstream))%value(size(upstream(size(upstream))%value, 1), :) = &
         delta
      do layer = size(model%lstm_layers), 1, -1
         units = model%lstm_layers(layer)%units
         allocate(hidden_delta(units), source=0.0_dp)
         allocate(cell_delta(units), source=0.0_dp)
         do time = size(lstm_cache(layer)%input, 1), 1, -1
            hidden_delta = hidden_delta + upstream(layer)%value(time, :)
            cell_tanh = activation_value(lstm_cache(layer)%cell(time, :), &
               model%lstm_activation)
            output_delta = hidden_delta*cell_tanh* &
               activation_derivative_from_output( &
               lstm_cache(layer)%output_gate(time, :), &
               model%recurrent_activation)
            cell_total = cell_delta + hidden_delta* &
               lstm_cache(layer)%output_gate(time, :)*activation_derivative( &
               lstm_cache(layer)%cell(time, :), model%lstm_activation)
            forget_delta = cell_total*lstm_cache(layer)%cell(time - 1, :)* &
               activation_derivative_from_output( &
               lstm_cache(layer)%forget_gate(time, :), &
               model%recurrent_activation)
            input_gate_delta = cell_total*lstm_cache(layer)%candidate(time, :)* &
               activation_derivative_from_output( &
               lstm_cache(layer)%input_gate(time, :), &
               model%recurrent_activation)
            candidate_delta = cell_total*lstm_cache(layer)%input_gate(time, :)* &
               activation_derivative_from_output( &
               lstm_cache(layer)%candidate(time, :), model%lstm_activation)
            gate_delta = [input_gate_delta, forget_delta, candidate_delta, &
               output_delta]
            lstm_gradient(layer)%input_weights = &
               lstm_gradient(layer)%input_weights + outer_product(gate_delta, &
               lstm_cache(layer)%input(time, :))
            lstm_gradient(layer)%recurrent_weights = &
               lstm_gradient(layer)%recurrent_weights + &
               outer_product(gate_delta, &
               lstm_cache(layer)%hidden(time - 1, :))
            lstm_gradient(layer)%bias = lstm_gradient(layer)%bias + gate_delta
            input_delta = matmul( &
               transpose(model%lstm_layers(layer)%input_weights), gate_delta)* &
               lstm_cache(layer)%multiplier(time, :)
            if (layer > 1) upstream(layer - 1)%value(time, :) = &
               upstream(layer - 1)%value(time, :) + input_delta
            hidden_delta = matmul(transpose( &
               model%lstm_layers(layer)%recurrent_weights), gate_delta)
            cell_delta = cell_total*lstm_cache(layer)%forget_gate(time, :)
         end do
         deallocate(hidden_delta, cell_delta)
      end do
   end subroutine backward_sample

   pure subroutine allocate_gradients(model, lstm_gradient, dense_gradient)
      !! Allocate zero gradient arrays matching a recurrent model.
      type(recurrent_lstm_model_t), intent(in) :: model !! Recurrent model.
      type(lstm_gradient_t), allocatable, intent(out) :: lstm_gradient(:) !! LSTM gradients.
      type(dense_gradient_t), allocatable, intent(out) :: dense_gradient(:) !! Dense gradients.
      integer :: layer

      allocate(lstm_gradient(size(model%lstm_layers)))
      do layer = 1, size(lstm_gradient)
         allocate(lstm_gradient(layer)%input_weights, &
            mold=model%lstm_layers(layer)%input_weights)
         allocate(lstm_gradient(layer)%recurrent_weights, &
            mold=model%lstm_layers(layer)%recurrent_weights)
         allocate(lstm_gradient(layer)%bias, mold=model%lstm_layers(layer)%bias)
      end do
      allocate(dense_gradient(size(model%dense_layers)))
      do layer = 1, size(dense_gradient)
         allocate(dense_gradient(layer)%weights, &
            mold=model%dense_layers(layer)%weights)
         allocate(dense_gradient(layer)%bias, mold=model%dense_layers(layer)%bias)
      end do
      call zero_gradients(lstm_gradient, dense_gradient)
   end subroutine allocate_gradients

   pure subroutine zero_gradients(lstm_gradient, dense_gradient)
      !! Reset all accumulated gradient arrays to zero.
      type(lstm_gradient_t), intent(inout) :: lstm_gradient(:) !! LSTM gradients.
      type(dense_gradient_t), intent(inout) :: dense_gradient(:) !! Dense gradients.
      integer :: layer

      do layer = 1, size(lstm_gradient)
         lstm_gradient(layer)%input_weights = 0.0_dp
         lstm_gradient(layer)%recurrent_weights = 0.0_dp
         lstm_gradient(layer)%bias = 0.0_dp
      end do
      do layer = 1, size(dense_gradient)
         dense_gradient(layer)%weights = 0.0_dp
         dense_gradient(layer)%bias = 0.0_dp
      end do
   end subroutine zero_gradients

   pure subroutine update_parameters(model, lstm_gradient, dense_gradient, &
      batch_count, update)
      !! Apply one SGD, Adam, or RMSprop parameter update.
      type(recurrent_lstm_model_t), intent(inout) :: model !! Model parameters.
      type(lstm_gradient_t), intent(in) :: lstm_gradient(:) !! LSTM gradients.
      type(dense_gradient_t), intent(in) :: dense_gradient(:) !! Dense gradients.
      integer, intent(in) :: batch_count !! Samples accumulated in the batch.
      integer, intent(in) :: update !! One-based optimizer update number.
      integer :: layer

      do layer = 1, size(model%lstm_layers)
         call update_array(model%lstm_layers(layer)%input_weights, &
            model%lstm_layers(layer)%input_first_moment, &
            model%lstm_layers(layer)%input_second_moment, &
            lstm_gradient(layer)%input_weights/real(batch_count, dp), model, update)
         call update_array(model%lstm_layers(layer)%recurrent_weights, &
            model%lstm_layers(layer)%recurrent_first_moment, &
            model%lstm_layers(layer)%recurrent_second_moment, &
            lstm_gradient(layer)%recurrent_weights/real(batch_count, dp), &
            model, update)
         call update_vector(model%lstm_layers(layer)%bias, &
            model%lstm_layers(layer)%bias_first_moment, &
            model%lstm_layers(layer)%bias_second_moment, &
            lstm_gradient(layer)%bias/real(batch_count, dp), model, update)
      end do
      do layer = 1, size(model%dense_layers)
         call update_array(model%dense_layers(layer)%weights, &
            model%dense_layers(layer)%weight_first_moment, &
            model%dense_layers(layer)%weight_second_moment, &
            dense_gradient(layer)%weights/real(batch_count, dp), model, update)
         call update_vector(model%dense_layers(layer)%bias, &
            model%dense_layers(layer)%bias_first_moment, &
            model%dense_layers(layer)%bias_second_moment, &
            dense_gradient(layer)%bias/real(batch_count, dp), model, update)
      end do
   end subroutine update_parameters

   pure subroutine update_array(parameter, first_moment, second_moment, &
      gradient, model, update)
      !! Update one matrix parameter and its optimizer state.
      real(dp), intent(inout) :: parameter(:, :) !! Parameter matrix.
      real(dp), intent(inout) :: first_moment(:, :) !! First-moment state.
      real(dp), intent(inout) :: second_moment(:, :) !! Second-moment state.
      real(dp), intent(in) :: gradient(:, :) !! Mean batch gradient.
      type(recurrent_lstm_model_t), intent(in) :: model !! Optimizer settings.
      integer, intent(in) :: update !! One-based update number.

      select case (model%optimizer)
      case (recurrent_optimizer_adam)
         first_moment = 0.9_dp*first_moment + 0.1_dp*gradient
         second_moment = 0.999_dp*second_moment + 0.001_dp*gradient**2
         parameter = parameter - model%learning_rate* &
            (first_moment/(1.0_dp - 0.9_dp**update))/ &
            (sqrt(second_moment/(1.0_dp - 0.999_dp**update)) + 1.0e-8_dp)
      case (recurrent_optimizer_rmsprop)
         second_moment = 0.9_dp*second_moment + 0.1_dp*gradient**2
         parameter = parameter - model%learning_rate*gradient/ &
            (sqrt(second_moment) + 1.0e-7_dp)
      case default
         parameter = parameter - model%learning_rate*gradient
      end select
   end subroutine update_array

   pure subroutine update_vector(parameter, first_moment, second_moment, &
      gradient, model, update)
      !! Update one vector parameter and its optimizer state.
      real(dp), intent(inout) :: parameter(:) !! Parameter vector.
      real(dp), intent(inout) :: first_moment(:) !! First-moment state.
      real(dp), intent(inout) :: second_moment(:) !! Second-moment state.
      real(dp), intent(in) :: gradient(:) !! Mean batch gradient.
      type(recurrent_lstm_model_t), intent(in) :: model !! Optimizer settings.
      integer, intent(in) :: update !! One-based update number.

      select case (model%optimizer)
      case (recurrent_optimizer_adam)
         first_moment = 0.9_dp*first_moment + 0.1_dp*gradient
         second_moment = 0.999_dp*second_moment + 0.001_dp*gradient**2
         parameter = parameter - model%learning_rate* &
            (first_moment/(1.0_dp - 0.9_dp**update))/ &
            (sqrt(second_moment/(1.0_dp - 0.999_dp**update)) + 1.0e-8_dp)
      case (recurrent_optimizer_rmsprop)
         second_moment = 0.9_dp*second_moment + 0.1_dp*gradient**2
         parameter = parameter - model%learning_rate*gradient/ &
            (sqrt(second_moment) + 1.0e-7_dp)
      case default
         parameter = parameter - model%learning_rate*gradient
      end select
   end subroutine update_vector

   pure real(dp) function evaluate_loss(model, sequences, targets) result(loss)
      !! Evaluate mean prediction loss without dropout.
      type(recurrent_lstm_model_t), intent(in) :: model !! Recurrent model.
      real(dp), intent(in) :: sequences(:, :, :) !! Validation sequences.
      real(dp), intent(in) :: targets(:, :) !! Validation targets.
      real(dp), allocatable :: prediction(:, :)
      integer :: sample

      prediction = recurrent_lstm_predict(model, sequences)
      loss = 0.0_dp
      do sample = 1, size(targets, 1)
         loss = loss + sample_loss(prediction(sample, :), targets(sample, :), &
            model%loss)
      end do
      loss = loss/real(size(targets, 1), dp)
   end function evaluate_loss

   pure real(dp) function sample_loss(prediction, target, loss_code) result(loss)
      !! Evaluate MSE or MAE for one output vector.
      real(dp), intent(in) :: prediction(:) !! Predicted values.
      real(dp), intent(in) :: target(:) !! Desired values.
      integer, intent(in) :: loss_code !! Loss-function code.

      if (loss_code == recurrent_loss_mae) then
         loss = sum(abs(prediction - target))/real(size(target), dp)
      else
         loss = 0.5_dp*sum((prediction - target)**2)/real(size(target), dp)
      end if
   end function sample_loss

   pure subroutine zero_states(model, hidden_state, cell_state, slots)
      !! Allocate zero hidden and cell states for every recurrent layer.
      type(recurrent_lstm_model_t), intent(in) :: model !! Recurrent model.
      type(vector_t), allocatable, intent(out) :: hidden_state(:) !! Hidden states.
      type(vector_t), allocatable, intent(out) :: cell_state(:) !! Cell states.
      integer, intent(in), optional :: slots !! Independent batch-position states.
      integer :: layer, slot, selected_slots, index

      selected_slots = 1
      if (present(slots)) selected_slots = slots
      allocate(hidden_state(size(model%lstm_layers)*selected_slots))
      allocate(cell_state(size(model%lstm_layers)*selected_slots))
      do slot = 1, selected_slots
         do layer = 1, size(model%lstm_layers)
            index = (slot - 1)*size(model%lstm_layers) + layer
            allocate(hidden_state(index)%value( &
               model%lstm_layers(layer)%units), source=0.0_dp)
            allocate(cell_state(index)%value( &
               model%lstm_layers(layer)%units), source=0.0_dp)
         end do
      end do
   end subroutine zero_states

   pure subroutine clear_states(hidden_state, cell_state)
      !! Reset allocated recurrent states to zero.
      type(vector_t), intent(inout) :: hidden_state(:) !! Hidden states.
      type(vector_t), intent(inout) :: cell_state(:) !! Cell states.
      integer :: layer

      do layer = 1, size(hidden_state)
         hidden_state(layer)%value = 0.0_dp
         cell_state(layer)%value = 0.0_dp
      end do
   end subroutine clear_states

   pure subroutine dropout_multiplier(multiplier, dropout_rate, seed, sample, &
      layer)
      !! Construct deterministic inverted input-dropout multipliers.
      real(dp), intent(inout) :: multiplier(:, :) !! Time-by-feature multipliers.
      real(dp), intent(in) :: dropout_rate !! Dropout fraction.
      integer, intent(in) :: seed !! Deterministic seed.
      integer, intent(in) :: sample !! Sample index.
      integer, intent(in) :: layer !! Layer index.
      real(dp) :: uniform
      integer :: feature

      do feature = 1, size(multiplier, 2)
         uniform = 0.5_dp + 0.5_dp*sin(real(seed + 1009*sample + &
            9176*layer + 7919*feature, dp)*0.6180339887498949_dp)
         if (uniform < dropout_rate) then
            multiplier(:, feature) = 0.0_dp
         else
            multiplier(:, feature) = 1.0_dp/(1.0_dp - dropout_rate)
         end if
      end do
   end subroutine dropout_multiplier

   pure function outer_product(first, second) result(product)
      !! Return the outer product of two vectors.
      real(dp), intent(in) :: first(:) !! Row dimension values.
      real(dp), intent(in) :: second(:) !! Column dimension values.
      real(dp) :: product(size(first), size(second))

      product = spread(first, 2, size(second))*spread(second, 1, size(first))
   end function outer_product

   pure elemental real(dp) function sigmoid(value) result(transformed)
      !! Evaluate a numerically stable logistic sigmoid.
      real(dp), intent(in) :: value !! Input value.

      if (value >= 0.0_dp) then
         transformed = 1.0_dp/(1.0_dp + exp(-value))
      else
         transformed = exp(value)/(1.0_dp + exp(value))
      end if
   end function sigmoid

   pure function activation_value(value, activation) result(transformed)
      !! Apply a dense-layer activation to a vector.
      real(dp), intent(in) :: value(:) !! Preactivation values.
      integer, intent(in) :: activation !! Activation code.
      real(dp) :: transformed(size(value))

      select case (activation)
      case (recurrent_activation_sigmoid)
         transformed = sigmoid(value)
      case (recurrent_activation_relu)
         transformed = max(0.0_dp, value)
      case (recurrent_activation_tanh)
         transformed = tanh(value)
      case default
         transformed = value
      end select
   end function activation_value

   pure function activation_derivative(value, activation) result(derivative)
      !! Evaluate a dense-layer activation derivative.
      real(dp), intent(in) :: value(:) !! Preactivation values.
      integer, intent(in) :: activation !! Activation code.
      real(dp) :: derivative(size(value)), transformed(size(value))

      select case (activation)
      case (recurrent_activation_sigmoid)
         transformed = sigmoid(value)
         derivative = transformed*(1.0_dp - transformed)
      case (recurrent_activation_relu)
         derivative = merge(1.0_dp, 0.0_dp, value > 0.0_dp)
      case (recurrent_activation_tanh)
         derivative = 1.0_dp - tanh(value)**2
      case default
         derivative = 1.0_dp
      end select
   end function activation_derivative

   pure function activation_derivative_from_output(value, activation) &
      result(derivative)
      !! Evaluate an activation derivative from already activated values.
      real(dp), intent(in) :: value(:) !! Activated values.
      integer, intent(in) :: activation !! Activation code.
      real(dp) :: derivative(size(value))

      select case (activation)
      case (recurrent_activation_sigmoid)
         derivative = value*(1.0_dp - value)
      case (recurrent_activation_relu)
         derivative = merge(1.0_dp, 0.0_dp, value > 0.0_dp)
      case (recurrent_activation_tanh)
         derivative = 1.0_dp - value**2
      case default
         derivative = 1.0_dp
      end select
   end function activation_derivative_from_output

   pure real(dp) function deterministic_weight(seed, row, column) result(value)
      !! Return one deterministic centered initialization value.
      integer, intent(in) :: seed !! Initialization seed.
      integer, intent(in) :: row !! Matrix row.
      integer, intent(in) :: column !! Matrix column.

      value = sin(real(seed + 104729*row + 13007*column, dp)* &
         0.0174532925199433_dp)
   end function deterministic_weight

end module recurrent_network_mod
