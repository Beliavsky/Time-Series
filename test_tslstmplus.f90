! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the R TSLSTMplus package translation.
program test_tslstmplus
   use kind_mod, only: dp
   use tslstmplus_mod
   use recurrent_network_mod, only: recurrent_lstm_model_t, &
      recurrent_lstm_predict, recurrent_activation_linear
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   implicit none

   type(tslstmplus_prepared_t) :: sequences, flattened
   type(tslstmplus_model_t) :: model, flat_model, stateful_model
   type(recurrent_lstm_model_t) :: fixed_network
   real(dp) :: series(100), regressors(100, 2), future_regressors(6, 2)
   real(dp), allocatable :: data(:, :), prediction(:), forecast(:)
   real(dp) :: fixed_sequence(1, 3, 1)
   real(dp), allocatable :: fixed_prediction(:, :)
   integer :: time, scratch_unit, display_size

   series(1) = 0.2_dp
   regressors(1, :) = [0.1_dp, -0.2_dp]
   do time = 2, size(series)
      regressors(time, 1) = sin(0.13_dp*real(time, dp))
      regressors(time, 2) = cos(0.07_dp*real(time, dp))
      series(time) = 0.62_dp*series(time - 1) + &
         0.28_dp*regressors(time, 1) - 0.17_dp*regressors(time, 2) + &
         0.03_dp*sin(0.91_dp*real(time, dp))
   end do
   data = tslstmplus_prepare_data(series, regressors, 3, [2])
   sequences = tslstmplus_prepare_sequences(series, regressors, 3, [2], .true.)
   flattened = tslstmplus_prepare_sequences(series, regressors, 3, [2], .false.)
   call check(all(shape(data) == [97, 10]) .and. sequences%info == 0 .and. &
      all(shape(sequences%sequences) == [97, 3, 3]) .and. &
      all(shape(flattened%sequences) == [97, 1, 9]), &
      'lagged data and tensor preparation')

   fixed_network%input_count = 1
   fixed_network%output_count = 1
   allocate(fixed_network%lstm_layers(1))
   fixed_network%lstm_layers(1)%inputs = 1
   fixed_network%lstm_layers(1)%units = 1
   fixed_network%lstm_layers(1)%input_weights = reshape( &
      [0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp], [4, 1])
   fixed_network%lstm_layers(1)%recurrent_weights = reshape( &
      [0.05_dp, 0.07_dp, -0.02_dp, 0.1_dp], [4, 1])
   fixed_network%lstm_layers(1)%bias = [0.1_dp, 0.2_dp, -0.1_dp, 0.05_dp]
   allocate(fixed_network%dense_layers(1))
   fixed_network%dense_layers(1)%weights = reshape([1.5_dp], [1, 1])
   fixed_network%dense_layers(1)%bias = [-0.2_dp]
   fixed_network%dense_layers(1)%activation = recurrent_activation_linear
   fixed_sequence(1, :, 1) = [0.4_dp, -0.2_dp, 0.7_dp]
   fixed_prediction = recurrent_lstm_predict(fixed_network, fixed_sequence)
   call check(abs(fixed_prediction(1, 1) + 0.1813494593932795_dp) < &
      1.0e-14_dp, 'TensorFlow-order fixed-weight LSTM recurrence')

   model = tslstmplus_fit(series, [5, 3], regressors, response_lag=3, &
      xreg_lags=[2], dense_units=[3], dropout_rate=0.05_dp, epochs=60, &
      learning_rate=0.004_dp, optimizer=tslstmplus_optimizer_adam, &
      scale_output=tslstmplus_scale_standard, &
      scale_input=tslstmplus_scale_minmax, batch_size=10, &
      validation_split=0.0_dp, patience=0, seed=29)
   call check(model%info == 0 .and. model%network%epochs_trained == 60 .and. &
      ieee_is_nan(model%fitted(1)) .and. &
      all(ieee_is_finite(model%fitted(4:))) .and. &
      model%network%training_loss(60) < model%network%training_loss(1), &
      'stacked Adam LSTM fit')
   prediction = tslstmplus_predict(model, series, regressors)
   call check(size(prediction) == 97 .and. all(ieee_is_finite(prediction)) .and. &
      maxval(abs(prediction - model%fitted(4:))) < 1.0e-12_dp, &
      'rolling LSTM prediction')
   do time = 1, size(future_regressors, 1)
      future_regressors(time, 1) = sin(0.13_dp*real(100 + time, dp))
      future_regressors(time, 2) = cos(0.07_dp*real(100 + time, dp))
   end do
   forecast = tslstmplus_forecast(model, series, 6, regressors, &
      future_regressors)
   call check(size(forecast) == 6 .and. all(ieee_is_finite(forecast)), &
      'recursive forecast with future regressors')

   flat_model = tslstmplus_fit(series, [4], response_lag=3, epochs=15, &
      learning_rate=0.002_dp, optimizer=tslstmplus_optimizer_rmsprop, &
      scale_output=tslstmplus_scale_minmax, validation_split=0.2_dp, &
      patience=3, lags_as_sequences=.false., seed=31)
   call check(flat_model%info == 0 .and. &
      flat_model%network%epochs_trained >= 1 .and. &
      all(ieee_is_finite(flat_model%fitted(4:))), &
      'flattened-lag RMSprop LSTM with early stopping')

   stateful_model = tslstmplus_fit(series, [3], response_lag=2, epochs=5, &
      learning_rate=0.001_dp, optimizer=tslstmplus_optimizer_sgd, &
      batch_size=5, validation_split=0.0_dp, patience=0, stateful=.true., &
      seed=37)
   call check(stateful_model%info == 0 .and. stateful_model%stateful .and. &
      stateful_model%lags_as_sequences .and. &
      all(ieee_is_finite(stateful_model%fitted(3:))), &
      'stateful SGD LSTM fit')

   open(newunit=scratch_unit, status='scratch', action='readwrite')
   call display(model, scratch_unit)
   inquire(unit=scratch_unit, size=display_size)
   close(scratch_unit)
   call check(display_size > 0, 'TSLSTMplus display method')

   print '(a)', 'TSLSTMplus tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//message
         error stop 1
      end if
   end subroutine check

end program test_tslstmplus
