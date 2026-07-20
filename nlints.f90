! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R NlinTS package.
module nlints_mod
   !! Nonlinear causality, neural VAR, and information-theory algorithms.
   use kind_mod, only: dp
   use special_functions_mod, only: regularized_beta
   use stats_mod, only: ols_fit
   use urca_mod, only: adf_result_t, adf_test
   use utils_mod, only: quiet_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: nlints_activation_linear = 0
   integer, parameter, public :: nlints_activation_sigmoid = 1
   integer, parameter, public :: nlints_activation_relu = 2
   integer, parameter, public :: nlints_activation_tanh = 3
   integer, parameter, public :: nlints_optimizer_sgd = 1
   integer, parameter, public :: nlints_optimizer_adam = 2

   type :: real_vector_t
      !! Internal vector container for differently sized neural layers.
      real(dp), allocatable :: value(:)
   end type real_vector_t

   type :: real_matrix_t
      !! Internal matrix container for differently sized neural layers.
      real(dp), allocatable :: value(:, :)
   end type real_matrix_t

   type, public :: nlints_dense_layer_t
      !! One fully connected neural-network layer.
      real(dp), allocatable :: weights(:, :)
      real(dp), allocatable :: first_moment(:, :)
      real(dp), allocatable :: second_moment(:, :)
      integer :: activation = nlints_activation_linear
   end type nlints_dense_layer_t

   type, public :: nlints_varmlp_t
      !! Min-max-scaled multivariate autoregressive neural network.
      type(nlints_dense_layer_t), allocatable :: layers(:)
      real(dp), allocatable :: minimum(:)
      real(dp), allocatable :: maximum(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: rss(:)
      integer, allocatable :: hidden_counts(:)
      integer, allocatable :: activations(:)
      integer :: lag = 0
      integer :: variables = 0
      integer :: iterations = 0
      integer :: optimizer = nlints_optimizer_sgd
      real(dp) :: learning_rate = 0.01_dp
      logical :: bias = .true.
      integer :: info = 0
   end type nlints_varmlp_t

   type, public :: nlints_causality_test_t
      !! Linear or neural Granger non-causality test result.
      real(dp) :: gci = 0.0_dp
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: numerator_df = 0
      integer :: denominator_df = 0
      integer :: lag = 0
      integer :: info = 0
   end type nlints_causality_test_t

   type, public :: nlints_neural_causality_t
      !! Neural Granger result and its restricted and unrestricted networks.
      type(nlints_causality_test_t) :: test
      type(nlints_varmlp_t) :: univariate
      type(nlints_varmlp_t) :: bivariate
   end type nlints_neural_causality_t

   interface display
      module procedure display_nlints_varmlp
      module procedure display_nlints_causality
      module procedure display_nlints_neural_causality
   end interface display

   public :: nlints_entropy_discrete, nlints_joint_entropy_discrete
   public :: nlints_mutual_information_discrete
   public :: nlints_multivariate_information_discrete
   public :: nlints_transfer_entropy_discrete
   public :: nlints_entropy_continuous, nlints_joint_entropy_continuous
   public :: nlints_mutual_information_continuous
   public :: nlints_transfer_entropy_continuous
   public :: nlints_varmlp_fit, nlints_varmlp_update
   public :: nlints_varmlp_predict, nlints_varmlp_forecast
   public :: nlints_varmlp_forecast_table
   public :: nlints_varmlp_save, nlints_varmlp_load
   public :: nlints_granger_test, nlints_neural_granger_test
   public :: nlints_df_test, display

contains

   pure real(dp) function nlints_entropy_discrete(values, logarithm) &
      result(entropy)
      !! Estimate Shannon entropy from integer category counts.
      integer, intent(in) :: values(:) !! Integer-valued observations.
      character(len=*), intent(in), optional :: logarithm !! Logarithm base name.
      logical, allocatable :: counted(:)
      real(dp) :: probability
      integer :: observation, matches

      entropy = 0.0_dp
      if (size(values) < 1) return
      allocate(counted(size(values)), source=.false.)
      do observation = 1, size(values)
         if (counted(observation)) cycle
         matches = count(values == values(observation))
         counted = counted .or. values == values(observation)
         probability = real(matches, dp)/real(size(values), dp)
         entropy = entropy - probability*selected_log(probability, logarithm)
      end do
   end function nlints_entropy_discrete

   pure real(dp) function nlints_joint_entropy_discrete(values, logarithm) &
      result(entropy)
      !! Estimate joint Shannon entropy from integer observation rows.
      integer, intent(in) :: values(:, :) !! Observation-by-variable categories.
      character(len=*), intent(in), optional :: logarithm !! Logarithm base name.
      logical, allocatable :: counted(:)
      real(dp) :: probability
      integer :: observation, other, matches

      entropy = 0.0_dp
      if (size(values, 1) < 1 .or. size(values, 2) < 1) return
      allocate(counted(size(values, 1)), source=.false.)
      do observation = 1, size(values, 1)
         if (counted(observation)) cycle
         matches = 0
         do other = observation, size(values, 1)
            if (all(values(other, :) == values(observation, :))) then
               matches = matches + 1
               counted(other) = .true.
            end if
         end do
         probability = real(matches, dp)/real(size(values, 1), dp)
         entropy = entropy - probability*selected_log(probability, logarithm)
      end do
   end function nlints_joint_entropy_discrete

   pure real(dp) function nlints_mutual_information_discrete(first, second, &
      logarithm, normalize) result(information)
      !! Estimate discrete bivariate mutual information.
      integer, intent(in) :: first(:) !! First categorical variable.
      integer, intent(in) :: second(:) !! Second categorical variable.
      character(len=*), intent(in), optional :: logarithm !! Logarithm base name.
      logical, intent(in), optional :: normalize !! Divide by the larger marginal entropy.
      integer, allocatable :: joint(:, :)
      real(dp) :: first_entropy, second_entropy, denominator
      logical :: use_normalization

      information = 0.0_dp
      if (size(first) < 1 .or. size(first) /= size(second)) return
      allocate(joint(size(first), 2))
      joint(:, 1) = first
      joint(:, 2) = second
      first_entropy = nlints_entropy_discrete(first, logarithm)
      second_entropy = nlints_entropy_discrete(second, logarithm)
      information = first_entropy + second_entropy - &
         nlints_joint_entropy_discrete(joint, logarithm)
      use_normalization = .false.
      if (present(normalize)) use_normalization = normalize
      denominator = max(first_entropy, second_entropy)
      if (use_normalization .and. denominator > 0.0_dp) &
         information = information/denominator
   end function nlints_mutual_information_discrete

   pure real(dp) function nlints_multivariate_information_discrete(values, &
      logarithm, normalize) result(information)
      !! Estimate multivariate interaction information by inclusion-exclusion.
      integer, intent(in) :: values(:, :) !! Observation-by-variable categories.
      character(len=*), intent(in), optional :: logarithm !! Logarithm base name.
      logical, intent(in), optional :: normalize !! Divide by the largest marginal entropy.
      integer, allocatable :: subset(:, :)
      integer :: mask, variable, variables, columns, column
      real(dp) :: entropy, largest
      logical :: use_normalization

      information = 0.0_dp
      variables = size(values, 2)
      if (size(values, 1) < 1 .or. variables < 1 .or. variables > 20) return
      largest = 0.0_dp
      do variable = 1, variables
         largest = max(largest, nlints_entropy_discrete(values(:, variable), &
            logarithm))
      end do
      do mask = 1, 2**variables - 1
         columns = popcnt(mask)
         allocate(subset(size(values, 1), columns))
         column = 0
         do variable = 1, variables
            if (.not. btest(mask, variable - 1)) cycle
            column = column + 1
            subset(:, column) = values(:, variable)
         end do
         entropy = nlints_joint_entropy_discrete(subset, logarithm)
         if (mod(columns, 2) == 1) then
            information = information + entropy
         else
            information = information - entropy
         end if
         deallocate(subset)
      end do
      use_normalization = .false.
      if (present(normalize)) use_normalization = normalize
      if (use_normalization .and. largest > 0.0_dp) &
         information = information/largest
   end function nlints_multivariate_information_discrete

   pure real(dp) function nlints_transfer_entropy_discrete(target, source, &
      target_lag, source_lag, logarithm, normalize) result(entropy)
      !! Estimate discrete transfer entropy from source to target.
      integer, intent(in) :: target(:) !! Target categorical time series.
      integer, intent(in) :: source(:) !! Source categorical time series.
      integer, intent(in), optional :: target_lag !! Target history order.
      integer, intent(in), optional :: source_lag !! Source history order.
      character(len=*), intent(in), optional :: logarithm !! Logarithm base name.
      logical, intent(in), optional :: normalize !! Divide by target range entropy.
      integer, allocatable :: target_past(:, :), target_present(:, :)
      integer, allocatable :: source_past(:, :), past_joint(:, :), full_joint(:, :)
      integer :: p, q, maximum_lag, observations, row, lag, time
      real(dp) :: denominator
      logical :: use_normalization

      entropy = 0.0_dp
      p = 1
      q = 1
      if (present(target_lag)) p = target_lag
      if (present(source_lag)) q = source_lag
      if (size(target) /= size(source) .or. p < 1 .or. q < 1) return
      maximum_lag = max(p, q)
      observations = size(target) - maximum_lag
      if (observations < 1) return
      allocate(target_past(observations, p))
      allocate(target_present(observations, p + 1))
      allocate(source_past(observations, q))
      do row = 1, observations
         time = maximum_lag + row
         target_present(row, 1) = target(time)
         do lag = 1, p
            target_past(row, lag) = target(time - lag)
            target_present(row, lag + 1) = target(time - lag)
         end do
         do lag = 1, q
            source_past(row, lag) = source(time - lag)
         end do
      end do
      allocate(past_joint(observations, p + q))
      allocate(full_joint(observations, p + q + 1))
      past_joint(:, :p) = target_past
      past_joint(:, p + 1:) = source_past
      full_joint(:, :p + 1) = target_present
      full_joint(:, p + 2:) = source_past
      denominator = nlints_joint_entropy_discrete(target_present, logarithm) - &
         nlints_joint_entropy_discrete(target_past, logarithm)
      entropy = denominator - nlints_joint_entropy_discrete(full_joint, &
         logarithm) + nlints_joint_entropy_discrete(past_joint, logarithm)
      use_normalization = .false.
      if (present(normalize)) use_normalization = normalize
      if (use_normalization) then
         denominator = selected_log(real(maxval(target) - minval(target), dp), &
            logarithm)
         if (denominator /= 0.0_dp) entropy = entropy/denominator
      end if
   end function nlints_transfer_entropy_discrete

   pure real(dp) function nlints_entropy_continuous(values, neighbors, &
      logarithm) result(entropy)
      !! Estimate continuous entropy with the Kozachenko nearest-neighbor formula.
      real(dp), intent(in) :: values(:) !! Continuous observations.
      integer, intent(in), optional :: neighbors !! Neighbor order.
      character(len=*), intent(in), optional :: logarithm !! Logarithm base name.
      real(dp), allocatable :: matrix(:, :)
      integer :: k

      entropy = quiet_nan()
      k = 3
      if (present(neighbors)) k = neighbors
      if (size(values) <= k .or. k < 1 .or. &
         .not. all(ieee_is_finite(values))) return
      allocate(matrix(size(values), 1))
      matrix(:, 1) = values
      entropy = nlints_joint_entropy_continuous(matrix, k, logarithm)
   end function nlints_entropy_continuous

   pure real(dp) function nlints_joint_entropy_continuous(values, neighbors, &
      logarithm) result(entropy)
      !! Estimate multivariate continuous entropy with maximum-norm distances.
      real(dp), intent(in) :: values(:, :) !! Observation-by-variable values.
      integer, intent(in), optional :: neighbors !! Neighbor order.
      character(len=*), intent(in), optional :: logarithm !! Logarithm base name.
      real(dp), allocatable :: distances(:)
      integer :: k, observation

      entropy = quiet_nan()
      k = 3
      if (present(neighbors)) k = neighbors
      if (size(values, 1) <= k .or. size(values, 2) < 1 .or. k < 1 .or. &
         .not. all(ieee_is_finite(values))) return
      distances = kth_neighbor_distances(values, k)
      if (any(distances <= 0.0_dp)) return
      entropy = digamma_value(real(size(values, 1), dp)) - &
         digamma_value(real(k, dp))
      do observation = 1, size(values, 1)
         entropy = entropy + real(size(values, 2), dp)* &
            selected_log(2.0_dp*distances(observation), logarithm)/ &
            real(size(values, 1), dp)
      end do
   end function nlints_joint_entropy_continuous

   pure real(dp) function nlints_mutual_information_continuous(first, second, &
      neighbors, algorithm, normalize) result(information)
      !! Estimate continuous mutual information with a KSG estimator.
      real(dp), intent(in) :: first(:) !! First continuous variable.
      real(dp), intent(in) :: second(:) !! Second continuous variable.
      integer, intent(in), optional :: neighbors !! Neighbor order.
      integer, intent(in), optional :: algorithm !! One for KSG-1 or two for KSG-2.
      logical, intent(in), optional :: normalize !! Divide by estimated joint entropy.
      real(dp), allocatable :: joint(:, :), radius(:), marginal_radius(:)
      integer, allocatable :: first_count(:), second_count(:)
      integer :: k, method, observation
      real(dp) :: total, joint_entropy
      logical :: use_normalization

      information = quiet_nan()
      k = 3
      method = 1
      if (present(neighbors)) k = neighbors
      if (present(algorithm)) method = algorithm
      if (size(first) /= size(second) .or. size(first) <= k .or. k < 1 .or. &
         (method /= 1 .and. method /= 2) .or. &
         .not. all(ieee_is_finite(first)) .or. &
         .not. all(ieee_is_finite(second))) return
      allocate(joint(size(first), 2))
      joint(:, 1) = first
      joint(:, 2) = second
      radius = kth_neighbor_distances(joint, k)
      if (method == 1) then
         first_count = neighbor_counts_vector(first, radius, .false.)
         second_count = neighbor_counts_vector(second, radius, .false.)
         total = 0.0_dp
         do observation = 1, size(first)
            total = total + digamma_value(real(first_count(observation) + 1, dp)) + &
               digamma_value(real(second_count(observation) + 1, dp))
         end do
         information = digamma_value(real(k, dp)) + &
            digamma_value(real(size(first), dp)) - total/real(size(first), dp)
      else
         allocate(marginal_radius(size(first)))
         marginal_radius = kth_neighbor_distances(reshape(first, &
            [size(first), 1]), k)
         first_count = neighbor_counts_vector(first, marginal_radius, .true.)
         marginal_radius = kth_neighbor_distances(reshape(second, &
            [size(second), 1]), k)
         second_count = neighbor_counts_vector(second, marginal_radius, .true.)
         total = 0.0_dp
         do observation = 1, size(first)
            total = total + digamma_value(real(first_count(observation), dp)) + &
               digamma_value(real(second_count(observation), dp))
         end do
         information = digamma_value(real(k, dp)) - 1.0_dp/real(k, dp) + &
            digamma_value(real(size(first), dp)) - total/real(size(first), dp)
      end if
      use_normalization = .false.
      if (present(normalize)) use_normalization = normalize
      if (use_normalization) then
         joint_entropy = sum(selected_log(2.0_dp*radius, 'log2'))/ &
            real(size(first), dp) + digamma_value(real(size(first), dp)) - &
            digamma_value(real(k, dp))
         if (joint_entropy /= 0.0_dp) information = information/joint_entropy
      end if
   end function nlints_mutual_information_continuous

   pure real(dp) function nlints_transfer_entropy_continuous(target, source, &
      target_lag, source_lag, neighbors, normalize) result(entropy)
      !! Estimate continuous transfer entropy with the package KSG formula.
      real(dp), intent(in) :: target(:) !! Target continuous time series.
      real(dp), intent(in) :: source(:) !! Source continuous time series.
      integer, intent(in), optional :: target_lag !! Target history order.
      integer, intent(in), optional :: source_lag !! Source history order.
      integer, intent(in), optional :: neighbors !! Neighbor order.
      logical, intent(in), optional :: normalize !! Apply the package entropy normalization.
      real(dp), allocatable :: target_past(:, :), target_present(:, :)
      real(dp), allocatable :: source_past(:, :), past_joint(:, :), full_joint(:, :)
      real(dp), allocatable :: radius(:)
      integer, allocatable :: target_count(:), past_count(:), present_count(:)
      integer :: p, q, k, maximum_lag, observations, row, lag, time
      real(dp) :: total, denominator, maximum_entropy
      logical :: use_normalization

      entropy = quiet_nan()
      p = 1
      q = 1
      k = 3
      if (present(target_lag)) p = target_lag
      if (present(source_lag)) q = source_lag
      if (present(neighbors)) k = neighbors
      maximum_lag = max(p, q)
      observations = size(target) - maximum_lag
      if (size(target) /= size(source) .or. p < 1 .or. q < 1 .or. &
         observations <= k .or. k < 1 .or. &
         .not. all(ieee_is_finite(target)) .or. &
         .not. all(ieee_is_finite(source))) return
      allocate(target_past(observations, p))
      allocate(target_present(observations, p + 1))
      allocate(source_past(observations, q))
      do row = 1, observations
         time = maximum_lag + row
         target_present(row, 1) = target(time)
         do lag = 1, p
            target_past(row, lag) = target(time - lag)
            target_present(row, lag + 1) = target(time - lag)
         end do
         do lag = 1, q
            source_past(row, lag) = source(time - lag)
         end do
      end do
      allocate(past_joint(observations, p + q))
      allocate(full_joint(observations, p + q + 1))
      past_joint(:, :p) = target_past
      past_joint(:, p + 1:) = source_past
      full_joint(:, :p + 1) = target_present
      full_joint(:, p + 2:) = source_past
      radius = kth_neighbor_distances(full_joint, k)
      target_count = neighbor_counts_matrix(target_past, radius, .false.)
      past_count = neighbor_counts_matrix(past_joint, radius, .false.)
      present_count = neighbor_counts_matrix(target_present, radius, .false.)
      total = 0.0_dp
      do row = 1, observations
         total = total + digamma_value(real(target_count(row) + 1, dp)) - &
            digamma_value(real(past_count(row) + 1, dp)) - &
            digamma_value(real(present_count(row) + 1, dp))
      end do
      entropy = digamma_value(real(k, dp)) + total/real(observations, dp)
      use_normalization = .false.
      if (present(normalize)) use_normalization = normalize
      if (use_normalization) then
         maximum_entropy = log(maxval(target) - minval(target))
         denominator = 0.0_dp
         do row = 1, observations
            denominator = denominator + log(2.0_dp*radius(row)) + &
               digamma_value(real(past_count(row) + 1, dp))
         end do
         denominator = maximum_entropy - denominator/real(observations, dp) + &
            digamma_value(real(k, dp))
         if (denominator /= 0.0_dp) entropy = entropy/denominator
      end if
   end function nlints_transfer_entropy_continuous

   pure function nlints_varmlp_fit(data, lag, hidden_counts, iterations, &
      learning_rate, optimizer, batch_size, bias, activations, seed) result(model)
      !! Fit the NlinTS multivariate autoregressive multilayer perceptron.
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable training data.
      integer, intent(in) :: lag !! Positive autoregressive order.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      integer, intent(in), optional :: iterations !! Training epochs.
      real(dp), intent(in), optional :: learning_rate !! SGD or Adam learning rate.
      integer, intent(in), optional :: optimizer !! SGD or Adam optimizer code.
      integer, intent(in), optional :: batch_size !! Mini-batch size.
      logical, intent(in), optional :: bias !! Include a bias in every layer.
      integer, intent(in), optional :: activations(:) !! Hidden and output activation codes.
      integer, intent(in), optional :: seed !! Deterministic weight initialization seed.
      type(nlints_varmlp_t) :: model
      real(dp), allocatable :: scaled(:, :), predictors(:, :), response(:, :)
      real(dp), allocatable :: prediction(:, :)
      integer :: epochs, batch, selected_seed, layer, input_count, output_count

      model%lag = lag
      model%variables = size(data, 2)
      model%hidden_counts = hidden_counts
      epochs = 50
      batch = 10
      selected_seed = 5
      if (present(iterations)) epochs = iterations
      if (present(batch_size)) batch = batch_size
      if (present(seed)) selected_seed = seed
      model%learning_rate = 0.01_dp
      if (present(learning_rate)) model%learning_rate = learning_rate
      model%optimizer = nlints_optimizer_sgd
      if (present(optimizer)) model%optimizer = optimizer
      model%bias = .true.
      if (present(bias)) model%bias = bias
      if (size(data, 1) <= lag .or. size(data, 2) < 1 .or. lag < 1 .or. &
         size(hidden_counts) < 1 .or. any(hidden_counts < 1) .or. epochs < 1 .or. &
         batch < 1 .or. model%learning_rate <= 0.0_dp .or. &
         (model%optimizer /= nlints_optimizer_sgd .and. &
         model%optimizer /= nlints_optimizer_adam) .or. &
         .not. all(ieee_is_finite(data))) then
         model%info = 1
         return
      end if
      allocate(model%activations(size(hidden_counts) + 1))
      model%activations(:size(hidden_counts)) = nlints_activation_relu
      model%activations(size(model%activations)) = nlints_activation_sigmoid
      if (present(activations)) then
         if (size(activations) /= size(model%activations) .or. &
            any(activations < nlints_activation_linear) .or. &
            any(activations > nlints_activation_tanh)) then
            model%info = 2
            return
         end if
         model%activations = activations
      end if
      call minmax_scale(data, scaled, model%minimum, model%maximum)
      call build_lagged_design(scaled, lag, predictors, response)
      allocate(model%layers(size(hidden_counts) + 1))
      input_count = size(predictors, 2)
      do layer = 1, size(model%layers)
         if (layer <= size(hidden_counts)) then
            output_count = hidden_counts(layer)
         else
            output_count = size(data, 2)
         end if
         call initialize_layer(model%layers(layer), input_count, output_count, &
            model%bias, model%activations(layer), selected_seed + 104729*layer)
         input_count = output_count
      end do
      call train_network(model, predictors, response, epochs, batch)
      prediction = network_predict(model, predictors)
      allocate(model%fitted(size(data, 1), size(data, 2)))
      allocate(model%residuals(size(data, 1), size(data, 2)))
      model%fitted = quiet_nan()
      model%residuals = quiet_nan()
      model%fitted(lag + 1:, :) = undo_minmax(prediction, model%minimum, &
         model%maximum)
      model%residuals(lag + 1:, :) = data(lag + 1:, :) - &
         model%fitted(lag + 1:, :)
      allocate(model%rss(size(data, 2)))
      do output_count = 1, size(data, 2)
         model%rss(output_count) = sum((prediction(:, output_count) - &
            response(:, output_count))**2)
      end do
      model%iterations = epochs
   end function nlints_varmlp_fit

   pure function nlints_varmlp_update(model, data, iterations, batch_size) &
      result(updated)
      !! Continue fitting an existing VARNN on additional data.
      type(nlints_varmlp_t), intent(in) :: model !! Existing fitted network.
      real(dp), intent(in) :: data(:, :) !! New observation-by-variable data.
      integer, intent(in), optional :: iterations !! Additional epochs.
      integer, intent(in), optional :: batch_size !! Mini-batch size.
      type(nlints_varmlp_t) :: updated
      real(dp), allocatable :: scaled(:, :), predictors(:, :), response(:, :)
      real(dp), allocatable :: prediction(:, :)
      integer :: epochs, batch, variable

      updated = model
      epochs = 50
      batch = 10
      if (present(iterations)) epochs = iterations
      if (present(batch_size)) batch = batch_size
      if (model%info /= 0 .or. size(data, 2) /= model%variables .or. &
         size(data, 1) <= model%lag .or. epochs < 1 .or. batch < 1 .or. &
         .not. all(ieee_is_finite(data))) then
         updated%info = 1
         return
      end if
      scaled = apply_minmax(data, model%minimum, model%maximum)
      call build_lagged_design(scaled, model%lag, predictors, response)
      call train_network(updated, predictors, response, epochs, batch)
      prediction = network_predict(updated, predictors)
      if (allocated(updated%fitted)) deallocate(updated%fitted)
      if (allocated(updated%residuals)) deallocate(updated%residuals)
      allocate(updated%fitted(size(data, 1), size(data, 2)))
      allocate(updated%residuals(size(data, 1), size(data, 2)))
      updated%fitted = quiet_nan()
      updated%residuals = quiet_nan()
      updated%fitted(model%lag + 1:, :) = undo_minmax(prediction, &
         model%minimum, model%maximum)
      updated%residuals(model%lag + 1:, :) = data(model%lag + 1:, :) - &
         updated%fitted(model%lag + 1:, :)
      if (allocated(updated%rss)) deallocate(updated%rss)
      allocate(updated%rss(model%variables))
      do variable = 1, model%variables
         updated%rss(variable) = sum((prediction(:, variable) - &
            response(:, variable))**2)
      end do
      updated%iterations = model%iterations + epochs
   end function nlints_varmlp_update

   pure function nlints_varmlp_predict(model, data) result(prediction)
      !! Return aligned one-step fitted values for supplied observations.
      type(nlints_varmlp_t), intent(in) :: model !! Fitted VARNN model.
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable histories.
      real(dp), allocatable :: prediction(:, :)
      real(dp), allocatable :: scaled(:, :), predictors(:, :), response(:, :)

      if (model%info /= 0 .or. size(data, 2) /= model%variables .or. &
         size(data, 1) <= model%lag .or. .not. all(ieee_is_finite(data))) then
         allocate(prediction(0, 0))
         return
      end if
      scaled = apply_minmax(data, model%minimum, model%maximum)
      call build_lagged_design(scaled, model%lag, predictors, response)
      prediction = undo_minmax(network_predict(model, predictors), &
         model%minimum, model%maximum)
   end function nlints_varmlp_predict

   pure function nlints_varmlp_forecast(model, data) result(forecast)
      !! Forecast the next observation recursively from the final lag block.
      type(nlints_varmlp_t), intent(in) :: model !! Fitted VARNN model.
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable history.
      real(dp), allocatable :: forecast(:)
      real(dp), allocatable :: scaled(:, :), predictor(:, :), normalized(:, :)
      real(dp), allocatable :: restored(:, :)
      integer :: variable, lag, column

      if (model%info /= 0 .or. size(data, 2) /= model%variables .or. &
         size(data, 1) < model%lag .or. .not. all(ieee_is_finite(data))) then
         allocate(forecast(0))
         return
      end if
      scaled = apply_minmax(data, model%minimum, model%maximum)
      allocate(predictor(1, model%variables*model%lag))
      column = 0
      do variable = 1, model%variables
         do lag = 1, model%lag
            column = column + 1
            predictor(1, column) = scaled(size(data, 1) + 1 - lag, variable)
         end do
      end do
      normalized = network_predict(model, predictor)
      restored = undo_minmax(normalized, model%minimum, model%maximum)
      forecast = restored(1, :)
   end function nlints_varmlp_forecast

   pure function nlints_varmlp_forecast_table(model, data) result(forecast)
      !! Return rolling predictions followed by the next-period forecast.
      type(nlints_varmlp_t), intent(in) :: model !! Fitted VARNN model.
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable history.
      real(dp), allocatable :: forecast(:, :)
      real(dp), allocatable :: scaled(:, :), predictor(:, :), normalized(:, :)
      real(dp), allocatable :: local_minimum(:), local_maximum(:)
      integer :: rows, variable, lag, column, row, time

      if (model%info /= 0 .or. size(data, 2) /= model%variables .or. &
         size(data, 1) < model%lag .or. .not. all(ieee_is_finite(data))) then
         allocate(forecast(0, 0))
         return
      end if
      call minmax_scale(data, scaled, local_minimum, local_maximum)
      rows = size(data, 1) - model%lag + 1
      allocate(predictor(rows, model%variables*model%lag))
      do row = 1, rows
         time = model%lag + row
         column = 0
         do variable = 1, model%variables
            do lag = 1, model%lag
               column = column + 1
               predictor(row, column) = scaled(time - lag, variable)
            end do
         end do
      end do
      normalized = network_predict(model, predictor)
      forecast = undo_minmax(normalized, local_minimum, local_maximum)
   end function nlints_varmlp_forecast_table

   subroutine nlints_varmlp_save(model, filename, info)
      !! Save a complete VARNN model in a versioned formatted text file.
      type(nlints_varmlp_t), intent(in) :: model !! Fitted VARNN model.
      character(len=*), intent(in) :: filename !! Destination filename.
      integer, intent(out) :: info !! Status code; zero indicates success.
      integer :: unit, layer, row, status

      info = 1
      if (model%info /= 0 .or. .not. allocated(model%layers) .or. &
         .not. allocated(model%minimum) .or. .not. allocated(model%maximum) .or. &
         .not. allocated(model%rss)) return
      open(newunit=unit, file=filename, status='replace', action='write', &
         iostat=status)
      if (status /= 0) then
         info = status
         return
      end if
      write(unit, '(a)', iostat=status) 'NLINTS_VARMLP 1'
      if (status == 0) write(unit, *, iostat=status) model%lag, &
         model%variables, model%iterations, model%optimizer, &
         model%learning_rate, model%bias
      if (status == 0) write(unit, *, iostat=status) size(model%hidden_counts)
      if (status == 0) write(unit, *, iostat=status) model%hidden_counts
      if (status == 0) write(unit, *, iostat=status) model%activations
      if (status == 0) write(unit, *, iostat=status) model%minimum
      if (status == 0) write(unit, *, iostat=status) model%maximum
      if (status == 0) write(unit, *, iostat=status) model%rss
      if (status == 0) write(unit, *, iostat=status) size(model%layers)
      do layer = 1, size(model%layers)
         if (status /= 0) exit
         write(unit, *, iostat=status) model%layers(layer)%activation, &
            size(model%layers(layer)%weights, 1), &
            size(model%layers(layer)%weights, 2)
         do row = 1, size(model%layers(layer)%weights, 1)
            if (status == 0) write(unit, *, iostat=status) &
               model%layers(layer)%weights(row, :)
         end do
         do row = 1, size(model%layers(layer)%first_moment, 1)
            if (status == 0) write(unit, *, iostat=status) &
               model%layers(layer)%first_moment(row, :)
         end do
         do row = 1, size(model%layers(layer)%second_moment, 1)
            if (status == 0) write(unit, *, iostat=status) &
               model%layers(layer)%second_moment(row, :)
         end do
      end do
      close(unit, iostat=info)
      if (status /= 0) info = status
   end subroutine nlints_varmlp_save

   function nlints_varmlp_load(filename) result(model)
      !! Load a VARNN model written by `nlints_varmlp_save`.
      character(len=*), intent(in) :: filename !! Source filename.
      type(nlints_varmlp_t) :: model
      character(len=32) :: header
      integer :: unit, status, hidden_layers, layer_count, layer, row
      integer :: outputs, inputs

      model%info = 1
      open(newunit=unit, file=filename, status='old', action='read', iostat=status)
      if (status /= 0) return
      read(unit, '(a)', iostat=status) header
      if (status /= 0 .or. trim(header) /= 'NLINTS_VARMLP 1') then
         close(unit)
         model%info = 2
         return
      end if
      read(unit, *, iostat=status) model%lag, model%variables, model%iterations, &
         model%optimizer, model%learning_rate, model%bias
      if (status == 0) read(unit, *, iostat=status) hidden_layers
      if (status /= 0 .or. hidden_layers < 1 .or. model%lag < 1 .or. &
         model%variables < 1) then
         close(unit)
         model%info = 3
         return
      end if
      allocate(model%hidden_counts(hidden_layers))
      allocate(model%activations(hidden_layers + 1))
      allocate(model%minimum(model%variables), model%maximum(model%variables))
      allocate(model%rss(model%variables))
      read(unit, *, iostat=status) model%hidden_counts
      if (status == 0) read(unit, *, iostat=status) model%activations
      if (status == 0) read(unit, *, iostat=status) model%minimum
      if (status == 0) read(unit, *, iostat=status) model%maximum
      if (status == 0) read(unit, *, iostat=status) model%rss
      if (status == 0) read(unit, *, iostat=status) layer_count
      if (status /= 0 .or. layer_count /= hidden_layers + 1) then
         close(unit)
         model%info = 3
         return
      end if
      allocate(model%layers(layer_count))
      do layer = 1, layer_count
         read(unit, *, iostat=status) model%layers(layer)%activation, outputs, inputs
         if (status /= 0 .or. outputs < 1 .or. inputs < 1) exit
         allocate(model%layers(layer)%weights(outputs, inputs))
         allocate(model%layers(layer)%first_moment(outputs, inputs))
         allocate(model%layers(layer)%second_moment(outputs, inputs))
         do row = 1, outputs
            read(unit, *, iostat=status) model%layers(layer)%weights(row, :)
            if (status /= 0) exit
         end do
         do row = 1, outputs
            read(unit, *, iostat=status) model%layers(layer)%first_moment(row, :)
            if (status /= 0) exit
         end do
         do row = 1, outputs
            read(unit, *, iostat=status) model%layers(layer)%second_moment(row, :)
            if (status /= 0) exit
         end do
         if (status /= 0) exit
      end do
      close(unit)
      if (status /= 0) then
         model%info = 4
         return
      end if
      allocate(model%fitted(0, model%variables))
      allocate(model%residuals(0, model%variables))
      model%info = 0
   end function nlints_varmlp_load

   pure function nlints_granger_test(target, source, lag, difference) result(out)
      !! Perform the package's classical pairwise Granger F test.
      real(dp), intent(in) :: target(:) !! Target time series.
      real(dp), intent(in) :: source(:) !! Candidate causing time series.
      integer, intent(in) :: lag !! Positive common lag order.
      logical, intent(in), optional :: difference !! Difference inferred integrated series.
      type(nlints_causality_test_t) :: out
      real(dp), allocatable :: restricted(:, :), unrestricted(:, :), response(:)
      real(dp), allocatable :: beta(:), standard_errors(:), residuals(:)
      real(dp), allocatable :: working_target(:), working_source(:)
      real(dp) :: restricted_rss, unrestricted_rss
      integer :: observations, row, column, info, target_order, source_order
      integer :: maximum_order
      logical :: use_differencing

      out%lag = lag
      if (size(target) /= size(source) .or. size(target) <= 2*lag + 1 .or. &
         lag < 1 .or. .not. all(ieee_is_finite(target)) .or. &
         .not. all(ieee_is_finite(source))) then
         out%info = 1
         return
      end if
      use_differencing = .false.
      if (present(difference)) use_differencing = difference
      if (use_differencing) then
         target_order = nlints_integration_order(target, lag)
         source_order = nlints_integration_order(source, lag)
         maximum_order = max(target_order, source_order)
         working_target = aligned_difference(target, target_order, maximum_order)
         working_source = aligned_difference(source, source_order, maximum_order)
      else
         working_target = target
         working_source = source
      end if
      if (size(working_target) <= 2*lag + 1) then
         out%info = 1
         return
      end if
      observations = size(working_target) - lag
      allocate(response(observations))
      allocate(restricted(observations, lag + 1))
      allocate(unrestricted(observations, 2*lag + 1))
      restricted(:, 1) = 1.0_dp
      unrestricted(:, 1) = 1.0_dp
      do row = 1, observations
         response(row) = working_target(lag + row)
         do column = 1, lag
            restricted(row, column + 1) = &
               working_target(lag + row - column)
            unrestricted(row, column + 1) = restricted(row, column + 1)
            unrestricted(row, lag + column + 1) = &
               working_source(lag + row - column)
         end do
      end do
      call ols_fit(restricted, response, beta, standard_errors, residuals, &
         restricted_rss, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      call ols_fit(unrestricted, response, beta, standard_errors, residuals, &
         unrestricted_rss, info)
      if (info /= 0 .or. unrestricted_rss <= 0.0_dp) then
         out%info = 2
         return
      end if
      out%gci = log(restricted_rss/unrestricted_rss)
      out%numerator_df = lag
      out%denominator_df = observations - 2*lag - 1
      out%statistic = ((restricted_rss - unrestricted_rss)/real(lag, dp))/ &
         (unrestricted_rss/real(out%denominator_df, dp))
      out%p_value = f_upper_probability(out%statistic, out%numerator_df, &
         out%denominator_df)
   end function nlints_granger_test

   pure function nlints_neural_granger_test(target, source, lag, &
      univariate_hidden, bivariate_hidden, iterations, learning_rate, &
      optimizer, batch_size, bias, seed) result(out)
      !! Compare restricted and bivariate VARNNs for nonlinear Granger causality.
      real(dp), intent(in) :: target(:) !! Target time series.
      real(dp), intent(in) :: source(:) !! Candidate causing time series.
      integer, intent(in) :: lag !! Positive common lag order.
      integer, intent(in) :: univariate_hidden(:) !! Restricted hidden-layer sizes.
      integer, intent(in) :: bivariate_hidden(:) !! Bivariate hidden-layer sizes.
      integer, intent(in), optional :: iterations !! Training epochs.
      real(dp), intent(in), optional :: learning_rate !! Network learning rate.
      integer, intent(in), optional :: optimizer !! SGD or Adam optimizer code.
      integer, intent(in), optional :: batch_size !! Mini-batch size.
      logical, intent(in), optional :: bias !! Include network biases.
      integer, intent(in), optional :: seed !! Deterministic initialization seed.
      type(nlints_neural_causality_t) :: out
      real(dp), allocatable :: univariate_data(:, :), bivariate_data(:, :)
      real(dp) :: restricted_rss, unrestricted_rss
      integer :: restricted_parameters, unrestricted_parameters, observations

      if (size(target) /= size(source) .or. size(target) <= lag + 2 .or. &
         lag < 1) then
         out%test%info = 1
         return
      end if
      allocate(univariate_data(size(target), 1))
      allocate(bivariate_data(size(target), 2))
      univariate_data(:, 1) = target
      bivariate_data(:, 1) = target
      bivariate_data(:, 2) = source
      out%univariate = nlints_varmlp_fit(univariate_data, lag, &
         univariate_hidden, iterations, learning_rate, optimizer, batch_size, &
         bias, seed=seed)
      out%bivariate = nlints_varmlp_fit(bivariate_data, lag, bivariate_hidden, &
         iterations, learning_rate, optimizer, batch_size, bias, seed=seed)
      if (out%univariate%info /= 0 .or. out%bivariate%info /= 0) then
         out%test%info = 2
         return
      end if
      restricted_rss = out%univariate%rss(1)
      unrestricted_rss = out%bivariate%rss(1)
      out%test%lag = lag
      if (unrestricted_rss >= restricted_rss) then
         out%test%gci = 0.0_dp
      else
         out%test%gci = log(restricted_rss/unrestricted_rss)
      end if
      restricted_parameters = causality_parameter_count(univariate_hidden, &
         lag + 1, out%univariate%bias)
      unrestricted_parameters = causality_parameter_count(bivariate_hidden, &
         2*lag + 1, out%bivariate%bias)
      observations = size(target) - lag
      unrestricted_parameters = min(unrestricted_parameters, observations)
      out%test%numerator_df = unrestricted_parameters - restricted_parameters
      out%test%denominator_df = observations - unrestricted_parameters
      if (out%test%numerator_df <= 0 .or. out%test%denominator_df <= 0 .or. &
         unrestricted_rss <= 0.0_dp) then
         out%test%info = 3
         out%test%statistic = quiet_nan()
         out%test%p_value = quiet_nan()
         return
      end if
      out%test%statistic = (restricted_rss - unrestricted_rss)/ &
         unrestricted_rss*real(out%test%denominator_df, dp)/ &
         real(out%test%numerator_df, dp)
      out%test%p_value = f_upper_probability(out%test%statistic, &
         out%test%numerator_df, out%test%denominator_df)
   end function nlints_neural_granger_test

   pure function nlints_df_test(series, lag) result(out)
      !! Reuse the shared trend ADF test for NlinTS df.test compatibility.
      real(dp), intent(in) :: series(:) !! Univariate observations.
      integer, intent(in) :: lag !! Augmentation lag order.
      type(adf_result_t) :: out

      out = adf_test(series, model='trend', lags=lag)
   end function nlints_df_test

   pure integer function nlints_integration_order(series, lag) result(order)
      !! Infer integration order using repeated NlinTS trend ADF decisions.
      real(dp), intent(in) :: series(:) !! Univariate observations.
      integer, intent(in) :: lag !! Fixed ADF augmentation order.
      type(adf_result_t) :: test
      real(dp), allocatable :: working(:)
      real(dp) :: critical_value

      order = 0
      working = series
      do while (size(working) > lag + 3)
         test = nlints_df_test(working, lag)
         if (test%info /= 0 .or. .not. allocated(test%statistic)) exit
         critical_value = nlints_adf_five_percent(size(working) - 1)
         if (test%statistic(1) < critical_value) exit
         working = working(2:) - working(:size(working) - 1)
         order = order + 1
      end do
   end function nlints_integration_order

   pure function aligned_difference(series, order, maximum_order) result(values)
      !! Difference by one order and trim to a common maximum-order sample.
      real(dp), intent(in) :: series(:) !! Original observations.
      integer, intent(in) :: order !! Differencing order for this series.
      integer, intent(in) :: maximum_order !! Largest order across paired series.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: working(:)
      integer :: difference_order, retained

      working = series
      do difference_order = 1, order
         working = working(2:) - working(:size(working) - 1)
      end do
      retained = size(series) - maximum_order
      values = working(size(working) - retained + 1:)
   end function aligned_difference

   pure real(dp) function nlints_adf_five_percent(observations) result(critical)
      !! Return the NlinTS finite-sample five-percent trend ADF critical value.
      integer, intent(in) :: observations !! Number of differenced observations.

      if (observations <= 25) then
         critical = -3.60_dp
      else if (observations <= 50) then
         critical = -3.50_dp
      else if (observations <= 100) then
         critical = -3.45_dp
      else if (observations <= 250) then
         critical = -3.43_dp
      else if (observations <= 500) then
         critical = -3.42_dp
      else
         critical = -3.41_dp
      end if
   end function nlints_adf_five_percent

   subroutine display_nlints_varmlp(model, unit, print_obs)
      !! Display a fitted VARNN and optionally its aligned observations.
      type(nlints_varmlp_t), intent(in) :: model !! Fitted VARNN model.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Print fitted values and residuals.
      integer :: destination, layer
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'NlinTS VARNN fit'
      write(destination, '(a,i0)') '  lag: ', model%lag
      write(destination, '(a,i0)') '  variables: ', model%variables
      write(destination, '(a,i0)') '  iterations: ', model%iterations
      write(destination, '(a,*(i0,1x))') '  hidden layers: ', model%hidden_counts
      write(destination, '(a,*(es12.4,1x))') '  normalized RSS: ', model%rss
      do layer = 1, size(model%layers)
         write(destination, '(a,i0,a,i0,a,i0)') '  layer ', layer, ': ', &
            size(model%layers(layer)%weights, 2) - merge(1, 0, model%bias), &
            ' -> ', size(model%layers(layer)%weights, 1)
      end do
      if (show_observations) then
         write(destination, '(a)') '  fitted values:'
         do layer = 1, size(model%fitted, 1)
            write(destination, '(*(es14.6,1x))') model%fitted(layer, :)
         end do
         write(destination, '(a)') '  residuals:'
         do layer = 1, size(model%residuals, 1)
            write(destination, '(*(es14.6,1x))') model%residuals(layer, :)
         end do
      end if
   end subroutine display_nlints_varmlp

   subroutine display_nlints_causality(test, unit)
      !! Display a linear or nonlinear Granger test result.
      type(nlints_causality_test_t), intent(in) :: test !! Causality test result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'NlinTS Granger causality test'
      write(destination, '(a,i0)') '  lag: ', test%lag
      write(destination, '(a,es14.6)') '  causality index: ', test%gci
      write(destination, '(a,es14.6)') '  F statistic: ', test%statistic
      write(destination, '(a,i0,a,i0)') '  degrees of freedom: ', &
         test%numerator_df, ', ', test%denominator_df
      write(destination, '(a,es14.6)') '  p-value: ', test%p_value
   end subroutine display_nlints_causality

   subroutine display_nlints_neural_causality(result, unit)
      !! Display a neural Granger test without printing its training samples.
      type(nlints_neural_causality_t), intent(in) :: result !! Neural causality result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.

      call display_nlints_causality(result%test, unit)
   end subroutine display_nlints_neural_causality

   pure subroutine minmax_scale(data, scaled, minimum, maximum)
      !! Scale every data column to the unit interval.
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable data.
      real(dp), allocatable, intent(out) :: scaled(:, :) !! Unit-interval data.
      real(dp), allocatable, intent(out) :: minimum(:) !! Column minima.
      real(dp), allocatable, intent(out) :: maximum(:) !! Column maxima.
      integer :: variable

      allocate(scaled, source=data)
      allocate(minimum(size(data, 2)), maximum(size(data, 2)))
      do variable = 1, size(data, 2)
         minimum(variable) = minval(data(:, variable))
         maximum(variable) = maxval(data(:, variable))
         if (maximum(variable) > minimum(variable)) then
            scaled(:, variable) = (data(:, variable) - minimum(variable))/ &
               (maximum(variable) - minimum(variable))
         else
            scaled(:, variable) = 0.0_dp
         end if
      end do
   end subroutine minmax_scale

   pure function apply_minmax(data, minimum, maximum) result(scaled)
      !! Apply retained columnwise min-max scaling.
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable data.
      real(dp), intent(in) :: minimum(:) !! Retained column minima.
      real(dp), intent(in) :: maximum(:) !! Retained column maxima.
      real(dp) :: scaled(size(data, 1), size(data, 2))
      integer :: variable

      scaled = data
      do variable = 1, size(data, 2)
         if (maximum(variable) > minimum(variable)) then
            scaled(:, variable) = (data(:, variable) - minimum(variable))/ &
               (maximum(variable) - minimum(variable))
         else
            scaled(:, variable) = 0.0_dp
         end if
      end do
   end function apply_minmax

   pure function undo_minmax(data, minimum, maximum) result(values)
      !! Restore original units after columnwise min-max scaling.
      real(dp), intent(in) :: data(:, :) !! Unit-interval values.
      real(dp), intent(in) :: minimum(:) !! Retained column minima.
      real(dp), intent(in) :: maximum(:) !! Retained column maxima.
      real(dp) :: values(size(data, 1), size(data, 2))
      integer :: variable

      do variable = 1, size(data, 2)
         values(:, variable) = minimum(variable) + data(:, variable)* &
            (maximum(variable) - minimum(variable))
      end do
   end function undo_minmax

   pure subroutine build_lagged_design(data, lag_order, predictors, response)
      !! Construct the NlinTS variable-major autoregressive design.
      real(dp), intent(in) :: data(:, :) !! Scaled observation matrix.
      integer, intent(in) :: lag_order !! Autoregressive order.
      real(dp), allocatable, intent(out) :: predictors(:, :) !! Lagged predictors.
      real(dp), allocatable, intent(out) :: response(:, :) !! Current responses.
      integer :: observations, variable, lag, column

      observations = size(data, 1) - lag_order
      allocate(predictors(observations, size(data, 2)*lag_order))
      allocate(response(observations, size(data, 2)))
      response = data(lag_order + 1:, :)
      column = 0
      do variable = 1, size(data, 2)
         do lag = 1, lag_order
            column = column + 1
            predictors(:, column) = data(lag_order + 1 - lag: &
               size(data, 1) - lag, variable)
         end do
      end do
   end subroutine build_lagged_design

   pure subroutine initialize_layer(layer, inputs, outputs, bias, activation, seed)
      !! Allocate deterministic normally shaped starting weights for one layer.
      type(nlints_dense_layer_t), intent(out) :: layer !! Initialized dense layer.
      integer, intent(in) :: inputs !! Number of incoming activations.
      integer, intent(in) :: outputs !! Number of outgoing activations.
      logical, intent(in) :: bias !! Include a leading bias weight.
      integer, intent(in) :: activation !! Activation-function code.
      integer, intent(in) :: seed !! Deterministic initialization seed.
      real(dp) :: scale
      integer :: output, input

      allocate(layer%weights(outputs, inputs + merge(1, 0, bias)))
      allocate(layer%first_moment, mold=layer%weights)
      allocate(layer%second_moment, mold=layer%weights)
      scale = 2.0_dp/real(inputs + outputs + merge(1, 0, bias), dp)
      do output = 1, outputs
         do input = 1, size(layer%weights, 2)
            layer%weights(output, input) = scale*sin(real(seed + &
               7919*output + 104729*input, dp)*0.0174532925199433_dp)
         end do
      end do
      layer%first_moment = 0.0_dp
      layer%second_moment = 0.0_dp
      layer%activation = activation
   end subroutine initialize_layer

   pure subroutine train_network(model, predictors, response, epochs, batch_size)
      !! Train all dense layers by mini-batch backpropagation.
      type(nlints_varmlp_t), intent(inout) :: model !! Network updated in place.
      real(dp), intent(in) :: predictors(:, :) !! Scaled predictor rows.
      real(dp), intent(in) :: response(:, :) !! Scaled response rows.
      integer, intent(in) :: epochs !! Number of passes over the data.
      integer, intent(in) :: batch_size !! Requested mini-batch size.
      type(real_vector_t), allocatable :: activation(:), preactivation(:), delta(:)
      type(real_matrix_t), allocatable :: gradient(:)
      real(dp), allocatable :: input(:), next_weights(:, :)
      real(dp) :: beta1_correction, beta2_correction, gradient_value
      integer :: epoch, observation, layer, batch_count, update, offset
      integer :: row

      allocate(activation(0:size(model%layers)))
      allocate(preactivation(size(model%layers)))
      allocate(delta(size(model%layers)))
      allocate(gradient(size(model%layers)))
      do layer = 1, size(model%layers)
         allocate(gradient(layer)%value, mold=model%layers(layer)%weights)
         gradient(layer)%value = 0.0_dp
      end do
      update = 0
      batch_count = 0
      do epoch = 1, epochs
         do observation = 1, size(predictors, 1)
            row = 1 + mod(observation - 1 + epoch - 1, size(predictors, 1))
            activation(0)%value = predictors(row, :)
            do layer = 1, size(model%layers)
               input = activation(layer - 1)%value
               if (model%bias) input = [1.0_dp, input]
               preactivation(layer)%value = matmul( &
                  model%layers(layer)%weights, input)
               activation(layer)%value = activation_value_vector( &
                  preactivation(layer)%value, model%layers(layer)%activation)
            end do
            layer = size(model%layers)
            delta(layer)%value = (activation(layer)%value - response(row, :))* &
               activation_derivative(preactivation(layer)%value, &
               model%layers(layer)%activation)
            do layer = size(model%layers) - 1, 1, -1
               offset = merge(2, 1, model%bias)
               next_weights = model%layers(layer + 1)%weights(:, offset:)
               delta(layer)%value = matmul(transpose(next_weights), &
                  delta(layer + 1)%value)*activation_derivative( &
                  preactivation(layer)%value, model%layers(layer)%activation)
            end do
            do layer = 1, size(model%layers)
               input = activation(layer - 1)%value
               if (model%bias) input = [1.0_dp, input]
               gradient(layer)%value = gradient(layer)%value + &
                  spread(delta(layer)%value, 2, size(input))* &
                  spread(input, 1, size(delta(layer)%value))
            end do
            batch_count = batch_count + 1
            if (batch_count < batch_size .and. &
               observation < size(predictors, 1)) cycle
            update = update + 1
            beta1_correction = 1.0_dp - 0.9_dp**update
            beta2_correction = 1.0_dp - 0.999_dp**update
            do layer = 1, size(model%layers)
               gradient(layer)%value = gradient(layer)%value/real(batch_count, dp)
               if (model%optimizer == nlints_optimizer_adam) then
                  model%layers(layer)%first_moment = &
                     0.9_dp*model%layers(layer)%first_moment + &
                     0.1_dp*gradient(layer)%value
                  model%layers(layer)%second_moment = &
                     0.999_dp*model%layers(layer)%second_moment + &
                     0.001_dp*gradient(layer)%value**2
                  model%layers(layer)%weights = model%layers(layer)%weights - &
                     model%learning_rate*(model%layers(layer)%first_moment/ &
                     beta1_correction)/(sqrt(model%layers(layer)%second_moment/ &
                     beta2_correction) + 1.0e-8_dp)
               else
                  model%layers(layer)%weights = model%layers(layer)%weights - &
                     model%learning_rate*gradient(layer)%value
               end if
               gradient(layer)%value = 0.0_dp
            end do
            batch_count = 0
         end do
      end do
   end subroutine train_network

   pure function network_predict(model, predictors) result(prediction)
      !! Evaluate all rows with a fitted NlinTS dense network.
      type(nlints_varmlp_t), intent(in) :: model !! Fitted neural model.
      real(dp), intent(in) :: predictors(:, :) !! Scaled predictor rows.
      real(dp), allocatable :: prediction(:, :)
      real(dp), allocatable :: activation(:, :), input(:, :)
      integer :: layer

      activation = predictors
      do layer = 1, size(model%layers)
         if (model%bias) then
            allocate(input(size(activation, 1), size(activation, 2) + 1))
            input(:, 1) = 1.0_dp
            input(:, 2:) = activation
         else
            input = activation
         end if
         activation = activation_value(matmul(input, &
            transpose(model%layers(layer)%weights)), &
            model%layers(layer)%activation)
         deallocate(input)
      end do
      prediction = activation
   end function network_predict

   pure elemental real(dp) function activation_scalar(value, activation) &
      result(output)
      !! Apply one supported dense-layer activation.
      real(dp), intent(in) :: value !! Preactivation value.
      integer, intent(in) :: activation !! Activation-function code.

      select case (activation)
      case (nlints_activation_sigmoid)
         if (value >= 0.0_dp) then
            output = 1.0_dp/(1.0_dp + exp(-value))
         else
            output = exp(value)/(1.0_dp + exp(value))
         end if
      case (nlints_activation_relu)
         output = max(0.0_dp, value)
      case (nlints_activation_tanh)
         output = tanh(value)
      case default
         output = value
      end select
   end function activation_scalar

   pure function activation_value(value, activation) result(output)
      !! Apply one activation code to an array.
      real(dp), intent(in) :: value(:, :) !! Preactivation matrix.
      integer, intent(in) :: activation !! Activation-function code.
      real(dp) :: output(size(value, 1), size(value, 2))

      output = activation_scalar(value, activation)
   end function activation_value

   pure function activation_value_vector(value, activation) result(output)
      !! Apply one activation code to a vector.
      real(dp), intent(in) :: value(:) !! Preactivation vector.
      integer, intent(in) :: activation !! Activation-function code.
      real(dp) :: output(size(value))

      output = activation_scalar(value, activation)
   end function activation_value_vector

   pure function activation_derivative(value, activation) result(derivative)
      !! Evaluate an activation derivative from preactivation values.
      real(dp), intent(in) :: value(:) !! Preactivation values.
      integer, intent(in) :: activation !! Activation-function code.
      real(dp) :: derivative(size(value)), activated(size(value))

      select case (activation)
      case (nlints_activation_sigmoid)
         activated = activation_scalar(value, activation)
         derivative = activated*(1.0_dp - activated)
      case (nlints_activation_relu)
         derivative = merge(1.0_dp, 0.0_dp, value > 0.0_dp)
      case (nlints_activation_tanh)
         derivative = 1.0_dp - tanh(value)**2
      case default
         derivative = 1.0_dp
      end select
   end function activation_derivative

   pure integer function causality_parameter_count(hidden_counts, first_width, &
      bias) result(count_parameters)
      !! Reproduce NlinTS's neural-test parameter-count convention.
      integer, intent(in) :: hidden_counts(:) !! Hidden-layer sizes.
      integer, intent(in) :: first_width !! Package first-layer width count.
      logical, intent(in) :: bias !! Whether layer biases are enabled.
      integer :: layer, previous

      count_parameters = 0
      previous = first_width
      do layer = 1, size(hidden_counts)
         count_parameters = count_parameters + hidden_counts(layer)*previous + &
            merge(1, 0, bias)
         previous = hidden_counts(layer)
      end do
      count_parameters = count_parameters + previous
   end function causality_parameter_count

   pure function kth_neighbor_distances(values, neighbors) result(distance)
      !! Return each row's maximum-norm distance to its kth neighbor.
      real(dp), intent(in) :: values(:, :) !! Observation points by coordinate.
      integer, intent(in) :: neighbors !! Positive neighbor order.
      real(dp) :: distance(size(values, 1))
      real(dp), allocatable :: row_distance(:)
      integer :: observation, other

      allocate(row_distance(size(values, 1)))
      do observation = 1, size(values, 1)
         do other = 1, size(values, 1)
            row_distance(other) = maxval(abs(values(observation, :) - &
               values(other, :)))
         end do
         call sort_in_place(row_distance)
         distance(observation) = row_distance(neighbors + 1)
      end do
   end function kth_neighbor_distances

   pure function neighbor_counts_vector(values, radius, inclusive) result(counts)
      !! Count scalar neighbors within observation-specific radii.
      real(dp), intent(in) :: values(:) !! Scalar observations.
      real(dp), intent(in) :: radius(:) !! Observation-specific radii.
      logical, intent(in) :: inclusive !! Include points on each radius.
      integer :: counts(size(values))
      integer :: observation, other

      counts = 0
      do observation = 1, size(values)
         do other = 1, size(values)
            if (other == observation) cycle
            if (inclusive) then
               if (abs(values(observation) - values(other)) <= &
                  radius(observation)) counts(observation) = counts(observation) + 1
            else
               if (abs(values(observation) - values(other)) < &
                  radius(observation)) counts(observation) = counts(observation) + 1
            end if
         end do
      end do
   end function neighbor_counts_vector

   pure function neighbor_counts_matrix(values, radius, inclusive) result(counts)
      !! Count maximum-norm neighbors within observation-specific radii.
      real(dp), intent(in) :: values(:, :) !! Observation points by coordinate.
      real(dp), intent(in) :: radius(:) !! Observation-specific radii.
      logical, intent(in) :: inclusive !! Include points on each radius.
      integer :: counts(size(values, 1))
      real(dp) :: distance
      integer :: observation, other

      counts = 0
      do observation = 1, size(values, 1)
         do other = 1, size(values, 1)
            if (other == observation) cycle
            distance = maxval(abs(values(observation, :) - values(other, :)))
            if (inclusive) then
               if (distance <= radius(observation)) &
                  counts(observation) = counts(observation) + 1
            else
               if (distance < radius(observation)) &
                  counts(observation) = counts(observation) + 1
            end if
         end do
      end do
   end function neighbor_counts_matrix

   pure subroutine sort_in_place(values)
      !! Sort a short real vector in ascending order by insertion sort.
      real(dp), intent(inout) :: values(:) !! Values replaced by sorted values.
      real(dp) :: item
      integer :: position, previous

      do position = 2, size(values)
         item = values(position)
         previous = position - 1
         do while (previous >= 1)
            if (values(previous) <= item) exit
            values(previous + 1) = values(previous)
            previous = previous - 1
         end do
         values(previous + 1) = item
      end do
   end subroutine sort_in_place

   pure elemental real(dp) function selected_log(value, logarithm) result(log_value)
      !! Evaluate the package's named logarithm base.
      real(dp), intent(in) :: value !! Positive logarithm argument.
      character(len=*), intent(in), optional :: logarithm !! Base name.
      character(len=5) :: selected

      if (value <= 0.0_dp) then
         log_value = 0.0_dp
         return
      end if
      selected = 'log2'
      if (present(logarithm)) selected = adjustl(logarithm)
      select case (trim(selected))
      case ('loge')
         log_value = log(value)
      case ('log10')
         log_value = log10(value)
      case default
         log_value = log(value)/log(2.0_dp)
      end select
   end function selected_log

   pure real(dp) function digamma_value(value) result(digamma)
      !! Evaluate the positive-argument digamma approximation used by NlinTS.
      real(dp), intent(in) :: value !! Positive argument.
      real(dp) :: shifted, inverse_square, correction

      shifted = value
      digamma = 0.0_dp
      do while (shifted <= 5.0_dp)
         digamma = digamma - 1.0_dp/shifted
         shifted = shifted + 1.0_dp
      end do
      inverse_square = 1.0_dp/(shifted*shifted)
      correction = inverse_square*(-1.0_dp/12.0_dp + &
         inverse_square*(1.0_dp/120.0_dp + &
         inverse_square*(-1.0_dp/252.0_dp + &
         inverse_square*(1.0_dp/240.0_dp + &
         inverse_square*(-1.0_dp/132.0_dp + &
         inverse_square*(691.0_dp/32760.0_dp + &
         inverse_square*(-1.0_dp/12.0_dp + &
         inverse_square*3617.0_dp/8160.0_dp)))))))
      digamma = digamma + log(shifted) - 0.5_dp/shifted + correction
   end function digamma_value

   pure real(dp) function f_upper_probability(statistic, numerator_df, &
      denominator_df) result(probability)
      !! Return an upper-tail F probability using the shared beta ratio.
      real(dp), intent(in) :: statistic !! F statistic.
      integer, intent(in) :: numerator_df !! Numerator degrees of freedom.
      integer, intent(in) :: denominator_df !! Denominator degrees of freedom.
      real(dp) :: argument

      if (statistic <= 0.0_dp) then
         probability = 1.0_dp
      else if (numerator_df < 1 .or. denominator_df < 1) then
         probability = quiet_nan()
      else
         argument = real(denominator_df, dp)/(real(denominator_df, dp) + &
            real(numerator_df, dp)*statistic)
         probability = regularized_beta(argument, &
            0.5_dp*real(denominator_df, dp), 0.5_dp*real(numerator_df, dp))
      end if
   end function f_upper_probability

end module nlints_mod
