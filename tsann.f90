! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R TSANN package.
module tsann_mod
   !! Automatic lag-order and hidden-width selection for neural autoregressions.
   use kind_mod, only: dp
   use nnfor_mod, only: nnfor_mlp_model_t, nnfor_forecast_t, nnfor_mlp, &
      nnfor_mlp_forecast, nnfor_combine_median
   use random_mod, only: set_random_seed
   use time_series_stats_mod, only: acf_values
   use utils_mod, only: quiet_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: tsann_select_validation = 1
   integer, parameter, public :: tsann_select_test = 2

   type, public :: tsann_trace_t
      !! Accuracy measurements for every lag-order and hidden-width candidate.
      integer, allocatable :: lag_order(:)
      integer, allocatable :: hidden_size(:)
      real(dp), allocatable :: train_rmse(:)
      real(dp), allocatable :: validation_rmse(:)
      real(dp), allocatable :: test_rmse(:)
   end type tsann_trace_t

   type, public :: tsann_fit_t
      !! Selected TSANN neural autoregression and chronological holdout results.
      type(nnfor_mlp_model_t) :: model
      type(tsann_trace_t) :: trace
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: predicted(:)
      real(dp) :: train_rmse = huge(1.0_dp)
      real(dp) :: validation_rmse = huge(1.0_dp)
      real(dp) :: test_rmse = huge(1.0_dp)
      integer :: selected_lag = 0
      integer :: selected_hidden = 0
      integer :: maximum_lag = 0
      integer :: training_observations = 0
      integer :: validation_observations = 0
      integer :: test_observations = 0
      integer :: selection_metric = tsann_select_validation
      integer :: info = 0
   end type tsann_fit_t

   interface display
      module procedure display_tsann_fit
   end interface display

   public :: tsann_maximum_lag, tsann_auto_fit, display

contains

   pure function tsann_maximum_lag(series, threshold, lag_max) result(order)
      !! Derive the TSANN maximum lag by counting ACF values above a threshold.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      real(dp), intent(in), optional :: threshold !! Positive ACF cutoff; defaults to 0.05.
      integer, intent(in), optional :: lag_max !! Largest ACF lag considered.
      integer :: order
      real(dp), allocatable :: autocorrelation(:)
      real(dp) :: cutoff
      integer :: maximum

      order = 0
      cutoff = 0.05_dp
      if (present(threshold)) cutoff = threshold
      if (size(series) < 2 .or. cutoff < -1.0_dp .or. cutoff > 1.0_dp .or. &
         .not. all(ieee_is_finite(series))) return
      maximum = min(size(series) - 1, &
         max(1, floor(10.0_dp*log10(real(size(series), dp)))))
      if (present(lag_max)) maximum = min(size(series) - 1, lag_max)
      if (maximum < 1) return
      autocorrelation = acf_values(series, maximum)
      order = count(autocorrelation > cutoff)
      order = max(1, min(order, maximum))
   end function tsann_maximum_lag

   function tsann_auto_fit(series, minimum_hidden, maximum_hidden, &
      split_ratio, maximum_lag, acf_threshold, validation_fraction, &
      repetitions, combination, max_iterations, tolerance, decay, seed, &
      selection_metric) result(out)
      !! Search TSANN lag and hidden sizes and refit on training plus validation data.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      integer, intent(in) :: minimum_hidden !! Smallest hidden-layer width.
      integer, intent(in) :: maximum_hidden !! Largest hidden-layer width.
      real(dp), intent(in) :: split_ratio !! Fraction reserved for training plus validation.
      integer, intent(in), optional :: maximum_lag !! Largest autoregressive order searched.
      real(dp), intent(in), optional :: acf_threshold !! ACF cutoff used to derive the lag limit.
      real(dp), intent(in), optional :: validation_fraction !! Fraction of pre-test data used for validation.
      integer, intent(in), optional :: repetitions !! Neural-network ensemble size.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations per network.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! Neural-network weight decay.
      integer, intent(in), optional :: seed !! Shared random-number seed.
      integer, intent(in), optional :: selection_metric !! Validation or compatibility test criterion.
      type(tsann_fit_t) :: out
      type(nnfor_mlp_model_t) :: train_model, pretest_model
      type(nnfor_forecast_t) :: validation_forecast, test_forecast
      integer, allocatable :: lags(:)
      real(dp) :: validation_share, score, best_score
      integer :: observations, pretest_count, train_count, validation_count
      integer :: test_count, lag_order, hidden, candidate, candidate_count
      integer :: lag_index

      observations = size(series)
      validation_share = 0.2_dp
      if (present(validation_fraction)) validation_share = validation_fraction
      out%selection_metric = tsann_select_validation
      if (present(selection_metric)) out%selection_metric = selection_metric
      if (observations < 10 .or. minimum_hidden < 1 .or. &
         maximum_hidden < minimum_hidden .or. split_ratio <= 0.0_dp .or. &
         split_ratio >= 1.0_dp .or. validation_share <= 0.0_dp .or. &
         validation_share >= 1.0_dp .or. &
         (out%selection_metric /= tsann_select_validation .and. &
         out%selection_metric /= tsann_select_test) .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      pretest_count = nint(real(observations, dp)*split_ratio)
      pretest_count = max(2, min(observations - 1, pretest_count))
      train_count = nint(real(pretest_count, dp)*(1.0_dp - validation_share))
      train_count = max(1, min(pretest_count - 1, train_count))
      validation_count = pretest_count - train_count
      test_count = observations - pretest_count
      out%maximum_lag = tsann_maximum_lag(series(:pretest_count), &
         acf_threshold)
      if (present(maximum_lag)) out%maximum_lag = maximum_lag
      out%maximum_lag = min(out%maximum_lag, train_count - 3)
      if (out%maximum_lag < 1 .or. validation_count < 1 .or. test_count < 1) then
         out%info = 2
         return
      end if
      out%training_observations = train_count
      out%validation_observations = validation_count
      out%test_observations = test_count
      out%series = series
      candidate_count = out%maximum_lag*(maximum_hidden - minimum_hidden + 1)
      allocate(out%trace%lag_order(candidate_count))
      allocate(out%trace%hidden_size(candidate_count))
      allocate(out%trace%train_rmse(candidate_count), source=huge(1.0_dp))
      allocate(out%trace%validation_rmse(candidate_count), source=huge(1.0_dp))
      allocate(out%trace%test_rmse(candidate_count), source=huge(1.0_dp))
      if (present(seed)) call set_random_seed(seed)
      best_score = huge(1.0_dp)
      candidate = 0
      do hidden = minimum_hidden, maximum_hidden
         do lag_order = 1, out%maximum_lag
            candidate = candidate + 1
            out%trace%lag_order(candidate) = lag_order
            out%trace%hidden_size(candidate) = hidden
            lags = [(lag_index, lag_index = 1, lag_order)]
            train_model = nnfor_mlp(series(:train_count), lags, hidden, &
               repetitions, combination, max_iterations, tolerance, decay)
            if (train_model%info /= 0) cycle
            out%trace%train_rmse(candidate) = sqrt(train_model%mse)
            validation_forecast = nnfor_mlp_forecast(train_model, &
               validation_count)
            if (validation_forecast%info /= 0) cycle
            out%trace%validation_rmse(candidate) = root_mean_square_error( &
               series(train_count + 1:pretest_count), &
               validation_forecast%mean)
            pretest_model = nnfor_mlp(series(:pretest_count), lags, hidden, &
               repetitions, combination, max_iterations, tolerance, decay)
            if (pretest_model%info /= 0) cycle
            test_forecast = nnfor_mlp_forecast(pretest_model, test_count)
            if (test_forecast%info /= 0) cycle
            out%trace%test_rmse(candidate) = root_mean_square_error( &
               series(pretest_count + 1:), test_forecast%mean)
            score = out%trace%validation_rmse(candidate)
            if (out%selection_metric == tsann_select_test) &
               score = out%trace%test_rmse(candidate)
            if (score < best_score) then
               best_score = score
               out%selected_lag = lag_order
               out%selected_hidden = hidden
               out%train_rmse = out%trace%train_rmse(candidate)
               out%validation_rmse = out%trace%validation_rmse(candidate)
               out%test_rmse = out%trace%test_rmse(candidate)
               out%model = pretest_model
               out%predicted = test_forecast%mean
            end if
         end do
      end do
      if (out%selected_lag == 0) then
         out%info = 3
         return
      end if
      allocate(out%fitted(observations), source=quiet_nan())
      out%fitted(out%selected_lag + 1:pretest_count) = out%model%fitted
   end function tsann_auto_fit

   pure function root_mean_square_error(actual, predicted) result(value)
      !! Return root mean square error for equally sized finite vectors.
      real(dp), intent(in) :: actual(:) !! Observed values.
      real(dp), intent(in) :: predicted(:) !! Predicted values.
      real(dp) :: value

      value = huge(1.0_dp)
      if (size(actual) < 1 .or. size(actual) /= size(predicted)) return
      if (.not. all(ieee_is_finite(actual)) .or. &
         .not. all(ieee_is_finite(predicted))) return
      value = sqrt(sum((actual - predicted)**2)/real(size(actual), dp))
   end function root_mean_square_error

   subroutine display_tsann_fit(fit, unit, print_obs)
      !! Display a selected TSANN model and optionally its holdout predictions.
      type(tsann_fit_t), intent(in) :: fit !! Automatic TSANN result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Print fitted and test observations.
      integer :: destination, observation, first_test
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Automatic TSANN fit'
      write(destination, '(a, i0)') '  status: ', fit%info
      write(destination, '(a, i0)') '  selected lag: ', fit%selected_lag
      write(destination, '(a, i0)') '  selected hidden size: ', &
         fit%selected_hidden
      write(destination, '(a, es14.6)') '  training RMSE: ', fit%train_rmse
      write(destination, '(a, es14.6)') '  validation RMSE: ', &
         fit%validation_rmse
      write(destination, '(a, es14.6)') '  test RMSE: ', fit%test_rmse
      if (.not. show_observations .or. fit%info /= 0) return
      write(destination, '(a)') '  index       observed         fitted       predicted'
      first_test = fit%training_observations + fit%validation_observations + 1
      do observation = 1, size(fit%series)
         if (observation < first_test) then
            write(destination, '(i7, 3(1x, es14.6))') observation, &
               fit%series(observation), fit%fitted(observation), quiet_nan()
         else
            write(destination, '(i7, 3(1x, es14.6))') observation, &
               fit%series(observation), quiet_nan(), &
               fit%predicted(observation - first_test + 1)
         end if
      end do
   end subroutine display_tsann_fit

end module tsann_mod
