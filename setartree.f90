! SPDX-License-Identifier: MIT
! SPDX-FileComment: Algorithms translated from the R setartree package.
module setartree_mod
   !! Pooled SETAR trees and subsampled randomized SETAR forests.
   use kind_mod, only: dp
   use linalg_mod, only: symmetric_pseudoinverse
   use stats_mod, only: normal_quantile
   use special_functions_mod, only: regularized_beta
   use random_mod, only: random_uniform
   use utils_mod, only: real_vector_t
   implicit none
   private

   integer, parameter, public :: setartree_stop_both = 1
   integer, parameter, public :: setartree_stop_linearity = 2
   integer, parameter, public :: setartree_stop_error = 3

   type, public :: setartree_category_levels_t
      !! Training levels retained for one integer-coded categorical predictor.
      integer, allocatable :: value(:)
   end type setartree_category_levels_t

   type, public :: setartree_model_t
      !! Flat binary tree with a pooled linear regression in every node.
      real(dp), allocatable :: coefficient(:, :)
      real(dp), allocatable :: threshold(:)
      real(dp), allocatable :: residual_standard_deviation(:)
      real(dp), allocatable :: effective_sample_size(:)
      integer, allocatable :: split_feature(:)
      integer, allocatable :: left_child(:)
      integer, allocatable :: right_child(:)
      integer, allocatable :: sample_count(:)
      logical, allocatable :: terminal(:)
      integer :: feature_count = 0
      integer :: numerical_feature_count = 0
      type(setartree_category_levels_t), allocatable :: categorical_levels(:)
      integer :: node_count = 0
      integer :: depth = 0
      integer :: threshold_count = 0
      integer :: stopping_criterion = setartree_stop_both
      real(dp) :: significance = 0.0_dp
      real(dp) :: significance_divider = 0.0_dp
      real(dp) :: error_threshold = 0.0_dp
      logical :: mean_normalization = .false.
      logical :: window_normalization = .false.
      integer :: info = 0
   end type setartree_model_t

   type, public :: setartree_prediction_t
      !! Point predictions, normal intervals, and terminal-node diagnostics.
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
      real(dp), allocatable :: probability(:)
      real(dp), allocatable :: standard_deviation(:)
      real(dp), allocatable :: effective_leaf_size(:)
      integer, allocatable :: leaf_size(:)
      integer :: info = 0
   end type setartree_prediction_t

   type, public :: setartree_forecast_t
      !! Recursive forecasts by series, horizon, and interval probability.
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: lower(:, :, :)
      real(dp), allocatable :: upper(:, :, :)
      real(dp), allocatable :: probability(:)
      integer :: info = 0
   end type setartree_forecast_t

   type, public :: setarforest_model_t
      !! Collection of subsampled and optionally randomized SETAR trees.
      type(setartree_model_t), allocatable :: tree(:)
      integer :: feature_count = 0
      integer :: info = 0
   end type setarforest_model_t

   public :: setartree_fit, setartree_fit_series
   public :: setartree_fit_categorical, setartree_predict_categorical
   public :: setartree_predict, setartree_forecast
   public :: setarforest_fit, setarforest_fit_series
   public :: setarforest_fit_categorical, setarforest_predict_categorical
   public :: setarforest_predict, setarforest_forecast

   interface setartree_fit_series
      module procedure setartree_fit_series_matrix
      module procedure setartree_fit_series_ragged
   end interface setartree_fit_series

   interface setarforest_fit_series
      module procedure setarforest_fit_series_matrix
      module procedure setarforest_fit_series_ragged
   end interface setarforest_fit_series

   interface setartree_forecast
      module procedure setartree_forecast_matrix
      module procedure setartree_forecast_ragged
   end interface setartree_forecast

   interface setarforest_forecast
      module procedure setarforest_forecast_matrix
      module procedure setarforest_forecast_ragged
   end interface setarforest_forecast

contains

   pure function setartree_fit(x, y, maximum_depth, threshold_count, &
      significance, significance_divider, error_threshold, stopping_criterion, &
      weights) result(out)
      !! Fit a pooled regression tree using lag-specific threshold splits.
      real(dp), intent(in) :: x(:, :) !! Predictor or embedded lag matrix.
      real(dp), intent(in) :: y(:) !! Response values corresponding to rows of x.
      integer, intent(in), optional :: maximum_depth !! Maximum number of splitting levels.
      integer, intent(in), optional :: threshold_count !! Equally spaced thresholds per feature.
      real(dp), intent(in), optional :: significance !! Root-node F-test significance level.
      real(dp), intent(in), optional :: significance_divider !! Per-level significance divisor.
      real(dp), intent(in), optional :: error_threshold !! Minimum relative SSE reduction.
      integer, intent(in), optional :: stopping_criterion !! Both, linearity, or error stopping code.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights.
      type(setartree_model_t) :: out
      integer, allocatable :: root_indices(:)
      real(dp), allocatable :: case_weight(:)
      integer :: depth_limit, grid_count, criterion, capacity, row
      real(dp) :: alpha, alpha_divider, improvement_limit

      depth_limit = 1000
      if (present(maximum_depth)) depth_limit = maximum_depth
      grid_count = 15
      if (present(threshold_count)) grid_count = threshold_count
      alpha = 0.05_dp
      if (present(significance)) alpha = significance
      alpha_divider = 2.0_dp
      if (present(significance_divider)) alpha_divider = significance_divider
      improvement_limit = 0.03_dp
      if (present(error_threshold)) improvement_limit = error_threshold
      criterion = setartree_stop_both
      if (present(stopping_criterion)) criterion = stopping_criterion
      allocate(case_weight(size(y)), source=1.0_dp)
      if (present(weights)) then
         if (size(weights) /= size(y) .or. any(weights < 0.0_dp) .or. &
            sum(weights) <= 0.0_dp) then
            out%info = 1
            return
         end if
         case_weight = weights
      end if
      if (size(x, 1) /= size(y) .or. size(x, 1) < 2 .or. &
         size(x, 2) < 1 .or. depth_limit < 0 .or. grid_count < 3 .or. &
         alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. alpha_divider <= 0.0_dp .or. &
         improvement_limit < 0.0_dp .or. criterion < setartree_stop_both .or. &
         criterion > setartree_stop_error .or. &
         sum(case_weight) <= real(size(x, 2), dp)) then
         out%info = 1
         return
      end if
      capacity = max(1, 2*size(y) - 1)
      allocate(out%coefficient(size(x, 2), capacity), source=0.0_dp)
      allocate(out%threshold(capacity), source=0.0_dp)
      allocate(out%residual_standard_deviation(capacity), source=0.0_dp)
      allocate(out%effective_sample_size(capacity), source=0.0_dp)
      allocate(out%split_feature(capacity), source=0)
      allocate(out%left_child(capacity), source=0)
      allocate(out%right_child(capacity), source=0)
      allocate(out%sample_count(capacity), source=0)
      allocate(out%terminal(capacity), source=.true.)
      out%feature_count = size(x, 2)
      out%numerical_feature_count = size(x, 2)
      out%threshold_count = grid_count
      out%stopping_criterion = criterion
      out%significance = alpha
      out%significance_divider = alpha_divider
      out%error_threshold = improvement_limit
      root_indices = [(row, row=1, size(y))]
      out%node_count = 1
      call grow_node(1, root_indices, 0)

   contains

      recursive pure subroutine grow_node(node, indices, level)
         !! Fit one node and recursively create accepted child splits.
         integer, intent(in) :: node !! Flat-array node index.
         integer, intent(in) :: indices(:) !! Training rows assigned to the node.
         integer, intent(in) :: level !! Zero-based node depth.
         real(dp), allocatable :: beta(:), residual(:)
         integer, allocatable :: left_indices(:), right_indices(:)
         real(dp) :: rss, child_rss, split_value, p_value, improvement
         real(dp) :: level_alpha
         integer :: info, feature, left, right, parameters
         logical :: linearity_ok, improvement_ok, accept

         parameters = size(x, 2)
         call weighted_ols_fit(x(indices, :), y(indices), &
            case_weight(indices), beta, residual, rss, info)
         out%sample_count(node) = size(indices)
         out%effective_sample_size(node) = sum(case_weight(indices))
         out%depth = max(out%depth, level)
         if (info /= 0) then
            if (node == 1) out%info = 2
            return
         end if
         out%coefficient(:, node) = beta
         if (out%effective_sample_size(node) > real(parameters, dp)) &
            out%residual_standard_deviation(node) = sqrt(max(0.0_dp, rss/ &
            (out%effective_sample_size(node) - real(parameters, dp))))
         if (level >= depth_limit .or. &
            out%effective_sample_size(node) <= real(2*parameters + 2, dp)) return
         call best_split(x, y, case_weight, indices, grid_count, feature, &
            split_value, child_rss, left_indices, right_indices)
         if (feature == 0 .or. rss <= tiny(1.0_dp)) return
         improvement = (rss - child_rss)/rss
         level_alpha = alpha/alpha_divider**real(level, dp)
         p_value = split_f_p_value(rss, child_rss, &
            out%effective_sample_size(node), parameters)
         linearity_ok = p_value <= level_alpha
         improvement_ok = improvement >= improvement_limit
         select case (criterion)
         case (setartree_stop_linearity)
            accept = linearity_ok
         case (setartree_stop_error)
            accept = improvement_ok
         case default
            accept = linearity_ok .and. improvement_ok
         end select
         if (.not. accept .or. out%node_count + 2 > capacity) return
         left = out%node_count + 1
         right = out%node_count + 2
         out%node_count = right
         out%terminal(node) = .false.
         out%split_feature(node) = feature
         out%threshold(node) = split_value
         out%left_child(node) = left
         out%right_child(node) = right
         call grow_node(left, left_indices, level + 1)
         call grow_node(right, right_indices, level + 1)
      end subroutine grow_node

   end function setartree_fit

   pure function setartree_fit_categorical(x, categories, y, maximum_depth, &
      threshold_count, significance, significance_divider, error_threshold, &
      stopping_criterion, weights) result(out)
      !! Fit a tree after reference-level encoding of integer categories.
      real(dp), intent(in) :: x(:, :) !! Numerical predictor columns.
      integer, intent(in) :: categories(:, :) !! Integer-coded categorical columns.
      real(dp), intent(in) :: y(:) !! Response values corresponding to input rows.
      integer, intent(in), optional :: maximum_depth !! Maximum number of splitting levels.
      integer, intent(in), optional :: threshold_count !! Equally spaced thresholds per feature.
      real(dp), intent(in), optional :: significance !! Root-node F-test significance level.
      real(dp), intent(in), optional :: significance_divider !! Per-level significance divisor.
      real(dp), intent(in), optional :: error_threshold !! Minimum relative SSE reduction.
      integer, intent(in), optional :: stopping_criterion !! Both, linearity, or error stopping code.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights.
      type(setartree_model_t) :: out
      type(setartree_category_levels_t), allocatable :: levels(:)
      real(dp), allocatable :: design(:, :)
      integer :: info

      call fit_categorical_design(x, categories, design, levels, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      out = setartree_fit(design, y, maximum_depth, threshold_count, &
         significance, significance_divider, error_threshold, &
         stopping_criterion, weights)
      if (out%info == 0) then
         out%numerical_feature_count = size(x, 2)
         out%categorical_levels = levels
      end if
   end function setartree_fit_categorical

   pure function setartree_predict_categorical(model, x, categories, &
      probabilities) result(out)
      !! Predict using the categorical levels retained during tree fitting.
      type(setartree_model_t), intent(in) :: model !! Categorically encoded fitted tree.
      real(dp), intent(in) :: x(:, :) !! Numerical predictor columns.
      integer, intent(in) :: categories(:, :) !! Integer-coded categorical columns.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_prediction_t) :: out
      real(dp), allocatable :: design(:, :)
      integer :: info

      call apply_categorical_design(model, x, categories, design, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      out = setartree_predict(model, design, probabilities)
   end function setartree_predict_categorical

   pure function setartree_fit_series_matrix(series, order, maximum_depth, &
      threshold_count, significance, significance_divider, error_threshold, &
      stopping_criterion, mean_normalization, window_normalization) result(out)
      !! Fit one global SETAR tree to pooled lag windows from multiple series.
      real(dp), intent(in) :: series(:, :) !! Observations by time and series.
      integer, intent(in) :: order !! Positive autoregressive window length.
      integer, intent(in), optional :: maximum_depth !! Maximum splitting depth.
      integer, intent(in), optional :: threshold_count !! Threshold grid size per lag.
      real(dp), intent(in), optional :: significance !! Root F-test significance level.
      real(dp), intent(in), optional :: significance_divider !! Per-level significance divisor.
      real(dp), intent(in), optional :: error_threshold !! Minimum relative SSE reduction.
      integer, intent(in), optional :: stopping_criterion !! Both, linearity, or error stopping code.
      logical, intent(in), optional :: mean_normalization !! Center each source series.
      logical, intent(in), optional :: window_normalization !! Center each embedded lag window.
      type(setartree_model_t) :: out
      real(dp), allocatable :: x(:, :), y(:)
      logical :: normalize_mean, normalize_window
      integer :: info

      normalize_mean = .false.
      if (present(mean_normalization)) normalize_mean = mean_normalization
      normalize_window = .false.
      if (present(window_normalization)) normalize_window = window_normalization
      call embed_series(series, order, normalize_mean, normalize_window, x, y, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      out = setartree_fit(x, y, maximum_depth, threshold_count, significance, &
         significance_divider, error_threshold, stopping_criterion)
      out%mean_normalization = normalize_mean
      out%window_normalization = normalize_window
   end function setartree_fit_series_matrix

   pure function setartree_fit_series_ragged(series, order, maximum_depth, &
      threshold_count, significance, significance_divider, error_threshold, &
      stopping_criterion, mean_normalization, window_normalization) result(out)
      !! Fit one global SETAR tree to unequal-length pooled time series.
      type(real_vector_t), intent(in) :: series(:) !! Ragged collection of source series.
      integer, intent(in) :: order !! Positive autoregressive window length.
      integer, intent(in), optional :: maximum_depth !! Maximum splitting depth.
      integer, intent(in), optional :: threshold_count !! Threshold grid size per lag.
      real(dp), intent(in), optional :: significance !! Root F-test significance level.
      real(dp), intent(in), optional :: significance_divider !! Per-level significance divisor.
      real(dp), intent(in), optional :: error_threshold !! Minimum relative SSE reduction.
      integer, intent(in), optional :: stopping_criterion !! Both, linearity, or error stopping code.
      logical, intent(in), optional :: mean_normalization !! Center each source series.
      logical, intent(in), optional :: window_normalization !! Center each embedded lag window.
      type(setartree_model_t) :: out
      real(dp), allocatable :: x(:, :), y(:)
      logical :: normalize_mean, normalize_window
      integer :: info

      normalize_mean = .false.
      if (present(mean_normalization)) normalize_mean = mean_normalization
      normalize_window = .false.
      if (present(window_normalization)) normalize_window = window_normalization
      call embed_ragged_series(series, order, normalize_mean, &
         normalize_window, x, y, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      out = setartree_fit(x, y, maximum_depth, threshold_count, significance, &
         significance_divider, error_threshold, stopping_criterion)
      out%mean_normalization = normalize_mean
      out%window_normalization = normalize_window
   end function setartree_fit_series_ragged

   pure function setartree_predict(model, x, probabilities) result(out)
      !! Predict rows and construct normal intervals from terminal residual scales.
      type(setartree_model_t), intent(in) :: model !! Fitted SETAR tree.
      real(dp), intent(in) :: x(:, :) !! Predictor rows in training-column order.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_prediction_t) :: out
      real(dp) :: multiplier
      integer :: row, interval, leaf

      if (model%info /= 0 .or. model%node_count < 1 .or. &
         size(x, 2) /= model%feature_count) then
         out%info = 1
         return
      end if
      if (present(probabilities)) then
         if (any(probabilities <= 0.0_dp) .or. any(probabilities >= 1.0_dp)) then
            out%info = 2
            return
         end if
         out%probability = probabilities
      else
         out%probability = [0.8_dp, 0.95_dp]
      end if
      allocate(out%value(size(x, 1)), out%standard_deviation(size(x, 1)))
      allocate(out%effective_leaf_size(size(x, 1)))
      allocate(out%leaf_size(size(x, 1)))
      allocate(out%lower(size(x, 1), size(out%probability)))
      allocate(out%upper(size(x, 1), size(out%probability)))
      do row = 1, size(x, 1)
         leaf = terminal_node(model, x(row, :))
         out%value(row) = dot_product(x(row, :), model%coefficient(:, leaf))
         out%standard_deviation(row) = &
            model%residual_standard_deviation(leaf)
         out%effective_leaf_size(row) = model%effective_sample_size(leaf)
         out%leaf_size(row) = model%sample_count(leaf)
         do interval = 1, size(out%probability)
            multiplier = normal_quantile( &
               0.5_dp*(1.0_dp + out%probability(interval)))
            out%lower(row, interval) = out%value(row) - &
               multiplier*out%standard_deviation(row)
            out%upper(row, interval) = out%value(row) + &
               multiplier*out%standard_deviation(row)
         end do
      end do
   end function setartree_predict

   pure function setartree_forecast_matrix(model, series, horizon, probabilities) &
      result(out)
      !! Recursively forecast one or more series with a fitted global tree.
      type(setartree_model_t), intent(in) :: model !! Tree fitted to lag windows.
      real(dp), intent(in) :: series(:, :) !! Historical observations by time and series.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_forecast_t) :: out

      call forecast_tree_collection([model], series, horizon, probabilities, out)
   end function setartree_forecast_matrix

   pure function setartree_forecast_ragged(model, series, horizon, &
      probabilities) result(out)
      !! Recursively forecast an unequal-length collection with one global tree.
      type(setartree_model_t), intent(in) :: model !! Tree fitted to lag windows.
      type(real_vector_t), intent(in) :: series(:) !! Ragged historical series.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_forecast_t) :: out

      call forecast_ragged_collection([model], series, horizon, probabilities, out)
   end function setartree_forecast_ragged

   function setarforest_fit(x, y, tree_count, bagging_fraction, &
      randomized_parameters, maximum_depth, threshold_count, significance, &
      significance_divider, error_threshold, stopping_criterion, &
      randomize_significance, randomize_significance_divider, &
      randomize_error_threshold, weights) result(out)
      !! Fit a forest of subsampled SETAR trees with optional hyperparameter draws.
      real(dp), intent(in) :: x(:, :) !! Predictor or embedded lag matrix.
      real(dp), intent(in) :: y(:) !! Response values corresponding to rows of x.
      integer, intent(in), optional :: tree_count !! Number of trees; default ten.
      real(dp), intent(in), optional :: bagging_fraction !! Fraction sampled without replacement.
      logical, intent(in), optional :: randomized_parameters !! Randomize stopping hyperparameters.
      integer, intent(in), optional :: maximum_depth !! Maximum splitting depth per tree.
      integer, intent(in), optional :: threshold_count !! Threshold grid size per feature.
      real(dp), intent(in), optional :: significance !! Fixed F-test significance level.
      real(dp), intent(in), optional :: significance_divider !! Fixed per-level divisor.
      real(dp), intent(in), optional :: error_threshold !! Fixed minimum SSE reduction.
      integer, intent(in), optional :: stopping_criterion !! Both, linearity, or error stopping code.
      logical, intent(in), optional :: randomize_significance !! Randomize each tree's F-test level.
      logical, intent(in), optional :: randomize_significance_divider !! Randomize each tree's level divisor.
      logical, intent(in), optional :: randomize_error_threshold !! Randomize each tree's SSE threshold.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights.
      type(setarforest_model_t) :: out
      integer, allocatable :: indices(:)
      real(dp), allocatable :: case_weight(:)
      real(dp) :: fraction, alpha, divider, improvement
      integer :: trees, sample_size, tree
      logical :: randomize, randomize_alpha, randomize_divider
      logical :: randomize_improvement

      trees = 10
      if (present(tree_count)) trees = tree_count
      fraction = 0.8_dp
      if (present(bagging_fraction)) fraction = bagging_fraction
      randomize = .true.
      if (present(randomized_parameters)) randomize = randomized_parameters
      randomize_alpha = randomize
      if (present(randomize_significance)) &
         randomize_alpha = randomize_significance
      randomize_divider = randomize
      if (present(randomize_significance_divider)) &
         randomize_divider = randomize_significance_divider
      randomize_improvement = randomize
      if (present(randomize_error_threshold)) &
         randomize_improvement = randomize_error_threshold
      if (trees < 1 .or. fraction <= 0.0_dp .or. fraction > 1.0_dp .or. &
         size(x, 1) /= size(y) .or. size(x, 2) < 1) then
         out%info = 1
         return
      end if
      allocate(case_weight(size(y)), source=1.0_dp)
      if (present(weights)) then
         if (size(weights) /= size(y) .or. any(weights < 0.0_dp) .or. &
            sum(weights) <= 0.0_dp) then
            out%info = 1
            return
         end if
         case_weight = weights
      end if
      sample_size = nint(fraction*real(size(y), dp))
      if (sample_size <= 2*size(x, 2) + 2) then
         out%info = 2
         return
      end if
      allocate(out%tree(trees))
      do tree = 1, trees
         call sample_without_replacement(size(y), sample_size, indices)
         alpha = 0.05_dp
         if (present(significance)) alpha = significance
         divider = 2.0_dp
         if (present(significance_divider)) divider = significance_divider
         improvement = 0.03_dp
         if (present(error_threshold)) improvement = error_threshold
         if (randomize_alpha) alpha = 0.01_dp + 0.09_dp*random_uniform()
         if (randomize_divider) &
            divider = real(2 + int(9.0_dp*random_uniform()), dp)
         if (randomize_improvement) &
            improvement = 0.001_dp + 0.049_dp*random_uniform()
         out%tree(tree) = setartree_fit(x(indices, :), y(indices), &
            maximum_depth, threshold_count, alpha, divider, improvement, &
            stopping_criterion, case_weight(indices))
         if (out%tree(tree)%info /= 0) then
            out%info = 100 + tree
            return
         end if
      end do
      out%feature_count = size(x, 2)
   end function setarforest_fit

   function setarforest_fit_categorical(x, categories, y, tree_count, &
      bagging_fraction, randomized_parameters, maximum_depth, threshold_count, &
      significance, significance_divider, error_threshold, stopping_criterion, &
      randomize_significance, randomize_significance_divider, &
      randomize_error_threshold, weights) result(out)
      !! Fit a forest after reference-level encoding of integer categories.
      real(dp), intent(in) :: x(:, :) !! Numerical predictor columns.
      integer, intent(in) :: categories(:, :) !! Integer-coded categorical columns.
      real(dp), intent(in) :: y(:) !! Response values corresponding to input rows.
      integer, intent(in), optional :: tree_count !! Number of trees; default ten.
      real(dp), intent(in), optional :: bagging_fraction !! Fraction sampled without replacement.
      logical, intent(in), optional :: randomized_parameters !! Randomize stopping hyperparameters.
      integer, intent(in), optional :: maximum_depth !! Maximum splitting depth per tree.
      integer, intent(in), optional :: threshold_count !! Threshold grid size per feature.
      real(dp), intent(in), optional :: significance !! Fixed F-test significance level.
      real(dp), intent(in), optional :: significance_divider !! Fixed per-level divisor.
      real(dp), intent(in), optional :: error_threshold !! Fixed minimum SSE reduction.
      integer, intent(in), optional :: stopping_criterion !! Both, linearity, or error stopping code.
      logical, intent(in), optional :: randomize_significance !! Randomize each tree's F-test level.
      logical, intent(in), optional :: randomize_significance_divider !! Randomize each tree's level divisor.
      logical, intent(in), optional :: randomize_error_threshold !! Randomize each tree's SSE threshold.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights.
      type(setarforest_model_t) :: out
      type(setartree_category_levels_t), allocatable :: levels(:)
      real(dp), allocatable :: design(:, :)
      integer :: info, tree

      call fit_categorical_design(x, categories, design, levels, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      out = setarforest_fit(design, y, tree_count, bagging_fraction, &
         randomized_parameters, maximum_depth, threshold_count, significance, &
         significance_divider, error_threshold, stopping_criterion, &
         randomize_significance, randomize_significance_divider, &
         randomize_error_threshold, weights)
      if (out%info == 0) then
         do tree = 1, size(out%tree)
            out%tree(tree)%numerical_feature_count = size(x, 2)
            out%tree(tree)%categorical_levels = levels
         end do
      end if
   end function setarforest_fit_categorical

   pure function setarforest_predict_categorical(model, x, categories, &
      probabilities) result(out)
      !! Predict with a forest using its retained categorical levels.
      type(setarforest_model_t), intent(in) :: model !! Categorically encoded fitted forest.
      real(dp), intent(in) :: x(:, :) !! Numerical predictor columns.
      integer, intent(in) :: categories(:, :) !! Integer-coded categorical columns.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_prediction_t) :: out
      real(dp), allocatable :: design(:, :)
      integer :: info

      if (model%info /= 0 .or. .not. allocated(model%tree) .or. &
         size(model%tree) < 1) then
         out%info = 1
         return
      end if
      call apply_categorical_design(model%tree(1), x, categories, design, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      out = setarforest_predict(model, design, probabilities)
   end function setarforest_predict_categorical

   function setarforest_fit_series_matrix(series, order, tree_count, &
      bagging_fraction, randomized_parameters, maximum_depth, threshold_count, &
      significance, significance_divider, error_threshold, stopping_criterion, &
      mean_normalization, window_normalization, randomize_significance, &
      randomize_significance_divider, randomize_error_threshold) result(out)
      !! Fit a SETAR forest to pooled lag windows from multiple series.
      real(dp), intent(in) :: series(:, :) !! Observations by time and series.
      integer, intent(in) :: order !! Positive autoregressive window length.
      integer, intent(in), optional :: tree_count !! Number of trees.
      real(dp), intent(in), optional :: bagging_fraction !! Fraction sampled per tree.
      logical, intent(in), optional :: randomized_parameters !! Randomize stopping hyperparameters.
      integer, intent(in), optional :: maximum_depth !! Maximum splitting depth per tree.
      integer, intent(in), optional :: threshold_count !! Threshold grid size per lag.
      real(dp), intent(in), optional :: significance !! Fixed F-test significance level.
      real(dp), intent(in), optional :: significance_divider !! Fixed per-level divisor.
      real(dp), intent(in), optional :: error_threshold !! Fixed minimum SSE reduction.
      integer, intent(in), optional :: stopping_criterion !! Both, linearity, or error stopping code.
      logical, intent(in), optional :: mean_normalization !! Center each source series.
      logical, intent(in), optional :: window_normalization !! Center each lag window.
      logical, intent(in), optional :: randomize_significance !! Randomize each tree's F-test level.
      logical, intent(in), optional :: randomize_significance_divider !! Randomize each tree's level divisor.
      logical, intent(in), optional :: randomize_error_threshold !! Randomize each tree's SSE threshold.
      type(setarforest_model_t) :: out
      real(dp), allocatable :: x(:, :), y(:)
      logical :: normalize_mean, normalize_window
      integer :: tree, info

      normalize_mean = .false.
      if (present(mean_normalization)) normalize_mean = mean_normalization
      normalize_window = .false.
      if (present(window_normalization)) normalize_window = window_normalization
      call embed_series(series, order, normalize_mean, normalize_window, x, y, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      out = setarforest_fit(x, y, tree_count, bagging_fraction, &
         randomized_parameters, maximum_depth, threshold_count, significance, &
         significance_divider, error_threshold, stopping_criterion, &
         randomize_significance, randomize_significance_divider, &
         randomize_error_threshold)
      if (out%info == 0) then
         do tree = 1, size(out%tree)
            out%tree(tree)%mean_normalization = normalize_mean
            out%tree(tree)%window_normalization = normalize_window
         end do
      end if
   end function setarforest_fit_series_matrix

   function setarforest_fit_series_ragged(series, order, tree_count, &
      bagging_fraction, randomized_parameters, maximum_depth, threshold_count, &
      significance, significance_divider, error_threshold, stopping_criterion, &
      mean_normalization, window_normalization, randomize_significance, &
      randomize_significance_divider, randomize_error_threshold) result(out)
      !! Fit a SETAR forest to unequal-length pooled time series.
      type(real_vector_t), intent(in) :: series(:) !! Ragged collection of source series.
      integer, intent(in) :: order !! Positive autoregressive window length.
      integer, intent(in), optional :: tree_count !! Number of trees.
      real(dp), intent(in), optional :: bagging_fraction !! Fraction sampled per tree.
      logical, intent(in), optional :: randomized_parameters !! Randomize stopping hyperparameters.
      integer, intent(in), optional :: maximum_depth !! Maximum splitting depth per tree.
      integer, intent(in), optional :: threshold_count !! Threshold grid size per lag.
      real(dp), intent(in), optional :: significance !! Fixed F-test significance level.
      real(dp), intent(in), optional :: significance_divider !! Fixed per-level divisor.
      real(dp), intent(in), optional :: error_threshold !! Fixed minimum SSE reduction.
      integer, intent(in), optional :: stopping_criterion !! Both, linearity, or error stopping code.
      logical, intent(in), optional :: mean_normalization !! Center each source series.
      logical, intent(in), optional :: window_normalization !! Center each lag window.
      logical, intent(in), optional :: randomize_significance !! Randomize each tree's F-test level.
      logical, intent(in), optional :: randomize_significance_divider !! Randomize each tree's level divisor.
      logical, intent(in), optional :: randomize_error_threshold !! Randomize each tree's SSE threshold.
      type(setarforest_model_t) :: out
      real(dp), allocatable :: x(:, :), y(:)
      logical :: normalize_mean, normalize_window
      integer :: tree, info

      normalize_mean = .false.
      if (present(mean_normalization)) normalize_mean = mean_normalization
      normalize_window = .false.
      if (present(window_normalization)) normalize_window = window_normalization
      call embed_ragged_series(series, order, normalize_mean, &
         normalize_window, x, y, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      out = setarforest_fit(x, y, tree_count, bagging_fraction, &
         randomized_parameters, maximum_depth, threshold_count, significance, &
         significance_divider, error_threshold, stopping_criterion, &
         randomize_significance, randomize_significance_divider, &
         randomize_error_threshold)
      if (out%info == 0) then
         do tree = 1, size(out%tree)
            out%tree(tree)%mean_normalization = normalize_mean
            out%tree(tree)%window_normalization = normalize_window
         end do
      end if
   end function setarforest_fit_series_ragged

   pure function setarforest_predict(model, x, probabilities) result(out)
      !! Average tree predictions and pool their within-leaf variances.
      type(setarforest_model_t), intent(in) :: model !! Fitted SETAR forest.
      real(dp), intent(in) :: x(:, :) !! Predictor rows in training-column order.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_prediction_t) :: out

      call predict_tree_collection(model%tree, x, probabilities, out)
   end function setarforest_predict

   pure function setarforest_forecast_matrix(model, series, horizon, &
      probabilities) &
      result(out)
      !! Recursively forecast series using averaged forest predictions.
      type(setarforest_model_t), intent(in) :: model !! Fitted SETAR forest.
      real(dp), intent(in) :: series(:, :) !! Historical observations by time and series.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_forecast_t) :: out

      if (model%info /= 0 .or. .not. allocated(model%tree)) then
         out%info = 1
         return
      end if
      call forecast_tree_collection(model%tree, series, horizon, probabilities, out)
   end function setarforest_forecast_matrix

   pure function setarforest_forecast_ragged(model, series, horizon, &
      probabilities) result(out)
      !! Recursively forecast unequal-length series with an averaged forest.
      type(setarforest_model_t), intent(in) :: model !! Fitted SETAR forest.
      type(real_vector_t), intent(in) :: series(:) !! Ragged historical series.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_forecast_t) :: out

      if (model%info /= 0 .or. .not. allocated(model%tree)) then
         out%info = 1
         return
      end if
      call forecast_ragged_collection(model%tree, series, horizon, &
         probabilities, out)
   end function setarforest_forecast_ragged

   pure subroutine fit_categorical_design(x, categories, design, levels, info)
      !! Learn first-seen category levels and construct reference-level indicators.
      real(dp), intent(in) :: x(:, :) !! Numerical predictor columns.
      integer, intent(in) :: categories(:, :) !! Integer-coded categorical columns.
      real(dp), allocatable, intent(out) :: design(:, :) !! Expanded numerical design matrix.
      type(setartree_category_levels_t), allocatable, intent(out) :: levels(:) !! Retained levels by column.
      integer, intent(out) :: info !! Status code; zero indicates success.
      integer, allocatable :: found(:)
      integer :: column, row, level, count, offset, feature_count

      info = 0
      if (size(x, 1) /= size(categories, 1) .or. size(x, 1) < 1 .or. &
         size(x, 2) + size(categories, 2) < 1) then
         info = 1
         return
      end if
      allocate(levels(size(categories, 2)))
      feature_count = size(x, 2)
      do column = 1, size(categories, 2)
         allocate(found(size(categories, 1)))
         count = 0
         do row = 1, size(categories, 1)
            if (count == 0 .or. &
               .not. any(found(:count) == categories(row, column))) then
               count = count + 1
               found(count) = categories(row, column)
            end if
         end do
         levels(column)%value = found(:count)
         feature_count = feature_count + max(0, count - 1)
         deallocate(found)
      end do
      if (feature_count < 1) then
         info = 1
         return
      end if
      allocate(design(size(x, 1), feature_count), source=0.0_dp)
      if (size(x, 2) > 0) design(:, :size(x, 2)) = x
      offset = size(x, 2)
      do column = 1, size(categories, 2)
         do level = 1, size(levels(column)%value) - 1
            design(:, offset + level) = merge(1.0_dp, 0.0_dp, &
               categories(:, column) == levels(column)%value(level))
         end do
         offset = offset + max(0, size(levels(column)%value) - 1)
      end do
   end subroutine fit_categorical_design

   pure subroutine apply_categorical_design(model, x, categories, design, info)
      !! Apply retained category levels and reject values absent during fitting.
      type(setartree_model_t), intent(in) :: model !! Fitted tree containing category metadata.
      real(dp), intent(in) :: x(:, :) !! Numerical predictor columns.
      integer, intent(in) :: categories(:, :) !! Integer-coded categorical columns.
      real(dp), allocatable, intent(out) :: design(:, :) !! Expanded numerical design matrix.
      integer, intent(out) :: info !! Zero, invalid dimensions, or unseen-level status.
      integer :: column, row, level, offset

      info = 0
      if (.not. allocated(model%categorical_levels) .or. &
         size(x, 1) /= size(categories, 1) .or. &
         size(x, 2) /= model%numerical_feature_count .or. &
         size(categories, 2) /= size(model%categorical_levels)) then
         info = 1
         return
      end if
      do column = 1, size(categories, 2)
         do row = 1, size(categories, 1)
            if (.not. any(model%categorical_levels(column)%value == &
               categories(row, column))) then
               info = 2
               return
            end if
         end do
      end do
      allocate(design(size(x, 1), model%feature_count), source=0.0_dp)
      if (size(x, 2) > 0) design(:, :size(x, 2)) = x
      offset = size(x, 2)
      do column = 1, size(categories, 2)
         do level = 1, size(model%categorical_levels(column)%value) - 1
            design(:, offset + level) = merge(1.0_dp, 0.0_dp, &
               categories(:, column) == &
               model%categorical_levels(column)%value(level))
         end do
         offset = offset + &
            max(0, size(model%categorical_levels(column)%value) - 1)
      end do
   end subroutine apply_categorical_design

   pure subroutine weighted_ols_fit(x, y, weights, beta, residual, rss, info)
      !! Fit zero-intercept weighted least squares by dense normal equations.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix.
      real(dp), intent(in) :: y(:) !! Response vector.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights.
      real(dp), allocatable, intent(out) :: beta(:) !! Fitted coefficients.
      real(dp), allocatable, intent(out) :: residual(:) !! Unweighted raw residuals.
      real(dp), intent(out) :: rss !! Weighted residual sum of squares.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: xtx(:, :), inverse(:, :), weighted_x(:, :)

      info = 0
      rss = huge(1.0_dp)
      if (size(x, 1) /= size(y) .or. size(weights) /= size(y) .or. &
         any(weights < 0.0_dp) .or. sum(weights) <= real(size(x, 2), dp)) then
         info = 1
         return
      end if
      weighted_x = x*spread(weights, 2, size(x, 2))
      xtx = matmul(transpose(x), weighted_x)
      allocate(inverse(size(xtx, 1), size(xtx, 2)))
      call symmetric_pseudoinverse(xtx, inverse, info)
      if (info /= 0) return
      beta = matmul(inverse, matmul(transpose(x), weights*y))
      residual = y - matmul(x, beta)
      rss = sum(weights*residual**2)
   end subroutine weighted_ols_fit

   pure subroutine best_split(x, y, weights, indices, threshold_count, feature, &
      threshold, rss, left_indices, right_indices)
      !! Find the best split by cumulatively updating binned sufficient statistics.
      real(dp), intent(in) :: x(:, :) !! Full predictor matrix.
      real(dp), intent(in) :: y(:) !! Full response vector.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights.
      integer, intent(in) :: indices(:) !! Rows assigned to the current node.
      integer, intent(in) :: threshold_count !! Number of grid points per feature.
      integer, intent(out) :: feature !! Selected feature, or zero if none is valid.
      real(dp), intent(out) :: threshold !! Selected split threshold.
      real(dp), intent(out) :: rss !! Combined child residual sum of squares.
      integer, allocatable, intent(out) :: left_indices(:) !! Rows below the threshold.
      integer, allocatable, intent(out) :: right_indices(:) !! Rows at or above the threshold.
      integer, allocatable :: candidate_left(:), candidate_right(:)
      real(dp), allocatable :: bin_xtx(:, :, :), bin_xty(:, :)
      real(dp), allocatable :: bin_yty(:), bin_count(:), grid_value(:)
      real(dp), allocatable :: left_xtx(:, :), right_xtx(:, :)
      real(dp), allocatable :: left_xty(:), right_xty(:)
      real(dp), allocatable :: inverse(:, :), beta(:), row(:)
      real(dp) :: left_yty, right_yty, left_count, right_count
      real(dp) :: minimum, maximum, candidate, left_rss, right_rss, cost
      integer :: column, grid, bin, observation, info, parameters

      feature = 0
      threshold = 0.0_dp
      rss = huge(1.0_dp)
      parameters = size(x, 2)
      allocate(bin_xtx(parameters, parameters, threshold_count + 1))
      allocate(bin_xty(parameters, threshold_count + 1))
      allocate(bin_yty(threshold_count + 1), bin_count(threshold_count + 1))
      allocate(grid_value(threshold_count), left_xtx(parameters, parameters))
      allocate(right_xtx(parameters, parameters), left_xty(parameters))
      allocate(right_xty(parameters), row(parameters))
      do column = 1, parameters
         minimum = minval(x(indices, column))
         maximum = maxval(x(indices, column))
         if (maximum <= minimum) cycle
         do grid = 1, threshold_count
            grid_value(grid) = minimum + real(grid - 1, dp)* &
               (maximum - minimum)/real(threshold_count - 1, dp)
         end do
         bin_xtx = 0.0_dp
         bin_xty = 0.0_dp
         bin_yty = 0.0_dp
         bin_count = 0.0_dp
         do observation = 1, size(indices)
            row = x(indices(observation), :)
            bin = threshold_count + 1
            do grid = 1, threshold_count - 1
               if (x(indices(observation), column) < grid_value(grid + 1)) then
                  bin = grid + 1
                  exit
               end if
            end do
            bin_xtx(:, :, bin) = bin_xtx(:, :, bin) + &
               weights(indices(observation))*spread(row, 2, parameters)* &
               spread(row, 1, parameters)
            bin_xty(:, bin) = bin_xty(:, bin) + &
               weights(indices(observation))*row*y(indices(observation))
            bin_yty(bin) = bin_yty(bin) + weights(indices(observation))* &
               y(indices(observation))**2
            bin_count(bin) = bin_count(bin) + weights(indices(observation))
         end do
         left_xtx = 0.0_dp
         left_xty = 0.0_dp
         left_yty = 0.0_dp
         left_count = 0.0_dp
         right_xtx = sum(bin_xtx, dim=3)
         right_xty = sum(bin_xty, dim=2)
         right_yty = sum(bin_yty)
         right_count = sum(bin_count)
         do grid = 1, threshold_count
            left_xtx = left_xtx + bin_xtx(:, :, grid)
            left_xty = left_xty + bin_xty(:, grid)
            left_yty = left_yty + bin_yty(grid)
            left_count = left_count + bin_count(grid)
            right_xtx = right_xtx - bin_xtx(:, :, grid)
            right_xty = right_xty - bin_xty(:, grid)
            right_yty = right_yty - bin_yty(grid)
            right_count = right_count - bin_count(grid)
            if (left_count <= real(parameters, dp) .or. &
               right_count <= real(parameters, dp)) cycle
            if (.not. allocated(inverse)) &
               allocate(inverse(parameters, parameters))
            call symmetric_pseudoinverse(left_xtx, inverse, info)
            if (info /= 0) cycle
            beta = matmul(inverse, left_xty)
            left_rss = left_yty - 2.0_dp*dot_product(beta, left_xty) + &
               dot_product(beta, matmul(left_xtx, beta))
            call symmetric_pseudoinverse(right_xtx, inverse, info)
            if (info /= 0) cycle
            beta = matmul(inverse, right_xty)
            right_rss = right_yty - 2.0_dp*dot_product(beta, right_xty) + &
               dot_product(beta, matmul(right_xtx, beta))
            cost = left_rss + right_rss
            if (cost <= rss) then
               feature = column
               threshold = grid_value(grid)
               rss = cost
               candidate = grid_value(grid)
               candidate_left = pack(indices, x(indices, column) < candidate)
               candidate_right = pack(indices, x(indices, column) >= candidate)
               left_indices = candidate_left
               right_indices = candidate_right
            end if
         end do
      end do
      if (feature == 0) then
         allocate(left_indices(0), right_indices(0))
      end if
   end subroutine best_split

   pure real(dp) function split_f_p_value(parent_rss, child_rss, observations, &
      predictors) result(p_value)
      !! Return the package's nested-regression F-test upper-tail probability.
      real(dp), intent(in) :: parent_rss !! Parent-node residual sum of squares.
      real(dp), intent(in) :: child_rss !! Combined child residual sum of squares.
      real(dp), intent(in) :: observations !! Effective parent-node sample size.
      integer, intent(in) :: predictors !! Number of regression predictors.
      real(dp) :: statistic, first_df, second_df, beta_argument

      first_df = real(predictors + 1, dp)
      second_df = observations - real(2*predictors + 2, dp)
      if (child_rss <= 0.0_dp .or. second_df <= 0.0_dp .or. &
         parent_rss <= child_rss) then
         p_value = merge(0.0_dp, 1.0_dp, parent_rss > child_rss)
         return
      end if
      statistic = ((parent_rss - child_rss)/first_df)/ &
         (child_rss/second_df)
      beta_argument = second_df/(second_df + first_df*statistic)
      p_value = regularized_beta(beta_argument, 0.5_dp*second_df, &
         0.5_dp*first_df)
   end function split_f_p_value

   pure integer function terminal_node(model, row) result(node)
      !! Traverse a flat SETAR tree to the terminal node for one row.
      type(setartree_model_t), intent(in) :: model !! Fitted SETAR tree.
      real(dp), intent(in) :: row(:) !! Predictor row.

      node = 1
      do while (.not. model%terminal(node))
         if (row(model%split_feature(node)) < model%threshold(node)) then
            node = model%left_child(node)
         else
            node = model%right_child(node)
         end if
      end do
   end function terminal_node

   pure subroutine embed_series(series, order, mean_normalization, &
      window_normalization, x, y, info)
      !! Pool embedded lag windows from equal-length input series.
      real(dp), intent(in) :: series(:, :) !! Observations by time and series.
      integer, intent(in) :: order !! Positive lag-window length.
      logical, intent(in) :: mean_normalization !! Center each series first.
      logical, intent(in) :: window_normalization !! Center each embedded window.
      real(dp), allocatable, intent(out) :: x(:, :) !! Pooled lag predictors.
      real(dp), allocatable, intent(out) :: y(:) !! Pooled responses.
      integer, intent(out) :: info !! Status code; zero indicates success.
      type(real_vector_t), allocatable :: ragged(:)
      integer :: source

      allocate(ragged(size(series, 2)))
      do source = 1, size(series, 2)
         ragged(source)%values = series(:, source)
      end do
      call embed_ragged_series(ragged, order, mean_normalization, &
         window_normalization, x, y, info)
   end subroutine embed_series

   pure subroutine embed_ragged_series(series, order, mean_normalization, &
      window_normalization, x, y, info)
      !! Pool embedded lag windows from unequal-length input series.
      type(real_vector_t), intent(in) :: series(:) !! Ragged source-series collection.
      integer, intent(in) :: order !! Positive lag-window length.
      logical, intent(in) :: mean_normalization !! Center each series first.
      logical, intent(in) :: window_normalization !! Center each embedded window.
      real(dp), allocatable, intent(out) :: x(:, :) !! Pooled lag predictors.
      real(dp), allocatable, intent(out) :: y(:) !! Pooled responses.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: values(:)
      real(dp) :: center, series_mean
      integer :: total, source, time, lag, row, observations

      info = 0
      total = 0
      if (order < 1 .or. size(series) < 1) then
         info = 1
         allocate(x(0, 0), y(0))
         return
      end if
      do source = 1, size(series)
         if (.not. allocated(series(source)%values) .or. &
            size(series(source)%values) <= order) then
            info = 2
            allocate(x(0, 0), y(0))
            return
         end if
         total = total + size(series(source)%values) - order
      end do
      allocate(x(total, order), y(total))
      row = 0
      do source = 1, size(series)
         values = series(source)%values
         observations = size(values)
         if (mean_normalization) then
            series_mean = sum(values)/real(observations, dp)
            values = values - series_mean
         end if
         do time = order + 1, observations
            row = row + 1
            y(row) = values(time)
            do lag = 1, order
               x(row, lag) = values(time - lag)
            end do
            if (window_normalization) then
               center = sum(x(row, :))/real(order, dp)
               x(row, :) = x(row, :) - center
               y(row) = y(row) - center
            end if
         end do
      end do
   end subroutine embed_ragged_series

   subroutine sample_without_replacement(population, sample_size, sample)
      !! Draw an unordered simple random sample using a partial Fisher-Yates shuffle.
      integer, intent(in) :: population !! Number of available row indices.
      integer, intent(in) :: sample_size !! Number of indices drawn.
      integer, allocatable, intent(out) :: sample(:) !! Selected row indices.
      integer, allocatable :: index(:)
      integer :: position, selected, held

      index = [(position, position=1, population)]
      do position = 1, sample_size
         selected = position + int(random_uniform()*real(population - position + 1, dp))
         selected = min(selected, population)
         held = index(position)
         index(position) = index(selected)
         index(selected) = held
      end do
      sample = index(:sample_size)
   end subroutine sample_without_replacement

   pure subroutine predict_tree_collection(tree, x, probabilities, out)
      !! Average tree predictions and pool within-leaf residual variances.
      type(setartree_model_t), intent(in) :: tree(:) !! Fitted tree collection.
      real(dp), intent(in) :: x(:, :) !! Predictor rows.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_prediction_t), intent(out) :: out !! Averaged prediction result.
      type(setartree_prediction_t) :: one
      real(dp), allocatable :: variance_numerator(:), degrees(:)
      real(dp) :: multiplier
      integer :: member, interval

      if (size(tree) < 1) then
         out%info = 1
         return
      end if
      one = setartree_predict(tree(1), x, probabilities)
      if (one%info /= 0) then
         out%info = one%info
         return
      end if
      out = one
      variance_numerator = max(0.0_dp, one%effective_leaf_size - 1.0_dp)* &
         one%standard_deviation**2
      degrees = max(0.0_dp, one%effective_leaf_size - 1.0_dp)
      do member = 2, size(tree)
         one = setartree_predict(tree(member), x, probabilities)
         if (one%info /= 0) then
            out%info = 100 + member
            return
         end if
         out%value = out%value + one%value
         variance_numerator = variance_numerator + &
            max(0.0_dp, one%effective_leaf_size - 1.0_dp)* &
            one%standard_deviation**2
         degrees = degrees + max(0.0_dp, one%effective_leaf_size - 1.0_dp)
      end do
      out%value = out%value/real(size(tree), dp)
      out%standard_deviation = sqrt(variance_numerator/max(1.0_dp, degrees))
      out%effective_leaf_size = degrees + 1.0_dp
      out%leaf_size = nint(degrees) + 1
      do interval = 1, size(out%probability)
         multiplier = normal_quantile( &
            0.5_dp*(1.0_dp + out%probability(interval)))
         out%lower(:, interval) = out%value - &
            multiplier*out%standard_deviation
         out%upper(:, interval) = out%value + &
            multiplier*out%standard_deviation
      end do
   end subroutine predict_tree_collection

   pure subroutine forecast_tree_collection(tree, series, horizon, &
      probabilities, out)
      !! Recursively forecast using one tree or an averaged tree collection.
      type(setartree_model_t), intent(in) :: tree(:) !! Fitted tree collection.
      real(dp), intent(in) :: series(:, :) !! Historical observations by time and series.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_forecast_t), intent(out) :: out !! Recursive forecasts and intervals.
      type(real_vector_t), allocatable :: ragged(:)
      integer :: source

      allocate(ragged(size(series, 2)))
      do source = 1, size(series, 2)
         ragged(source)%values = series(:, source)
      end do
      call forecast_ragged_collection(tree, ragged, horizon, probabilities, out)
   end subroutine forecast_tree_collection

   pure subroutine forecast_ragged_collection(tree, series, horizon, &
      probabilities, out)
      !! Recursively forecast unequal-length series with a tree collection.
      type(setartree_model_t), intent(in) :: tree(:) !! Fitted tree collection.
      type(real_vector_t), intent(in) :: series(:) !! Ragged historical series.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: probabilities(:) !! Central interval probabilities.
      type(setartree_forecast_t), intent(out) :: out !! Recursive forecasts and intervals.
      type(setartree_prediction_t) :: prediction
      real(dp), allocatable :: state(:, :), predictor(:, :), series_mean(:)
      real(dp) :: center
      integer :: order, count, time, source, lag, interval, observations

      if (size(tree) < 1 .or. horizon < 1) then
         out%info = 1
         return
      end if
      order = tree(1)%feature_count
      count = size(series)
      if (count < 1) then
         out%info = 2
         return
      end if
      do source = 1, count
         if (.not. allocated(series(source)%values) .or. &
            size(series(source)%values) < order) then
            out%info = 2
            return
         end if
      end do
      allocate(state(order + horizon, count), predictor(count, order))
      allocate(series_mean(count), source=0.0_dp)
      do source = 1, count
         observations = size(series(source)%values)
         if (tree(1)%mean_normalization) series_mean(source) = &
            sum(series(source)%values)/real(observations, dp)
         state(:order, source) = &
            series(source)%values(observations - order + 1:) - &
            series_mean(source)
      end do
      do time = 1, horizon
         do source = 1, count
            do lag = 1, order
               predictor(source, lag) = &
                  state(order + time - lag, source)
            end do
            if (tree(1)%window_normalization) then
               center = sum(predictor(source, :))/real(order, dp)
               predictor(source, :) = predictor(source, :) - center
            end if
         end do
         call predict_tree_collection(tree, predictor, probabilities, prediction)
         if (prediction%info /= 0) then
            out%info = 100 + time
            return
         end if
         if (time == 1) then
            out%probability = prediction%probability
            allocate(out%mean(count, horizon))
            allocate(out%lower(count, horizon, size(out%probability)))
            allocate(out%upper(count, horizon, size(out%probability)))
         end if
         do source = 1, count
            center = 0.0_dp
            if (tree(1)%window_normalization) center = &
               sum(state(time:order + time - 1, source))/real(order, dp)
            state(order + time, source) = &
               prediction%value(source) + center
            out%mean(source, time) = state(order + time, source) + &
               series_mean(source)
            do interval = 1, size(out%probability)
               out%lower(source, time, interval) = &
                  prediction%lower(source, interval) + center + &
                  series_mean(source)
               out%upper(source, time, interval) = &
                  prediction%upper(source, interval) + center + &
                  series_mean(source)
            end do
         end do
      end do
   end subroutine forecast_ragged_collection

end module setartree_mod
