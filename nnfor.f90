! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R nnfor package.
module nnfor_mod
   !! Neural-network ensembles for univariate time-series forecasting.
   use kind_mod, only: dp
   use linalg_mod, only: symmetric_pseudoinverse
   use gmdh_mod, only: gmdh_ridge_fit_t, gmdh_ridge_fit
   use neural_network_mod, only: neural_network_t, neural_network_fit, &
      neural_network_predict, neural_network_parameter_count
   use random_mod, only: random_uniform
   use special_functions_mod, only: regularized_beta, regularized_gamma_q
   use stats_mod, only: covariance, median, quantile, sorted, &
      standard_deviation
   use utils_mod, only: real_vector_t, integer_vector_t, quiet_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: nnfor_estimator_least_squares = 1
   integer, parameter, public :: nnfor_estimator_ridge = 2
   integer, parameter, public :: nnfor_estimator_stepwise = 3
   integer, parameter, public :: nnfor_estimator_lasso = 4
   integer, parameter, public :: nnfor_combine_mean = 1
   integer, parameter, public :: nnfor_combine_median = 2
   integer, parameter, public :: nnfor_combine_mode = 3
   integer, parameter, public :: nnfor_seasonal_none = 0
   integer, parameter, public :: nnfor_seasonal_binary = 1
   integer, parameter, public :: nnfor_seasonal_trigonometric = 2
   integer, parameter, public :: nnfor_hidden_terminal = 1
   integer, parameter, public :: nnfor_hidden_holdout = 2
   integer, parameter, public :: nnfor_hidden_cross_validation = 3
   ! Canova-Hansen code and tables follow uroot 2.1-3, as used by nnfor.
   integer, parameter, public :: nnfor_ch_trigonometric = 1
   integer, parameter, public :: nnfor_ch_dummy = 2
   real(dp), parameter :: ch_response_probabilities(19) = [ &
      0.900_dp, 0.905_dp, 0.910_dp, 0.915_dp, 0.920_dp, 0.925_dp, &
      0.930_dp, 0.935_dp, 0.940_dp, 0.945_dp, 0.950_dp, 0.955_dp, &
      0.960_dp, 0.965_dp, 0.970_dp, 0.975_dp, 0.980_dp, 0.985_dp, &
      0.990_dp]
   real(dp), parameter :: ch_joint_coefficients(7, 19) = reshape([ &
      0.225883431594717_dp, 40.785912780134_dp, -1223.10813512854_dp, &
      -9.49538996613502_dp, 241.869094038328_dp, 0.213206946773136_dp, &
      -0.00094320183958978_dp, &
      0.236596242525799_dp, 40.4675149902581_dp, -1223.09283412158_dp, &
      -9.57540801813603_dp, 244.138726781229_dp, 0.213904558560353_dp, &
      -0.000952779716644947_dp, &
      0.247973922126445_dp, 40.1038364301467_dp, -1221.7016669941_dp, &
      -9.65500231274856_dp, 246.372666496303_dp, 0.214601888453645_dp, &
      -0.000962372474829586_dp, &
      0.259986350173644_dp, 39.692464011608_dp, -1219.0836304927_dp, &
      -9.73554823797049_dp, 248.616545049371_dp, 0.215325184223994_dp, &
      -0.000972418184555302_dp, &
      0.272562750091551_dp, 39.283360159758_dp, -1217.8970791641_dp, &
      -9.82040814361997_dp, 251.008639900973_dp, 0.2160922393077_dp, &
      -0.000983154762963007_dp, &
      0.286104337773487_dp, 38.7762639228097_dp, -1213.26184585229_dp, &
      -9.90527747280634_dp, 253.344733975778_dp, 0.216876792389571_dp, &
      -0.000994271131518386_dp, &
      0.300285491168037_dp, 38.2819144125574_dp, -1210.53239094133_dp, &
      -9.99620184597908_dp, 255.886272082389_dp, 0.217720596696401_dp, &
      -0.00100624314944497_dp, &
      0.315697528102994_dp, 37.6946472903104_dp, -1205.75251283096_dp, &
      -10.088015334406_dp, 258.437785308615_dp, 0.218586128322271_dp, &
      -0.00101867719750197_dp, &
      0.332437711770686_dp, 36.9980063320313_dp, -1197.63136968853_dp, &
      -10.1832603228807_dp, 261.0407672875_dp, 0.219486820722537_dp, &
      -0.00103166550391064_dp, &
      0.350494229843421_dp, 36.2396721518201_dp, -1188.74765469053_dp, &
      -10.2800983545934_dp, 263.684495368906_dp, 0.220426166130944_dp, &
      -0.00104541041522041_dp, &
      0.3702961729347_dp, 35.3814993726306_dp, -1178.12965548839_dp, &
      -10.3813239153467_dp, 266.445699639894_dp, 0.221405709557784_dp, &
      -0.00105986751629488_dp, &
      0.391996227047414_dp, 34.3923625934811_dp, -1165.06804966711_dp, &
      -10.4888022893712_dp, 269.360871246441_dp, 0.222491166831091_dp, &
      -0.00107613366559805_dp, &
      0.416281776255618_dp, 33.2128546357423_dp, -1147.58111741195_dp, &
      -10.6037109905073_dp, 272.45255829158_dp, 0.223669916824912_dp, &
      -0.00109405527991665_dp, &
      0.443885585281462_dp, 31.7964506996133_dp, -1124.1019083257_dp, &
      -10.7258281790835_dp, 275.689783771037_dp, 0.224926773511144_dp, &
      -0.00111333084511144_dp, &
      0.475850715844896_dp, 30.0495012015624_dp, -1091.85131756088_dp, &
      -10.8599473411407_dp, 279.173479343288_dp, 0.226308484907728_dp, &
      -0.00113485866790207_dp, &
      0.513589227193938_dp, 27.895033873042_dp, -1051.01580209343_dp, &
      -11.0066456564546_dp, 282.966264187653_dp, 0.227839822438738_dp, &
      -0.00115893912845004_dp, &
      0.559033256531683_dp, 25.2325036467302_dp, -999.397907921899_dp, &
      -11.1775800075243_dp, 287.365962114776_dp, 0.229661544959428_dp, &
      -0.00118829287377202_dp, &
      0.61822259161555_dp, 21.5203391933771_dp, -921.032879412685_dp, &
      -11.3824641723439_dp, 292.522550935924_dp, 0.231756196698891_dp, &
      -0.00122247845484498_dp, &
      0.700480583101547_dp, 15.9908659933099_dp, -796.986602362939_dp, &
      -11.6554227482472_dp, 299.235670490034_dp, 0.234590735499554_dp, &
      -0.00126984486026105_dp], [7, 19])
   real(dp), parameter :: ch_raw_probabilities(6) = [ &
      0.20_dp, 0.10_dp, 0.075_dp, 0.05_dp, 0.025_dp, 0.01_dp]
   real(dp), parameter :: ch_raw_critical_values(6, 12) = reshape([ &
      0.243_dp, 0.353_dp, 0.398_dp, 0.470_dp, 0.593_dp, 0.748_dp, &
      0.469_dp, 0.610_dp, 0.670_dp, 0.749_dp, 0.898_dp, 1.070_dp, &
      0.679_dp, 0.846_dp, 0.913_dp, 1.010_dp, 1.160_dp, 1.350_dp, &
      0.883_dp, 1.070_dp, 1.140_dp, 1.240_dp, 1.390_dp, 1.600_dp, &
      1.080_dp, 1.280_dp, 1.360_dp, 1.470_dp, 1.630_dp, 1.880_dp, &
      1.280_dp, 1.490_dp, 1.580_dp, 1.680_dp, 1.890_dp, 2.120_dp, &
      1.460_dp, 1.690_dp, 1.780_dp, 1.900_dp, 2.100_dp, 2.350_dp, &
      1.660_dp, 1.890_dp, 1.990_dp, 2.110_dp, 2.330_dp, 2.590_dp, &
      1.850_dp, 2.100_dp, 2.190_dp, 2.320_dp, 2.550_dp, 2.820_dp, &
      2.030_dp, 2.290_dp, 2.400_dp, 2.540_dp, 2.760_dp, 3.050_dp, &
      2.220_dp, 2.490_dp, 2.600_dp, 2.750_dp, 2.990_dp, 3.270_dp, &
      2.410_dp, 2.690_dp, 2.810_dp, 2.960_dp, 3.180_dp, 3.510_dp], [6, 12])

   type, public :: nnfor_scale_t
      !! Parameters of an affine min-max transformation.
      real(dp) :: target_minimum = -1.0_dp
      real(dp) :: target_maximum = 1.0_dp
      real(dp) :: original_minimum = 0.0_dp
      real(dp) :: original_maximum = 1.0_dp
      integer :: info = 0
   end type nnfor_scale_t

   type, public :: nnfor_scaled_t
      !! Scaled values and the transformation needed to reverse them.
      real(dp), allocatable :: values(:)
      type(nnfor_scale_t) :: scale
      integer :: info = 0
   end type nnfor_scaled_t

   type, public :: nnfor_elm_layer_t
      !! Bias-augmented incoming weights for one ELM hidden layer.
      real(dp), allocatable :: input_weights(:, :)
   end type nnfor_elm_layer_t

   type, public :: nnfor_elm_member_t
      !! Weights for one fast extreme-learning-machine ensemble member.
      real(dp), allocatable :: input_weights(:, :)
      type(nnfor_elm_layer_t), allocatable :: layers(:)
      integer, allocatable :: hidden_counts(:)
      real(dp), allocatable :: output_weights(:)
      real(dp), allocatable :: direct_weights(:)
      real(dp) :: output_bias = 0.0_dp
      real(dp) :: lambda = 0.0_dp
      integer :: hidden_count = 0
      integer :: info = 0
   end type nnfor_elm_member_t

   type, public :: nnfor_output_fit_t
      !! Selected ELM output-layer regression coefficients.
      real(dp), allocatable :: coefficients(:)
      logical, allocatable :: active(:)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: validation_mse = huge(1.0_dp)
      real(dp) :: rss = huge(1.0_dp)
      integer :: info = 0
   end type nnfor_output_fit_t

   type, public :: nnfor_elm_fast_t
      !! Fast ELM ensemble fitted to a regression matrix.
      type(nnfor_elm_member_t), allocatable :: members(:)
      type(nnfor_scale_t), allocatable :: predictor_scales(:)
      type(nnfor_scale_t) :: response_scale
      real(dp), allocatable :: fitted_all(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: mse = huge(1.0_dp)
      integer :: predictor_count = 0
      integer :: repetitions = 0
      integer :: estimator = nnfor_estimator_least_squares
      integer :: combination = nnfor_combine_median
      integer :: orthogonalized_layer_count = 0
      logical :: direct = .false.
      logical :: scaled = .true.
      logical :: orthogonal = .false.
      integer :: info = 0
   end type nnfor_elm_fast_t

   type, public :: nnfor_preprocessing_t
      !! Aligned differencing, lag, exogenous, and seasonal network inputs.
      type(real_vector_t), allocatable :: difference_levels(:)
      type(integer_vector_t), allocatable :: exogenous_lags(:)
      integer, allocatable :: lags(:)
      integer, allocatable :: difference_lags(:)
      integer, allocatable :: periods(:)
      integer, allocatable :: seasonal_types(:)
      real(dp), allocatable :: response(:)
      real(dp), allocatable :: original_response(:)
      real(dp), allocatable :: predictors(:, :)
      integer :: maximum_lag = 0
      integer :: difference_offset = 0
      integer :: period = 1
      integer :: seasonal_type = nnfor_seasonal_none
      integer :: exogenous_count = 0
      integer :: start_index = 0
      integer :: info = 0
   end type nnfor_preprocessing_t

   type, public :: nnfor_lag_selection_t
      !! Backward-AIC selection result for candidate response lags.
      integer, allocatable :: selected_lags(:)
      logical, allocatable :: active(:)
      real(dp) :: aic = huge(1.0_dp)
      integer :: info = 0
   end type nnfor_lag_selection_t

   type, public :: nnfor_elm_model_t
      !! Lagged univariate ELM forecasting model.
      type(nnfor_elm_fast_t) :: network
      type(nnfor_preprocessing_t) :: preprocessing
      integer, allocatable :: lags(:)
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: mse = huge(1.0_dp)
      integer :: info = 0
      logical :: extended_preprocessing = .false.
   end type nnfor_elm_model_t

   type, public :: nnfor_mlp_model_t
      !! Ensemble of multilayer perceptrons fitted to lagged observations.
      type(neural_network_t), allocatable :: members(:)
      type(nnfor_preprocessing_t) :: preprocessing
      type(nnfor_scale_t) :: response_scale
      integer, allocatable :: lags(:)
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: fitted_all(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      integer, allocatable :: hidden_counts(:)
      real(dp) :: mse = huge(1.0_dp)
      integer :: hidden_count = 0
      integer :: repetitions = 0
      integer :: combination = nnfor_combine_median
      integer :: info = 0
      logical :: extended_preprocessing = .false.
   end type nnfor_mlp_model_t

   type, public :: nnfor_hidden_selection_t
      !! Validation errors and selected one-hidden-layer MLP size.
      real(dp), allocatable :: mse(:)
      integer :: selected = 0
      integer :: maximum = 0
      integer :: repetitions = 0
      integer :: info = 0
   end type nnfor_hidden_selection_t

   type, public :: nnfor_elm_hidden_selection_t
      !! Significant hidden-unit counts and automatic ELM hidden size.
      integer, allocatable :: significant_count(:)
      integer :: selected = 0
      integer :: candidate_hidden = 0
      integer :: repetitions = 0
      real(dp) :: alpha = 0.05_dp
      integer :: info = 0
   end type nnfor_elm_hidden_selection_t

   type, public :: nnfor_seasonality_t
      !! Friedman seasonality diagnostic and centered moving average.
      real(dp), allocatable :: moving_average(:)
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: period = 1
      integer :: cycles = 0
      integer :: info = 0
      logical :: seasonal = .false.
      logical :: multiplicative = .false.
   end type nnfor_seasonality_t

   type, public :: nnfor_canova_hansen_t
      !! Joint Canova-Hansen seasonal-stability test result.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: critical_value = huge(1.0_dp)
      real(dp) :: alpha = 0.05_dp
      integer :: period = 1
      integer :: seasonal_type = nnfor_ch_trigonometric
      integer :: newey_west_order = 0
      integer :: degrees_of_freedom = 0
      integer :: info = 0
      logical :: difference_required = .false.
   end type nnfor_canova_hansen_t

   type, public :: nnfor_mseason_t
      !! Correlation test for multiplicative seasonal magnitude.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: alpha = 0.05_dp
      integer :: period = 1
      integer :: seasonal_count = 1
      integer :: info = 0
      logical :: multiplicative = .false.
   end type nnfor_mseason_t

   type, public :: nnfor_trend_t
      !! AIC comparison of level-only and additive-trend exponential smoothing.
      real(dp) :: level_aic = huge(1.0_dp)
      real(dp) :: trend_aic = huge(1.0_dp)
      integer :: info = 0
      logical :: trending = .false.
   end type nnfor_trend_t

   type, public :: nnfor_difference_selection_t
      !! Automatically selected ordinary and seasonal difference lags.
      integer, allocatable :: difference_lags(:)
      type(nnfor_seasonality_t) :: seasonality
      type(nnfor_mseason_t) :: multiplicative_seasonality
      type(nnfor_canova_hansen_t) :: canova_hansen
      type(nnfor_trend_t) :: trend
      integer :: info = 0
   end type nnfor_difference_selection_t

   type, public :: nnfor_forecast_t
      !! Ensemble point forecasts and individual network paths.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: all_mean(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      integer :: info = 0
   end type nnfor_forecast_t

   type, public :: nnfor_elm_auto_t
      !! Automatically specified ELM fit and its selection diagnostics.
      type(nnfor_elm_model_t) :: model
      type(nnfor_difference_selection_t) :: difference_selection
      type(nnfor_lag_selection_t) :: lag_selection
      type(nnfor_elm_hidden_selection_t) :: hidden_selection
      integer, allocatable :: candidate_lags(:)
      integer, allocatable :: difference_lags(:)
      integer, allocatable :: periods(:)
      integer, allocatable :: seasonal_types(:)
      integer :: seasonal_type = nnfor_seasonal_none
      integer :: info = 0
      logical :: automatic_differences = .true.
   end type nnfor_elm_auto_t

   type, public :: nnfor_mlp_auto_t
      !! Automatically specified MLP fit and its selection diagnostics.
      type(nnfor_mlp_model_t) :: model
      type(nnfor_difference_selection_t) :: difference_selection
      type(nnfor_lag_selection_t) :: lag_selection
      type(nnfor_hidden_selection_t) :: hidden_selection
      integer, allocatable :: candidate_lags(:)
      integer, allocatable :: difference_lags(:)
      integer, allocatable :: periods(:)
      integer, allocatable :: seasonal_types(:)
      integer :: seasonal_type = nnfor_seasonal_none
      integer :: info = 0
      logical :: automatic_differences = .true.
   end type nnfor_mlp_auto_t

   type :: least_squares_fit_t
      !! Internal least-squares fit for a selected design subset.
      real(dp), allocatable :: coefficients(:)
      real(dp) :: rss = huge(1.0_dp)
      integer :: info = 0
   end type least_squares_fit_t

   type :: significance_result_t
      !! Internal coefficient significance calculation.
      real(dp), allocatable :: p_value(:)
      integer :: info = 0
   end type significance_result_t

   type :: automatic_specification_t
      !! Internal automatic preprocessing and lag-selection result.
      type(nnfor_preprocessing_t) :: preprocessing
      type(nnfor_difference_selection_t) :: difference_selection
      type(nnfor_lag_selection_t) :: lag_selection
      integer, allocatable :: candidate_lags(:)
      integer, allocatable :: selected_lags(:)
      integer, allocatable :: difference_lags(:)
      integer, allocatable :: periods(:)
      integer, allocatable :: seasonal_types(:)
      integer :: period = 1
      integer :: seasonal_type = nnfor_seasonal_none
      integer :: info = 0
      logical :: automatic_differences = .true.
   end type automatic_specification_t

   interface display
      module procedure display_nnfor_elm_fast
      module procedure display_nnfor_elm_model
      module procedure display_nnfor_mlp_model
      module procedure display_nnfor_forecast
      module procedure display_nnfor_mseason
      module procedure display_nnfor_canova_hansen
      module procedure display_nnfor_elm_auto
      module procedure display_nnfor_mlp_auto
   end interface display

   public :: nnfor_linscale, nnfor_apply_scale, nnfor_fast_sigmoid
   public :: integer_vector_t
   public :: nnfor_difference, nnfor_preprocess
   public :: nnfor_select_lags
   public :: nnfor_elm_fast_from_weights, nnfor_elm_fast
   public :: nnfor_elm_fast_layers_from_weights, nnfor_elm_fast_layers
   public :: nnfor_elm_fast_predict
   public :: nnfor_elm_from_weights, nnfor_elm, nnfor_elm_forecast
   public :: nnfor_elm_layers_from_weights, nnfor_elm_layers
   public :: nnfor_elm_preprocessed_from_weights
   public :: nnfor_elm_preprocessed, nnfor_elm_preprocessed_forecast
   public :: nnfor_elm_preprocessed_layers_from_weights
   public :: nnfor_elm_preprocessed_layers
   public :: nnfor_mlp_from_initial, nnfor_mlp, nnfor_mlp_forecast
   public :: nnfor_mlp_layers_from_initial, nnfor_mlp_layers
   public :: nnfor_mlp_preprocessed_from_initial
   public :: nnfor_mlp_preprocessed, nnfor_mlp_preprocessed_forecast
   public :: nnfor_mlp_preprocessed_layers_from_initial
   public :: nnfor_mlp_preprocessed_layers
   public :: nnfor_select_hidden_count, nnfor_select_hidden_count_folds
   public :: nnfor_select_hidden_count_random
   public :: nnfor_select_elm_hidden_from_weights, nnfor_select_elm_hidden
   public :: nnfor_elm_auto, nnfor_elm_auto_forecast
   public :: nnfor_mlp_auto, nnfor_mlp_auto_forecast
   public :: nnfor_elm_thief, nnfor_mlp_thief
   public :: nnfor_elm_refit, nnfor_elm_retrain
   public :: nnfor_mlp_refit, nnfor_mlp_retrain
   public :: nnfor_season_check, nnfor_mseason_test, nnfor_trend_check
   public :: nnfor_canova_hansen, nnfor_select_differences
   public :: nnfor_ridge_output_fit, nnfor_lasso_output_fit
   public :: nnfor_stepwise_output_fit, nnfor_kde_mode, nnfor_combine
   public :: display, display_nnfor_elm_fast, display_nnfor_elm_model
   public :: display_nnfor_mlp_model, display_nnfor_forecast
   public :: display_nnfor_mseason, display_nnfor_canova_hansen

contains

   pure function nnfor_linscale(values, target_minimum, target_maximum) &
      result(out)
      !! Fit and apply nnfor's affine min-max scaling transformation.
      real(dp), intent(in) :: values(:) !! Values used to estimate and apply scaling.
      real(dp), intent(in), optional :: target_minimum !! Target lower bound; defaults to -1.
      real(dp), intent(in), optional :: target_maximum !! Target upper bound; defaults to 1.
      type(nnfor_scaled_t) :: out
      real(dp) :: width

      if (present(target_minimum)) out%scale%target_minimum = target_minimum
      if (present(target_maximum)) out%scale%target_maximum = target_maximum
      if (size(values) < 1 .or. &
         out%scale%target_maximum <= out%scale%target_minimum .or. &
         .not. all(ieee_is_finite(values))) then
         out%info = 1
         out%scale%info = 1
         return
      end if
      out%scale%original_minimum = minval(values)
      out%scale%original_maximum = maxval(values)
      width = out%scale%original_maximum - out%scale%original_minimum
      allocate(out%values(size(values)))
      if (width <= epsilon(1.0_dp)*max(1.0_dp, &
         abs(out%scale%original_maximum))) then
         out%values = 0.5_dp*(out%scale%target_minimum + &
            out%scale%target_maximum)
      else
         out%values = (out%scale%target_maximum - &
            out%scale%target_minimum)*(values - &
            out%scale%original_minimum)/width + out%scale%target_minimum
      end if
   end function nnfor_linscale

   pure function nnfor_apply_scale(values, scale, reverse) result(out)
      !! Apply a stored min-max scale in the forward or reverse direction.
      real(dp), intent(in) :: values(:) !! Values to transform.
      type(nnfor_scale_t), intent(in) :: scale !! Stored affine transformation.
      logical, intent(in), optional :: reverse !! Reverse to original units when true.
      real(dp), allocatable :: out(:)
      real(dp) :: source_width, destination_width
      logical :: inverse

      inverse = .false.
      if (present(reverse)) inverse = reverse
      allocate(out(size(values)))
      if (inverse) then
         source_width = scale%target_maximum - scale%target_minimum
         destination_width = scale%original_maximum - scale%original_minimum
         if (abs(source_width) <= tiny(1.0_dp)) then
            out = scale%original_minimum
         else
            out = destination_width*(values - scale%target_minimum)/ &
               source_width + scale%original_minimum
         end if
      else
         source_width = scale%original_maximum - scale%original_minimum
         destination_width = scale%target_maximum - scale%target_minimum
         if (abs(source_width) <= epsilon(1.0_dp)*max(1.0_dp, &
            abs(scale%original_maximum))) then
            out = 0.5_dp*(scale%target_minimum + scale%target_maximum)
         else
            out = destination_width*(values - scale%original_minimum)/ &
               source_width + scale%target_minimum
         end if
      end if
   end function nnfor_apply_scale

   pure elemental real(dp) function nnfor_fast_sigmoid(value) result(out)
      !! Evaluate nnfor's symmetric fast sigmoid x/(1+abs(x)).
      real(dp), intent(in) :: value !! Linear hidden-unit activation.

      out = value/(1.0_dp + abs(value))
   end function nnfor_fast_sigmoid

   pure function nnfor_difference(series, difference_lags) result(values)
      !! Apply sequential ordinary or seasonal differences at supplied lags.
      real(dp), intent(in) :: series(:) !! Undifferenced time-series observations.
      integer, intent(in) :: difference_lags(:) !! Positive sequential difference lags.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: next(:)
      integer :: difference, lag

      values = series
      do difference = 1, size(difference_lags)
         lag = difference_lags(difference)
         if (lag < 1 .or. size(values) <= lag) then
            allocate(next(0))
            call move_alloc(next, values)
            return
         end if
         next = values(lag + 1:) - values(:size(values) - lag)
         call move_alloc(next, values)
      end do
   end function nnfor_difference

   pure function nnfor_preprocess(series, lags, difference_lags, period, &
      seasonal_type, exogenous, exogenous_lags, periods, seasonal_types) &
      result(out)
      !! Build nnfor response and lagged deterministic/exogenous network inputs.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive lags of the differenced response.
      integer, intent(in), optional :: difference_lags(:) !! Sequential positive difference lags.
      integer, intent(in), optional :: period !! Seasonal period for deterministic inputs.
      integer, intent(in), optional :: seasonal_type !! None, binary, or trigonometric seasonality.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time exogenous regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Nonnegative lags by regressor.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Input type for each seasonal period.
      type(nnfor_preprocessing_t) :: out
      real(dp), allocatable :: next(:)
      integer :: differences, rows, predictor_count, seasonal_columns
      integer :: difference, lag, row, column, variable, time, phase
      integer :: seasonality

      out%period = 1
      if (present(period)) out%period = period
      out%seasonal_type = nnfor_seasonal_none
      if (present(seasonal_type)) out%seasonal_type = seasonal_type
      if (present(periods)) then
         out%periods = periods
      else
         out%periods = [out%period]
      end if
      if (present(seasonal_types)) then
         out%seasonal_types = seasonal_types
      else
         allocate(out%seasonal_types(size(out%periods)), &
            source=out%seasonal_type)
      end if
      if (size(out%periods) > 0) out%period = maxval(out%periods)
      if (size(out%periods) == 1) then
         out%seasonal_type = out%seasonal_types(1)
      else
         out%seasonal_type = nnfor_seasonal_none
      end if
      if (present(difference_lags)) then
         out%difference_lags = difference_lags
      else
         allocate(out%difference_lags(0))
      end if
      out%lags = lags
      differences = size(out%difference_lags)
      if (size(series) < 3 .or. any(lags < 1) .or. &
         any(out%difference_lags < 1) .or. &
         sum(out%difference_lags) >= size(series) - 2 .or. &
         size(out%periods) < 1 .or. &
         size(out%seasonal_types) /= size(out%periods) .or. &
         any(out%periods < 1) .or. &
         any(out%seasonal_types < nnfor_seasonal_none) .or. &
         any(out%seasonal_types > nnfor_seasonal_trigonometric) .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      do seasonality = 1, size(out%periods)
         if (out%seasonal_types(seasonality) /= nnfor_seasonal_none .and. &
            out%periods(seasonality) < 2) then
            out%info = 1
            return
         end if
         if (seasonality > 1) then
            if (any(out%periods(:seasonality - 1) == &
               out%periods(seasonality))) then
               out%info = 1
               return
            end if
         end if
      end do
      out%exogenous_count = 0
      if (present(exogenous)) out%exogenous_count = size(exogenous, 2)
      if (out%exogenous_count > 0) then
         if (.not. present(exogenous_lags)) then
            out%info = 1
            return
         end if
         if (size(exogenous_lags) /= out%exogenous_count .or. &
            size(exogenous, 1) < size(series) .or. &
            .not. all(ieee_is_finite(exogenous))) then
            out%info = 1
            return
         end if
         out%exogenous_lags = exogenous_lags
         do variable = 1, out%exogenous_count
            if (.not. allocated(out%exogenous_lags(variable)%values)) then
               out%info = 1
               return
            end if
            if (any(out%exogenous_lags(variable)%values < 0)) then
               out%info = 1
               return
            end if
         end do
      else
         allocate(out%exogenous_lags(0))
      end if
      allocate(out%difference_levels(differences + 1))
      out%difference_levels(1)%values = series
      do difference = 1, differences
         lag = out%difference_lags(difference)
         next = out%difference_levels(difference)%values(lag + 1:) - &
            out%difference_levels(difference)%values( &
            :size(out%difference_levels(difference)%values) - lag)
         out%difference_levels(difference + 1)%values = next
      end do
      out%difference_offset = sum(out%difference_lags)
      out%maximum_lag = 0
      if (size(lags) > 0) out%maximum_lag = maxval(lags)
      do variable = 1, out%exogenous_count
         if (size(out%exogenous_lags(variable)%values) > 0) &
            out%maximum_lag = max(out%maximum_lag, &
               maxval(out%exogenous_lags(variable)%values))
      end do
      rows = size(out%difference_levels(differences + 1)%values) - &
         out%maximum_lag
      if (rows < 3) then
         out%info = 2
         return
      end if
      seasonal_columns = 0
      do seasonality = 1, size(out%periods)
         if (out%seasonal_types(seasonality) == nnfor_seasonal_binary) &
            seasonal_columns = seasonal_columns + out%periods(seasonality) - 1
         if (out%seasonal_types(seasonality) == &
            nnfor_seasonal_trigonometric) &
            seasonal_columns = seasonal_columns + 2
      end do
      predictor_count = size(lags) + seasonal_columns
      do variable = 1, out%exogenous_count
         predictor_count = predictor_count + &
            size(out%exogenous_lags(variable)%values)
      end do
      if (predictor_count < 1) then
         out%info = 2
         return
      end if
      allocate(out%response(rows), out%original_response(rows))
      allocate(out%predictors(rows, predictor_count), source=0.0_dp)
      out%response = out%difference_levels(differences + 1)%values( &
         out%maximum_lag + 1:)
      out%start_index = out%difference_offset + out%maximum_lag + 1
      out%original_response = series(out%start_index:)
      do row = 1, rows
         column = 0
         do lag = 1, size(lags)
            column = column + 1
            out%predictors(row, column) = &
               out%difference_levels(differences + 1)%values( &
               out%maximum_lag + row - lags(lag))
         end do
         time = out%difference_offset + out%maximum_lag + row
         do variable = 1, out%exogenous_count
            do lag = 1, size(out%exogenous_lags(variable)%values)
               column = column + 1
               out%predictors(row, column) = exogenous(time - &
                  out%exogenous_lags(variable)%values(lag), variable)
            end do
         end do
         do seasonality = 1, size(out%periods)
            phase = modulo(time - 1, out%periods(seasonality)) + 1
            if (out%seasonal_types(seasonality) == nnfor_seasonal_binary) then
               do lag = 1, out%periods(seasonality) - 1
                  column = column + 1
                  out%predictors(row, column) = &
                     merge(1.0_dp, 0.0_dp, phase == lag)
               end do
            else if (out%seasonal_types(seasonality) == &
               nnfor_seasonal_trigonometric) then
               column = column + 1
               out%predictors(row, column) = sin(2.0_dp*acos(-1.0_dp)* &
                  real(time, dp)/real(out%periods(seasonality), dp))
               column = column + 1
               out%predictors(row, column) = cos(2.0_dp*acos(-1.0_dp)* &
                  real(time, dp)/real(out%periods(seasonality), dp))
            end if
         end do
      end do
   end function nnfor_preprocess

   pure function nnfor_select_lags(preprocessing, keep) result(out)
      !! Select response lags by backward AIC while retaining other inputs.
      type(nnfor_preprocessing_t), intent(in) :: preprocessing !! Candidate network inputs.
      logical, intent(in), optional :: keep(:) !! Response-lag mask that cannot be removed.
      type(nnfor_lag_selection_t) :: out
      type(least_squares_fit_t) :: current, candidate_fit, best_fit
      real(dp), allocatable :: design(:, :)
      logical, allocatable :: active(:), candidate_active(:), best_active(:)
      logical, allocatable :: forced(:)
      real(dp) :: current_aic, candidate_aic, best_aic
      integer :: observations, candidate

      if (preprocessing%info /= 0) then
         out%info = 1
         return
      end if
      if (.not. allocated(preprocessing%predictors) .or. &
         .not. allocated(preprocessing%response) .or. &
         .not. allocated(preprocessing%lags)) then
         out%info = 1
         return
      end if
      if (present(keep)) then
         if (size(keep) /= size(preprocessing%lags)) then
            out%info = 1
            return
         end if
         forced = keep
      else
         allocate(forced(size(preprocessing%lags)), source=.false.)
      end if
      observations = size(preprocessing%response)
      allocate(design(observations, size(preprocessing%predictors, 2) + 1))
      design(:, 1) = 1.0_dp
      design(:, 2:) = preprocessing%predictors
      allocate(active(size(design, 2)), source=.true.)
      current = fit_selected_design(design, preprocessing%response, active)
      if (current%info /= 0) then
         out%info = 2
         return
      end if
      current_aic = regression_aic(current%rss, observations, count(active))
      do
         best_aic = current_aic
         best_active = active
         best_fit = current
         do candidate = 1, size(preprocessing%lags)
            if (.not. active(candidate + 1) .or. forced(candidate)) cycle
            candidate_active = active
            candidate_active(candidate + 1) = .false.
            if (count(candidate_active(2:size(preprocessing%lags) + 1)) == 0 .and. &
               size(preprocessing%predictors, 2) == size(preprocessing%lags)) &
               cycle
            candidate_fit = fit_selected_design(design, preprocessing%response, &
               candidate_active)
            if (candidate_fit%info /= 0) cycle
            candidate_aic = regression_aic(candidate_fit%rss, observations, &
               count(candidate_active))
            if (candidate_aic < best_aic - 1.0e-10_dp) then
               best_aic = candidate_aic
               best_active = candidate_active
               best_fit = candidate_fit
            end if
         end do
         if (all(best_active .eqv. active)) exit
         active = best_active
         current = best_fit
         current_aic = best_aic
      end do
      out%active = active(2:size(preprocessing%lags) + 1)
      out%selected_lags = pack(preprocessing%lags, out%active)
      out%aic = current_aic
   end function nnfor_select_lags

   pure function nnfor_elm_fast_from_weights(response, predictors, &
      input_weights, estimator, combination, direct, scale_data, lambdas, &
      validation_weight) result(model)
      !! Fit a reproducible fast ELM ensemble from supplied hidden-layer weights.
      real(dp), intent(in) :: response(:) !! Regression target.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-predictor matrix.
      real(dp), intent(in) :: input_weights(:, :, :) !! Bias-and-input weights by repetition.
      integer, intent(in), optional :: estimator !! Least-squares or ridge output estimator.
      integer, intent(in), optional :: combination !! Mean or median ensemble operator.
      logical, intent(in), optional :: direct !! Include direct linear input-output connections.
      logical, intent(in), optional :: scale_data !! Scale inputs and target to [-0.8,0.8].
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge penalties.
      real(dp), intent(in), optional :: validation_weight !! Ridge-training fraction.
      type(nnfor_elm_fast_t) :: model
      type(nnfor_scaled_t) :: scaled
      type(nnfor_output_fit_t) :: output_fit
      real(dp), allocatable :: x(:, :), y(:), augmented(:, :), hidden(:, :)
      real(dp), allocatable :: design(:, :), output_design(:, :)
      real(dp), allocatable :: selected_lambdas(:), prediction(:)
      real(dp) :: weight
      integer :: observations, predictors_count, hidden_count, repetitions
      integer :: member, predictor, columns, offset

      observations = size(predictors, 1)
      predictors_count = size(predictors, 2)
      hidden_count = size(input_weights, 2)
      repetitions = size(input_weights, 3)
      model%estimator = nnfor_estimator_least_squares
      if (present(estimator)) model%estimator = estimator
      model%combination = nnfor_combine_median
      if (present(combination)) model%combination = combination
      if (present(direct)) model%direct = direct
      model%scaled = .true.
      if (present(scale_data)) model%scaled = scale_data
      weight = 0.7_dp
      if (present(validation_weight)) weight = validation_weight
      if (present(lambdas)) then
         selected_lambdas = lambdas
      else
         selected_lambdas = [0.0_dp, 0.01_dp, 0.1_dp, 1.0_dp, 10.0_dp]
      end if
      if (size(response) /= observations .or. observations < 3 .or. &
         predictors_count < 1 .or. hidden_count < 1 .or. repetitions < 1 .or. &
         size(input_weights, 1) /= predictors_count + 1 .or. &
         (model%estimator < nnfor_estimator_least_squares .or. &
         model%estimator > nnfor_estimator_lasso) .or. &
         (model%combination /= nnfor_combine_mean .and. &
         model%combination /= nnfor_combine_median .and. &
         model%combination /= nnfor_combine_mode) .or. &
         .not. all(ieee_is_finite(response)) .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         .not. all(ieee_is_finite(input_weights))) then
         model%info = 1
         return
      end if
      model%predictor_count = predictors_count
      model%repetitions = repetitions
      allocate(x(observations, predictors_count), y(observations))
      allocate(model%predictor_scales(predictors_count))
      if (model%scaled) then
         scaled = nnfor_linscale(response, -0.8_dp, 0.8_dp)
         y = scaled%values
         model%response_scale = scaled%scale
         do predictor = 1, predictors_count
            scaled = nnfor_linscale(predictors(:, predictor), -0.8_dp, 0.8_dp)
            x(:, predictor) = scaled%values
            model%predictor_scales(predictor) = scaled%scale
         end do
      else
         y = response
         x = predictors
         model%response_scale%original_minimum = 0.0_dp
         model%response_scale%original_maximum = 1.0_dp
         do predictor = 1, predictors_count
            model%predictor_scales(predictor)%original_minimum = 0.0_dp
            model%predictor_scales(predictor)%original_maximum = 1.0_dp
         end do
      end if
      allocate(augmented(observations, predictors_count + 1))
      augmented(:, 1) = 1.0_dp
      augmented(:, 2:) = x
      columns = hidden_count + merge(predictors_count, 0, model%direct)
      allocate(model%members(repetitions), model%fitted_all(observations, repetitions))
      do member = 1, repetitions
         hidden = nnfor_fast_sigmoid(matmul(augmented, input_weights(:, :, member)))
         allocate(design(observations, columns))
         design(:, :hidden_count) = hidden
         if (model%direct) design(:, hidden_count + 1:) = x
         allocate(output_design(observations, columns + 1))
         output_design(:, 1) = 1.0_dp
         output_design(:, 2:) = design
         select case (model%estimator)
         case (nnfor_estimator_least_squares)
            output_fit = nnfor_ridge_output_fit(output_design, y, [0.0_dp], &
               weight)
         case (nnfor_estimator_ridge)
            output_fit = nnfor_ridge_output_fit(output_design, y, &
               selected_lambdas, weight)
         case (nnfor_estimator_stepwise)
            output_fit = nnfor_stepwise_output_fit(output_design, y)
         case (nnfor_estimator_lasso)
            output_fit = nnfor_lasso_output_fit(output_design, y, &
               selected_lambdas, weight)
         end select
         if (output_fit%info /= 0) then
            model%info = 2
            return
         end if
         model%members(member)%input_weights = input_weights(:, :, member)
         model%members(member)%hidden_count = hidden_count
         model%members(member)%output_bias = output_fit%coefficients(1)
         model%members(member)%output_weights = &
            output_fit%coefficients(2:hidden_count + 1)
         model%members(member)%lambda = output_fit%lambda
         offset = hidden_count + 1
         if (model%direct) model%members(member)%direct_weights = &
            output_fit%coefficients(offset + 1:offset + predictors_count)
         prediction = model%members(member)%output_bias + &
            matmul(hidden, model%members(member)%output_weights)
         if (model%direct) prediction = prediction + &
            matmul(x, model%members(member)%direct_weights)
         if (model%scaled) prediction = nnfor_apply_scale(prediction, &
            model%response_scale, .true.)
         model%fitted_all(:, member) = prediction
         deallocate(design, output_design)
      end do
      model%fitted = nnfor_combine(model%fitted_all, model%combination)
      model%residuals = response - model%fitted
      model%mse = sum(model%residuals**2)/real(observations, dp)
   end function nnfor_elm_fast_from_weights

   pure function nnfor_elm_fast_layers_from_weights(response, predictors, &
      layer_weights, estimator, combination, direct, scale_data, lambdas, &
      validation_weight) result(model)
      !! Fit a multilayer ELM ensemble from supplied hidden-layer weights.
      real(dp), intent(in) :: response(:) !! Regression target.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-predictor matrix.
      type(nnfor_elm_layer_t), intent(in) :: layer_weights(:, :) !! Layers by repetition.
      integer, intent(in), optional :: estimator !! Output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      logical, intent(in), optional :: scale_data !! Scale inputs and target.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate regularization penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      type(nnfor_elm_fast_t) :: model
      type(nnfor_scaled_t) :: scaled
      type(nnfor_output_fit_t) :: output_fit
      real(dp), allocatable :: x(:, :), y(:), hidden(:, :), design(:, :)
      real(dp), allocatable :: output_design(:, :), selected_lambdas(:)
      real(dp), allocatable :: prediction(:)
      real(dp) :: weight
      integer :: observations, predictors_count, layers_count, repetitions
      integer :: member, layer, predictor, previous, hidden_count, columns, offset

      observations = size(predictors, 1)
      predictors_count = size(predictors, 2)
      layers_count = size(layer_weights, 1)
      repetitions = size(layer_weights, 2)
      model%estimator = nnfor_estimator_least_squares
      if (present(estimator)) model%estimator = estimator
      model%combination = nnfor_combine_median
      if (present(combination)) model%combination = combination
      if (present(direct)) model%direct = direct
      model%scaled = .true.
      if (present(scale_data)) model%scaled = scale_data
      weight = 0.7_dp
      if (present(validation_weight)) weight = validation_weight
      if (present(lambdas)) then
         selected_lambdas = lambdas
      else
         selected_lambdas = [0.0_dp, 0.01_dp, 0.1_dp, 1.0_dp, 10.0_dp]
      end if
      if (size(response) /= observations .or. observations < 3 .or. &
         predictors_count < 1 .or. layers_count < 1 .or. repetitions < 1 .or. &
         (model%estimator < nnfor_estimator_least_squares .or. &
         model%estimator > nnfor_estimator_lasso) .or. &
         (model%combination /= nnfor_combine_mean .and. &
         model%combination /= nnfor_combine_median .and. &
         model%combination /= nnfor_combine_mode) .or. &
         .not. all(ieee_is_finite(response)) .or. &
         .not. all(ieee_is_finite(predictors))) then
         model%info = 1
         return
      end if
      previous = predictors_count
      do layer = 1, layers_count
         do member = 1, repetitions
            if (.not. allocated(layer_weights(layer, member)%input_weights)) then
               model%info = 1
               return
            end if
            if (size(layer_weights(layer, member)%input_weights, 1) /= &
               previous + 1 .or. &
               size(layer_weights(layer, member)%input_weights, 2) < 1 .or. &
               .not. all(ieee_is_finite( &
               layer_weights(layer, member)%input_weights))) then
               model%info = 1
               return
            end if
            if (member > 1 .and. size(layer_weights(layer, member)%input_weights, &
               2) /= size(layer_weights(layer, 1)%input_weights, 2)) then
               model%info = 1
               return
            end if
         end do
         previous = size(layer_weights(layer, 1)%input_weights, 2)
      end do
      hidden_count = previous
      model%predictor_count = predictors_count
      model%repetitions = repetitions
      allocate(x(observations, predictors_count), y(observations))
      allocate(model%predictor_scales(predictors_count))
      if (model%scaled) then
         scaled = nnfor_linscale(response, -0.8_dp, 0.8_dp)
         y = scaled%values
         model%response_scale = scaled%scale
         do predictor = 1, predictors_count
            scaled = nnfor_linscale(predictors(:, predictor), -0.8_dp, 0.8_dp)
            x(:, predictor) = scaled%values
            model%predictor_scales(predictor) = scaled%scale
         end do
      else
         y = response
         x = predictors
         model%response_scale%original_minimum = 0.0_dp
         model%response_scale%original_maximum = 1.0_dp
         do predictor = 1, predictors_count
            model%predictor_scales(predictor)%original_minimum = 0.0_dp
            model%predictor_scales(predictor)%original_maximum = 1.0_dp
         end do
      end if
      columns = hidden_count + merge(predictors_count, 0, model%direct)
      allocate(model%members(repetitions), model%fitted_all(observations, repetitions))
      do member = 1, repetitions
         hidden = evaluate_elm_layers(layer_weights(:, member), x)
         allocate(design(observations, columns))
         design(:, :hidden_count) = hidden
         if (model%direct) design(:, hidden_count + 1:) = x
         allocate(output_design(observations, columns + 1))
         output_design(:, 1) = 1.0_dp
         output_design(:, 2:) = design
         select case (model%estimator)
         case (nnfor_estimator_least_squares)
            output_fit = nnfor_ridge_output_fit(output_design, y, [0.0_dp], weight)
         case (nnfor_estimator_ridge)
            output_fit = nnfor_ridge_output_fit(output_design, y, &
               selected_lambdas, weight)
         case (nnfor_estimator_stepwise)
            output_fit = nnfor_stepwise_output_fit(output_design, y)
         case (nnfor_estimator_lasso)
            output_fit = nnfor_lasso_output_fit(output_design, y, &
               selected_lambdas, weight)
         end select
         if (output_fit%info /= 0) then
            model%info = 2
            return
         end if
         model%members(member)%layers = layer_weights(:, member)
         model%members(member)%input_weights = &
            layer_weights(1, member)%input_weights
         allocate(model%members(member)%hidden_counts(layers_count))
         do layer = 1, layers_count
            model%members(member)%hidden_counts(layer) = &
               size(layer_weights(layer, member)%input_weights, 2)
         end do
         model%members(member)%hidden_count = hidden_count
         model%members(member)%output_bias = output_fit%coefficients(1)
         model%members(member)%output_weights = &
            output_fit%coefficients(2:hidden_count + 1)
         model%members(member)%lambda = output_fit%lambda
         offset = hidden_count + 1
         if (model%direct) model%members(member)%direct_weights = &
            output_fit%coefficients(offset + 1:offset + predictors_count)
         prediction = model%members(member)%output_bias + &
            matmul(hidden, model%members(member)%output_weights)
         if (model%direct) prediction = prediction + &
            matmul(x, model%members(member)%direct_weights)
         if (model%scaled) prediction = nnfor_apply_scale(prediction, &
            model%response_scale, .true.)
         model%fitted_all(:, member) = prediction
         deallocate(design, output_design)
      end do
      model%fitted = nnfor_combine(model%fitted_all, model%combination)
      model%residuals = response - model%fitted
      model%mse = sum(model%residuals**2)/real(observations, dp)
   end function nnfor_elm_fast_layers_from_weights

   function nnfor_elm_fast_layers(response, predictors, hidden_counts, &
      repetitions, estimator, combination, direct, scale_data, lambdas, &
      validation_weight, orthogonal) result(model)
      !! Fit a multilayer ELM ensemble using shared uniform random weights.
      real(dp), intent(in) :: response(:) !! Regression target.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-predictor matrix.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: estimator !! Output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      logical, intent(in), optional :: scale_data !! Scale inputs and target.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate regularization penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      logical, intent(in), optional :: orthogonal !! Orthogonalize feasible random hidden weights.
      type(nnfor_elm_fast_t) :: model
      type(nnfor_elm_layer_t), allocatable :: weights(:, :)
      real(dp) :: bound
      integer :: reps, layer, member, row, column, previous, orthogonalized
      logical :: selected_orthogonal, success, layer_success

      reps = 20
      if (present(repetitions)) reps = repetitions
      selected_orthogonal = .false.
      if (present(orthogonal)) selected_orthogonal = orthogonal
      if (size(hidden_counts) < 1 .or. any(hidden_counts < 1) .or. &
         reps < 1 .or. size(predictors, 2) < 1) then
         model%info = 1
         return
      end if
      allocate(weights(size(hidden_counts), reps))
      previous = size(predictors, 2)
      orthogonalized = 0
      do layer = 1, size(hidden_counts)
         bound = 1.0_dp/sqrt(real(previous, dp))
         do member = 1, reps
            allocate(weights(layer, member)%input_weights( &
               previous + 1, hidden_counts(layer)))
            do column = 1, hidden_counts(layer)
               do row = 1, previous + 1
                  weights(layer, member)%input_weights(row, column) = &
                     bound*(2.0_dp*random_uniform() - 1.0_dp)
               end do
            end do
         end do
         if (selected_orthogonal .and. previous + 1 >= hidden_counts(layer)) then
            layer_success = .true.
            do member = 1, reps
               call orthogonalize_columns( &
                  weights(layer, member)%input_weights, success)
               layer_success = layer_success .and. success
            end do
            if (layer_success) orthogonalized = orthogonalized + 1
         end if
         previous = hidden_counts(layer)
      end do
      model = nnfor_elm_fast_layers_from_weights(response, predictors, weights, &
         estimator, combination, direct, scale_data, lambdas, validation_weight)
      model%orthogonalized_layer_count = orthogonalized
      model%orthogonal = orthogonalized > 0
   end function nnfor_elm_fast_layers

   function nnfor_elm_fast(response, predictors, hidden_count, repetitions, &
      estimator, combination, direct, scale_data, lambdas, &
      validation_weight, orthogonal) result(model)
      !! Fit a fast ELM ensemble using shared uniform random weights.
      real(dp), intent(in) :: response(:) !! Regression target.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-predictor matrix.
      integer, intent(in), optional :: hidden_count !! Hidden neurons per repetition.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: estimator !! Least-squares or ridge output estimator.
      integer, intent(in), optional :: combination !! Mean or median ensemble operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      logical, intent(in), optional :: scale_data !! Scale inputs and target.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge penalties.
      real(dp), intent(in), optional :: validation_weight !! Ridge-training fraction.
      logical, intent(in), optional :: orthogonal !! Orthogonalize feasible random hidden weights.
      type(nnfor_elm_fast_t) :: model
      real(dp), allocatable :: weights(:, :, :)
      real(dp) :: bound
      integer :: hidden, reps, row, column, member
      logical :: selected_orthogonal, success, all_success

      hidden = 100
      if (present(hidden_count)) hidden = hidden_count
      reps = 20
      if (present(repetitions)) reps = repetitions
      selected_orthogonal = .false.
      if (present(orthogonal)) selected_orthogonal = orthogonal
      all_success = .true.
      if (hidden < 1 .or. reps < 1 .or. size(predictors, 2) < 1) then
         model%info = 1
         return
      end if
      allocate(weights(size(predictors, 2) + 1, hidden, reps))
      bound = 1.0_dp/sqrt(real(size(predictors, 2), dp))
      do member = 1, reps
         do column = 1, hidden
            do row = 1, size(predictors, 2) + 1
               weights(row, column, member) = &
                  bound*(2.0_dp*random_uniform() - 1.0_dp)
            end do
         end do
         if (selected_orthogonal .and. size(predictors, 2) + 1 >= hidden) then
            call orthogonalize_columns(weights(:, :, member), success)
            all_success = all_success .and. success
         end if
      end do
      model = nnfor_elm_fast_from_weights(response, predictors, weights, &
         estimator, combination, direct, scale_data, lambdas, validation_weight)
      if (selected_orthogonal .and. size(predictors, 2) + 1 >= hidden .and. &
         all_success) then
         model%orthogonalized_layer_count = 1
         model%orthogonal = .true.
      end if
   end function nnfor_elm_fast

   pure function nnfor_elm_fast_predict(model, predictors) result(out)
      !! Predict a regression matrix with every ELM ensemble member.
      type(nnfor_elm_fast_t), intent(in) :: model !! Fitted fast ELM ensemble.
      real(dp), intent(in) :: predictors(:, :) !! New observation-by-predictor matrix.
      type(nnfor_forecast_t) :: out
      real(dp), allocatable :: prediction(:)
      integer :: observations, member

      observations = size(predictors, 1)
      if (model%info /= 0 .or. size(predictors, 2) /= model%predictor_count .or. &
         .not. all(ieee_is_finite(predictors))) then
         out%info = 1
         return
      end if
      allocate(out%all_mean(observations, model%repetitions))
      do member = 1, model%repetitions
         prediction = predict_elm_member(model, member, predictors)
         out%all_mean(:, member) = prediction
      end do
      out%mean = nnfor_combine(out%all_mean, model%combination)
   end function nnfor_elm_fast_predict

   pure function nnfor_elm_from_weights(series, lags, input_weights, &
      estimator, combination, direct, lambdas, validation_weight) result(model)
      !! Fit an autoregressive ELM from supplied hidden-layer weights.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive autoregressive lags.
      real(dp), intent(in) :: input_weights(:, :, :) !! ELM hidden-layer weights.
      integer, intent(in), optional :: estimator !! Least-squares or ridge output estimator.
      integer, intent(in), optional :: combination !! Mean or median ensemble operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge penalties.
      real(dp), intent(in), optional :: validation_weight !! Ridge-training fraction.
      type(nnfor_elm_model_t) :: model
      real(dp), allocatable :: predictors(:, :), response(:)
      integer :: maximum_lag, row, lag

      if (size(lags) < 1 .or. any(lags < 1) .or. &
         size(series) <= maxval(lags) + 2 .or. &
         .not. all(ieee_is_finite(series))) then
         model%info = 1
         return
      end if
      maximum_lag = maxval(lags)
      allocate(response(size(series) - maximum_lag))
      allocate(predictors(size(response), size(lags)))
      response = series(maximum_lag + 1:)
      do lag = 1, size(lags)
         do row = 1, size(response)
            predictors(row, lag) = &
               series(maximum_lag + row - lags(lag))
         end do
      end do
      model%network = nnfor_elm_fast_from_weights(response, predictors, &
         input_weights, estimator, combination, direct, scale_data=.true., &
         lambdas=lambdas, validation_weight=validation_weight)
      model%lags = lags
      model%series = series
      model%fitted = model%network%fitted
      model%residuals = response - model%fitted
      model%mse = sum(model%residuals**2)/real(size(response), dp)
      model%info = model%network%info
   end function nnfor_elm_from_weights

   function nnfor_elm(series, lags, hidden_count, repetitions, estimator, &
      combination, direct, lambdas, validation_weight, orthogonal) result(model)
      !! Fit an autoregressive ELM using shared uniform random weights.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive autoregressive lags.
      integer, intent(in), optional :: hidden_count !! Hidden neurons per repetition.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: estimator !! Least-squares or ridge output estimator.
      integer, intent(in), optional :: combination !! Mean or median ensemble operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge penalties.
      real(dp), intent(in), optional :: validation_weight !! Ridge-training fraction.
      logical, intent(in), optional :: orthogonal !! Orthogonalize feasible random hidden weights.
      type(nnfor_elm_model_t) :: model
      real(dp), allocatable :: weights(:, :, :)
      real(dp) :: bound
      integer :: hidden, reps, row, column, member
      logical :: selected_orthogonal, success, all_success

      hidden = 5
      if (present(hidden_count)) hidden = hidden_count
      reps = 20
      if (present(repetitions)) reps = repetitions
      selected_orthogonal = .false.
      if (present(orthogonal)) selected_orthogonal = orthogonal
      all_success = .true.
      if (size(lags) < 1 .or. hidden < 1 .or. reps < 1) then
         model%info = 1
         return
      end if
      allocate(weights(size(lags) + 1, hidden, reps))
      bound = 1.0_dp/sqrt(real(size(lags), dp))
      do member = 1, reps
         do column = 1, hidden
            do row = 1, size(lags) + 1
               weights(row, column, member) = &
                  bound*(2.0_dp*random_uniform() - 1.0_dp)
            end do
         end do
         if (selected_orthogonal .and. size(lags) + 1 >= hidden) then
            call orthogonalize_columns(weights(:, :, member), success)
            all_success = all_success .and. success
         end if
      end do
      model = nnfor_elm_from_weights(series, lags, weights, estimator, &
         combination, direct, lambdas, validation_weight)
      if (selected_orthogonal .and. size(lags) + 1 >= hidden .and. &
         all_success) then
         model%network%orthogonalized_layer_count = 1
         model%network%orthogonal = .true.
      end if
   end function nnfor_elm

   pure function nnfor_elm_layers_from_weights(series, lags, layer_weights, &
      estimator, combination, direct, lambdas, validation_weight) result(model)
      !! Fit an autoregressive multilayer ELM from supplied hidden weights.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive autoregressive lags.
      type(nnfor_elm_layer_t), intent(in) :: layer_weights(:, :) !! Layers by repetition.
      integer, intent(in), optional :: estimator !! Output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate regularization penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      type(nnfor_elm_model_t) :: model
      real(dp), allocatable :: predictors(:, :), response(:)
      integer :: maximum_lag, row, lag

      if (size(lags) < 1 .or. any(lags < 1) .or. &
         size(series) <= maxval(lags) + 2 .or. &
         .not. all(ieee_is_finite(series))) then
         model%info = 1
         return
      end if
      maximum_lag = maxval(lags)
      allocate(response(size(series) - maximum_lag))
      allocate(predictors(size(response), size(lags)))
      response = series(maximum_lag + 1:)
      do lag = 1, size(lags)
         do row = 1, size(response)
            predictors(row, lag) = series(maximum_lag + row - lags(lag))
         end do
      end do
      model%network = nnfor_elm_fast_layers_from_weights(response, predictors, &
         layer_weights, estimator, combination, direct, scale_data=.true., &
         lambdas=lambdas, validation_weight=validation_weight)
      model%lags = lags
      model%series = series
      model%fitted = model%network%fitted
      model%residuals = response - model%fitted
      model%mse = sum(model%residuals**2)/real(size(response), dp)
      model%info = model%network%info
   end function nnfor_elm_layers_from_weights

   function nnfor_elm_layers(series, lags, hidden_counts, repetitions, &
      estimator, combination, direct, lambdas, validation_weight, orthogonal) &
      result(model)
      !! Fit an autoregressive multilayer ELM using shared random weights.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive autoregressive lags.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: estimator !! Output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate regularization penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      logical, intent(in), optional :: orthogonal !! Orthogonalize feasible random hidden weights.
      type(nnfor_elm_model_t) :: model
      type(nnfor_elm_layer_t), allocatable :: weights(:, :)
      real(dp) :: bound
      integer :: reps, layer, member, row, column, previous, orthogonalized
      logical :: selected_orthogonal, success, layer_success

      reps = 20
      if (present(repetitions)) reps = repetitions
      selected_orthogonal = .false.
      if (present(orthogonal)) selected_orthogonal = orthogonal
      if (size(lags) < 1 .or. size(hidden_counts) < 1 .or. &
         any(hidden_counts < 1) .or. reps < 1) then
         model%info = 1
         return
      end if
      allocate(weights(size(hidden_counts), reps))
      previous = size(lags)
      orthogonalized = 0
      do layer = 1, size(hidden_counts)
         bound = 1.0_dp/sqrt(real(previous, dp))
         do member = 1, reps
            allocate(weights(layer, member)%input_weights( &
               previous + 1, hidden_counts(layer)))
            do column = 1, hidden_counts(layer)
               do row = 1, previous + 1
                  weights(layer, member)%input_weights(row, column) = &
                     bound*(2.0_dp*random_uniform() - 1.0_dp)
               end do
            end do
         end do
         if (selected_orthogonal .and. previous + 1 >= hidden_counts(layer)) then
            layer_success = .true.
            do member = 1, reps
               call orthogonalize_columns( &
                  weights(layer, member)%input_weights, success)
               layer_success = layer_success .and. success
            end do
            if (layer_success) orthogonalized = orthogonalized + 1
         end if
         previous = hidden_counts(layer)
      end do
      model = nnfor_elm_layers_from_weights(series, lags, weights, estimator, &
         combination, direct, lambdas, validation_weight)
      model%network%orthogonalized_layer_count = orthogonalized
      model%network%orthogonal = orthogonalized > 0
   end function nnfor_elm_layers

   pure function nnfor_elm_forecast(model, horizon, series) result(out)
      !! Recursively forecast an autoregressive ELM ensemble.
      type(nnfor_elm_model_t), intent(in) :: model !! Fitted autoregressive ELM.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: series(:) !! Optional replacement history.
      type(nnfor_forecast_t) :: out
      real(dp), allocatable :: history(:, :), predictors(:, :), prediction(:)
      integer :: member, step, lag, observations

      if (model%info /= 0 .or. horizon < 1) then
         out%info = 1
         return
      end if
      observations = size(model%series)
      if (present(series)) observations = size(series)
      if (observations < maxval(model%lags)) then
         out%info = 1
         return
      end if
      allocate(history(observations + horizon, model%network%repetitions))
      if (present(series)) then
         history(:observations, :) = spread(series, 2, model%network%repetitions)
      else
         history(:observations, :) = spread(model%series, 2, &
            model%network%repetitions)
      end if
      allocate(out%all_mean(horizon, model%network%repetitions))
      allocate(predictors(1, size(model%lags)))
      do member = 1, model%network%repetitions
         do step = 1, horizon
            do lag = 1, size(model%lags)
               predictors(1, lag) = history(observations + step - &
                  model%lags(lag), member)
            end do
            prediction = predict_elm_member(model%network, member, predictors)
            history(observations + step, member) = prediction(1)
            out%all_mean(step, member) = prediction(1)
         end do
      end do
      out%mean = nnfor_combine(out%all_mean, model%network%combination)
      out%fitted = model%fitted
      out%residuals = model%residuals
   end function nnfor_elm_forecast

   pure function nnfor_elm_preprocessed_from_weights(series, lags, &
      input_weights, difference_lags, period, seasonal_type, exogenous, &
      exogenous_lags, estimator, combination, direct, lambdas, &
      validation_weight, periods, seasonal_types) result(model)
      !! Fit a differenced, seasonal, or exogenous autoregressive ELM.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive lags of the final differenced series.
      real(dp), intent(in) :: input_weights(:, :, :) !! ELM hidden-layer weights.
      integer, intent(in), optional :: difference_lags(:) !! Sequential positive difference lags.
      integer, intent(in), optional :: period !! Seasonal period for deterministic inputs.
      integer, intent(in), optional :: seasonal_type !! None, binary, or trigonometric seasonality.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time exogenous regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Nonnegative lags by regressor.
      integer, intent(in), optional :: estimator !! ELM output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge or lasso penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Input type for each seasonal period.
      type(nnfor_elm_model_t) :: model
      real(dp), allocatable :: restored_all(:, :)

      model%preprocessing = nnfor_preprocess(series, lags, difference_lags, &
         period, seasonal_type, exogenous, exogenous_lags, periods, &
         seasonal_types)
      if (model%preprocessing%info /= 0) then
         model%info = model%preprocessing%info
         return
      end if
      model%network = nnfor_elm_fast_from_weights( &
         model%preprocessing%response, model%preprocessing%predictors, &
         input_weights, estimator, combination, direct, scale_data=.true., &
         lambdas=lambdas, validation_weight=validation_weight)
      if (model%network%info /= 0) then
         model%info = model%network%info
         return
      end if
      model%lags = lags
      model%series = series
      restored_all = restore_fitted_values(model%preprocessing, &
         model%network%fitted_all)
      model%fitted = nnfor_combine(restored_all, model%network%combination)
      model%residuals = model%preprocessing%original_response - model%fitted
      model%mse = sum(model%residuals**2)/real(size(model%residuals), dp)
      model%extended_preprocessing = .true.
   end function nnfor_elm_preprocessed_from_weights

   function nnfor_elm_preprocessed(series, lags, hidden_count, repetitions, &
      difference_lags, period, seasonal_type, exogenous, exogenous_lags, &
      estimator, combination, direct, lambdas, validation_weight, periods, &
      seasonal_types, orthogonal) result(model)
      !! Fit a preprocessed autoregressive ELM using shared random weights.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive lags of the final differenced series.
      integer, intent(in), optional :: hidden_count !! Hidden neurons per repetition.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: difference_lags(:) !! Sequential positive difference lags.
      integer, intent(in), optional :: period !! Seasonal period for deterministic inputs.
      integer, intent(in), optional :: seasonal_type !! None, binary, or trigonometric seasonality.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time exogenous regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Nonnegative lags by regressor.
      integer, intent(in), optional :: estimator !! ELM output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge or lasso penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Input type for each seasonal period.
      logical, intent(in), optional :: orthogonal !! Orthogonalize feasible random hidden weights.
      type(nnfor_elm_model_t) :: model
      type(nnfor_preprocessing_t) :: preprocessing
      real(dp), allocatable :: weights(:, :, :)
      real(dp) :: bound
      integer :: hidden, reps, predictors_count, row, column, member
      logical :: selected_orthogonal, success, all_success

      hidden = 5
      if (present(hidden_count)) hidden = hidden_count
      reps = 20
      if (present(repetitions)) reps = repetitions
      selected_orthogonal = .false.
      if (present(orthogonal)) selected_orthogonal = orthogonal
      all_success = .true.
      preprocessing = nnfor_preprocess(series, lags, difference_lags, period, &
         seasonal_type, exogenous, exogenous_lags, periods, seasonal_types)
      if (preprocessing%info /= 0 .or. hidden < 1 .or. reps < 1) then
         model%info = merge(preprocessing%info, 1, preprocessing%info /= 0)
         return
      end if
      predictors_count = size(preprocessing%predictors, 2)
      allocate(weights(predictors_count + 1, hidden, reps))
      bound = 1.0_dp/sqrt(real(predictors_count, dp))
      do member = 1, reps
         do column = 1, hidden
            do row = 1, predictors_count + 1
               weights(row, column, member) = &
                  bound*(2.0_dp*random_uniform() - 1.0_dp)
            end do
         end do
         if (selected_orthogonal .and. predictors_count + 1 >= hidden) then
            call orthogonalize_columns(weights(:, :, member), success)
            all_success = all_success .and. success
         end if
      end do
      model = nnfor_elm_preprocessed_from_weights(series, lags, weights, &
         difference_lags, period, seasonal_type, exogenous, exogenous_lags, &
         estimator, combination, direct, lambdas, validation_weight, periods, &
         seasonal_types)
      if (selected_orthogonal .and. predictors_count + 1 >= hidden .and. &
         all_success) then
         model%network%orthogonalized_layer_count = 1
         model%network%orthogonal = .true.
      end if
   end function nnfor_elm_preprocessed

   pure function nnfor_elm_preprocessed_layers_from_weights(series, lags, &
      layer_weights, difference_lags, period, seasonal_type, exogenous, &
      exogenous_lags, estimator, combination, direct, lambdas, &
      validation_weight, periods, seasonal_types) result(model)
      !! Fit a preprocessed multilayer ELM from supplied hidden weights.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive lags of the differenced series.
      type(nnfor_elm_layer_t), intent(in) :: layer_weights(:, :) !! Layers by repetition.
      integer, intent(in), optional :: difference_lags(:) !! Sequential difference lags.
      integer, intent(in), optional :: period !! Scalar seasonal period.
      integer, intent(in), optional :: seasonal_type !! Scalar seasonal input type.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Lags by regressor.
      integer, intent(in), optional :: estimator !! Output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate regularization penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Type for each seasonal period.
      type(nnfor_elm_model_t) :: model
      real(dp), allocatable :: restored_all(:, :)

      model%preprocessing = nnfor_preprocess(series, lags, difference_lags, &
         period, seasonal_type, exogenous, exogenous_lags, periods, &
         seasonal_types)
      if (model%preprocessing%info /= 0) then
         model%info = model%preprocessing%info
         return
      end if
      model%network = nnfor_elm_fast_layers_from_weights( &
         model%preprocessing%response, model%preprocessing%predictors, &
         layer_weights, estimator, combination, direct, scale_data=.true., &
         lambdas=lambdas, validation_weight=validation_weight)
      if (model%network%info /= 0) then
         model%info = model%network%info
         return
      end if
      model%lags = lags
      model%series = series
      restored_all = restore_fitted_values(model%preprocessing, &
         model%network%fitted_all)
      model%fitted = nnfor_combine(restored_all, model%network%combination)
      model%residuals = model%preprocessing%original_response - model%fitted
      model%mse = sum(model%residuals**2)/real(size(model%residuals), dp)
      model%extended_preprocessing = .true.
   end function nnfor_elm_preprocessed_layers_from_weights

   function nnfor_elm_preprocessed_layers(series, lags, hidden_counts, &
      repetitions, difference_lags, period, seasonal_type, exogenous, &
      exogenous_lags, estimator, combination, direct, lambdas, &
      validation_weight, periods, seasonal_types, orthogonal) result(model)
      !! Fit a preprocessed multilayer ELM using shared random weights.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive lags of the differenced series.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: difference_lags(:) !! Sequential difference lags.
      integer, intent(in), optional :: period !! Scalar seasonal period.
      integer, intent(in), optional :: seasonal_type !! Scalar seasonal input type.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Lags by regressor.
      integer, intent(in), optional :: estimator !! Output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate regularization penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Type for each seasonal period.
      logical, intent(in), optional :: orthogonal !! Orthogonalize feasible random hidden weights.
      type(nnfor_elm_model_t) :: model
      type(nnfor_preprocessing_t) :: preprocessing
      type(nnfor_elm_layer_t), allocatable :: weights(:, :)
      real(dp) :: bound
      integer :: reps, layer, member, row, column, previous, orthogonalized
      logical :: selected_orthogonal, success, layer_success

      reps = 20
      if (present(repetitions)) reps = repetitions
      selected_orthogonal = .false.
      if (present(orthogonal)) selected_orthogonal = orthogonal
      preprocessing = nnfor_preprocess(series, lags, difference_lags, period, &
         seasonal_type, exogenous, exogenous_lags, periods, seasonal_types)
      if (preprocessing%info /= 0 .or. size(hidden_counts) < 1 .or. &
         any(hidden_counts < 1) .or. reps < 1) then
         model%info = merge(preprocessing%info, 1, preprocessing%info /= 0)
         return
      end if
      allocate(weights(size(hidden_counts), reps))
      previous = size(preprocessing%predictors, 2)
      orthogonalized = 0
      do layer = 1, size(hidden_counts)
         bound = 1.0_dp/sqrt(real(previous, dp))
         do member = 1, reps
            allocate(weights(layer, member)%input_weights( &
               previous + 1, hidden_counts(layer)))
            do column = 1, hidden_counts(layer)
               do row = 1, previous + 1
                  weights(layer, member)%input_weights(row, column) = &
                     bound*(2.0_dp*random_uniform() - 1.0_dp)
               end do
            end do
         end do
         if (selected_orthogonal .and. previous + 1 >= hidden_counts(layer)) then
            layer_success = .true.
            do member = 1, reps
               call orthogonalize_columns( &
                  weights(layer, member)%input_weights, success)
               layer_success = layer_success .and. success
            end do
            if (layer_success) orthogonalized = orthogonalized + 1
         end if
         previous = hidden_counts(layer)
      end do
      model = nnfor_elm_preprocessed_layers_from_weights(series, lags, weights, &
         difference_lags, period, seasonal_type, exogenous, exogenous_lags, &
         estimator, combination, direct, lambdas, validation_weight, periods, &
         seasonal_types)
      model%network%orthogonalized_layer_count = orthogonalized
      model%network%orthogonal = orthogonalized > 0
   end function nnfor_elm_preprocessed_layers

   pure function nnfor_elm_preprocessed_forecast(model, horizon, exogenous) &
      result(out)
      !! Recursively forecast and inverse-difference a preprocessed ELM ensemble.
      type(nnfor_elm_model_t), intent(in) :: model !! Fitted preprocessed ELM model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: exogenous(:, :) !! Exogenous history and future values.
      type(nnfor_forecast_t) :: out
      real(dp), allocatable :: future(:, :, :), predictors(:, :), prediction(:)
      real(dp) :: base
      integer :: differences, deepest, member, step, lag, variable
      integer :: column, time, stage, source_step, source_index

      if (model%info /= 0 .or. .not. model%extended_preprocessing .or. &
         horizon < 1) then
         out%info = 1
         return
      end if
      if (model%preprocessing%exogenous_count > 0) then
         if (.not. present(exogenous)) then
            out%info = 1
            return
         end if
         if (size(exogenous, 1) < size(model%series) + horizon .or. &
            size(exogenous, 2) /= model%preprocessing%exogenous_count .or. &
            .not. all(ieee_is_finite(exogenous))) then
            out%info = 1
            return
         end if
      end if
      differences = size(model%preprocessing%difference_lags)
      deepest = differences + 1
      allocate(future(horizon, deepest, model%network%repetitions), &
         source=0.0_dp)
      allocate(out%all_mean(horizon, model%network%repetitions))
      allocate(predictors(1, size(model%preprocessing%predictors, 2)))
      do member = 1, model%network%repetitions
         do step = 1, horizon
            column = 0
            do lag = 1, size(model%preprocessing%lags)
               column = column + 1
               source_step = step - model%preprocessing%lags(lag)
               if (source_step <= 0) then
                  source_index = size(model%preprocessing% &
                     difference_levels(deepest)%values) + source_step
                  predictors(1, column) = model%preprocessing% &
                     difference_levels(deepest)%values(source_index)
               else
                  predictors(1, column) = future(source_step, deepest, member)
               end if
            end do
            time = size(model%series) + step
            do variable = 1, model%preprocessing%exogenous_count
               do lag = 1, size(model%preprocessing% &
                  exogenous_lags(variable)%values)
                  column = column + 1
                  predictors(1, column) = exogenous(time - &
                     model%preprocessing%exogenous_lags(variable)%values(lag), &
                     variable)
               end do
            end do
            call append_seasonal_predictors(model%preprocessing, time, &
               predictors, column)
            prediction = predict_elm_member(model%network, member, predictors)
            future(step, deepest, member) = prediction(1)
            do stage = differences, 1, -1
               lag = model%preprocessing%difference_lags(stage)
               source_step = step - lag
               if (source_step <= 0) then
                  source_index = size(model%preprocessing% &
                     difference_levels(stage)%values) + source_step
                  base = model%preprocessing%difference_levels(stage)% &
                     values(source_index)
               else
                  base = future(source_step, stage, member)
               end if
               future(step, stage, member) = future(step, stage + 1, member) + &
                  base
            end do
            out%all_mean(step, member) = future(step, 1, member)
         end do
      end do
      out%mean = nnfor_combine(out%all_mean, model%network%combination)
      out%fitted = model%fitted
      out%residuals = model%residuals
   end function nnfor_elm_preprocessed_forecast

   pure function nnfor_mlp_from_initial(series, lags, hidden_count, &
      initial_parameters, combination, max_iterations, tolerance, decay) &
      result(model)
      !! Fit a reproducible lagged MLP ensemble from supplied initial parameters.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive autoregressive lags.
      integer, intent(in) :: hidden_count !! Hidden neurons in each network.
      real(dp), intent(in) :: initial_parameters(:, :) !! Parameter vectors by repetition.
      integer, intent(in), optional :: combination !! Mean or median ensemble operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      type(nnfor_mlp_model_t) :: model
      type(nnfor_scaled_t) :: scaled
      real(dp), allocatable :: predictors(:, :), response(:), response_matrix(:, :)
      real(dp), allocatable :: prediction(:, :)
      integer :: maximum_lag, row, lag, member, parameter_count

      model%combination = nnfor_combine_median
      if (present(combination)) model%combination = combination
      parameter_count = neural_network_parameter_count(size(lags), &
         hidden_count, 1)
      if (size(lags) < 1 .or. any(lags < 1) .or. hidden_count < 1 .or. &
         size(initial_parameters, 1) /= parameter_count .or. &
         size(initial_parameters, 2) < 1 .or. &
         size(series) <= maxval(lags) + 2 .or. &
         (model%combination /= nnfor_combine_mean .and. &
         model%combination /= nnfor_combine_median .and. &
         model%combination /= nnfor_combine_mode) .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(initial_parameters))) then
         model%info = 1
         return
      end if
      maximum_lag = maxval(lags)
      allocate(response(size(series) - maximum_lag))
      allocate(predictors(size(response), size(lags)))
      response = series(maximum_lag + 1:)
      do lag = 1, size(lags)
         do row = 1, size(response)
            predictors(row, lag) = series(maximum_lag + row - lags(lag))
         end do
      end do
      scaled = nnfor_linscale(response, -0.8_dp, 0.8_dp)
      model%response_scale = scaled%scale
      allocate(response_matrix(size(response), 1))
      response_matrix(:, 1) = scaled%values
      model%hidden_count = hidden_count
      model%hidden_counts = [hidden_count]
      model%repetitions = size(initial_parameters, 2)
      model%lags = lags
      model%series = series
      allocate(model%members(model%repetitions))
      allocate(model%fitted_all(size(response), model%repetitions))
      do member = 1, model%repetitions
         model%members(member) = neural_network_fit(predictors, response_matrix, &
            hidden_count, max_iterations, tolerance, decay, &
            initial_parameters(:, member))
         if (model%members(member)%info /= 0) then
            model%info = 2
            return
         end if
         prediction = neural_network_predict(model%members(member), predictors)
         model%fitted_all(:, member) = nnfor_apply_scale(prediction(:, 1), &
            model%response_scale, .true.)
      end do
      model%fitted = nnfor_combine(model%fitted_all, model%combination)
      model%residuals = response - model%fitted
      model%mse = sum(model%residuals**2)/real(size(response), dp)
   end function nnfor_mlp_from_initial

   function nnfor_mlp(series, lags, hidden_count, repetitions, combination, &
      max_iterations, tolerance, decay) result(model)
      !! Fit a lagged MLP ensemble using shared random initial parameters.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive autoregressive lags.
      integer, intent(in), optional :: hidden_count !! Hidden neurons per network.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: combination !! Mean or median ensemble operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      type(nnfor_mlp_model_t) :: model
      real(dp), allocatable :: initial(:, :)
      integer :: hidden, reps, parameters, member, parameter

      hidden = 5
      if (present(hidden_count)) hidden = hidden_count
      reps = 20
      if (present(repetitions)) reps = repetitions
      if (hidden < 1 .or. reps < 1 .or. size(lags) < 1) then
         model%info = 1
         return
      end if
      parameters = neural_network_parameter_count(size(lags), hidden, 1)
      allocate(initial(parameters, reps))
      do member = 1, reps
         do parameter = 1, parameters
            initial(parameter, member) = 0.3_dp*(2.0_dp*random_uniform() - 1.0_dp)
         end do
      end do
      model = nnfor_mlp_from_initial(series, lags, hidden, initial, &
         combination, max_iterations, tolerance, decay)
   end function nnfor_mlp

   pure function nnfor_mlp_layers_from_initial(series, lags, hidden_counts, &
      initial_parameters, combination, max_iterations, tolerance, decay) &
      result(model)
      !! Fit a reproducible multilayer MLP ensemble from supplied initial parameters.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive autoregressive lags.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      real(dp), intent(in) :: initial_parameters(:, :) !! Parameter vectors by repetition.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      type(nnfor_mlp_model_t) :: model
      type(nnfor_scaled_t) :: scaled
      real(dp), allocatable :: predictors(:, :), response(:), response_matrix(:, :)
      real(dp), allocatable :: prediction(:, :)
      integer :: maximum_lag, row, lag, member, parameter_count

      model%combination = nnfor_combine_median
      if (present(combination)) model%combination = combination
      parameter_count = neural_network_parameter_count(size(lags), &
         hidden_counts, 1)
      if (size(lags) < 1 .or. any(lags < 1) .or. &
         size(hidden_counts) < 1 .or. any(hidden_counts < 1) .or. &
         size(initial_parameters, 1) /= parameter_count .or. &
         size(initial_parameters, 2) < 1 .or. &
         size(series) <= maxval(lags) + 2 .or. &
         (model%combination /= nnfor_combine_mean .and. &
         model%combination /= nnfor_combine_median .and. &
         model%combination /= nnfor_combine_mode) .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(initial_parameters))) then
         model%info = 1
         return
      end if
      maximum_lag = maxval(lags)
      allocate(response(size(series) - maximum_lag))
      allocate(predictors(size(response), size(lags)))
      response = series(maximum_lag + 1:)
      do lag = 1, size(lags)
         do row = 1, size(response)
            predictors(row, lag) = series(maximum_lag + row - lags(lag))
         end do
      end do
      scaled = nnfor_linscale(response, -0.8_dp, 0.8_dp)
      model%response_scale = scaled%scale
      allocate(response_matrix(size(response), 1))
      response_matrix(:, 1) = scaled%values
      model%hidden_count = hidden_counts(1)
      model%hidden_counts = hidden_counts
      model%repetitions = size(initial_parameters, 2)
      model%lags = lags
      model%series = series
      allocate(model%members(model%repetitions))
      allocate(model%fitted_all(size(response), model%repetitions))
      do member = 1, model%repetitions
         model%members(member) = neural_network_fit(predictors, response_matrix, &
            hidden_counts, max_iterations, tolerance, decay, &
            initial_parameters(:, member))
         if (model%members(member)%info /= 0) then
            model%info = 2
            return
         end if
         prediction = neural_network_predict(model%members(member), predictors)
         model%fitted_all(:, member) = nnfor_apply_scale(prediction(:, 1), &
            model%response_scale, .true.)
      end do
      model%fitted = nnfor_combine(model%fitted_all, model%combination)
      model%residuals = response - model%fitted
      model%mse = sum(model%residuals**2)/real(size(response), dp)
   end function nnfor_mlp_layers_from_initial

   function nnfor_mlp_layers(series, lags, hidden_counts, repetitions, &
      combination, max_iterations, tolerance, decay) result(model)
      !! Fit a multilayer MLP ensemble using shared random initial parameters.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive autoregressive lags.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      type(nnfor_mlp_model_t) :: model
      real(dp), allocatable :: initial(:, :)
      integer :: reps, parameters, member, parameter

      reps = 20
      if (present(repetitions)) reps = repetitions
      if (size(hidden_counts) < 1 .or. any(hidden_counts < 1) .or. &
         reps < 1 .or. size(lags) < 1) then
         model%info = 1
         return
      end if
      parameters = neural_network_parameter_count(size(lags), hidden_counts, 1)
      allocate(initial(parameters, reps))
      do member = 1, reps
         do parameter = 1, parameters
            initial(parameter, member) = 0.3_dp*(2.0_dp*random_uniform() - 1.0_dp)
         end do
      end do
      model = nnfor_mlp_layers_from_initial(series, lags, hidden_counts, initial, &
         combination, max_iterations, tolerance, decay)
   end function nnfor_mlp_layers

   pure function nnfor_mlp_forecast(model, horizon, series) result(out)
      !! Recursively forecast every member of a lagged MLP ensemble.
      type(nnfor_mlp_model_t), intent(in) :: model !! Fitted MLP ensemble.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: series(:) !! Optional replacement history.
      type(nnfor_forecast_t) :: out
      real(dp), allocatable :: history(:, :), predictors(:, :), prediction(:, :)
      real(dp), allocatable :: restored(:)
      integer :: member, step, lag, observations

      if (model%info /= 0 .or. horizon < 1) then
         out%info = 1
         return
      end if
      observations = size(model%series)
      if (present(series)) observations = size(series)
      if (observations < maxval(model%lags)) then
         out%info = 1
         return
      end if
      allocate(history(observations + horizon, model%repetitions))
      if (present(series)) then
         history(:observations, :) = spread(series, 2, model%repetitions)
      else
         history(:observations, :) = spread(model%series, 2, model%repetitions)
      end if
      allocate(out%all_mean(horizon, model%repetitions))
      allocate(predictors(1, size(model%lags)))
      do member = 1, model%repetitions
         do step = 1, horizon
            do lag = 1, size(model%lags)
               predictors(1, lag) = history(observations + step - &
                  model%lags(lag), member)
            end do
            prediction = neural_network_predict(model%members(member), predictors)
            restored = nnfor_apply_scale(prediction(:, 1), &
               model%response_scale, .true.)
            history(observations + step, member) = restored(1)
            out%all_mean(step, member) = restored(1)
         end do
      end do
      out%mean = nnfor_combine(out%all_mean, model%combination)
      out%fitted = model%fitted
      out%residuals = model%residuals
   end function nnfor_mlp_forecast

   pure function nnfor_mlp_preprocessed_from_initial(series, lags, &
      hidden_count, initial_parameters, difference_lags, period, &
      seasonal_type, exogenous, exogenous_lags, combination, max_iterations, &
      tolerance, decay, periods, seasonal_types) result(model)
      !! Fit a reproducible preprocessed MLP ensemble from supplied starts.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive lags of the final differenced series.
      integer, intent(in) :: hidden_count !! Hidden neurons in each network.
      real(dp), intent(in) :: initial_parameters(:, :) !! Parameter vectors by repetition.
      integer, intent(in), optional :: difference_lags(:) !! Sequential positive difference lags.
      integer, intent(in), optional :: period !! Seasonal period for deterministic inputs.
      integer, intent(in), optional :: seasonal_type !! None, binary, or trigonometric seasonality.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time exogenous regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Nonnegative lags by regressor.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Input type for each seasonal period.
      type(nnfor_mlp_model_t) :: model
      type(nnfor_scaled_t) :: scaled
      real(dp), allocatable :: response_matrix(:, :), prediction(:, :)
      real(dp), allocatable :: deepest_all(:, :), restored_all(:, :)
      integer :: member, parameter_count

      model%combination = nnfor_combine_median
      if (present(combination)) model%combination = combination
      model%preprocessing = nnfor_preprocess(series, lags, difference_lags, &
         period, seasonal_type, exogenous, exogenous_lags, periods, &
         seasonal_types)
      if (model%preprocessing%info /= 0) then
         model%info = model%preprocessing%info
         return
      end if
      parameter_count = neural_network_parameter_count( &
         size(model%preprocessing%predictors, 2), hidden_count, 1)
      if (hidden_count < 1 .or. size(initial_parameters, 1) /= parameter_count .or. &
         size(initial_parameters, 2) < 1 .or. &
         (model%combination /= nnfor_combine_mean .and. &
         model%combination /= nnfor_combine_median .and. &
         model%combination /= nnfor_combine_mode) .or. &
         .not. all(ieee_is_finite(initial_parameters))) then
         model%info = 1
         return
      end if
      scaled = nnfor_linscale(model%preprocessing%response, -0.8_dp, 0.8_dp)
      model%response_scale = scaled%scale
      allocate(response_matrix(size(scaled%values), 1))
      response_matrix(:, 1) = scaled%values
      model%hidden_count = hidden_count
      model%hidden_counts = [hidden_count]
      model%repetitions = size(initial_parameters, 2)
      model%lags = lags
      model%series = series
      allocate(model%members(model%repetitions))
      allocate(deepest_all(size(scaled%values), model%repetitions))
      do member = 1, model%repetitions
         model%members(member) = neural_network_fit( &
            model%preprocessing%predictors, response_matrix, hidden_count, &
            max_iterations, tolerance, decay, initial_parameters(:, member))
         if (model%members(member)%info /= 0) then
            model%info = 2
            return
         end if
         prediction = neural_network_predict(model%members(member), &
            model%preprocessing%predictors)
         deepest_all(:, member) = nnfor_apply_scale(prediction(:, 1), &
            model%response_scale, .true.)
      end do
      restored_all = restore_fitted_values(model%preprocessing, deepest_all)
      model%fitted_all = restored_all
      model%fitted = nnfor_combine(restored_all, model%combination)
      model%residuals = model%preprocessing%original_response - model%fitted
      model%mse = sum(model%residuals**2)/real(size(model%residuals), dp)
      model%extended_preprocessing = .true.
   end function nnfor_mlp_preprocessed_from_initial

   function nnfor_mlp_preprocessed(series, lags, hidden_count, repetitions, &
      difference_lags, period, seasonal_type, exogenous, exogenous_lags, &
      combination, max_iterations, tolerance, decay, periods, seasonal_types) &
      result(model)
      !! Fit a preprocessed MLP ensemble using shared random initial parameters.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive lags of the final differenced series.
      integer, intent(in), optional :: hidden_count !! Hidden neurons per network.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: difference_lags(:) !! Sequential positive difference lags.
      integer, intent(in), optional :: period !! Seasonal period for deterministic inputs.
      integer, intent(in), optional :: seasonal_type !! None, binary, or trigonometric seasonality.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time exogenous regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Nonnegative lags by regressor.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Input type for each seasonal period.
      type(nnfor_mlp_model_t) :: model
      type(nnfor_preprocessing_t) :: preprocessing
      real(dp), allocatable :: initial(:, :)
      integer :: hidden, reps, parameters, member, parameter

      hidden = 5
      if (present(hidden_count)) hidden = hidden_count
      reps = 20
      if (present(repetitions)) reps = repetitions
      preprocessing = nnfor_preprocess(series, lags, difference_lags, period, &
         seasonal_type, exogenous, exogenous_lags, periods, seasonal_types)
      if (preprocessing%info /= 0 .or. hidden < 1 .or. reps < 1) then
         model%info = merge(preprocessing%info, 1, preprocessing%info /= 0)
         return
      end if
      parameters = neural_network_parameter_count(size(preprocessing%predictors, 2), &
         hidden, 1)
      allocate(initial(parameters, reps))
      do member = 1, reps
         do parameter = 1, parameters
            initial(parameter, member) = &
               0.3_dp*(2.0_dp*random_uniform() - 1.0_dp)
         end do
      end do
      model = nnfor_mlp_preprocessed_from_initial(series, lags, hidden, initial, &
         difference_lags, period, seasonal_type, exogenous, exogenous_lags, &
         combination, max_iterations, tolerance, decay, periods, seasonal_types)
   end function nnfor_mlp_preprocessed

   pure function nnfor_mlp_preprocessed_layers_from_initial(series, lags, &
      hidden_counts, initial_parameters, difference_lags, period, &
      seasonal_type, exogenous, exogenous_lags, combination, max_iterations, &
      tolerance, decay, periods, seasonal_types) result(model)
      !! Fit a reproducible preprocessed multilayer MLP ensemble.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive lags of the differenced series.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      real(dp), intent(in) :: initial_parameters(:, :) !! Parameter vectors by repetition.
      integer, intent(in), optional :: difference_lags(:) !! Sequential difference lags.
      integer, intent(in), optional :: period !! Scalar seasonal period.
      integer, intent(in), optional :: seasonal_type !! Scalar seasonal input type.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Lags by regressor.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Type for each seasonal period.
      type(nnfor_mlp_model_t) :: model
      type(nnfor_scaled_t) :: scaled
      real(dp), allocatable :: response_matrix(:, :), prediction(:, :)
      real(dp), allocatable :: deepest_all(:, :), restored_all(:, :)
      integer :: member, parameter_count

      model%combination = nnfor_combine_median
      if (present(combination)) model%combination = combination
      model%preprocessing = nnfor_preprocess(series, lags, difference_lags, &
         period, seasonal_type, exogenous, exogenous_lags, periods, &
         seasonal_types)
      if (model%preprocessing%info /= 0) then
         model%info = model%preprocessing%info
         return
      end if
      parameter_count = neural_network_parameter_count( &
         size(model%preprocessing%predictors, 2), hidden_counts, 1)
      if (size(hidden_counts) < 1 .or. any(hidden_counts < 1) .or. &
         size(initial_parameters, 1) /= parameter_count .or. &
         size(initial_parameters, 2) < 1 .or. &
         (model%combination /= nnfor_combine_mean .and. &
         model%combination /= nnfor_combine_median .and. &
         model%combination /= nnfor_combine_mode) .or. &
         .not. all(ieee_is_finite(initial_parameters))) then
         model%info = 1
         return
      end if
      scaled = nnfor_linscale(model%preprocessing%response, -0.8_dp, 0.8_dp)
      model%response_scale = scaled%scale
      allocate(response_matrix(size(scaled%values), 1))
      response_matrix(:, 1) = scaled%values
      model%hidden_count = hidden_counts(1)
      model%hidden_counts = hidden_counts
      model%repetitions = size(initial_parameters, 2)
      model%lags = lags
      model%series = series
      allocate(model%members(model%repetitions))
      allocate(deepest_all(size(scaled%values), model%repetitions))
      do member = 1, model%repetitions
         model%members(member) = neural_network_fit( &
            model%preprocessing%predictors, response_matrix, hidden_counts, &
            max_iterations, tolerance, decay, initial_parameters(:, member))
         if (model%members(member)%info /= 0) then
            model%info = 2
            return
         end if
         prediction = neural_network_predict(model%members(member), &
            model%preprocessing%predictors)
         deepest_all(:, member) = nnfor_apply_scale(prediction(:, 1), &
            model%response_scale, .true.)
      end do
      restored_all = restore_fitted_values(model%preprocessing, deepest_all)
      model%fitted_all = restored_all
      model%fitted = nnfor_combine(restored_all, model%combination)
      model%residuals = model%preprocessing%original_response - model%fitted
      model%mse = sum(model%residuals**2)/real(size(model%residuals), dp)
      model%extended_preprocessing = .true.
   end function nnfor_mlp_preprocessed_layers_from_initial

   function nnfor_mlp_preprocessed_layers(series, lags, hidden_counts, &
      repetitions, difference_lags, period, seasonal_type, exogenous, &
      exogenous_lags, combination, max_iterations, tolerance, decay, periods, &
      seasonal_types) result(model)
      !! Fit a preprocessed multilayer MLP ensemble from shared random starts.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive lags of the differenced series.
      integer, intent(in) :: hidden_counts(:) !! Units in successive hidden layers.
      integer, intent(in), optional :: repetitions !! Number of ensemble members.
      integer, intent(in), optional :: difference_lags(:) !! Sequential difference lags.
      integer, intent(in), optional :: period !! Scalar seasonal period.
      integer, intent(in), optional :: seasonal_type !! Scalar seasonal input type.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Lags by regressor.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Type for each seasonal period.
      type(nnfor_mlp_model_t) :: model
      type(nnfor_preprocessing_t) :: preprocessing
      real(dp), allocatable :: initial(:, :)
      integer :: reps, parameters, member, parameter

      reps = 20
      if (present(repetitions)) reps = repetitions
      preprocessing = nnfor_preprocess(series, lags, difference_lags, period, &
         seasonal_type, exogenous, exogenous_lags, periods, seasonal_types)
      if (preprocessing%info /= 0 .or. size(hidden_counts) < 1 .or. &
         any(hidden_counts < 1) .or. reps < 1) then
         model%info = merge(preprocessing%info, 1, preprocessing%info /= 0)
         return
      end if
      parameters = neural_network_parameter_count( &
         size(preprocessing%predictors, 2), hidden_counts, 1)
      allocate(initial(parameters, reps))
      do member = 1, reps
         do parameter = 1, parameters
            initial(parameter, member) = &
               0.3_dp*(2.0_dp*random_uniform() - 1.0_dp)
         end do
      end do
      model = nnfor_mlp_preprocessed_layers_from_initial(series, lags, &
         hidden_counts, initial, difference_lags, period, seasonal_type, &
         exogenous, exogenous_lags, combination, max_iterations, tolerance, &
         decay, periods, seasonal_types)
   end function nnfor_mlp_preprocessed_layers

   pure function nnfor_mlp_preprocessed_forecast(model, horizon, exogenous) &
      result(out)
      !! Recursively forecast and inverse-difference a preprocessed MLP ensemble.
      type(nnfor_mlp_model_t), intent(in) :: model !! Fitted preprocessed MLP model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: exogenous(:, :) !! Exogenous history and future values.
      type(nnfor_forecast_t) :: out
      real(dp), allocatable :: future(:, :, :), predictors(:, :)
      real(dp), allocatable :: prediction(:, :), restored(:)
      real(dp) :: base
      integer :: differences, deepest, member, step, lag, variable
      integer :: column, time, stage, source_step, source_index

      if (model%info /= 0 .or. .not. model%extended_preprocessing .or. &
         horizon < 1) then
         out%info = 1
         return
      end if
      if (model%preprocessing%exogenous_count > 0) then
         if (.not. present(exogenous)) then
            out%info = 1
            return
         end if
         if (size(exogenous, 1) < size(model%series) + horizon .or. &
            size(exogenous, 2) /= model%preprocessing%exogenous_count .or. &
            .not. all(ieee_is_finite(exogenous))) then
            out%info = 1
            return
         end if
      end if
      differences = size(model%preprocessing%difference_lags)
      deepest = differences + 1
      allocate(future(horizon, deepest, model%repetitions), source=0.0_dp)
      allocate(out%all_mean(horizon, model%repetitions))
      allocate(predictors(1, size(model%preprocessing%predictors, 2)))
      do member = 1, model%repetitions
         do step = 1, horizon
            column = 0
            do lag = 1, size(model%preprocessing%lags)
               column = column + 1
               source_step = step - model%preprocessing%lags(lag)
               if (source_step <= 0) then
                  source_index = size(model%preprocessing% &
                     difference_levels(deepest)%values) + source_step
                  predictors(1, column) = model%preprocessing% &
                     difference_levels(deepest)%values(source_index)
               else
                  predictors(1, column) = future(source_step, deepest, member)
               end if
            end do
            time = size(model%series) + step
            do variable = 1, model%preprocessing%exogenous_count
               do lag = 1, size(model%preprocessing% &
                  exogenous_lags(variable)%values)
                  column = column + 1
                  predictors(1, column) = exogenous(time - &
                     model%preprocessing%exogenous_lags(variable)%values(lag), &
                     variable)
               end do
            end do
            call append_seasonal_predictors(model%preprocessing, time, &
               predictors, column)
            prediction = neural_network_predict(model%members(member), predictors)
            restored = nnfor_apply_scale(prediction(:, 1), &
               model%response_scale, .true.)
            future(step, deepest, member) = restored(1)
            do stage = differences, 1, -1
               lag = model%preprocessing%difference_lags(stage)
               source_step = step - lag
               if (source_step <= 0) then
                  source_index = size(model%preprocessing% &
                     difference_levels(stage)%values) + source_step
                  base = model%preprocessing%difference_levels(stage)% &
                     values(source_index)
               else
                  base = future(source_step, stage, member)
               end if
               future(step, stage, member) = future(step, stage + 1, member) + &
                  base
            end do
            out%all_mean(step, member) = future(step, 1, member)
         end do
      end do
      out%mean = nnfor_combine(out%all_mean, model%combination)
      out%fitted = model%fitted
      out%residuals = model%residuals
   end function nnfor_mlp_preprocessed_forecast

   pure function nnfor_elm_refit(model, series, exogenous) result(out)
      !! Reuse fixed ELM weights on a new compatible sample and recompute diagnostics.
      type(nnfor_elm_model_t), intent(in) :: model !! Existing ELM model specification.
      real(dp), intent(in) :: series(:) !! New univariate time-series sample.
      real(dp), intent(in), optional :: exogenous(:, :) !! New exogenous-regressor sample.
      type(nnfor_elm_model_t) :: out
      type(nnfor_forecast_t) :: prediction
      real(dp), allocatable :: response(:), predictors(:, :), restored_all(:, :)
      integer :: info

      if (model%info /= 0 .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      out = model
      out%series = series
      if (model%extended_preprocessing) then
         out%preprocessing = nnfor_preprocess(series, model%preprocessing%lags, &
            model%preprocessing%difference_lags, exogenous=exogenous, &
            exogenous_lags=model%preprocessing%exogenous_lags, &
            periods=model%preprocessing%periods, &
            seasonal_types=model%preprocessing%seasonal_types)
         if (out%preprocessing%info /= 0) then
            out%info = out%preprocessing%info
            return
         end if
         response = out%preprocessing%response
         predictors = out%preprocessing%predictors
      else
         call lagged_regression_inputs(series, model%lags, response, predictors, &
            info)
         if (info /= 0) then
            out%info = info
            return
         end if
      end if
      call update_elm_scales(out%network, response, predictors)
      prediction = nnfor_elm_fast_predict(out%network, predictors)
      if (prediction%info /= 0) then
         out%info = 2
         return
      end if
      if (model%extended_preprocessing) then
         restored_all = restore_fitted_values(out%preprocessing, &
            prediction%all_mean)
         out%fitted = nnfor_combine(restored_all, out%network%combination)
         out%residuals = out%preprocessing%original_response - out%fitted
      else
         out%fitted = prediction%mean
         out%residuals = response - out%fitted
      end if
      out%mse = sum(out%residuals**2)/real(size(out%residuals), dp)
      out%info = 0
   end function nnfor_elm_refit

   function nnfor_elm_retrain(model, series, exogenous, lambdas, &
      validation_weight) result(out)
      !! Retrain an ELM ensemble while preserving an existing model specification.
      type(nnfor_elm_model_t), intent(in) :: model !! Existing ELM model specification.
      real(dp), intent(in) :: series(:) !! New univariate time-series sample.
      real(dp), intent(in), optional :: exogenous(:, :) !! New exogenous-regressor sample.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate output penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      type(nnfor_elm_model_t) :: out
      integer, allocatable :: hidden_counts(:)
      integer :: repetitions

      if (model%info /= 0 .or. .not. allocated(model%network%members) .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      repetitions = model%network%repetitions
      if (allocated(model%network%members(1)%hidden_counts)) then
         hidden_counts = model%network%members(1)%hidden_counts
      else
         hidden_counts = [model%network%members(1)%hidden_count]
      end if
      if (model%extended_preprocessing) then
         if (size(hidden_counts) > 1) then
            out = nnfor_elm_preprocessed_layers(series, model%preprocessing%lags, &
               hidden_counts, repetitions, model%preprocessing%difference_lags, &
               exogenous=exogenous, &
               exogenous_lags=model%preprocessing%exogenous_lags, &
               estimator=model%network%estimator, &
               combination=model%network%combination, direct=model%network%direct, &
               lambdas=lambdas, validation_weight=validation_weight, &
               periods=model%preprocessing%periods, &
               seasonal_types=model%preprocessing%seasonal_types, &
               orthogonal=model%network%orthogonal)
         else
            out = nnfor_elm_preprocessed(series, model%preprocessing%lags, &
               hidden_counts(1), repetitions, &
               model%preprocessing%difference_lags, exogenous=exogenous, &
               exogenous_lags=model%preprocessing%exogenous_lags, &
               estimator=model%network%estimator, &
               combination=model%network%combination, direct=model%network%direct, &
               lambdas=lambdas, validation_weight=validation_weight, &
               periods=model%preprocessing%periods, &
               seasonal_types=model%preprocessing%seasonal_types, &
               orthogonal=model%network%orthogonal)
         end if
      else if (size(hidden_counts) > 1) then
         out = nnfor_elm_layers(series, model%lags, hidden_counts, repetitions, &
            model%network%estimator, model%network%combination, &
            model%network%direct, lambdas, validation_weight, &
            model%network%orthogonal)
      else
         out = nnfor_elm(series, model%lags, hidden_counts(1), repetitions, &
            model%network%estimator, model%network%combination, &
            model%network%direct, lambdas, validation_weight, &
            model%network%orthogonal)
      end if
   end function nnfor_elm_retrain

   pure function nnfor_mlp_refit(model, series, exogenous) result(out)
      !! Reuse fixed MLP weights on a new compatible sample and recompute diagnostics.
      type(nnfor_mlp_model_t), intent(in) :: model !! Existing MLP model specification.
      real(dp), intent(in) :: series(:) !! New univariate time-series sample.
      real(dp), intent(in), optional :: exogenous(:, :) !! New exogenous-regressor sample.
      type(nnfor_mlp_model_t) :: out
      type(nnfor_scaled_t) :: scaled
      real(dp), allocatable :: response(:), predictors(:, :)
      real(dp), allocatable :: prediction(:, :), deepest_all(:, :), restored_all(:, :)
      integer :: member, info

      if (model%info /= 0 .or. .not. allocated(model%members) .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      out = model
      out%series = series
      if (model%extended_preprocessing) then
         out%preprocessing = nnfor_preprocess(series, model%preprocessing%lags, &
            model%preprocessing%difference_lags, exogenous=exogenous, &
            exogenous_lags=model%preprocessing%exogenous_lags, &
            periods=model%preprocessing%periods, &
            seasonal_types=model%preprocessing%seasonal_types)
         if (out%preprocessing%info /= 0) then
            out%info = out%preprocessing%info
            return
         end if
         response = out%preprocessing%response
         predictors = out%preprocessing%predictors
      else
         call lagged_regression_inputs(series, model%lags, response, predictors, &
            info)
         if (info /= 0) then
            out%info = info
            return
         end if
      end if
      scaled = nnfor_linscale(response, -0.8_dp, 0.8_dp)
      out%response_scale = scaled%scale
      allocate(deepest_all(size(response), out%repetitions))
      do member = 1, out%repetitions
         call update_neural_network_scales(out%members(member), predictors)
         prediction = neural_network_predict(out%members(member), predictors)
         if (size(prediction, 1) /= size(response)) then
            out%info = 2
            return
         end if
         deepest_all(:, member) = nnfor_apply_scale(prediction(:, 1), &
            out%response_scale, .true.)
      end do
      if (model%extended_preprocessing) then
         restored_all = restore_fitted_values(out%preprocessing, deepest_all)
         out%fitted_all = restored_all
         out%fitted = nnfor_combine(restored_all, out%combination)
         out%residuals = out%preprocessing%original_response - out%fitted
      else
         out%fitted_all = deepest_all
         out%fitted = nnfor_combine(deepest_all, out%combination)
         out%residuals = response - out%fitted
      end if
      out%mse = sum(out%residuals**2)/real(size(out%residuals), dp)
      out%info = 0
   end function nnfor_mlp_refit

   function nnfor_mlp_retrain(model, series, exogenous, max_iterations, &
      tolerance, decay) result(out)
      !! Retrain an MLP ensemble while preserving an existing model specification.
      type(nnfor_mlp_model_t), intent(in) :: model !! Existing MLP model specification.
      real(dp), intent(in) :: series(:) !! New univariate time-series sample.
      real(dp), intent(in), optional :: exogenous(:, :) !! New exogenous-regressor sample.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! Replacement L2 weight decay.
      type(nnfor_mlp_model_t) :: out
      integer, allocatable :: hidden_counts(:)
      real(dp) :: selected_decay

      if (model%info /= 0 .or. .not. allocated(model%members) .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (allocated(model%hidden_counts)) then
         hidden_counts = model%hidden_counts
      else
         hidden_counts = [model%hidden_count]
      end if
      selected_decay = model%members(1)%decay
      if (present(decay)) selected_decay = decay
      if (model%extended_preprocessing) then
         if (size(hidden_counts) > 1) then
            out = nnfor_mlp_preprocessed_layers(series, &
               model%preprocessing%lags, hidden_counts, model%repetitions, &
               model%preprocessing%difference_lags, exogenous=exogenous, &
               exogenous_lags=model%preprocessing%exogenous_lags, &
               combination=model%combination, max_iterations=max_iterations, &
               tolerance=tolerance, decay=selected_decay, &
               periods=model%preprocessing%periods, &
               seasonal_types=model%preprocessing%seasonal_types)
         else
            out = nnfor_mlp_preprocessed(series, model%preprocessing%lags, &
               hidden_counts(1), model%repetitions, &
               model%preprocessing%difference_lags, exogenous=exogenous, &
               exogenous_lags=model%preprocessing%exogenous_lags, &
               combination=model%combination, max_iterations=max_iterations, &
               tolerance=tolerance, decay=selected_decay, &
               periods=model%preprocessing%periods, &
               seasonal_types=model%preprocessing%seasonal_types)
         end if
      else if (size(hidden_counts) > 1) then
         out = nnfor_mlp_layers(series, model%lags, hidden_counts, &
            model%repetitions, model%combination, max_iterations, tolerance, &
            selected_decay)
      else
         out = nnfor_mlp(series, model%lags, hidden_counts(1), &
            model%repetitions, model%combination, max_iterations, tolerance, &
            selected_decay)
      end if
   end function nnfor_mlp_retrain

   pure function nnfor_select_elm_hidden_from_weights(response, predictors, &
      input_weights, alpha) result(out)
      !! Select ELM hidden size from significant stepwise output coefficients.
      real(dp), intent(in) :: response(:) !! Network response.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-predictor matrix.
      real(dp), intent(in) :: input_weights(:, :, :) !! Bias and input weights by repetition.
      real(dp), intent(in), optional :: alpha !! Coefficient significance threshold.
      type(nnfor_elm_hidden_selection_t) :: out
      type(nnfor_output_fit_t) :: fit
      type(significance_result_t) :: significance
      real(dp), allocatable :: augmented(:, :), hidden(:, :), design(:, :)
      real(dp), allocatable :: counts(:)
      integer :: member

      if (present(alpha)) out%alpha = alpha
      if (size(response) /= size(predictors, 1) .or. size(response) < 4 .or. &
         size(predictors, 2) < 1 .or. &
         size(input_weights, 1) /= size(predictors, 2) + 1 .or. &
         size(input_weights, 2) < 1 .or. size(input_weights, 3) < 1 .or. &
         size(response) <= size(input_weights, 2) + 1 .or. &
         out%alpha <= 0.0_dp .or. out%alpha >= 1.0_dp .or. &
         .not. all(ieee_is_finite(response)) .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         .not. all(ieee_is_finite(input_weights))) then
         out%info = 1
         return
      end if
      out%candidate_hidden = size(input_weights, 2)
      out%repetitions = size(input_weights, 3)
      allocate(out%significant_count(out%repetitions))
      allocate(augmented(size(response), size(predictors, 2) + 1))
      augmented(:, 1) = 1.0_dp
      augmented(:, 2:) = predictors
      allocate(design(size(response), out%candidate_hidden + 1))
      design(:, 1) = 1.0_dp
      do member = 1, out%repetitions
         hidden = logistic_sigmoid(matmul(augmented, &
            input_weights(:, :, member)))
         design(:, 2:) = hidden
         fit = nnfor_stepwise_output_fit(design, response)
         if (fit%info /= 0) then
            out%info = 2
            return
         end if
         significance = selected_coefficient_significance(design, response, &
            fit%active, fit%coefficients)
         if (significance%info /= 0) then
            out%info = 2
            return
         end if
         out%significant_count(member) = count( &
            significance%p_value(2:) < out%alpha)
      end do
      counts = real(out%significant_count, dp)
      out%selected = max(1, nint(median(counts)))
   end function nnfor_select_elm_hidden_from_weights

   function nnfor_select_elm_hidden(response, predictors, repetitions, &
      maximum_hidden, alpha) result(out)
      !! Select ELM hidden size using shared-RNG candidate hidden weights.
      real(dp), intent(in) :: response(:) !! Network response.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-predictor matrix.
      integer, intent(in), optional :: repetitions !! Candidate ELM repetitions.
      integer, intent(in), optional :: maximum_hidden !! Candidate hidden-layer size.
      real(dp), intent(in), optional :: alpha !! Coefficient significance threshold.
      type(nnfor_elm_hidden_selection_t) :: out
      real(dp), allocatable :: input_weights(:, :, :)
      integer :: reps, hidden, row, column, member

      reps = 20
      if (present(repetitions)) reps = repetitions
      hidden = min(size(predictors, 2) + 2, size(response) - 2)
      if (present(maximum_hidden)) hidden = maximum_hidden
      if (size(response) /= size(predictors, 1) .or. &
         size(predictors, 2) < 1 .or. reps < 1 .or. hidden < 1 .or. &
         size(response) <= hidden + 1) then
         out%info = 1
         return
      end if
      allocate(input_weights(size(predictors, 2) + 1, hidden, reps))
      do member = 1, reps
         do column = 1, hidden
            do row = 1, size(input_weights, 1)
               input_weights(row, column, member) = &
                  4.0_dp*random_uniform() - 2.0_dp
            end do
         end do
      end do
      out = nnfor_select_elm_hidden_from_weights(response, predictors, &
         input_weights, alpha)
   end function nnfor_select_elm_hidden

   pure function nnfor_select_hidden_count(response, predictors, &
      maximum_hidden, repetitions, validation_weight, max_iterations, &
      tolerance, decay, combination) result(out)
      !! Select one hidden-layer size using a terminal validation sample.
      real(dp), intent(in) :: response(:) !! Network response.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-predictor matrix.
      integer, intent(in), optional :: maximum_hidden !! Largest hidden size evaluated.
      integer, intent(in), optional :: repetitions !! Networks fitted at each size.
      real(dp), intent(in), optional :: validation_weight !! Fraction used for training.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      type(nnfor_hidden_selection_t) :: out
      real(dp) :: weight
      integer, allocatable :: fold_ids(:)
      integer :: training

      weight = 0.8_dp
      if (present(validation_weight)) weight = validation_weight
      training = nint(weight*real(size(response), dp))
      if (training < 3 .or. training >= size(response)) then
         out%info = 1
         return
      end if
      allocate(fold_ids(size(response)), source=0)
      fold_ids(training + 1:) = 1
      out = nnfor_select_hidden_count_folds(response, predictors, fold_ids, &
         maximum_hidden, repetitions, max_iterations, tolerance, decay, &
         combination)
   end function nnfor_select_hidden_count

   pure function nnfor_select_hidden_count_folds(response, predictors, &
      fold_ids, maximum_hidden, repetitions, max_iterations, tolerance, &
      decay, combination) result(out)
      !! Select hidden size from supplied holdout or cross-validation folds.
      real(dp), intent(in) :: response(:) !! Network response.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-predictor matrix.
      integer, intent(in) :: fold_ids(:) !! Zero for fixed training rows or positive fold index.
      integer, intent(in), optional :: maximum_hidden !! Largest hidden size evaluated.
      integer, intent(in), optional :: repetitions !! Networks fitted for each fold and size.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      type(nnfor_hidden_selection_t) :: out
      type(neural_network_t) :: network
      real(dp), allocatable :: training_response(:, :), training_predictors(:, :)
      real(dp), allocatable :: test_response(:), test_predictors(:, :)
      real(dp), allocatable :: initial(:), member_predictions(:, :)
      real(dp), allocatable :: prediction(:, :), combined(:), fold_error(:)
      integer :: limit, reps, method, folds, fold, hidden, member, parameter
      integer :: parameters, successful, training_count, test_count

      reps = 3
      if (present(repetitions)) reps = repetitions
      limit = min(size(predictors, 2) + 2, size(response) - 3)
      if (present(maximum_hidden)) limit = maximum_hidden
      method = nnfor_combine_median
      if (present(combination)) method = combination
      folds = maxval(fold_ids)
      if (size(response) /= size(predictors, 1) .or. &
         size(response) /= size(fold_ids) .or. size(response) < 6 .or. &
         size(predictors, 2) < 1 .or. limit < 1 .or. reps < 1 .or. &
         folds < 1 .or. minval(fold_ids) < 0 .or. &
         (method < nnfor_combine_mean .or. method > nnfor_combine_mode) .or. &
         .not. all(ieee_is_finite(response)) .or. &
         .not. all(ieee_is_finite(predictors))) then
         out%info = 1
         return
      end if
      do fold = 1, folds
         if (count(fold_ids == fold) < 1 .or. count(fold_ids /= fold) < 3) then
            out%info = 1
            return
         end if
      end do
      out%maximum = limit
      out%repetitions = reps
      allocate(out%mse(limit), source=huge(1.0_dp))
      allocate(fold_error(folds))
      do hidden = 1, limit
         fold_error = huge(1.0_dp)
         do fold = 1, folds
            training_count = count(fold_ids /= fold)
            test_count = count(fold_ids == fold)
            allocate(training_response(training_count, 1))
            allocate(training_predictors(training_count, size(predictors, 2)))
            allocate(test_response(test_count))
            allocate(test_predictors(test_count, size(predictors, 2)))
            training_response(:, 1) = pack(response, fold_ids /= fold)
            test_response = pack(response, fold_ids == fold)
            call pack_predictor_rows(predictors, fold_ids /= fold, &
               training_predictors)
            call pack_predictor_rows(predictors, fold_ids == fold, &
               test_predictors)
            parameters = neural_network_parameter_count(size(predictors, 2), &
               hidden, 1)
            allocate(member_predictions(test_count, reps), source=0.0_dp)
            successful = 0
            do member = 1, reps
               allocate(initial(parameters))
               do parameter = 1, parameters
                  initial(parameter) = 0.2_dp*sin(0.6180339887498948_dp* &
                     real(parameter + 11*member + 17*hidden + 23*fold, dp))
               end do
               network = neural_network_fit(training_predictors, &
                  training_response, hidden, max_iterations, tolerance, decay, &
                  initial)
               deallocate(initial)
               if (network%info /= 0) cycle
               prediction = neural_network_predict(network, test_predictors)
               successful = successful + 1
               member_predictions(:, successful) = prediction(:, 1)
            end do
            if (successful > 0) then
               combined = nnfor_combine( &
                  member_predictions(:, :successful), method)
               fold_error(fold) = sum((test_response - combined)**2)/ &
                  real(test_count, dp)
            end if
            deallocate(training_response, training_predictors, test_response)
            deallocate(test_predictors, member_predictions)
         end do
         if (all(fold_error < huge(1.0_dp))) &
            out%mse(hidden) = sum(fold_error)/real(folds, dp)
      end do
      if (minval(out%mse) >= huge(1.0_dp)) then
         out%info = 2
         return
      end if
      out%selected = minloc(out%mse, dim=1)
   end function nnfor_select_hidden_count_folds

   function nnfor_select_hidden_count_random(response, predictors, method, &
      maximum_hidden, repetitions, fold_count, validation_fraction, &
      max_iterations, tolerance, decay, combination) result(out)
      !! Select hidden size using randomized holdout or cross-validation rows.
      real(dp), intent(in) :: response(:) !! Network response.
      real(dp), intent(in) :: predictors(:, :) !! Observation-by-predictor matrix.
      integer, intent(in), optional :: method !! Holdout or cross-validation selector.
      integer, intent(in), optional :: maximum_hidden !! Largest hidden size evaluated.
      integer, intent(in), optional :: repetitions !! Networks fitted for each fold and size.
      integer, intent(in), optional :: fold_count !! Number of cross-validation folds.
      real(dp), intent(in), optional :: validation_fraction !! Random holdout fraction.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      type(nnfor_hidden_selection_t) :: out
      integer, allocatable :: order(:), fold_ids(:)
      integer :: selection_method, folds, validation, i, j, temporary
      real(dp) :: fraction

      selection_method = nnfor_hidden_cross_validation
      if (present(method)) selection_method = method
      folds = 5
      if (present(fold_count)) folds = fold_count
      fraction = 0.2_dp
      if (present(validation_fraction)) fraction = validation_fraction
      if (size(response) /= size(predictors, 1) .or. size(response) < 6 .or. &
         (selection_method /= nnfor_hidden_holdout .and. &
         selection_method /= nnfor_hidden_cross_validation) .or. &
         (selection_method == nnfor_hidden_cross_validation .and. &
         (folds < 2 .or. folds >= size(response))) .or. &
         (selection_method == nnfor_hidden_holdout .and. &
         (fraction <= 0.0_dp .or. fraction >= 1.0_dp))) then
         out%info = 1
         return
      end if
      allocate(order(size(response)), fold_ids(size(response)), source=0)
      order = [(i, i=1, size(response))]
      do i = size(order), 2, -1
         j = 1 + int(random_uniform()*real(i, dp))
         temporary = order(i)
         order(i) = order(j)
         order(j) = temporary
      end do
      if (selection_method == nnfor_hidden_holdout) then
         validation = max(1, nint(fraction*real(size(response), dp)))
         if (size(response) - validation < 3) then
            out%info = 1
            return
         end if
         fold_ids(order(:validation)) = 1
      else
         do i = 1, size(order)
            fold_ids(order(i)) = 1 + mod(i - 1, folds)
         end do
      end if
      out = nnfor_select_hidden_count_folds(response, predictors, fold_ids, &
         maximum_hidden, repetitions, max_iterations, tolerance, decay, &
         combination)
   end function nnfor_select_hidden_count_random

   pure function nnfor_season_check(series, period, alpha, multiplicative, &
      moving_average) result(out)
      !! Test detrended seasonal positions using nnfor's Friedman diagnostic.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), intent(in), optional :: alpha !! Test significance level.
      logical, intent(in), optional :: multiplicative !! Request ratio detrending.
      real(dp), intent(in), optional :: moving_average(:) !! Supplied level or trend.
      type(nnfor_seasonality_t) :: out
      real(dp), allocatable :: detrended(:), rank_sums(:), row(:)
      real(dp) :: significance, tie_total, tie_count, denominator
      integer :: cycle, season, other, first, occurrences
      logical :: use_multiplicative

      significance = 0.05_dp
      if (present(alpha)) significance = alpha
      use_multiplicative = .true.
      if (present(multiplicative)) use_multiplicative = multiplicative
      out%period = period
      if (size(series) < 2*period .or. period < 2 .or. &
         significance <= 0.0_dp .or. significance >= 1.0_dp .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (present(moving_average)) then
         if (size(moving_average) /= size(series)) then
            out%info = 1
            return
         end if
      end if
      if (present(moving_average)) then
         if (.not. all(ieee_is_finite(moving_average))) then
            out%info = 1
            return
         end if
      end if
      if (present(moving_average)) then
         out%moving_average = moving_average
      else
         out%moving_average = centered_moving_average_fill(series, period)
      end if
      if (.not. all(ieee_is_finite(out%moving_average))) then
         out%info = 1
         return
      end if
      use_multiplicative = use_multiplicative .and. minval(series) > 0.0_dp &
         .and. minval(abs(out%moving_average)) > tiny(1.0_dp)
      out%multiplicative = use_multiplicative
      if (use_multiplicative) then
         detrended = series/out%moving_average
      else
         detrended = series - out%moving_average
      end if
      out%cycles = size(series)/period
      allocate(rank_sums(period), source=0.0_dp)
      allocate(row(period))
      tie_total = 0.0_dp
      do cycle = 1, out%cycles
         first = 1 + (cycle - 1)*period
         row = detrended(first:first + period - 1)
         do season = 1, period
            occurrences = count(row == row(season))
            rank_sums(season) = rank_sums(season) + 1.0_dp + &
               real(count(row < row(season)), dp) + &
               0.5_dp*real(occurrences - 1, dp)
            if (occurrences > 1) then
               if (season == 1) then
                  other = 0
               else
                  other = count(row(:season - 1) == row(season))
               end if
               if (other == 0) then
                  tie_count = real(occurrences, dp)
                  tie_total = tie_total + tie_count**3 - tie_count
               end if
            end if
         end do
      end do
      out%statistic = 12.0_dp*sum(rank_sums**2)/ &
         (real(out%cycles, dp)*real(period, dp)*real(period + 1, dp)) - &
         3.0_dp*real(out%cycles, dp)*real(period + 1, dp)
      denominator = real(out%cycles, dp)*(real(period, dp)**3 - &
         real(period, dp))
      if (denominator > 0.0_dp .and. tie_total < denominator) &
         out%statistic = out%statistic/(1.0_dp - tie_total/denominator)
      out%p_value = regularized_gamma_q(0.5_dp*real(period - 1, dp), &
         0.5_dp*out%statistic)
      out%seasonal = out%p_value <= significance
   end function nnfor_season_check

   pure function nnfor_mseason_test(series, period, moving_average, &
      seasonal_count, alpha) result(out)
      !! Test whether seasonal magnitude increases with the series level.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      integer, intent(in) :: period !! Seasonal period of at least two.
      real(dp), intent(in), optional :: moving_average(:) !! Supplied level or trend.
      integer, intent(in), optional :: seasonal_count !! Strong seasonal positions to test.
      real(dp), intent(in), optional :: alpha !! Decision significance level.
      type(nnfor_mseason_t) :: out
      real(dp), allocatable :: level(:), residual(:), magnitudes(:)
      real(dp), allocatable :: correlations(:), p_values(:), x(:), y(:)
      logical, allocatable :: selected(:)
      real(dp) :: direction
      integer :: count_seasons, season, rank, position, observations

      out%period = period
      out%alpha = 0.05_dp
      if (present(alpha)) out%alpha = alpha
      count_seasons = 1
      if (present(seasonal_count)) count_seasons = min(period, seasonal_count)
      out%seasonal_count = count_seasons
      if (period < 2 .or. count_seasons < 1 .or. size(series) < 2*period .or. &
         out%alpha <= 0.0_dp .or. out%alpha >= 1.0_dp .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (present(moving_average)) then
         if (size(moving_average) /= size(series) .or. &
            .not. all(ieee_is_finite(moving_average))) then
            out%info = 1
            return
         end if
         level = moving_average
      else
         level = centered_moving_average_fill(series, period)
      end if
      residual = series - level
      allocate(magnitudes(period))
      allocate(selected(period), source=.false.)
      do season = 1, period
         observations = 1 + (size(series) - season)/period
         magnitudes(season) = sum(residual(season::period))/ &
            real(observations, dp)
      end do
      allocate(correlations(count_seasons), p_values(count_seasons))
      do rank = 1, count_seasons
         position = 0
         do season = 1, period
            if (selected(season)) cycle
            if (position == 0) then
               position = season
            else if (abs(magnitudes(season)) > &
               abs(magnitudes(position))) then
               position = season
            end if
         end do
         selected(position) = .true.
         observations = 1 + (size(series) - position)/period
         if (observations < 3) then
            out%info = 2
            return
         end if
         direction = merge(-1.0_dp, 1.0_dp, magnitudes(position) < 0.0_dp)
         allocate(x(observations), y(observations))
         x = direction*residual(position::period)
         y = level(position::period)
         call pearson_greater_test(x, y, correlations(rank), p_values(rank))
         deallocate(x, y)
      end do
      out%statistic = median(correlations)
      out%p_value = median(p_values)
      if (out%statistic <= 0.0_dp) out%p_value = 1.0_dp
      out%multiplicative = out%p_value <= out%alpha
   end function nnfor_mseason_test

   pure function nnfor_trend_check(series) result(out)
      !! Compare level-only and additive-trend exponential-smoothing AICs.
      real(dp), intent(in) :: series(:) !! Level or trend observations.
      type(nnfor_trend_t) :: out
      real(dp) :: alpha, beta, level, trend, old_level, residual, rss
      real(dp) :: best_level_rss, best_trend_rss, observations
      integer :: alpha_index, beta_index, observation

      if (size(series) < 4 .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      best_level_rss = huge(1.0_dp)
      best_trend_rss = huge(1.0_dp)
      do alpha_index = 1, 19
         alpha = 0.05_dp*real(alpha_index, dp)
         level = series(1)
         rss = 0.0_dp
         do observation = 2, size(series)
            residual = series(observation) - level
            rss = rss + residual**2
            level = alpha*series(observation) + (1.0_dp - alpha)*level
         end do
         best_level_rss = min(best_level_rss, rss)
         do beta_index = 1, 19
            beta = 0.05_dp*real(beta_index, dp)
            level = series(1)
            trend = series(2) - series(1)
            rss = 0.0_dp
            do observation = 2, size(series)
               residual = series(observation) - level - trend
               rss = rss + residual**2
               old_level = level
               level = alpha*series(observation) + &
                  (1.0_dp - alpha)*(level + trend)
               trend = beta*(level - old_level) + (1.0_dp - beta)*trend
            end do
            best_trend_rss = min(best_trend_rss, rss)
         end do
      end do
      observations = real(size(series) - 1, dp)
      out%level_aic = observations*log(max(best_level_rss/observations, &
         tiny(1.0_dp))) + 4.0_dp
      out%trend_aic = observations*log(max(best_trend_rss/observations, &
         tiny(1.0_dp))) + 8.0_dp
      out%trending = out%trend_aic < out%level_aic
   end function nnfor_trend_check

   pure function nnfor_canova_hansen(series, period, seasonal_type, alpha, &
      newey_west_order) result(out)
      !! Compute the joint Canova-Hansen seasonal-stability statistic.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      integer, intent(in) :: period !! Seasonal period of at least two.
      integer, intent(in), optional :: seasonal_type !! Trigonometric or dummy formulation.
      real(dp), intent(in), optional :: alpha !! Decision significance in [0.01, 0.10].
      integer, intent(in), optional :: newey_west_order !! Nonnegative covariance bandwidth.
      type(nnfor_canova_hansen_t) :: out
      real(dp), allocatable :: seasonal(:, :), design(:, :), cross_product(:, :)
      real(dp), allocatable :: inverse(:, :), coefficients(:), residuals(:)
      real(dp), allocatable :: scaled_residuals(:, :), omega(:, :)
      real(dp), allocatable :: lag_cross(:, :), cumulative(:, :), cumulative_cross(:, :)
      real(dp), allocatable :: product(:, :)
      real(dp) :: pi_value, angle, weight
      integer :: observations, components, cycles, harmonic, column
      integer :: time, season, lag, status, index

      observations = size(series)
      out%period = period
      out%seasonal_type = nnfor_ch_trigonometric
      if (present(seasonal_type)) out%seasonal_type = seasonal_type
      out%alpha = 0.05_dp
      if (present(alpha)) out%alpha = max(0.01_dp, min(0.10_dp, alpha))
      if (period < 2 .or. observations <= period .or. &
         .not. all(ieee_is_finite(series)) .or. &
         (out%seasonal_type /= nnfor_ch_trigonometric .and. &
         out%seasonal_type /= nnfor_ch_dummy)) then
         out%info = 1
         return
      end if
      if (out%seasonal_type == nnfor_ch_dummy .and. period > 12) then
         out%info = 2
         return
      end if
      out%newey_west_order = nint(real(period, dp)* &
         (real(observations, dp)/100.0_dp)**0.25_dp)
      if (present(newey_west_order)) &
         out%newey_west_order = newey_west_order
      if (out%newey_west_order < 0 .or. &
         out%newey_west_order >= observations) then
         out%info = 1
         return
      end if
      pi_value = acos(-1.0_dp)
      if (out%seasonal_type == nnfor_ch_dummy) then
         components = period
         allocate(seasonal(observations, components), source=0.0_dp)
         do time = 1, observations
            season = 1 + mod(time - 1, period)
            seasonal(time, season) = 1.0_dp
         end do
         design = seasonal
      else
         components = period - 1
         allocate(seasonal(observations, components))
         column = 0
         cycles = period/2 - merge(1, 0, mod(period, 2) == 0)
         do harmonic = 1, cycles
            do time = 1, observations
               angle = 2.0_dp*pi_value*real(harmonic*time, dp)/ &
                  real(period, dp)
               seasonal(time, column + 1) = cos(angle)
               seasonal(time, column + 2) = sin(angle)
            end do
            column = column + 2
         end do
         if (mod(period, 2) == 0) then
            do time = 1, observations
               seasonal(time, components) = merge(1.0_dp, -1.0_dp, &
                  mod(time, 2) == 0)
            end do
         end if
         allocate(design(observations, components + 1))
         design(:, 1) = 1.0_dp
         design(:, 2:) = seasonal
      end if
      allocate(cross_product(size(design, 2), size(design, 2)))
      allocate(inverse(size(design, 2), size(design, 2)))
      cross_product = matmul(transpose(design), design)
      call symmetric_pseudoinverse(cross_product, inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      coefficients = matmul(inverse, matmul(transpose(design), series))
      residuals = series - matmul(design, coefficients)
      scaled_residuals = seasonal*spread(residuals, 2, components)
      omega = matmul(transpose(scaled_residuals), scaled_residuals)
      do lag = 1, out%newey_west_order
         weight = 1.0_dp - real(lag, dp)/ &
            real(out%newey_west_order + 1, dp)
         lag_cross = matmul(transpose(scaled_residuals(:observations - lag, :)), &
            scaled_residuals(lag + 1:, :))
         omega = omega + weight*(lag_cross + transpose(lag_cross))
      end do
      omega = omega/real(observations, dp)
      deallocate(inverse)
      allocate(inverse(components, components))
      call symmetric_pseudoinverse(omega, inverse, status)
      if (status /= 0) then
         out%info = 4
         return
      end if
      allocate(cumulative(observations, components))
      cumulative(1, :) = scaled_residuals(1, :)
      do time = 2, observations
         cumulative(time, :) = cumulative(time - 1, :) + &
            scaled_residuals(time, :)
      end do
      cumulative_cross = matmul(transpose(cumulative), cumulative)
      product = matmul(inverse, cumulative_cross)
      out%statistic = 0.0_dp
      do index = 1, components
         out%statistic = out%statistic + product(index, index)
      end do
      out%statistic = out%statistic/real(observations, dp)**2
      out%degrees_of_freedom = components
      out%critical_value = canova_hansen_critical_value(period, &
         out%seasonal_type, observations - size(design, 2), out%alpha)
      out%p_value = canova_hansen_p_value(out%statistic, period, &
         out%seasonal_type, observations - size(design, 2))
      out%difference_required = out%statistic > out%critical_value
   end function nnfor_canova_hansen

   pure function nnfor_select_differences(series, period, alpha, &
      multiplicative) result(out)
      !! Select ordinary and seasonal difference lags from nnfor diagnostics.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in), optional :: period !! Seasonal period; one disables seasonality.
      real(dp), intent(in), optional :: alpha !! Seasonality-test significance level.
      logical, intent(in), optional :: multiplicative !! Request ratio detrending.
      type(nnfor_difference_selection_t) :: out
      real(dp), allocatable :: detrended(:)
      integer :: seasonal_period, count_lags

      seasonal_period = 1
      if (present(period)) seasonal_period = period
      if (size(series) < 4 .or. seasonal_period < 1 .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (seasonal_period > 1) then
         out%seasonality = nnfor_season_check(series, seasonal_period, alpha, &
            multiplicative)
         if (out%seasonality%info /= 0) then
            out%info = 1
            return
         end if
         out%trend = nnfor_trend_check(out%seasonality%moving_average)
         if (out%seasonality%seasonal) then
            out%multiplicative_seasonality%period = seasonal_period
            out%multiplicative_seasonality%alpha = 0.05_dp
            if (present(alpha)) out%multiplicative_seasonality%alpha = alpha
            if (present(multiplicative)) then
               out%multiplicative_seasonality%multiplicative = &
                  out%seasonality%multiplicative
            else if (size(series)/seasonal_period < 3) then
               out%multiplicative_seasonality%multiplicative = .true.
               out%seasonality%multiplicative = &
                  minval(abs(out%seasonality%moving_average)) > tiny(1.0_dp)
            else
               out%multiplicative_seasonality = nnfor_mseason_test(series, &
                  seasonal_period, out%seasonality%moving_average, alpha=alpha)
               if (out%multiplicative_seasonality%info /= 0) then
                  out%info = 1
                  return
               end if
               out%seasonality%multiplicative = &
                  out%multiplicative_seasonality%multiplicative .and. &
                  minval(abs(out%seasonality%moving_average)) > tiny(1.0_dp)
            end if
            if (out%seasonality%multiplicative) then
               detrended = series/out%seasonality%moving_average
            else
               detrended = series - out%seasonality%moving_average
            end if
            out%canova_hansen = nnfor_canova_hansen(detrended, &
               seasonal_period, nnfor_ch_trigonometric, alpha)
            if (out%canova_hansen%info /= 0) then
               out%info = 1
               return
            end if
         else
            out%canova_hansen%period = seasonal_period
            out%canova_hansen%alpha = 0.05_dp
            if (present(alpha)) out%canova_hansen%alpha = &
               max(0.01_dp, min(0.10_dp, alpha))
         end if
      else
         out%trend = nnfor_trend_check(series)
      end if
      if (out%trend%info /= 0) then
         out%info = 1
         return
      end if
      count_lags = 0
      if (out%trend%trending) count_lags = count_lags + 1
      if (seasonal_period > 1 .and. &
         out%canova_hansen%difference_required) &
         count_lags = count_lags + 1
      allocate(out%difference_lags(count_lags))
      count_lags = 0
      if (out%trend%trending) then
         count_lags = count_lags + 1
         out%difference_lags(count_lags) = 1
      end if
      if (seasonal_period > 1 .and. &
         out%canova_hansen%difference_required) then
         count_lags = count_lags + 1
         out%difference_lags(count_lags) = seasonal_period
      end if
   end function nnfor_select_differences

   function nnfor_elm_auto(series, period, candidate_lags, keep, &
      difference_lags, seasonal_type, exogenous, exogenous_lags, repetitions, &
      selection_repetitions, estimator, combination, direct, lambdas, &
      validation_weight, alpha, periods, seasonal_types) result(out)
      !! Automatically specify preprocessing, lags, hidden size, and an ELM fit.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in), optional :: period !! Seasonal period.
      integer, intent(in), optional :: candidate_lags(:) !! Candidate response lags.
      logical, intent(in), optional :: keep(:) !! Candidate lags forced into the model.
      integer, intent(in), optional :: difference_lags(:) !! Supplied sequential differences.
      integer, intent(in), optional :: seasonal_type !! Supplied deterministic seasonal type.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time exogenous regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Lags by exogenous regressor.
      integer, intent(in), optional :: repetitions !! Final ELM ensemble size.
      integer, intent(in), optional :: selection_repetitions !! Hidden-selection repetitions.
      integer, intent(in), optional :: estimator !! ELM output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge or lasso penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      real(dp), intent(in), optional :: alpha !! Automatic-test significance level.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Supplied input types by period.
      type(nnfor_elm_auto_t) :: out
      type(automatic_specification_t) :: specification
      integer :: reps, selection_reps

      reps = 20
      if (present(repetitions)) reps = repetitions
      selection_reps = 20
      if (present(selection_repetitions)) &
         selection_reps = selection_repetitions
      specification = build_automatic_specification(series, period, &
         candidate_lags, keep, difference_lags, seasonal_type, exogenous, &
         exogenous_lags, alpha, periods, seasonal_types)
      call copy_automatic_specification(specification, out%candidate_lags, &
         out%difference_lags, out%difference_selection, out%lag_selection, &
         out%seasonal_type, out%automatic_differences, out%periods, &
         out%seasonal_types)
      if (specification%info /= 0 .or. reps < 1 .or. selection_reps < 1) then
         out%info = merge(specification%info, 1, specification%info /= 0)
         return
      end if
      out%hidden_selection = nnfor_select_elm_hidden( &
         specification%preprocessing%response, &
         specification%preprocessing%predictors, selection_reps, &
         alpha=alpha)
      if (out%hidden_selection%info /= 0) then
         out%info = 10 + out%hidden_selection%info
         return
      end if
      out%model = nnfor_elm_preprocessed(series, &
         specification%selected_lags, out%hidden_selection%selected, reps, &
         specification%difference_lags, specification%period, &
         specification%seasonal_type, exogenous, exogenous_lags, estimator, &
         combination, direct, lambdas, validation_weight, &
         specification%periods, specification%seasonal_types)
      if (out%model%info /= 0) out%info = 20 + out%model%info
   end function nnfor_elm_auto

   pure function nnfor_elm_auto_forecast(fit, horizon, exogenous) result(out)
      !! Forecast from an automatically specified ELM fit.
      type(nnfor_elm_auto_t), intent(in) :: fit !! Automatic ELM result.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: exogenous(:, :) !! Exogenous history and future values.
      type(nnfor_forecast_t) :: out

      if (fit%info /= 0) then
         out%info = 1
         return
      end if
      out = nnfor_elm_preprocessed_forecast(fit%model, horizon, exogenous)
   end function nnfor_elm_auto_forecast

   function nnfor_mlp_auto(series, period, candidate_lags, keep, &
      difference_lags, seasonal_type, exogenous, exogenous_lags, repetitions, &
      hidden_method, maximum_hidden, selection_repetitions, fold_count, &
      validation_fraction, combination, max_iterations, tolerance, decay, &
      alpha, periods, seasonal_types) result(out)
      !! Automatically specify preprocessing, lags, hidden size, and an MLP fit.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in), optional :: period !! Seasonal period.
      integer, intent(in), optional :: candidate_lags(:) !! Candidate response lags.
      logical, intent(in), optional :: keep(:) !! Candidate lags forced into the model.
      integer, intent(in), optional :: difference_lags(:) !! Supplied sequential differences.
      integer, intent(in), optional :: seasonal_type !! Supplied deterministic seasonal type.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time exogenous regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Lags by exogenous regressor.
      integer, intent(in), optional :: repetitions !! Final MLP ensemble size.
      integer, intent(in), optional :: hidden_method !! Terminal, holdout, or cross-validation selector.
      integer, intent(in), optional :: maximum_hidden !! Largest hidden size evaluated.
      integer, intent(in), optional :: selection_repetitions !! Networks fitted per candidate size.
      integer, intent(in), optional :: fold_count !! Number of cross-validation folds.
      real(dp), intent(in), optional :: validation_fraction !! Holdout fraction.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      real(dp), intent(in), optional :: alpha !! Automatic-test significance level.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Supplied input types by period.
      type(nnfor_mlp_auto_t) :: out
      type(automatic_specification_t) :: specification
      integer :: reps, selection_reps, method
      real(dp) :: fraction

      reps = 20
      if (present(repetitions)) reps = repetitions
      selection_reps = 3
      if (present(selection_repetitions)) &
         selection_reps = selection_repetitions
      method = nnfor_hidden_cross_validation
      if (present(hidden_method)) method = hidden_method
      fraction = 0.2_dp
      if (present(validation_fraction)) fraction = validation_fraction
      specification = build_automatic_specification(series, period, &
         candidate_lags, keep, difference_lags, seasonal_type, exogenous, &
         exogenous_lags, alpha, periods, seasonal_types)
      call copy_automatic_specification(specification, out%candidate_lags, &
         out%difference_lags, out%difference_selection, out%lag_selection, &
         out%seasonal_type, out%automatic_differences, out%periods, &
         out%seasonal_types)
      if (specification%info /= 0 .or. reps < 1 .or. selection_reps < 1) then
         out%info = merge(specification%info, 1, specification%info /= 0)
         return
      end if
      if (method == nnfor_hidden_terminal) then
         out%hidden_selection = nnfor_select_hidden_count( &
            specification%preprocessing%response, &
            specification%preprocessing%predictors, maximum_hidden, &
            selection_reps, 1.0_dp - fraction, max_iterations, tolerance, &
            decay, combination)
      else
         out%hidden_selection = nnfor_select_hidden_count_random( &
            specification%preprocessing%response, &
            specification%preprocessing%predictors, method, maximum_hidden, &
            selection_reps, fold_count, fraction, max_iterations, tolerance, &
            decay, combination)
      end if
      if (out%hidden_selection%info /= 0) then
         out%info = 10 + out%hidden_selection%info
         return
      end if
      out%model = nnfor_mlp_preprocessed(series, &
         specification%selected_lags, out%hidden_selection%selected, reps, &
         specification%difference_lags, specification%period, &
         specification%seasonal_type, exogenous, exogenous_lags, combination, &
         max_iterations, tolerance, decay, specification%periods, &
         specification%seasonal_types)
      if (out%model%info /= 0) out%info = 20 + out%model%info
   end function nnfor_mlp_auto

   pure function nnfor_mlp_auto_forecast(fit, horizon, exogenous) result(out)
      !! Forecast from an automatically specified MLP fit.
      type(nnfor_mlp_auto_t), intent(in) :: fit !! Automatic MLP result.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: exogenous(:, :) !! Exogenous history and future values.
      type(nnfor_forecast_t) :: out

      if (fit%info /= 0) then
         out%info = 1
         return
      end if
      out = nnfor_mlp_preprocessed_forecast(fit%model, horizon, exogenous)
   end function nnfor_mlp_auto_forecast

   function nnfor_elm_thief(series, horizon, period, candidate_lags, keep, &
      difference_lags, seasonal_type, exogenous, exogenous_lags, repetitions, &
      selection_repetitions, estimator, combination, direct, lambdas, &
      validation_weight, alpha, periods, seasonal_types) result(out)
      !! Fit and forecast an ELM for a temporal-hierarchy aggregation level.
      real(dp), intent(in) :: series(:) !! Aggregated univariate time series.
      integer, intent(in), optional :: horizon !! Forecast horizon; defaults to the seasonal period.
      integer, intent(in), optional :: period !! Primary seasonal period.
      integer, intent(in), optional :: candidate_lags(:) !! Candidate response lags.
      logical, intent(in), optional :: keep(:) !! Candidate lags forced into the model.
      integer, intent(in), optional :: difference_lags(:) !! Supplied sequential differences.
      integer, intent(in), optional :: seasonal_type !! Deterministic seasonal input type.
      real(dp), intent(in), optional :: exogenous(:, :) !! Historical and future regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Lags by regressor.
      integer, intent(in), optional :: repetitions !! Final ELM ensemble size.
      integer, intent(in), optional :: selection_repetitions !! Hidden-selection repetitions.
      integer, intent(in), optional :: estimator !! ELM output-layer estimator.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: lambdas(:) !! Candidate ridge or lasso penalties.
      real(dp), intent(in), optional :: validation_weight !! Penalty-training fraction.
      real(dp), intent(in), optional :: alpha !! Automatic-test significance level.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Input types by period.
      type(nnfor_forecast_t) :: out
      type(nnfor_elm_auto_t) :: fit
      integer :: steps

      steps = thief_default_horizon(horizon, period, periods)
      if (steps < 1) then
         out%info = 1
         return
      end if
      fit = nnfor_elm_auto(series, period, candidate_lags, keep, &
         difference_lags, seasonal_type, exogenous, exogenous_lags, &
         repetitions, selection_repetitions, estimator, combination, direct, &
         lambdas, validation_weight, alpha, periods, seasonal_types)
      if (fit%info /= 0) then
         out%info = fit%info
         return
      end if
      out = nnfor_elm_auto_forecast(fit, steps, exogenous)
      if (out%info == 0) call align_thief_forecast(out, size(series))
   end function nnfor_elm_thief

   function nnfor_mlp_thief(series, horizon, period, candidate_lags, keep, &
      difference_lags, seasonal_type, exogenous, exogenous_lags, repetitions, &
      hidden_method, maximum_hidden, selection_repetitions, fold_count, &
      validation_fraction, combination, max_iterations, tolerance, decay, &
      alpha, periods, seasonal_types) result(out)
      !! Fit and forecast an MLP for a temporal-hierarchy aggregation level.
      real(dp), intent(in) :: series(:) !! Aggregated univariate time series.
      integer, intent(in), optional :: horizon !! Forecast horizon; defaults to the seasonal period.
      integer, intent(in), optional :: period !! Primary seasonal period.
      integer, intent(in), optional :: candidate_lags(:) !! Candidate response lags.
      logical, intent(in), optional :: keep(:) !! Candidate lags forced into the model.
      integer, intent(in), optional :: difference_lags(:) !! Supplied sequential differences.
      integer, intent(in), optional :: seasonal_type !! Deterministic seasonal input type.
      real(dp), intent(in), optional :: exogenous(:, :) !! Historical and future regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Lags by regressor.
      integer, intent(in), optional :: repetitions !! Final MLP ensemble size.
      integer, intent(in), optional :: hidden_method !! Hidden-size selection method.
      integer, intent(in), optional :: maximum_hidden !! Largest hidden size evaluated.
      integer, intent(in), optional :: selection_repetitions !! Networks per hidden size.
      integer, intent(in), optional :: fold_count !! Number of cross-validation folds.
      real(dp), intent(in), optional :: validation_fraction !! Holdout fraction.
      integer, intent(in), optional :: combination !! Ensemble combination operator.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      real(dp), intent(in), optional :: decay !! L2 weight decay.
      real(dp), intent(in), optional :: alpha !! Automatic-test significance level.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Input types by period.
      type(nnfor_forecast_t) :: out
      type(nnfor_mlp_auto_t) :: fit
      integer :: steps

      steps = thief_default_horizon(horizon, period, periods)
      if (steps < 1) then
         out%info = 1
         return
      end if
      fit = nnfor_mlp_auto(series, period, candidate_lags, keep, &
         difference_lags, seasonal_type, exogenous, exogenous_lags, &
         repetitions, hidden_method, maximum_hidden, selection_repetitions, &
         fold_count, validation_fraction, combination, max_iterations, &
         tolerance, decay, alpha, periods, seasonal_types)
      if (fit%info /= 0) then
         out%info = fit%info
         return
      end if
      out = nnfor_mlp_auto_forecast(fit, steps, exogenous)
      if (out%info == 0) call align_thief_forecast(out, size(series))
   end function nnfor_mlp_thief

   pure function nnfor_combine(values, method) result(combined)
      !! Combine ensemble members rowwise by their mean or median.
      real(dp), intent(in) :: values(:, :) !! Observation-by-member values.
      integer, intent(in) :: method !! Mean or median combination selector.
      real(dp), allocatable :: combined(:)
      integer :: row

      allocate(combined(size(values, 1)))
      if (method == nnfor_combine_mean) then
         combined = sum(values, 2)/real(size(values, 2), dp)
      else if (method == nnfor_combine_median) then
         do row = 1, size(values, 1)
            combined(row) = median(values(row, :))
         end do
      else
         do row = 1, size(values, 1)
            combined(row) = nnfor_kde_mode(values(row, :))
         end do
      end if
   end function nnfor_combine

   pure function nnfor_ridge_output_fit(design, response, lambdas, &
      validation_weight) result(out)
      !! Select and fit an intercept-preserving ridge output regression.
      real(dp), intent(in) :: design(:, :) !! Output-layer design with intercept first.
      real(dp), intent(in) :: response(:) !! Output-layer response.
      real(dp), intent(in) :: lambdas(:) !! Candidate nonnegative penalties.
      real(dp), intent(in) :: validation_weight !! Fraction used for penalty training.
      type(nnfor_output_fit_t) :: out
      type(gmdh_ridge_fit_t) :: ridge

      ridge = gmdh_ridge_fit(design, response, lambdas, validation_weight)
      out%info = ridge%info
      if (ridge%info /= 0) return
      out%coefficients = ridge%coefficients
      out%lambda = ridge%lambda
      out%validation_mse = ridge%validation_mse
      allocate(out%active(size(out%coefficients)), source=.true.)
      out%rss = sum((response - matmul(design, out%coefficients))**2)
   end function nnfor_ridge_output_fit

   pure function nnfor_lasso_output_fit(design, response, lambdas, &
      validation_weight) result(out)
      !! Select a lasso penalty by holdout error and refit the output layer.
      real(dp), intent(in) :: design(:, :) !! Output-layer design with intercept first.
      real(dp), intent(in) :: response(:) !! Output-layer response.
      real(dp), intent(in) :: lambdas(:) !! Candidate nonnegative lasso penalties.
      real(dp), intent(in) :: validation_weight !! Fraction used for penalty training.
      type(nnfor_output_fit_t) :: out
      real(dp), allocatable :: coefficients(:), prediction(:)
      real(dp) :: cost
      integer :: observations, training, candidate, best

      observations = size(design, 1)
      training = nint(real(observations, dp)*validation_weight)
      if (size(response) /= observations .or. size(design, 2) < 1 .or. &
         size(lambdas) < 1 .or. any(lambdas < 0.0_dp) .or. training < 2 .or. &
         training >= observations .or. validation_weight <= 0.0_dp .or. &
         validation_weight >= 1.0_dp .or. &
         .not. all(ieee_is_finite(design)) .or. &
         .not. all(ieee_is_finite(response))) then
         out%info = 1
         return
      end if
      best = 0
      do candidate = 1, size(lambdas)
         coefficients = lasso_coordinate_descent(design(:training, :), &
            response(:training), lambdas(candidate))
         prediction = matmul(design(training + 1:, :), coefficients)
         cost = sum((response(training + 1:) - prediction)**2)/ &
            real(observations - training, dp)
         if (cost < out%validation_mse) then
            out%validation_mse = cost
            out%lambda = lambdas(candidate)
            best = candidate
         end if
      end do
      if (best == 0) then
         out%info = 2
         return
      end if
      out%coefficients = lasso_coordinate_descent(design, response, out%lambda)
      out%active = abs(out%coefficients) > &
         sqrt(epsilon(1.0_dp))*max(1.0_dp, maxval(abs(out%coefficients)))
      out%active(1) = .true.
      out%rss = sum((response - matmul(design, out%coefficients))**2)
   end function nnfor_lasso_output_fit

   pure function nnfor_stepwise_output_fit(design, response) result(out)
      !! Perform backward output-variable selection using Akaike's criterion.
      real(dp), intent(in) :: design(:, :) !! Output-layer design with intercept first.
      real(dp), intent(in) :: response(:) !! Output-layer response.
      type(nnfor_output_fit_t) :: out
      type(least_squares_fit_t) :: current, candidate_fit, best_fit
      logical, allocatable :: active(:), candidate_active(:), best_active(:)
      real(dp) :: current_aic, candidate_aic, best_aic
      integer :: variable, observations

      observations = size(design, 1)
      if (size(response) /= observations .or. observations < 2 .or. &
         size(design, 2) < 1 .or. .not. all(ieee_is_finite(design)) .or. &
         .not. all(ieee_is_finite(response))) then
         out%info = 1
         return
      end if
      allocate(active(size(design, 2)), source=.true.)
      current = fit_selected_design(design, response, active)
      if (current%info /= 0) then
         out%info = 2
         return
      end if
      current_aic = regression_aic(current%rss, observations, count(active))
      do while (count(active) > 1)
         best_aic = current_aic
         best_active = active
         best_fit = current
         do variable = 2, size(active)
            if (.not. active(variable)) cycle
            candidate_active = active
            candidate_active(variable) = .false.
            candidate_fit = fit_selected_design(design, response, &
               candidate_active)
            if (candidate_fit%info /= 0) cycle
            candidate_aic = regression_aic(candidate_fit%rss, observations, &
               count(candidate_active))
            if (candidate_aic < best_aic - 1.0e-10_dp) then
               best_aic = candidate_aic
               best_active = candidate_active
               best_fit = candidate_fit
            end if
         end do
         if (all(best_active .eqv. active)) exit
         active = best_active
         current = best_fit
         current_aic = best_aic
      end do
      out%coefficients = current%coefficients
      out%active = active
      out%rss = current%rss
   end function nnfor_stepwise_output_fit

   pure real(dp) function nnfor_kde_mode(values) result(mode)
      !! Estimate an ensemble mode using a Gaussian kernel-density grid.
      real(dp), intent(in) :: values(:) !! Ensemble-member values.
      integer, parameter :: grid_count = 512
      real(dp) :: ordered(size(values)), bandwidth, spread_value, lower, upper
      real(dp) :: point, density, best_density, quartile_range
      integer :: grid, member

      if (size(values) == 0) then
         mode = 0.0_dp
         return
      end if
      if (size(values) == 1 .or. maxval(values) == minval(values)) then
         mode = values(1)
         return
      end if
      ordered = sorted(values)
      quartile_range = quantile(ordered, 0.75_dp) - quantile(ordered, 0.25_dp)
      spread_value = min(standard_deviation(values), quartile_range/1.34_dp)
      if (spread_value <= tiny(1.0_dp)) spread_value = &
         standard_deviation(values)
      bandwidth = 0.9_dp*spread_value*real(size(values), dp)**(-0.2_dp)
      bandwidth = max(bandwidth, sqrt(epsilon(1.0_dp))* &
         max(1.0_dp, maxval(abs(values))))
      lower = minval(values) - 0.1_dp*(maxval(values) - minval(values))
      upper = maxval(values) + 0.1_dp*(maxval(values) - minval(values))
      best_density = -1.0_dp
      mode = lower
      do grid = 1, grid_count
         point = lower + real(grid - 1, dp)*(upper - lower)/ &
            real(grid_count - 1, dp)
         density = 0.0_dp
         do member = 1, size(values)
            density = density + exp(-0.5_dp*((point - values(member))/ &
               bandwidth)**2)
         end do
         if (density > best_density) then
            best_density = density
            mode = point
         end if
      end do
   end function nnfor_kde_mode

   pure function lasso_coordinate_descent(design, response, lambda) &
      result(coefficients)
      !! Fit a standardized lasso with an unpenalized intercept.
      real(dp), intent(in) :: design(:, :) !! Design matrix with intercept first.
      real(dp), intent(in) :: response(:) !! Regression response.
      real(dp), intent(in) :: lambda !! Nonnegative lasso penalty.
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: centered(:, :), means(:), scales(:), beta(:)
      real(dp), allocatable :: residual(:)
      real(dp) :: response_mean, denominator, rho, old, change
      integer :: observations, variables, variable, iteration

      observations = size(design, 1)
      variables = size(design, 2) - 1
      allocate(coefficients(variables + 1), source=0.0_dp)
      if (variables == 0) then
         coefficients(1) = sum(response)/real(observations, dp)
         return
      end if
      allocate(centered(observations, variables), means(variables))
      allocate(scales(variables), beta(variables), source=0.0_dp)
      response_mean = sum(response)/real(observations, dp)
      do variable = 1, variables
         means(variable) = sum(design(:, variable + 1))/ &
            real(observations, dp)
         centered(:, variable) = design(:, variable + 1) - means(variable)
         scales(variable) = sqrt(sum(centered(:, variable)**2)/ &
            real(observations, dp))
         if (scales(variable) > sqrt(epsilon(1.0_dp))) then
            centered(:, variable) = centered(:, variable)/scales(variable)
         else
            centered(:, variable) = 0.0_dp
            scales(variable) = 1.0_dp
         end if
      end do
      residual = response - response_mean
      do iteration = 1, 5000
         change = 0.0_dp
         do variable = 1, variables
            old = beta(variable)
            residual = residual + centered(:, variable)*old
            denominator = sum(centered(:, variable)**2)/real(observations, dp)
            rho = sum(centered(:, variable)*residual)/real(observations, dp)
            if (denominator > tiny(1.0_dp)) then
               beta(variable) = soft_threshold(rho, lambda)/denominator
            else
               beta(variable) = 0.0_dp
            end if
            residual = residual - centered(:, variable)*beta(variable)
            change = max(change, abs(beta(variable) - old))
         end do
         if (change <= 1.0e-10_dp*(1.0_dp + maxval(abs(beta)))) exit
      end do
      coefficients(2:) = beta/scales
      coefficients(1) = response_mean - dot_product(means, coefficients(2:))
   end function lasso_coordinate_descent

   pure elemental real(dp) function soft_threshold(value, threshold) result(out)
      !! Apply the scalar lasso proximal operator.
      real(dp), intent(in) :: value !! Unthresholded coordinate value.
      real(dp), intent(in) :: threshold !! Nonnegative shrinkage threshold.

      out = sign(max(abs(value) - threshold, 0.0_dp), value)
   end function soft_threshold

   pure function fit_selected_design(design, response, active) result(out)
      !! Fit least squares using the selected output-layer columns.
      real(dp), intent(in) :: design(:, :) !! Complete output-layer design.
      real(dp), intent(in) :: response(:) !! Output-layer response.
      logical, intent(in) :: active(:) !! Active design-column mask.
      type(least_squares_fit_t) :: out
      real(dp), allocatable :: selected(:, :), gram(:, :), inverse(:, :)
      real(dp), allocatable :: selected_coefficients(:)
      integer, allocatable :: indices(:)
      integer :: columns, column, position, status

      columns = count(active)
      allocate(out%coefficients(size(active)), source=0.0_dp)
      if (size(active) /= size(design, 2) .or. columns < 1) then
         out%info = 1
         return
      end if
      allocate(indices(columns), selected(size(design, 1), columns))
      position = 0
      do column = 1, size(active)
         if (.not. active(column)) cycle
         position = position + 1
         indices(position) = column
         selected(:, position) = design(:, column)
      end do
      allocate(gram(columns, columns), inverse(columns, columns))
      gram = matmul(transpose(selected), selected)
      call symmetric_pseudoinverse(gram, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      selected_coefficients = matmul(inverse, &
         matmul(transpose(selected), response))
      do position = 1, columns
         out%coefficients(indices(position)) = selected_coefficients(position)
      end do
      out%rss = sum((response - matmul(design, out%coefficients))**2)
   end function fit_selected_design

   pure real(dp) function regression_aic(rss, observations, parameters) &
      result(value)
      !! Evaluate Gaussian regression AIC up to an additive constant.
      real(dp), intent(in) :: rss !! Residual sum of squares.
      integer, intent(in) :: observations !! Number of fitted observations.
      integer, intent(in) :: parameters !! Number of active coefficients.

      value = real(observations, dp)*log(max(rss/real(observations, dp), &
         epsilon(1.0_dp))) + 2.0_dp*real(parameters, dp)
   end function regression_aic

   subroutine display_nnfor_elm_fast(model, unit, print_obs)
      !! Display a fast ELM ensemble and optionally its fitted observations.
      type(nnfor_elm_fast_t), intent(in) :: model !! Fast ELM ensemble to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted values and residuals.
      integer :: destination, member, observation
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Fast extreme-learning-machine ensemble'
      write(destination, '(a, i0)') 'Status: ', model%info
      write(destination, '(a, i0)') 'Predictors: ', model%predictor_count
      write(destination, '(a, i0)') 'Repetitions: ', model%repetitions
      write(destination, '(a, l1)') 'Direct connections: ', model%direct
      write(destination, '(a, l1)') 'Orthogonal initialization: ', &
         model%orthogonal
      if (model%orthogonal) write(destination, '(a, i0)') &
         'Orthogonalized hidden layers: ', model%orthogonalized_layer_count
      write(destination, '(a, es14.6)') 'MSE: ', model%mse
      if (allocated(model%members)) then
         do member = 1, size(model%members)
            if (allocated(model%members(member)%hidden_counts)) then
               write(destination, '(a, i0, a, *(i0, 1x))') 'Member ', member, &
                  ', hidden layers ', model%members(member)%hidden_counts
               write(destination, '(a, es14.6)') '  Lambda ', &
                  model%members(member)%lambda
            else
               write(destination, '(a, i0, a, i0, a, es14.6)') 'Member ', &
                  member, ', hidden nodes ', &
                  model%members(member)%hidden_count, ', lambda ', &
                  model%members(member)%lambda
            end if
         end do
      end if
      if (show_observations .and. allocated(model%fitted)) then
         write(destination, '(a)') 'Index, fitted value, residual:'
         do observation = 1, size(model%fitted)
            write(destination, '(i8, 2(1x, es14.6))') observation, &
               model%fitted(observation), model%residuals(observation)
         end do
      end if
   end subroutine display_nnfor_elm_fast

   subroutine display_nnfor_elm_model(model, unit, print_obs)
      !! Display a lagged ELM forecasting model.
      type(nnfor_elm_model_t), intent(in) :: model !! Lagged ELM model to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted observations.
      integer :: destination, observation, original_index
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Autoregressive ELM model'
      if (allocated(model%lags)) then
         write(destination, '(a, *(i0, 1x))') 'Lags: ', model%lags
      end if
      if (model%extended_preprocessing) then
         write(destination, '(a, *(i0, 1x))') 'Difference lags: ', &
            model%preprocessing%difference_lags
         write(destination, '(a, *(i0, 1x))') 'Seasonal periods: ', &
            model%preprocessing%periods
         write(destination, '(a, *(i0, 1x))') 'Seasonal input types: ', &
            model%preprocessing%seasonal_types
         write(destination, '(a, i0)') 'Exogenous variables: ', &
            model%preprocessing%exogenous_count
         write(destination, '(a, es14.6)') 'Original-scale MSE: ', model%mse
      end if
      call display_nnfor_elm_fast(model%network, destination)
      if (show_observations .and. allocated(model%fitted)) then
         write(destination, '(a)') 'Original index, fitted value, residual:'
         do observation = 1, size(model%fitted)
            if (model%extended_preprocessing) then
               original_index = model%preprocessing%start_index + &
                  observation - 1
            else
               original_index = maxval(model%lags) + observation
            end if
            write(destination, '(i8, 2(1x, es14.6))') original_index, &
               model%fitted(observation), model%residuals(observation)
         end do
      end if
   end subroutine display_nnfor_elm_model

   subroutine display_nnfor_mlp_model(model, unit, print_obs)
      !! Display a lagged MLP ensemble and optionally its fitted observations.
      type(nnfor_mlp_model_t), intent(in) :: model !! Lagged MLP model to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted observations.
      integer :: destination, observation, original_index
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Autoregressive MLP ensemble'
      write(destination, '(a, i0)') 'Status: ', model%info
      if (allocated(model%hidden_counts)) then
         write(destination, '(a, *(i0, 1x))') 'Hidden layers: ', &
            model%hidden_counts
      else
         write(destination, '(a, i0)') 'Hidden nodes: ', model%hidden_count
      end if
      write(destination, '(a, i0)') 'Repetitions: ', model%repetitions
      write(destination, '(a, es14.6)') 'MSE: ', model%mse
      if (allocated(model%lags)) &
         write(destination, '(a, *(i0, 1x))') 'Lags: ', model%lags
      if (model%extended_preprocessing) then
         write(destination, '(a, *(i0, 1x))') 'Difference lags: ', &
            model%preprocessing%difference_lags
         write(destination, '(a, *(i0, 1x))') 'Seasonal periods: ', &
            model%preprocessing%periods
         write(destination, '(a, *(i0, 1x))') 'Seasonal input types: ', &
            model%preprocessing%seasonal_types
         write(destination, '(a, i0)') 'Exogenous variables: ', &
            model%preprocessing%exogenous_count
      end if
      if (show_observations .and. allocated(model%fitted)) then
         if (model%extended_preprocessing) then
            write(destination, '(a)') &
               'Original index, fitted value, residual:'
         else
            write(destination, '(a)') 'Index, fitted value, residual:'
         end if
         do observation = 1, size(model%fitted)
            if (model%extended_preprocessing) then
               original_index = model%preprocessing%start_index + &
                  observation - 1
            else
               original_index = maxval(model%lags) + observation
            end if
            write(destination, '(i8, 2(1x, es14.6))') original_index, &
               model%fitted(observation), model%residuals(observation)
         end do
      end if
   end subroutine display_nnfor_mlp_model

   subroutine display_nnfor_elm_auto(fit, unit, print_obs)
      !! Display automatic ELM selections and the resulting fitted model.
      type(nnfor_elm_auto_t), intent(in) :: fit !! Automatic ELM result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted observations.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Automatically specified ELM model'
      write(destination, '(a, i0)') 'Status: ', fit%info
      if (allocated(fit%candidate_lags)) &
         write(destination, '(a, *(i0, 1x))') 'Candidate lags: ', &
            fit%candidate_lags
      if (allocated(fit%lag_selection%selected_lags)) &
         write(destination, '(a, *(i0, 1x))') 'Selected lags: ', &
            fit%lag_selection%selected_lags
      if (allocated(fit%difference_lags)) &
         write(destination, '(a, *(i0, 1x))') 'Difference lags: ', &
            fit%difference_lags
      if (fit%difference_selection%multiplicative_seasonality%period > 1) &
         call display_nnfor_mseason( &
            fit%difference_selection%multiplicative_seasonality, destination)
      if (fit%difference_selection%canova_hansen%period > 1) &
         call display_nnfor_canova_hansen( &
            fit%difference_selection%canova_hansen, destination)
      write(destination, '(a, i0)') 'Selected hidden nodes: ', &
         fit%hidden_selection%selected
      call display_nnfor_elm_model(fit%model, destination, print_obs)
   end subroutine display_nnfor_elm_auto

   subroutine display_nnfor_mlp_auto(fit, unit, print_obs)
      !! Display automatic MLP selections and the resulting fitted model.
      type(nnfor_mlp_auto_t), intent(in) :: fit !! Automatic MLP result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted observations.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Automatically specified MLP model'
      write(destination, '(a, i0)') 'Status: ', fit%info
      if (allocated(fit%candidate_lags)) &
         write(destination, '(a, *(i0, 1x))') 'Candidate lags: ', &
            fit%candidate_lags
      if (allocated(fit%lag_selection%selected_lags)) &
         write(destination, '(a, *(i0, 1x))') 'Selected lags: ', &
            fit%lag_selection%selected_lags
      if (allocated(fit%difference_lags)) &
         write(destination, '(a, *(i0, 1x))') 'Difference lags: ', &
            fit%difference_lags
      if (fit%difference_selection%multiplicative_seasonality%period > 1) &
         call display_nnfor_mseason( &
            fit%difference_selection%multiplicative_seasonality, destination)
      if (fit%difference_selection%canova_hansen%period > 1) &
         call display_nnfor_canova_hansen( &
            fit%difference_selection%canova_hansen, destination)
      write(destination, '(a, i0)') 'Selected hidden nodes: ', &
         fit%hidden_selection%selected
      call display_nnfor_mlp_model(fit%model, destination, print_obs)
   end subroutine display_nnfor_mlp_auto

   subroutine display_nnfor_mseason(result, unit)
      !! Display a multiplicative-seasonality correlation test result.
      type(nnfor_mseason_t), intent(in) :: result !! Test result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Multiplicative seasonality test'
      write(destination, '(a, i0)') 'Status: ', result%info
      write(destination, '(a, i0)') 'Seasonal period: ', result%period
      write(destination, '(a, es14.6)') 'Correlation: ', result%statistic
      write(destination, '(a, es14.6)') 'P-value: ', result%p_value
      write(destination, '(a, l1)') 'Multiplicative: ', result%multiplicative
   end subroutine display_nnfor_mseason

   subroutine display_nnfor_canova_hansen(result, unit)
      !! Display a joint Canova-Hansen seasonal-stability test result.
      type(nnfor_canova_hansen_t), intent(in) :: result !! Test result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Canova-Hansen seasonal-stability test'
      write(destination, '(a, i0)') 'Status: ', result%info
      write(destination, '(a, i0)') 'Seasonal period: ', result%period
      write(destination, '(a, i0)') 'Newey-West order: ', &
         result%newey_west_order
      write(destination, '(a, es14.6)') 'Joint statistic: ', result%statistic
      write(destination, '(a, es14.6)') 'Critical value: ', &
         result%critical_value
      write(destination, '(a, es14.6)') 'Diagnostic p-value: ', result%p_value
      write(destination, '(a, l1)') 'Seasonal difference required: ', &
         result%difference_required
   end subroutine display_nnfor_canova_hansen

   subroutine display_nnfor_forecast(forecast, unit)
      !! Display ensemble point forecasts and their member range.
      type(nnfor_forecast_t), intent(in) :: forecast !! Neural-network forecast to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination, step

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Neural-network ensemble forecast'
      write(destination, '(a, i0)') 'Status: ', forecast%info
      if (allocated(forecast%mean) .and. allocated(forecast%all_mean)) then
         write(destination, '(a)') 'Horizon, combined, minimum, maximum:'
         do step = 1, size(forecast%mean)
            write(destination, '(i8, 3(1x, es14.6))') step, &
               forecast%mean(step), minval(forecast%all_mean(step, :)), &
               maxval(forecast%all_mean(step, :))
         end do
      end if
   end subroutine display_nnfor_forecast

   pure function build_automatic_specification(series, period, candidate_lags, &
      keep, difference_lags, seasonal_type, exogenous, exogenous_lags, alpha, &
      periods, seasonal_types) result(out)
      !! Construct automatic differences, deterministic inputs, and selected lags.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in), optional :: period !! Seasonal period.
      integer, intent(in), optional :: candidate_lags(:) !! Candidate response lags.
      logical, intent(in), optional :: keep(:) !! Candidate lags forced into the model.
      integer, intent(in), optional :: difference_lags(:) !! Supplied sequential differences.
      integer, intent(in), optional :: seasonal_type !! Supplied deterministic seasonal type.
      real(dp), intent(in), optional :: exogenous(:, :) !! Original-time exogenous regressors.
      type(integer_vector_t), intent(in), optional :: exogenous_lags(:) !! Lags by exogenous regressor.
      real(dp), intent(in), optional :: alpha !! Automatic-test significance level.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.
      integer, intent(in), optional :: seasonal_types(:) !! Supplied input types by period.
      type(automatic_specification_t) :: out
      type(nnfor_seasonality_t) :: seasonal_diagnostic
      real(dp), allocatable :: differenced(:)
      integer :: maximum_lag, lag, seasonal_index

      if (present(periods)) then
         out%periods = periods
      else if (present(period)) then
         out%periods = [period]
      else
         out%periods = [1]
      end if
      if (size(out%periods) < 1) then
         out%info = 1
         return
      end if
      out%period = maxval(out%periods)
      if (any(out%periods < 1) .or. size(series) < 6 .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (present(difference_lags)) then
         out%difference_lags = difference_lags
         out%difference_selection%difference_lags = difference_lags
         out%automatic_differences = .false.
      else
         out%difference_selection = nnfor_select_differences(series, &
            out%period, alpha)
         if (out%difference_selection%info /= 0) then
            out%info = 2
            return
         end if
         out%difference_lags = &
            out%difference_selection%difference_lags
      end if
      differenced = nnfor_difference(series, out%difference_lags)
      if (size(differenced) < 4) then
         out%info = 3
         return
      end if
      if (present(candidate_lags)) then
         out%candidate_lags = candidate_lags
      else
         maximum_lag = merge(out%period, 4, out%period > 3)
         maximum_lag = min(maximum_lag, size(differenced) - 3)
         if (maximum_lag < 1) then
            out%info = 3
            return
         end if
         out%candidate_lags = [(lag, lag=1, maximum_lag)]
      end if
      if (present(seasonal_types)) then
         if (size(seasonal_types) /= size(out%periods)) then
            out%info = 4
            return
         end if
         out%seasonal_types = seasonal_types
      else if (present(seasonal_type)) then
         allocate(out%seasonal_types(size(out%periods)), &
            source=seasonal_type)
      else
         allocate(out%seasonal_types(size(out%periods)), &
            source=nnfor_seasonal_none)
         do seasonal_index = 1, size(out%periods)
            if (out%periods(seasonal_index) <= 1 .or. &
               any(out%difference_lags == out%periods(seasonal_index))) cycle
            if (out%automatic_differences .and. &
               out%periods(seasonal_index) == out%period) then
               seasonal_diagnostic = out%difference_selection%seasonality
            else
               seasonal_diagnostic = nnfor_season_check(series, &
                  out%periods(seasonal_index), alpha)
            end if
            if (seasonal_diagnostic%info /= 0) then
               out%info = 4
               return
            end if
            if (seasonal_diagnostic%seasonal) then
               if (size(out%periods) == 1 .and. &
                  out%periods(seasonal_index) <= 12) then
                  out%seasonal_types(seasonal_index) = nnfor_seasonal_binary
               else
                  out%seasonal_types(seasonal_index) = &
                     nnfor_seasonal_trigonometric
               end if
            end if
         end do
      end if
      out%seasonal_type = nnfor_seasonal_none
      if (size(out%seasonal_types) == 1) &
         out%seasonal_type = out%seasonal_types(1)
      out%preprocessing = nnfor_preprocess(series, out%candidate_lags, &
         out%difference_lags, out%period, out%seasonal_type, exogenous, &
         exogenous_lags, out%periods, out%seasonal_types)
      if (out%preprocessing%info /= 0) then
         out%info = 5
         return
      end if
      out%lag_selection = nnfor_select_lags(out%preprocessing, keep)
      if (out%lag_selection%info /= 0) then
         out%info = 6
         return
      end if
      out%selected_lags = out%lag_selection%selected_lags
      out%preprocessing = nnfor_preprocess(series, out%selected_lags, &
         out%difference_lags, out%period, out%seasonal_type, exogenous, &
         exogenous_lags, out%periods, out%seasonal_types)
      if (out%preprocessing%info /= 0) out%info = 7
   end function build_automatic_specification

   pure subroutine copy_automatic_specification(specification, candidate_lags, &
      difference_lags, difference_selection, lag_selection, seasonal_type, &
      automatic_differences, periods, seasonal_types)
      !! Copy public diagnostics from the internal automatic specification.
      type(automatic_specification_t), intent(in) :: specification !! Internal specification result.
      integer, allocatable, intent(out) :: candidate_lags(:) !! Candidate response lags.
      integer, allocatable, intent(out) :: difference_lags(:) !! Applied sequential differences.
      type(nnfor_difference_selection_t), intent(out) :: difference_selection !! Difference diagnostics.
      type(nnfor_lag_selection_t), intent(out) :: lag_selection !! Lag-selection diagnostics.
      integer, intent(out) :: seasonal_type !! Applied deterministic seasonal type.
      logical, intent(out) :: automatic_differences !! Whether differences were selected.
      integer, allocatable, intent(out) :: periods(:) !! Applied seasonal periods.
      integer, allocatable, intent(out) :: seasonal_types(:) !! Applied input types by period.

      if (allocated(specification%candidate_lags)) &
         candidate_lags = specification%candidate_lags
      if (allocated(specification%difference_lags)) &
         difference_lags = specification%difference_lags
      if (allocated(specification%periods)) periods = specification%periods
      if (allocated(specification%seasonal_types)) &
         seasonal_types = specification%seasonal_types
      difference_selection = specification%difference_selection
      lag_selection = specification%lag_selection
      seasonal_type = specification%seasonal_type
      automatic_differences = specification%automatic_differences
   end subroutine copy_automatic_specification

   pure subroutine pack_predictor_rows(predictors, mask, packed)
      !! Pack selected rows of a predictor matrix without flattening columns.
      real(dp), intent(in) :: predictors(:, :) !! Source predictor matrix.
      logical, intent(in) :: mask(:) !! Row-selection mask.
      real(dp), intent(out) :: packed(:, :) !! Packed predictor rows.
      integer :: predictor

      do predictor = 1, size(predictors, 2)
         packed(:, predictor) = pack(predictors(:, predictor), mask)
      end do
   end subroutine pack_predictor_rows

   pure subroutine append_seasonal_predictors(preprocessing, time, &
      predictors, column)
      !! Append every stored seasonal-period input for one forecast time.
      type(nnfor_preprocessing_t), intent(in) :: preprocessing !! Stored seasonal specification.
      integer, intent(in) :: time !! Original-series time index.
      real(dp), intent(inout) :: predictors(:, :) !! Forecast predictor row.
      integer, intent(inout) :: column !! Last populated predictor column.
      integer :: seasonality, phase, level

      do seasonality = 1, size(preprocessing%periods)
         phase = modulo(time - 1, preprocessing%periods(seasonality)) + 1
         if (preprocessing%seasonal_types(seasonality) == &
            nnfor_seasonal_binary) then
            do level = 1, preprocessing%periods(seasonality) - 1
               column = column + 1
               predictors(1, column) = &
                  merge(1.0_dp, 0.0_dp, phase == level)
            end do
         else if (preprocessing%seasonal_types(seasonality) == &
            nnfor_seasonal_trigonometric) then
            column = column + 1
            predictors(1, column) = sin(2.0_dp*acos(-1.0_dp)* &
               real(time, dp)/real(preprocessing%periods(seasonality), dp))
            column = column + 1
            predictors(1, column) = cos(2.0_dp*acos(-1.0_dp)* &
               real(time, dp)/real(preprocessing%periods(seasonality), dp))
         end if
      end do
   end subroutine append_seasonal_predictors

   pure elemental real(dp) function logistic_sigmoid(value) result(out)
      !! Evaluate a numerically stable logistic hidden-unit activation.
      real(dp), intent(in) :: value !! Hidden-unit linear activation.

      if (value >= 0.0_dp) then
         out = 1.0_dp/(1.0_dp + exp(-value))
      else
         out = exp(value)/(1.0_dp + exp(value))
      end if
   end function logistic_sigmoid

   pure function selected_coefficient_significance(design, response, active, &
      coefficients) result(out)
      !! Compute two-sided Student-t p-values for selected regression columns.
      real(dp), intent(in) :: design(:, :) !! Complete regression design.
      real(dp), intent(in) :: response(:) !! Regression response.
      logical, intent(in) :: active(:) !! Selected design-column mask.
      real(dp), intent(in) :: coefficients(:) !! Full coefficient vector.
      type(significance_result_t) :: out
      real(dp), allocatable :: selected(:, :), gram(:, :), inverse(:, :)
      real(dp), allocatable :: residuals(:)
      real(dp) :: sigma_squared, standard_error, statistic, degrees_real
      integer, allocatable :: indices(:)
      integer :: columns, position, column, degrees, status

      columns = count(active)
      degrees = size(response) - columns
      if (size(design, 1) /= size(response) .or. &
         size(design, 2) /= size(active) .or. &
         size(coefficients) /= size(active) .or. columns < 1 .or. &
         degrees < 1) then
         out%info = 1
         return
      end if
      allocate(out%p_value(size(active)), source=1.0_dp)
      allocate(indices(columns), selected(size(response), columns))
      position = 0
      do column = 1, size(active)
         if (.not. active(column)) cycle
         position = position + 1
         indices(position) = column
         selected(:, position) = design(:, column)
      end do
      gram = matmul(transpose(selected), selected)
      allocate(inverse(columns, columns))
      call symmetric_pseudoinverse(gram, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      residuals = response - matmul(design, coefficients)
      sigma_squared = sum(residuals**2)/real(degrees, dp)
      degrees_real = real(degrees, dp)
      do position = 1, columns
         standard_error = sqrt(max(0.0_dp, &
            sigma_squared*inverse(position, position)))
         if (standard_error <= tiny(1.0_dp)) then
            if (abs(coefficients(indices(position))) > tiny(1.0_dp)) &
               out%p_value(indices(position)) = 0.0_dp
         else
            statistic = coefficients(indices(position))/standard_error
            out%p_value(indices(position)) = regularized_beta( &
               degrees_real/(degrees_real + statistic**2), &
               0.5_dp*degrees_real, 0.5_dp)
         end if
      end do
   end function selected_coefficient_significance

   pure subroutine pearson_greater_test(first, second, correlation, p_value)
      !! Compute Pearson correlation and its one-sided Student test.
      real(dp), intent(in) :: first(:) !! First sample.
      real(dp), intent(in) :: second(:) !! Second sample.
      real(dp), intent(out) :: correlation !! Sample correlation.
      real(dp), intent(out) :: p_value !! Upper-tail null probability.
      real(dp) :: first_scale, second_scale, degrees, statistic, beta_value

      correlation = 0.0_dp
      p_value = 1.0_dp
      if (size(first) /= size(second) .or. size(first) < 3) return
      first_scale = standard_deviation(first)
      second_scale = standard_deviation(second)
      if (first_scale <= tiny(1.0_dp) .or. &
         second_scale <= tiny(1.0_dp)) return
      correlation = max(-1.0_dp, min(1.0_dp, &
         covariance(first, second)/(first_scale*second_scale)))
      if (correlation >= 1.0_dp - epsilon(1.0_dp)) then
         p_value = 0.0_dp
         return
      end if
      if (correlation <= -1.0_dp + epsilon(1.0_dp)) return
      degrees = real(size(first) - 2, dp)
      statistic = correlation*sqrt(degrees/ &
         max(tiny(1.0_dp), 1.0_dp - correlation**2))
      beta_value = regularized_beta(degrees/(degrees + statistic**2), &
         0.5_dp*degrees, 0.5_dp)
      if (statistic >= 0.0_dp) then
         p_value = 0.5_dp*beta_value
      else
         p_value = 1.0_dp - 0.5_dp*beta_value
      end if
   end subroutine pearson_greater_test

   pure integer function thief_default_horizon(horizon, period, periods) &
      result(steps)
      !! Select the THieF callback horizon from explicit or seasonal inputs.
      integer, intent(in), optional :: horizon !! Explicit horizon.
      integer, intent(in), optional :: period !! Primary seasonal period.
      integer, intent(in), optional :: periods(:) !! Multiple seasonal periods.

      steps = 1
      if (present(period)) steps = period
      if (present(periods)) then
         if (size(periods) > 0) steps = maxval(periods)
      end if
      if (present(horizon)) steps = horizon
   end function thief_default_horizon

   pure subroutine align_thief_forecast(forecast, observations)
      !! Pad callback fitted values and residuals to the complete input sample.
      type(nnfor_forecast_t), intent(inout) :: forecast !! Forecast to align.
      integer, intent(in) :: observations !! Original sample size.
      real(dp), allocatable :: aligned(:)
      integer :: first

      if (allocated(forecast%fitted)) then
         if (size(forecast%fitted) <= observations) then
            allocate(aligned(observations), source=quiet_nan())
            first = observations - size(forecast%fitted) + 1
            aligned(first:) = forecast%fitted
            call move_alloc(aligned, forecast%fitted)
         end if
      end if
      if (allocated(forecast%residuals)) then
         if (size(forecast%residuals) <= observations) then
            allocate(aligned(observations), source=quiet_nan())
            first = observations - size(forecast%residuals) + 1
            aligned(first:) = forecast%residuals
            call move_alloc(aligned, forecast%residuals)
         end if
      end if
   end subroutine align_thief_forecast

   pure function centered_moving_average_fill(series, period) result(smoothed)
      !! Compute nnfor's centered moving average and extend its endpoints.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: period !! Moving-average period.
      real(dp), allocatable :: smoothed(:), weights(:)
      integer :: width, left, first_valid, last_valid, observation

      width = period
      if (mod(period, 2) == 0) width = period + 1
      allocate(weights(width), source=1.0_dp/real(period, dp))
      if (mod(period, 2) == 0) then
         weights(1) = 0.5_dp/real(period, dp)
         weights(width) = weights(1)
      end if
      allocate(smoothed(size(series)))
      left = width/2
      first_valid = left + 1
      last_valid = size(series) - (width - left - 1)
      do observation = first_valid, last_valid
         smoothed(observation) = dot_product(weights, &
            series(observation - left:observation - left + width - 1))
      end do
      smoothed(:first_valid - 1) = smoothed(first_valid)
      smoothed(last_valid + 1:) = smoothed(last_valid)
   end function centered_moving_average_fill

   pure function restore_fitted_values(preprocessing, deepest_values) &
      result(restored)
      !! Reverse sequential differences for aligned one-step fitted values.
      type(nnfor_preprocessing_t), intent(in) :: preprocessing !! Stored differencing levels.
      real(dp), intent(in) :: deepest_values(:, :) !! Final-difference fits by member.
      real(dp), allocatable :: restored(:, :)
      real(dp) :: value
      integer :: row, member, stage, deepest_index, offset_after, base_index

      allocate(restored(size(deepest_values, 1), size(deepest_values, 2)))
      do member = 1, size(deepest_values, 2)
         do row = 1, size(deepest_values, 1)
            value = deepest_values(row, member)
            deepest_index = preprocessing%maximum_lag + row
            offset_after = 0
            do stage = size(preprocessing%difference_lags), 1, -1
               base_index = deepest_index + offset_after
               value = value + preprocessing%difference_levels(stage)% &
                  values(base_index)
               offset_after = offset_after + &
                  preprocessing%difference_lags(stage)
            end do
            restored(row, member) = value
         end do
      end do
   end function restore_fitted_values

   pure function canova_hansen_response_critical_values(period, &
      effective_observations) result(values)
      !! Evaluate joint trigonometric response-surface critical quantiles.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: effective_observations !! Regression residual degrees of freedom.
      real(dp) :: values(size(ch_response_probabilities))
      real(dp) :: regressors(7), sample
      integer :: probability

      sample = real(effective_observations, dp)
      regressors = [1.0_dp, 1.0_dp/sample, 1.0_dp/sample**2, &
         real(period, dp)/sample, real(period, dp)/sample**2, &
         real(period - 1, dp), real(period - 1, dp)**2]
      do probability = 1, size(values)
         values(probability) = dot_product( &
            ch_joint_coefficients(:, probability), regressors)
      end do
   end function canova_hansen_response_critical_values

   pure real(dp) function canova_hansen_critical_value(period, seasonal_type, &
      effective_observations, alpha) result(value)
      !! Interpolate the Canova-Hansen critical value at a significance level.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: seasonal_type !! Trigonometric or dummy formulation.
      integer, intent(in) :: effective_observations !! Regression residual degrees of freedom.
      real(dp), intent(in) :: alpha !! Significance level in [0.01, 0.10].
      real(dp) :: values(size(ch_response_probabilities)), target, fraction
      integer :: index

      if (seasonal_type == nnfor_ch_dummy) then
         do index = 1, size(ch_raw_probabilities) - 1
            if (alpha <= ch_raw_probabilities(index) .and. &
               alpha >= ch_raw_probabilities(index + 1)) then
               fraction = (alpha - ch_raw_probabilities(index))/ &
                  (ch_raw_probabilities(index + 1) - &
                  ch_raw_probabilities(index))
               value = ch_raw_critical_values(index, period) + fraction* &
                  (ch_raw_critical_values(index + 1, period) - &
                  ch_raw_critical_values(index, period))
               return
            end if
         end do
         value = ch_raw_critical_values(size(ch_raw_probabilities), period)
         return
      end if
      values = canova_hansen_response_critical_values(period, &
         effective_observations)
      target = 1.0_dp - alpha
      do index = 1, size(ch_response_probabilities) - 1
         if (target >= ch_response_probabilities(index) .and. &
            target <= ch_response_probabilities(index + 1)) then
            fraction = (target - ch_response_probabilities(index))/ &
               (ch_response_probabilities(index + 1) - &
               ch_response_probabilities(index))
            value = values(index) + fraction*(values(index + 1) - values(index))
            return
         end if
      end do
      value = values(size(values))
   end function canova_hansen_critical_value

   pure real(dp) function canova_hansen_p_value(statistic, period, &
      seasonal_type, effective_observations) result(value)
      !! Interpolate a diagnostic p-value over the supported critical range.
      real(dp), intent(in) :: statistic !! Joint Canova-Hansen statistic.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: seasonal_type !! Trigonometric or dummy formulation.
      integer, intent(in) :: effective_observations !! Regression residual degrees of freedom.
      real(dp) :: values(size(ch_response_probabilities)), fraction
      integer :: index

      if (seasonal_type == nnfor_ch_dummy) then
         if (statistic <= ch_raw_critical_values(1, period)) then
            value = ch_raw_probabilities(1)
            return
         end if
         do index = 1, size(ch_raw_probabilities) - 1
            if (statistic >= ch_raw_critical_values(index, period) .and. &
               statistic <= ch_raw_critical_values(index + 1, period)) then
               fraction = (statistic - ch_raw_critical_values(index, period))/ &
                  (ch_raw_critical_values(index + 1, period) - &
                  ch_raw_critical_values(index, period))
               value = ch_raw_probabilities(index) + fraction* &
                  (ch_raw_probabilities(index + 1) - ch_raw_probabilities(index))
               return
            end if
         end do
         value = ch_raw_probabilities(size(ch_raw_probabilities))
         return
      end if
      values = canova_hansen_response_critical_values(period, &
         effective_observations)
      if (statistic <= values(1)) then
         value = 1.0_dp - ch_response_probabilities(1)
         return
      end if
      do index = 1, size(values) - 1
         if (statistic >= values(index) .and. statistic <= values(index + 1)) then
            fraction = (statistic - values(index))/ &
               (values(index + 1) - values(index))
            value = 1.0_dp - (ch_response_probabilities(index) + fraction* &
               (ch_response_probabilities(index + 1) - &
               ch_response_probabilities(index)))
            return
         end if
      end do
      value = 1.0_dp - ch_response_probabilities(size(values))
   end function canova_hansen_p_value

   pure subroutine lagged_regression_inputs(series, lags, response, predictors, &
      info)
      !! Construct aligned autoregressive response and predictor arrays.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: lags(:) !! Positive autoregressive lags.
      real(dp), allocatable, intent(out) :: response(:) !! Aligned response values.
      real(dp), allocatable, intent(out) :: predictors(:, :) !! Lagged predictor rows.
      integer, intent(out) :: info !! Zero on success and one for invalid inputs.
      integer :: maximum_lag, row, lag

      info = 1
      if (size(lags) < 1 .or. any(lags < 1) .or. &
         size(series) <= maxval(lags) + 2 .or. &
         .not. all(ieee_is_finite(series))) return
      maximum_lag = maxval(lags)
      allocate(response(size(series) - maximum_lag))
      allocate(predictors(size(response), size(lags)))
      response = series(maximum_lag + 1:)
      do lag = 1, size(lags)
         do row = 1, size(response)
            predictors(row, lag) = series(maximum_lag + row - lags(lag))
         end do
      end do
      info = 0
   end subroutine lagged_regression_inputs

   pure subroutine update_elm_scales(model, response, predictors)
      !! Re-estimate ELM min-max transformations without changing weights.
      type(nnfor_elm_fast_t), intent(inout) :: model !! ELM receiving new scale parameters.
      real(dp), intent(in) :: response(:) !! New regression response.
      real(dp), intent(in) :: predictors(:, :) !! New regression predictors.
      type(nnfor_scaled_t) :: scaled
      integer :: predictor

      if (.not. model%scaled) return
      scaled = nnfor_linscale(response, -0.8_dp, 0.8_dp)
      model%response_scale = scaled%scale
      do predictor = 1, size(predictors, 2)
         scaled = nnfor_linscale(predictors(:, predictor), -0.8_dp, 0.8_dp)
         model%predictor_scales(predictor) = scaled%scale
      end do
   end subroutine update_elm_scales

   pure subroutine update_neural_network_scales(network, predictors)
      !! Re-estimate standardized-input parameters without changing MLP weights.
      type(neural_network_t), intent(inout) :: network !! Network receiving new scales.
      real(dp), intent(in) :: predictors(:, :) !! New regression predictors.
      real(dp) :: center, variance
      integer :: predictor, observations

      observations = size(predictors, 1)
      do predictor = 1, size(predictors, 2)
         center = sum(predictors(:, predictor))/real(observations, dp)
         variance = sum((predictors(:, predictor) - center)**2)/ &
            real(max(1, observations - 1), dp)
         network%input_mean(predictor) = center
         network%input_scale(predictor) = &
            sqrt(max(variance, epsilon(1.0_dp)))
      end do
   end subroutine update_neural_network_scales

   pure subroutine orthogonalize_columns(weights, success)
      !! Orthonormalize matrix columns with twice-modified Gram-Schmidt.
      real(dp), intent(inout) :: weights(:, :) !! Matrix replaced by orthonormal columns.
      logical, intent(out) :: success !! True when all columns retain numerical rank.
      real(dp), allocatable :: original(:, :), vector(:)
      real(dp) :: norm_value, threshold
      integer :: column, previous, pass

      success = .false.
      if (size(weights, 1) < size(weights, 2) .or. size(weights, 2) < 1) return
      original = weights
      do column = 1, size(weights, 2)
         vector = original(:, column)
         do pass = 1, 2
            do previous = 1, column - 1
               vector = vector - dot_product(weights(:, previous), vector)* &
                  weights(:, previous)
            end do
         end do
         norm_value = sqrt(dot_product(vector, vector))
         threshold = sqrt(epsilon(1.0_dp))*max(1.0_dp, &
            sqrt(dot_product(original(:, column), original(:, column))))
         if (norm_value <= threshold) then
            weights = original
            return
         end if
         weights(:, column) = vector/norm_value
      end do
      success = .true.
   end subroutine orthogonalize_columns

   pure function evaluate_elm_layers(layers, predictors) result(hidden)
      !! Propagate predictor rows through fixed sigmoid ELM hidden layers.
      type(nnfor_elm_layer_t), intent(in) :: layers(:) !! Successive hidden-layer weights.
      real(dp), intent(in) :: predictors(:, :) !! Predictor rows entering the first layer.
      real(dp), allocatable :: hidden(:, :)
      real(dp), allocatable :: augmented(:, :)
      integer :: layer

      hidden = predictors
      do layer = 1, size(layers)
         allocate(augmented(size(hidden, 1), size(hidden, 2) + 1))
         augmented(:, 1) = 1.0_dp
         augmented(:, 2:) = hidden
         hidden = nnfor_fast_sigmoid(matmul(augmented, &
            layers(layer)%input_weights))
         deallocate(augmented)
      end do
   end function evaluate_elm_layers

   pure function predict_elm_member(model, member, predictors) result(prediction)
      !! Predict one ELM member without applying the ensemble operator.
      type(nnfor_elm_fast_t), intent(in) :: model !! Fitted fast ELM ensemble.
      integer, intent(in) :: member !! Ensemble-member index.
      real(dp), intent(in) :: predictors(:, :) !! New predictor matrix.
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: x(:, :), augmented(:, :), hidden(:, :)
      integer :: predictor

      allocate(x(size(predictors, 1), model%predictor_count))
      do predictor = 1, model%predictor_count
         if (model%scaled) then
            x(:, predictor) = nnfor_apply_scale(predictors(:, predictor), &
               model%predictor_scales(predictor))
         else
            x(:, predictor) = predictors(:, predictor)
         end if
      end do
      allocate(augmented(size(predictors, 1), model%predictor_count + 1))
      augmented(:, 1) = 1.0_dp
      augmented(:, 2:) = x
      if (allocated(model%members(member)%layers)) then
         hidden = evaluate_elm_layers(model%members(member)%layers, x)
      else
         hidden = nnfor_fast_sigmoid(matmul(augmented, &
            model%members(member)%input_weights))
      end if
      prediction = model%members(member)%output_bias + &
         matmul(hidden, model%members(member)%output_weights)
      if (model%direct) prediction = prediction + &
         matmul(x, model%members(member)%direct_weights)
      if (model%scaled) prediction = nnfor_apply_scale(prediction, &
         model%response_scale, .true.)
   end function predict_elm_member

end module nnfor_mod
