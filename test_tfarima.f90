! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the tfarima translation.
program test_tfarima
   use kind_mod, only: dp
   use kfas_mod, only: ssm_model_t, kfs_filter_t, kfs_smoother_t, &
      make_local_level, make_local_linear_trend, kfs_filter, kfs_smooth
   use tfarima_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   type(tfarima_polynomial_division_t) :: division
   type(tfarima_transfer_t) :: transfer
   type(tfarima_forecast_t) :: forecast
   type(tfarima_ma_factor_t) :: ma_factor
   type(tfarima_lag_polynomial_t) :: lag_model, updated_lag_model
   type(tfarima_lag_polynomial_t) :: numerator_model, denominator_model
   type(tfarima_transfer_fit_t) :: transfer_fit
   type(tfarima_outlier_result_t) :: outliers
   type(tfarima_calendar_result_t) :: calendar
   type(tfarima_ucarima_component_t), allocatable :: components(:)
   type(tfarima_ucarima_model_t) :: ucarima
   type(tfarima_wk_filter_t) :: wk_filter
   type(tfarima_ucarima_decomposition_t) :: decomposition
   type(tfarima_bezout_t) :: bezout
   type(tfarima_partial_fraction_t) :: partial_fraction
   type(tfarima_ucarima_conversion_t) :: conversion
   type(tfarima_ucarima_ssm_t) :: ucarima_ssm
   type(tfarima_ucarima_smoothing_t) :: smoothing
   type(tfarima_ucarima_fit_t) :: ucarima_fit
   type(tfarima_transfer_spec_t), allocatable :: transfer_specs(:)
   type(tfarima_exact_transfer_fit_t) :: exact_transfer_fit
   type(tfarima_exact_transfer_forecast_t) :: exact_transfer_forecast
   type(tfarima_prewhitened_ccf_t) :: prewhitened
   type(tfarima_transfer_identification_t) :: identified_transfer
   type(tfarima_prewhitened_ccf_t), allocatable :: transfer_ccf(:)
   type(tfarima_transfer_selection_t) :: transfer_selection
   type(tfarima_transfer_diagnostics_t) :: transfer_diagnostics
   type(tfarima_transfer_signal_t) :: transfer_signal
   type(tfarima_transfer_simulation_t) :: transfer_simulation, repeated_simulation
   type(tfarima_exact_transfer_fit_t) :: simulation_fit
   type(tfarima_ucarima_component_t) :: transfer_noise
   type(tfarima_ucarima_component_t) :: algebra_model
   type(tfarima_arima_algebra_t) :: algebra
   type(tfarima_weights_t) :: psi_weights, pi_weights
   type(tfarima_leverrier_t) :: leverrier
   type(tfarima_ssm_form_t) :: ssm_form, switched_form
   type(tfarima_ssm_reduction_t) :: ssm_reduction
   type(tfarima_decomposition_basis_t) :: decomposition_basis
   type(tfarima_root_decomposition_t) :: root_decomposition
   type(tfarima_structural_ssm_t) :: structural_ssm
   type(tfarima_structural_initialization_t) :: structural_initialization
   type(tfarima_structural_filter_t) :: structural_filter, comparison_filter
   type(tfarima_structural_smoother_t) :: structural_smoother
   type(tfarima_structural_forecast_t) :: structural_forecast
   type(tfarima_band_cholesky_t) :: band_factor
   type(tfarima_reduced_likelihood_t) :: reduced_likelihood
   type(tfarima_reduced_likelihood_t) :: wrapper_likelihood
   type(kfs_filter_t) :: reference_filter
   type(kfs_smoother_t) :: reference_smoother
   type(ssm_model_t) :: kfas_model
   type(tfarima_ucarima_component_t), allocatable :: input_models(:)
   real(dp), allocatable :: values(:), matrix(:, :), future_matrix(:, :)
   real(dp), allocatable :: selection_matrix(:, :), selection_regressors(:, :)
   real(dp), allocatable :: simulation_innovations(:, :)
   real(dp), allocatable :: input(:), output(:), signal(:)
   real(dp), allocatable :: decomposition_series(:)
   real(dp) :: noise
   real(dp) :: expected_likelihood, expected_variance
   complex(dp) :: derivative
   integer :: i

   values = tfarima_polynomial_multiply([1.0_dp, -0.5_dp], &
      [1.0_dp, 0.0_dp, -0.2_dp])
   call check_close(values, [1.0_dp, -0.5_dp, -0.2_dp, 0.1_dp], 1.0e-14_dp, &
      'polynomial multiplication')
   division = tfarima_polynomial_divide(values, [1.0_dp, -0.5_dp])
   call check(division%info == 0, 'polynomial division status')
   call check_close(division%quotient, [1.0_dp, 0.0_dp, -0.2_dp], 1.0e-14_dp, &
      'polynomial quotient')
   call check_close(division%remainder, [0.0_dp], 1.0e-14_dp, &
      'polynomial remainder')
   values = tfarima_polynomial_gcd( &
      tfarima_polynomial_multiply([1.0_dp, -0.5_dp], [1.0_dp, 0.2_dp]), &
      tfarima_polynomial_multiply([1.0_dp, -0.5_dp], [1.0_dp, -0.3_dp]))
   call check_close(values, [1.0_dp, -0.5_dp], 1.0e-12_dp, 'polynomial GCD')
   bezout = tfarima_extended_polynomial_gcd([1.0_dp, -1.0_dp], &
      [1.0_dp, 1.0_dp], tolerance=1.0e-12_dp)
   call check(bezout%info == 0, 'extended polynomial GCD status')
   values = tfarima_polynomial_multiply([1.0_dp, -1.0_dp], bezout%u) + &
      tfarima_polynomial_multiply([1.0_dp, 1.0_dp], bezout%v)
   call check(abs(values(1) - bezout%gcd(1)) < 1.0e-12_dp .and. &
      maxval(abs(values(2:))) < 1.0e-12_dp, &
      'extended polynomial Bezout identity')
   matrix = reshape([1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp], [2, 2])
   partial_fraction = tfarima_partial_fractions([1.0_dp], matrix, [2, 2], &
      tolerance=1.0e-12_dp)
   call check(partial_fraction%info == 0, 'partial-fraction status')
   call check_close(partial_fraction%numerator(1, :), [0.5_dp, 0.5_dp], &
      1.0e-12_dp, 'partial-fraction numerators')
   call check_close(partial_fraction%reconstruction, [1.0_dp, 0.0_dp], &
      1.0e-12_dp, 'partial-fraction reconstruction')
   values = tfarima_polynomial_power([1.0_dp, -1.0_dp], 3)
   call check_close(values, [1.0_dp, -3.0_dp, 3.0_dp, -1.0_dp], 1.0e-14_dp, &
      'polynomial integer power')
   values = tfarima_polynomial_ratio([1.0_dp, 0.2_dp], [1.0_dp, -0.5_dp], 3)
   call check_close(values, [1.0_dp, 0.7_dp, 0.35_dp, 0.175_dp], 1.0e-14_dp, &
      'rational polynomial expansion')
   algebra_model = tfarima_ucarima_component([1.0_dp, -0.5_dp], [1.0_dp], &
      [1.0_dp, 0.2_dp], 2.0_dp)
   matrix = reshape([0.01_dp, 0.0_dp, 0.0_dp, 0.01_dp], [2, 2])
   psi_weights = tfarima_psi_weights(algebra_model, 3, &
      coefficient_covariance=matrix)
   call check(psi_weights%info == 0, 'PSI weight status')
   call check_close(psi_weights%weight, [1.0_dp, 0.7_dp, 0.35_dp, &
      0.175_dp], 1.0e-14_dp, 'PSI weights')
   call check_close(psi_weights%cumulative, [1.0_dp, 1.7_dp, 2.05_dp, &
      2.225_dp], 1.0e-14_dp, 'cumulative PSI weights')
   call check(psi_weights%standard_error(1) == 0.0_dp .and. &
      all(psi_weights%standard_error(2:) > 0.0_dp), &
      'PSI delta-method standard errors')
   call check(all(psi_weights%cumulative_standard_error >= 0.0_dp), &
      'cumulative PSI standard errors')
   pi_weights = tfarima_pi_weights(algebra_model, 3, &
      coefficient_covariance=matrix)
   call check(pi_weights%info == 0, 'PI weight status')
   call check_close(pi_weights%weight, [1.0_dp, -0.7_dp, 0.14_dp, &
      -0.028_dp], 1.0e-14_dp, 'PI weights')
   algebra_model = tfarima_ucarima_component([1.0_dp], [1.0_dp, -1.0_dp], &
      [1.0_dp], 1.0_dp)
   psi_weights = tfarima_psi_weights(algebra_model, 3, &
      include_difference=.true.)
   call check_close(psi_weights%weight, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      1.0e-14_dp, 'integrated PSI weights')
   pi_weights = tfarima_pi_weights(algebra_model, 3, &
      include_difference=.true.)
   call check_close(pi_weights%weight, [1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp], &
      1.0e-14_dp, 'integrated PI weights')
   matrix = reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [2, 2])
   leverrier = tfarima_leverrier_faddeev([1.0_dp, 0.0_dp], matrix)
   call check(leverrier%info == 0, 'Leverrier-Faddeev status')
   call check_close(leverrier%characteristic, [1.0_dp, -2.0_dp, 1.0_dp], &
      1.0e-14_dp, 'Leverrier-Faddeev characteristic polynomial')
   call check_close(leverrier%adjoint_loading(1, :), [1.0_dp, 0.0_dp], &
      1.0e-14_dp, 'Leverrier-Faddeev leading adjoint loading')
   call check_close(leverrier%adjoint_loading(2, :), [-1.0_dp, 1.0_dp], &
      1.0e-14_dp, 'Leverrier-Faddeev trailing adjoint loading')
   ssm_form = tfarima_ssm_form([1.0_dp], reshape([0.5_dp], [1, 1]), &
      reshape([2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], [2, 2]), .true.)
   switched_form = tfarima_switch_ssm_form(ssm_form)
   call check(switched_form%info == 0 .and. &
      .not. switched_form%state_noise_contemporaneous, &
      'state-space lag-form switch status')
   call check_close(switched_form%observation_loading, [0.5_dp], 1.0e-14_dp, &
      'state-space lag-form loading')
   call check_close(switched_form%disturbance_covariance(1, :), &
      [5.0_dp, 3.0_dp], 1.0e-14_dp, 'state-space lag-form covariance')
   switched_form = tfarima_switch_ssm_form(switched_form)
   call check_close(switched_form%observation_loading, [1.0_dp], 1.0e-14_dp, &
      'state-space form-switch loading round trip')
   call check_close(switched_form%disturbance_covariance(1, :), &
      [2.0_dp, 0.0_dp], 1.0e-14_dp, &
      'state-space form-switch covariance round trip')
   kfas_model = make_local_level([0.0_dp, 0.0_dp], 2.0_dp, 3.0_dp)
   ssm_reduction = tfarima_ssm_to_arima(kfas_model, tolerance=1.0e-10_dp)
   call check(ssm_reduction%info == 0 .and. ssm_reduction%factor%converged, &
      'local-level ARIMA reduction status')
   call check_close(ssm_reduction%characteristic, [1.0_dp, -1.0_dp], &
      1.0e-14_dp, 'local-level characteristic polynomial')
   call check_close(ssm_reduction%target_autocovariance, [7.0_dp, -2.0_dp], &
      1.0e-12_dp, 'local-level reduced numerator covariance')
   call check_close(ssm_reduction%model%ar_polynomial, [1.0_dp], 1.0e-12_dp, &
      'local-level reduced AR polynomial')
   call check_close(ssm_reduction%model%difference_polynomial, &
      [1.0_dp, -1.0_dp], 1.0e-12_dp, &
      'local-level reduced difference polynomial')
   kfas_model%transition = 0.5_dp
   kfas_model%h = 0.0_dp
   kfas_model%q = 2.0_dp
   ssm_reduction = tfarima_ssm_to_arima(kfas_model, tolerance=1.0e-10_dp)
   call check(ssm_reduction%info == 0, 'stationary-state ARIMA reduction status')
   call check_close(ssm_reduction%model%ar_polynomial, [1.0_dp, -0.5_dp], &
      1.0e-12_dp, 'stationary-state reduced AR polynomial')
   call check_close(ssm_reduction%model%difference_polynomial, [1.0_dp], &
      1.0e-12_dp, 'stationary-state reduced difference polynomial')
   call check_close(ssm_reduction%model%ma_polynomial, [1.0_dp], 1.0e-12_dp, &
      'stationary-state reduced MA polynomial')
   call check(abs(ssm_reduction%model%innovation_variance - 2.0_dp) < 1.0e-10_dp, &
      'stationary-state reduced innovation variance')
   kfas_model = make_local_linear_trend([0.0_dp, 0.0_dp], 1.0_dp, 2.0_dp, &
      3.0_dp)
   ssm_reduction = tfarima_ssm_to_arima(kfas_model)
   call check(ssm_reduction%info == 0, &
      'local-linear-trend ARIMA reduction status')
   call check_close(ssm_reduction%model%difference_polynomial, &
      [1.0_dp, -2.0_dp, 1.0_dp], 1.0e-6_dp, &
      'local-linear-trend reduced difference polynomial')
   decomposition_basis = tfarima_root_decomposition_basis( &
      [1.0_dp, -0.8_dp], [1.0_dp, -1.0_dp])
   call check(decomposition_basis%info == 0, &
      'root decomposition basis status')
   call check(all(shape(decomposition_basis%basis) == [3, 2]), &
      'root decomposition basis shape')
   call check_close(decomposition_basis%basis(:, 1), &
      [0.8_dp, 0.64_dp, 0.512_dp], 1.0e-12_dp, &
      'stationary exponential root basis')
   call check_close(decomposition_basis%basis(:, 2), &
      [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-12_dp, 'unit trend root basis')
   call check(sum(decomposition_basis%classification( &
      tfarima_component_exponential, :)) == 1.0_dp .and. &
      sum(decomposition_basis%classification( &
      tfarima_component_trend, :)) == 1.0_dp, &
      'root component classification')
   decomposition_basis = tfarima_root_decomposition_basis([1.0_dp], &
      [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp])
   call check(decomposition_basis%info == 0 .and. &
      sum(decomposition_basis%classification( &
      tfarima_component_trend, :)) == 1.0_dp .and. &
      sum(decomposition_basis%classification( &
      tfarima_component_seasonal, :)) == 3.0_dp, &
      'seasonal unit-root classification')
   decomposition_series = [1.0_dp, 1.4_dp, 1.1_dp, 1.8_dp, 1.6_dp, 2.2_dp, &
      1.9_dp, 2.5_dp]
   root_decomposition = tfarima_root_decompose(decomposition_series, [1.0_dp], &
      [1.0_dp, -1.0_dp], [1.0_dp, 0.2_dp], 1.0_dp, &
      method=tfarima_decomposition_forecast)
   call check(root_decomposition%info == 0, &
      'forecast root decomposition status')
   call check_close(root_decomposition%reconstruction, decomposition_series, &
      1.0e-10_dp, &
      'forecast root decomposition reconstruction')
   root_decomposition = tfarima_root_decompose(decomposition_series, [1.0_dp], &
      [1.0_dp, -1.0_dp], [1.0_dp, 0.2_dp], 1.0_dp, &
      method=tfarima_decomposition_backcast)
   call check(root_decomposition%info == 0, &
      'backcast root decomposition status')
   call check_close(root_decomposition%reconstruction, decomposition_series, &
      1.0e-10_dp, &
      'backcast root decomposition reconstruction')
   root_decomposition = tfarima_root_decompose(decomposition_series, [1.0_dp], &
      [1.0_dp, -1.0_dp], [1.0_dp, 0.2_dp], 1.0_dp)
   call check(root_decomposition%info == 0, 'mixed root decomposition status')
   call check_close(root_decomposition%reconstruction, decomposition_series, &
      1.0e-10_dp, &
      'mixed root decomposition reconstruction')
   call check_close(root_decomposition%seasonally_adjusted, decomposition_series, &
      1.0e-10_dp, &
      'nonseasonal adjusted-series identity')
   decomposition_series = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 11.0_dp, &
      21.0_dp, 31.0_dp, 41.0_dp, 12.0_dp, 22.0_dp, 32.0_dp, 42.0_dp]
   root_decomposition = tfarima_root_decompose(decomposition_series, [1.0_dp], &
      [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [1.0_dp], 1.0_dp)
   call check(root_decomposition%info == 0 .and. &
      maxval(abs(root_decomposition%component(:, &
      tfarima_component_seasonal))) > 1.0_dp, &
      'seasonal root decomposition status')
   call check_close(root_decomposition%reconstruction, decomposition_series, &
      1.0e-9_dp, 'seasonal root decomposition reconstruction')
   call check_close(root_decomposition%seasonally_adjusted, &
      decomposition_series - root_decomposition%component(:, &
      tfarima_component_seasonal), 1.0e-12_dp, &
      'additive root seasonal adjustment')
   root_decomposition = tfarima_root_decompose(decomposition_series, [1.0_dp], &
      [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [1.0_dp], 1.0_dp, &
      log_transform=.true.)
   call check(root_decomposition%info == 0, &
      'multiplicative root decomposition status')
   call check_close(root_decomposition%reconstruction, log(decomposition_series), &
      1.0e-9_dp, 'multiplicative root decomposition reconstruction')
   call check_close(root_decomposition%seasonally_adjusted, &
      exp(log(decomposition_series) - root_decomposition%component(:, &
      tfarima_component_seasonal)), 1.0e-12_dp, &
      'multiplicative root seasonal adjustment')
   algebra_model = tfarima_ucarima_component([1.0_dp], [1.0_dp, -1.0_dp], &
      [1.0_dp, 0.2_dp], 2.0_dp)
   structural_ssm = tfarima_arima_to_structural_ssm(algebra_model, &
      multiple_sources=.false., contemporaneous=.false.)
   call check(structural_ssm%info == 0 .and. structural_ssm%admissible .and. &
      .not. structural_ssm%multiple_sources, &
      'lagged SSOE structural conversion status')
   call check_close(structural_ssm%target_autocovariance, [2.08_dp, 0.4_dp], &
      1.0e-12_dp, 'SSOE target autocovariance')
   call check_close(structural_ssm%fitted_autocovariance, [2.08_dp, 0.4_dp], &
      1.0e-10_dp, 'SSOE fitted autocovariance')
   call check(structural_ssm%reconstruction_error < 1.0e-10_dp, &
      'SSOE covariance reconstruction')
   ssm_form = structural_ssm%form
   ssm_reduction = tfarima_reduce_ssm(ssm_form, tolerance=1.0e-10_dp)
   call check(ssm_reduction%info == 0, 'SSOE structural round-trip status')
   call check_close(ssm_reduction%target_autocovariance, [2.08_dp, 0.4_dp], &
      1.0e-10_dp, 'SSOE structural round-trip covariance')
   structural_ssm = tfarima_arima_to_structural_ssm(algebra_model, &
      multiple_sources=.false., contemporaneous=.true.)
   switched_form = tfarima_switch_ssm_form(structural_ssm%form)
   call check(switched_form%info == 0, &
      'contemporaneous SSOE form switch status')
   call check_close(switched_form%observation_loading, &
      ssm_form%observation_loading, 1.0e-12_dp, &
      'SSOE form-switch loading equivalence')
   call check_close(switched_form%disturbance_covariance(1, :), &
      ssm_form%disturbance_covariance(1, :), 1.0e-10_dp, &
      'SSOE form-switch covariance equivalence')
   algebra_model = tfarima_ucarima_component([1.0_dp], [1.0_dp, -1.0_dp], &
      [1.0_dp, -0.2_dp], 2.0_dp)
   structural_ssm = tfarima_arima_to_structural_ssm(algebra_model, &
      multiple_sources=.true., contemporaneous=.false.)
   call check(structural_ssm%info == 0 .and. structural_ssm%admissible .and. &
      .not. structural_ssm%used_nonnegative_fit, &
      'admissible MSOE structural conversion status')
   call check_close(structural_ssm%disturbance_variance, [0.4_dp, 1.28_dp], &
      1.0e-10_dp, 'admissible MSOE disturbance variances')
   call check(structural_ssm%reconstruction_error < 1.0e-10_dp, &
      'admissible MSOE covariance reconstruction')
   algebra_model = tfarima_ucarima_component([1.0_dp], [1.0_dp, -1.0_dp], &
      [1.0_dp, 0.2_dp], 2.0_dp)
   structural_ssm = tfarima_arima_to_structural_ssm(algebra_model, &
      multiple_sources=.true., contemporaneous=.false.)
   call check(structural_ssm%info == 0 .and. .not. structural_ssm%admissible .and. &
      structural_ssm%used_nonnegative_fit .and. &
      all(structural_ssm%disturbance_variance >= 0.0_dp), &
      'nonnegative MSOE structural fallback')
   structural_ssm = tfarima_arima_to_structural_ssm(algebra_model, &
      multiple_sources=.true., contemporaneous=.false., &
      grouping=reshape([1.0_dp, 1.0_dp], [2, 1]))
   call check(structural_ssm%info == 0 .and. &
      abs(structural_ssm%disturbance_variance(1) - &
      structural_ssm%disturbance_variance(2)) < 1.0e-12_dp, &
      'grouped MSOE disturbance variances')
   decomposition_series = [1.0_dp, 1.5_dp, 0.7_dp, 2.0_dp, 1.8_dp]
   ssm_form = tfarima_ssm_form([1.0_dp], reshape([1.0_dp], [1, 1]), &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp], [2, 2]), .false.)
   structural_filter = tfarima_structural_filter(decomposition_series, ssm_form, &
      initial_state=[0.0_dp], initial_covariance=reshape([2.0_dp], [1, 1]))
   call check(structural_filter%info == 0 .and. &
      ieee_is_finite(structural_filter%log_likelihood), &
      'joint-disturbance structural filter status')
   kfas_model = make_local_level(decomposition_series, 1.0_dp, 0.5_dp, &
      a1=0.0_dp, p1=2.0_dp)
   reference_filter = kfs_filter(kfas_model)
   call check_close(structural_filter%innovation, &
      reference_filter%innovation(:, 1), 1.0e-12_dp, &
      'structural and KFAS filter innovations')
   call check_close(structural_filter%innovation_variance, &
      [(reference_filter%innovation_cov(1, 1, i), &
      i=1,size(decomposition_series))], 1.0e-12_dp, &
      'structural and KFAS innovation variances')
   call check_close(structural_filter%observation_state(1, :), &
      reference_filter%a_pred(1, :), 1.0e-12_dp, &
      'structural and KFAS observation states')
   call check_close(structural_filter%filtered_state(1, &
      :size(decomposition_series) - 1), reference_filter%a_pred(1, 2:), &
      1.0e-12_dp, 'structural next-state filtering')
   structural_smoother = tfarima_structural_smooth(decomposition_series, ssm_form, &
      initial_state=[0.0_dp], initial_covariance=reshape([2.0_dp], [1, 1]))
   reference_smoother = kfs_smooth(kfas_model, reference_filter)
   call check(structural_smoother%info == 0, &
      'joint-disturbance structural smoother status')
   call check_close(structural_smoother%state(1, :), &
      reference_smoother%state(1, :), 1.0e-10_dp, &
      'structural and KFAS smoothed states')
   structural_initialization = tfarima_structural_initialize( &
      decomposition_series, ssm_form)
   call check(structural_initialization%info == 0 .and. &
      all(ieee_is_finite(structural_initialization%state)) .and. &
      all(ieee_is_finite(structural_initialization%covariance)), &
      'structural GLS initialization')
   structural_forecast = tfarima_structural_forecast(structural_filter, &
      ssm_form, 3)
   call check(structural_forecast%info == 0 .and. &
      all(structural_forecast%variance > 0.0_dp) .and. &
      all(shape(structural_forecast%state) == [1, 3]), &
      'structural forecast moments')
   comparison_filter = structural_filter
   matrix = reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [5, 1])
   structural_filter = tfarima_structural_filter( &
      exp(decomposition_series + 0.2_dp), ssm_form, initial_state=[0.0_dp], &
      initial_covariance=reshape([2.0_dp], [1, 1]), log_transform=.true., &
      regressors=matrix, regression_coefficients=[0.2_dp])
   call check_close(structural_filter%innovation, comparison_filter%innovation, &
      1.0e-12_dp, 'structural log-regression preprocessing')
   structural_forecast = tfarima_structural_forecast(structural_filter, &
      ssm_form, 2, future_regressors=reshape([1.0_dp, 1.0_dp], [2, 1]), &
      regression_coefficients=[0.2_dp], log_transform=.true.)
   call check(structural_forecast%info == 0 .and. &
      all(structural_forecast%mean > 0.0_dp) .and. &
      all(structural_forecast%variance > 0.0_dp), &
      'structural log-regression forecasts')
   ssm_form = tfarima_ssm_form([1.0_dp], reshape([1.0_dp], [1, 1]), &
      reshape([2.0_dp, 2.4_dp, 2.4_dp, 2.88_dp], [2, 2]), .false.)
   structural_filter = tfarima_structural_filter(decomposition_series, &
      form=ssm_form, initial_state=[0.0_dp], &
      initial_covariance=reshape([1.0_dp], [1, 1]))
   call check(structural_filter%info == 0 .and. &
      all(structural_filter%innovation_variance > 0.0_dp), &
      'correlated SSOE structural filtering')

   matrix = reshape([4.0_dp, 3.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, 0.5_dp], &
      [3, 2])
   band_factor = tfarima_band_cholesky(matrix)
   call check(band_factor%info == 0, 'band Cholesky status')
   call check_close(band_factor%factor(:, 1), &
      [2.0_dp, sqrt(2.75_dp), sqrt(2.0_dp - 0.25_dp/2.75_dp)], &
      1.0e-12_dp, 'band Cholesky diagonal')
   values = tfarima_band_forward_solve(band_factor%factor, &
      [2.0_dp, 0.5_dp + 2.0_dp*sqrt(2.75_dp), &
      1.0_dp/sqrt(2.75_dp) + &
      3.0_dp*sqrt(2.0_dp - 0.25_dp/2.75_dp)])
   call check_close(values, [1.0_dp, 2.0_dp, 3.0_dp], 1.0e-12_dp, &
      'band forward solution')

   decomposition_series = [1.0_dp, 0.2_dp, -0.3_dp, 0.7_dp, 0.1_dp]
   reduced_likelihood = tfarima_reduced_likelihood(decomposition_series, &
      [1.0_dp, -0.4_dp], reshape([1.0_dp], [1, 1]), &
      reshape([1.0_dp], [1, 1]))
   expected_variance = (0.84_dp + 0.04_dp + 0.1444_dp + &
      0.6724_dp + 0.0324_dp)/5.0_dp
   expected_likelihood = -2.5_dp*(1.0_dp + &
      log(2.0_dp*acos(-1.0_dp)) + log(expected_variance)) + &
      0.5_dp*log(0.84_dp)
   call check(reduced_likelihood%info == 0 .and. &
      reduced_likelihood%bandwidth == 0, 'exact reduced AR likelihood status')
   call check(abs(reduced_likelihood%innovation_variance - &
      expected_variance) < 1.0e-12_dp .and. &
      abs(reduced_likelihood%log_likelihood - expected_likelihood) < &
      1.0e-12_dp, 'exact reduced AR likelihood values')
   call check_close(reduced_likelihood%residuals, &
      [sqrt(0.84_dp), -0.2_dp, -0.38_dp, 0.82_dp, -0.18_dp], &
      1.0e-12_dp, 'exact reduced AR residuals')

   reduced_likelihood = tfarima_reduced_likelihood(decomposition_series, &
      [1.0_dp], reshape([1.0_dp, 0.5_dp, 0.2_dp, -0.1_dp], [2, 2]), &
      reshape([2.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2, 2]))
   call check(reduced_likelihood%info == 0 .and. &
      reduced_likelihood%bandwidth == 1, &
      'correlated-disturbance reduced MA status')
   call check_close(reduced_likelihood%covariance_band(:, 1), &
      [2.628_dp, 2.628_dp, 2.628_dp, 2.628_dp, 2.628_dp], &
      1.0e-12_dp, 'correlated-disturbance reduced MA variance')
   call check_close(reduced_likelihood%covariance_band(2:, 2), &
      [0.35_dp, 0.35_dp, 0.35_dp, 0.35_dp], 1.0e-12_dp, &
      'correlated-disturbance reduced MA covariance')

   allocate(components(2))
   components(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 2.0_dp)
   components(2) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 1.0_dp)
   ucarima = tfarima_build_ucarima(components)
   wrapper_likelihood = tfarima_ucarima_reduced_likelihood( &
      decomposition_series, ucarima)
   reduced_likelihood = tfarima_reduced_likelihood(decomposition_series, &
      [1.0_dp], reshape([1.0_dp, 1.0_dp], [2, 1]), &
      reshape([2.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))
   call check(wrapper_likelihood%info == 0 .and. &
      abs(wrapper_likelihood%log_likelihood - &
      reduced_likelihood%log_likelihood) < 1.0e-12_dp, &
      'UCARIMA reduced likelihood wrapper')
   deallocate(components)

   ssm_form = tfarima_ssm_form([1.0_dp], reshape([1.0_dp], [1, 1]), &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp], [2, 2]), .false.)
   wrapper_likelihood = tfarima_structural_reduced_likelihood( &
      [1.0_dp, 1.5_dp, 0.7_dp, 2.0_dp, 1.8_dp, 2.2_dp], ssm_form, &
      [1.0_dp, -1.0_dp])
   call check(wrapper_likelihood%info == 0 .and. &
      ieee_is_finite(wrapper_likelihood%log_likelihood) .and. &
      wrapper_likelihood%innovation_variance > 0.0_dp, &
      'structural reduced likelihood wrapper')
   derivative = tfarima_polynomial_derivative([1.0_dp, 2.0_dp, 3.0_dp], &
      cmplx(2.0_dp, 0.0_dp, dp), 1)
   call check(abs(derivative - cmplx(14.0_dp, 0.0_dp, dp)) < 1.0e-14_dp, &
      'polynomial derivative evaluation')

   transfer = tfarima_transfer([2.0_dp], [1.0_dp, -0.5_dp], delay=1)
   values = tfarima_filter([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], transfer)
   call check_close(values, [0.0_dp, 2.0_dp, 1.0_dp, 0.5_dp, 0.25_dp], &
      1.0e-14_dp, 'delayed rational filter')
   values = tfarima_impulse_response(transfer, 4)
   call check_close(values, [0.0_dp, 2.0_dp, 1.0_dp, 0.5_dp, 0.25_dp], &
      1.0e-14_dp, 'transfer impulse response')
   values = tfarima_impulse_response(transfer, 4, cumulative=.true.)
   call check_close(values, [0.0_dp, 2.0_dp, 3.0_dp, 3.5_dp, 3.75_dp], &
      1.0e-14_dp, 'transfer step response')

   values = tfarima_difference([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [1.0_dp, -1.0_dp])
   call check_close(values, [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-14_dp, &
      'lag-polynomial differencing')
   values = tfarima_intervention(5, 3, tfarima_intervention_ramp)
   call check_close(values, [0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], &
      1.0e-14_dp, 'ramp intervention')
   matrix = tfarima_seasonal_dummies(4, 4, 1, 1)
   call check(all(shape(matrix) == [4, 3]), 'seasonal dummy shape')
   call check_close(matrix(1, :), [-1.0_dp, -1.0_dp, -1.0_dp], 1.0e-14_dp, &
      'seasonal dummy reference row')
   call check_close(matrix(2, :), [1.0_dp, 0.0_dp, 0.0_dp], 1.0e-14_dp, &
      'seasonal dummy active row')
   matrix = tfarima_harmonic_regressors(2, 4, 1)
   call check_close(matrix(1, :), [0.0_dp, 1.0_dp, -1.0_dp], 1.0e-14_dp, &
      'seasonal harmonic first row')
   call check_close(matrix(2, :), [-1.0_dp, 0.0_dp, 1.0_dp], 1.0e-14_dp, &
      'seasonal harmonic second row')

   forecast = tfarima_arima_forecast([1.0_dp, 0.5_dp, 0.25_dp, 0.125_dp], &
      [1.0_dp, -0.5_dp], [1.0_dp], [1.0_dp], 1.0_dp, 2, &
      levels=[0.95_dp])
   call check(forecast%info == 0, 'stationary AR forecast status')
   call check_close(forecast%mean, [0.0625_dp, 0.03125_dp], 1.0e-14_dp, &
      'stationary AR point forecast')
   call check_close(forecast%variance, [1.0_dp, 1.25_dp], 1.0e-14_dp, &
      'stationary AR forecast variance')
   call check(all(forecast%lower(:, 1) < forecast%mean) .and. &
      all(forecast%upper(:, 1) > forecast%mean), 'normal forecast intervals')
   forecast = tfarima_arima_forecast([1.0_dp, 2.0_dp, 3.0_dp], [1.0_dp], &
      [1.0_dp, -1.0_dp], [1.0_dp], 1.0_dp, 2, mean_value=1.0_dp)
   call check_close(forecast%mean, [4.0_dp, 5.0_dp], 1.0e-14_dp, &
      'integrated point forecast')
   call check_close(forecast%variance, [1.0_dp, 2.0_dp], 1.0e-14_dp, &
      'integrated forecast variance')

   values = tfarima_autocovariance([1.0_dp, -0.5_dp], [1.0_dp], 2.0_dp, 2)
   call check_close(values, [8.0_dp/3.0_dp, 4.0_dp/3.0_dp, 2.0_dp/3.0_dp], &
      1.0e-12_dp, 'theoretical AR autocovariance')
   values = tfarima_partial_autocorrelation([1.0_dp, -0.5_dp], [1.0_dp], 3)
   call check_close(values, [0.5_dp, 0.0_dp, 0.0_dp], 1.0e-12_dp, &
      'theoretical AR partial autocorrelation')
   ma_factor = tfarima_autocovariance_to_ma([2.5_dp, 1.0_dp], &
      tolerance=1.0e-12_dp)
   call check(ma_factor%info == 0 .and. ma_factor%converged, &
      'autocovariance-to-MA convergence')
   call check_close(ma_factor%coefficients, &
      [sqrt(2.0_dp), 0.5_dp*sqrt(2.0_dp)], 1.0e-10_dp, &
      'autocovariance-to-MA coefficients')
   ma_factor = tfarima_cramer_wold_factor([2.5_dp, 1.0_dp], &
      tfarima_factor_roots, tolerance=1.0e-10_dp)
   call check(ma_factor%info == 0 .and. ma_factor%converged .and. &
      ma_factor%method == tfarima_factor_roots .and. &
      ma_factor%residual_norm < 1.0e-10_dp, &
      'root Cramer-Wold factorization')
   ma_factor = tfarima_cramer_wold_factor([2.5_dp, 1.0_dp], &
      tfarima_factor_bauer, tolerance=1.0e-10_dp)
   call check(ma_factor%info == 0 .and. ma_factor%converged .and. &
      ma_factor%residual_norm < 1.0e-8_dp, 'Bauer Cramer-Wold factorization')
   ma_factor = tfarima_cramer_wold_factor([2.5_dp, 1.0_dp], &
      tfarima_factor_laurie, tolerance=1.0e-10_dp)
   call check(ma_factor%info == 0 .and. ma_factor%converged .and. &
      ma_factor%residual_norm < 1.0e-8_dp, 'Laurie Cramer-Wold factorization')
   ma_factor = tfarima_cramer_wold_factor([2.5_dp, 1.0_dp], &
      tfarima_factor_wilson, tolerance=1.0e-10_dp)
   call check(ma_factor%info == 0 .and. ma_factor%converged .and. &
      ma_factor%residual_norm < 1.0e-8_dp, 'Wilson Cramer-Wold factorization')
   ma_factor = tfarima_cramer_wold_factor( &
      [4.8961_dp, -3.258_dp, 0.81_dp], tfarima_factor_best, &
      tolerance=1.0e-10_dp)
   call check(ma_factor%info == 0 .and. ma_factor%converged .and. &
      ma_factor%residual_norm < 1.0e-8_dp, &
      'best near-repeated-root Cramer-Wold factorization')
   values = tfarima_palindromic_to_wold([6.0_dp, -4.0_dp, 1.0_dp])
   values = tfarima_wold_to_palindromic(values)
   call check_close(values, [6.0_dp, -4.0_dp, 1.0_dp], 1.0e-12_dp, &
      'palindromic and Wold conversion round trip')

   allocate(components(2))
   components(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 1.0_dp)
   components(2) = tfarima_ucarima_component([1.0_dp, -0.5_dp], [1.0_dp], &
      [1.0_dp], 2.0_dp)
   ucarima = tfarima_build_ucarima(components, tolerance=1.0e-12_dp)
   call check(ucarima%info == 0 .and. ucarima%factor_converged, &
      'UCARIMA Cramer-Wold factorization status')
   call check(ucarima%factor_method >= tfarima_factor_newton .and. &
      ucarima%factor_method <= tfarima_factor_wilson .and. &
      ucarima%factor_residual_norm < 1.0e-8_dp, &
      'UCARIMA automatic factor selection diagnostics')
   call check_close(ucarima%denominator, [1.0_dp, -0.5_dp], 1.0e-14_dp, &
      'UCARIMA common denominator')
   call check_close(ucarima%numerator_autocovariance, [3.25_dp, -0.5_dp], &
      1.0e-14_dp, 'UCARIMA numerator autocovariance')
   call check_close(ucarima%component_numerator(:, 1), [1.0_dp, -0.5_dp], &
      1.0e-14_dp, 'UCARIMA lifted white-noise numerator')
   call check_close(ucarima%component_numerator(:, 2), [1.0_dp, 0.0_dp], &
      1.0e-14_dp, 'UCARIMA lifted AR numerator')
   call check(abs(ucarima%innovation_variance* &
      sum(ucarima%ma_polynomial**2) - 3.25_dp) < 1.0e-10_dp .and. &
      abs(ucarima%innovation_variance*product(ucarima%ma_polynomial) + &
      0.5_dp) < 1.0e-10_dp, 'UCARIMA aggregate MA moments')
   algebra = tfarima_combine_arima(components(1), components(2), &
      tolerance=1.0e-10_dp)
   call check(algebra%info == 0 .and. algebra%factor%converged, &
      'ARIMA spectral addition status')
   call check_close(algebra%target_autocovariance, [3.25_dp, -0.5_dp], &
      1.0e-12_dp, 'ARIMA spectral addition target')
   call check(algebra%reconstruction_error < 1.0e-8_dp, &
      'ARIMA spectral addition reconstruction')
   algebra_model = algebra%model
   components(1) = tfarima_ucarima_component([1.0_dp, -0.5_dp], [1.0_dp], &
      [1.0_dp], 1.0_dp)
   components(2) = tfarima_ucarima_component([1.0_dp, -0.5_dp], [1.0_dp], &
      [1.0_dp], 2.0_dp)
   algebra = tfarima_combine_arima(components(1), components(2), &
      tolerance=1.0e-10_dp)
   call check(algebra%info == 0, 'shared-denominator ARIMA addition status')
   call check_close(algebra%model%ar_polynomial, [1.0_dp, -0.5_dp], &
      1.0e-12_dp, 'shared-denominator ARIMA polynomial')
   call check_close(algebra%model%ma_polynomial, [1.0_dp], 1.0e-12_dp, &
      'shared-denominator aggregate MA polynomial')
   call check(abs(algebra%model%innovation_variance - 3.0_dp) < 1.0e-10_dp, &
      'shared-denominator aggregate variance')
   algebra = tfarima_combine_arima(algebra_model, components(2), &
      subtract_second=.true., tolerance=1.0e-8_dp)
   call check(algebra%info == 0, 'ARIMA spectral subtraction status')
   call check_close(algebra%model%ar_polynomial, [1.0_dp], 1.0e-7_dp, &
      'ARIMA subtraction factor cancellation')
   call check_close(algebra%model%ma_polynomial, [1.0_dp], 1.0e-7_dp, &
      'ARIMA subtraction white-noise MA')
   call check(abs(algebra%model%innovation_variance - 1.0_dp) < 1.0e-7_dp, &
      'ARIMA subtraction white-noise variance')
   components(1) = tfarima_ucarima_component([2.0_dp], [3.0_dp], &
      [4.0_dp], 5.0_dp)
   call check_close([components(1)%innovation_variance], [20.0_dp/9.0_dp], &
      1.0e-14_dp, 'UCARIMA component normalization')
   components(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 1.0_dp)
   components(2) = tfarima_ucarima_component([1.0_dp], [1.0_dp, -1.0_dp], &
      [1.0_dp], 2.0_dp)
   ucarima = tfarima_build_ucarima(components, tolerance=1.0e-12_dp)
   call check(ucarima%info == 0, 'UCARIMA integrated component status')
   call check_close(ucarima%denominator, [1.0_dp, -1.0_dp], 1.0e-14_dp, &
      'UCARIMA integrated common denominator')
   components(1) = tfarima_ucarima_component([1.0_dp, -0.5_dp], [1.0_dp], &
      [1.0_dp], 1.0_dp)
   components(2) = tfarima_ucarima_component([1.0_dp, -0.5_dp], [1.0_dp], &
      [1.0_dp], 2.0_dp)
   ucarima = tfarima_build_ucarima(components)
   call check(ucarima%info == 2, 'UCARIMA common-factor rejection')
   deallocate(components)

   allocate(components(2))
   components(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 1.0_dp)
   components(2) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 3.0_dp)
   ucarima = tfarima_build_ucarima(components, tolerance=1.0e-12_dp)
   wk_filter = tfarima_wiener_kolmogorov_filter(ucarima, 1, 4, &
      tolerance=1.0e-12_dp)
   call check(wk_filter%info == 0, 'white-noise WK filter status')
   call check_close(wk_filter%symmetric_numerator, [0.25_dp], 1.0e-14_dp, &
      'white-noise WK rational numerator')
   call check_close(wk_filter%weights, &
      [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp], 1.0e-14_dp, 'white-noise WK weights')
   values = [(sin(0.2_dp*real(i, dp)), i=1, 30)]
   decomposition = tfarima_ucarima_decompose(values, ucarima, max_lag=4, &
      tolerance=1.0e-12_dp)
   call check(decomposition%info == 0, 'white-noise UCARIMA decomposition status')
   call check_close(decomposition%component(:, 1), 0.25_dp*values, &
      1.0e-13_dp, 'white-noise first component estimate')
   call check_close(decomposition%component(:, 2), 0.75_dp*values, &
      1.0e-13_dp, 'white-noise second component estimate')
   call check_close(decomposition%reconstruction, values, 1.0e-13_dp, &
      'white-noise UCARIMA reconstruction')
   ucarima_ssm = tfarima_ucarima_state_space(values, ucarima, &
      tolerance=1.0e-12_dp)
   call check(ucarima_ssm%info == 0 .and. .not. ucarima_ssm%diffuse .and. &
      all(shape(ucarima_ssm%component_loading) == [2, 2]), &
      'stationary UCARIMA state-space construction')
   smoothing = tfarima_ucarima_smooth(values, ucarima, forecast_horizon=3, &
      levels=[0.95_dp], tolerance=1.0e-12_dp)
   call check(smoothing%info == 0, 'stationary UCARIMA smoothing status')
   call check_close(smoothing%component(:, 1), 0.25_dp*values, &
      1.0e-12_dp, 'state-space first white-noise component')
   call check_close(smoothing%component(:, 2), 0.75_dp*values, &
      1.0e-12_dp, 'state-space second white-noise component')
   call check_close(smoothing%variance(:, 1), &
      [(0.75_dp, i=1,size(values))], 1.0e-12_dp, &
      'state-space white-noise component variance')
   call check_close(smoothing%forecast(:, 1), [0.0_dp, 0.0_dp, 0.0_dp], &
      1.0e-12_dp, 'state-space white-noise forecast')
   call check_close(smoothing%forecast_variance(:, 1), &
      [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-12_dp, &
      'state-space white-noise forecast variance')
   call check(all(smoothing%lower(:, :, 1) < smoothing%forecast) .and. &
      all(smoothing%upper(:, :, 1) > smoothing%forecast), &
      'state-space component forecast intervals')
   components(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 1.0_dp)
   components(2) = tfarima_ucarima_component([1.0_dp, -0.5_dp], [1.0_dp], &
      [1.0_dp], 2.0_dp)
   ucarima = tfarima_build_ucarima(components, tolerance=1.0e-12_dp)
   wk_filter = tfarima_wiener_kolmogorov_filter(ucarima, 1, 30, &
      tolerance=1.0e-12_dp)
   call check(wk_filter%info == 0, 'dynamic first WK filter status')
   values = wk_filter%weights
   wk_filter = tfarima_wiener_kolmogorov_filter(ucarima, 2, 30, &
      tolerance=1.0e-12_dp)
   call check(wk_filter%info == 0, 'dynamic second WK filter status')
   values = values + wk_filter%weights
   call check(abs(values(31) - 1.0_dp) < 1.0e-10_dp .and. &
      maxval(abs(values(:30))) < 1.0e-10_dp .and. &
      maxval(abs(values(32:))) < 1.0e-10_dp, &
      'dynamic WK component filters sum to identity')
   matrix = reshape([1.0_dp, -0.5_dp], [2, 1])
   conversion = tfarima_arima_to_ucarima([1.0_dp, -0.5_dp], [1.0_dp], &
      ucarima%ma_polynomial, ucarima%innovation_variance, matrix, [2], &
      [.true.], canonical=.true., tolerance=1.0e-10_dp)
   call check(conversion%info == 0 .and. conversion%admissible, &
      'aggregate ARIMA to UCARIMA conversion status')
   call check(size(conversion%model%components) == 2 .and. &
      conversion%reconstruction_error < 1.0e-8_dp .and. &
      conversion%partial_fraction%reconstruction_error < 1.0e-8_dp, &
      'aggregate ARIMA to UCARIMA reconstruction')
   call check_close(conversion%model%denominator, ucarima%denominator, &
      1.0e-10_dp, 'converted UCARIMA denominator')
   components(1) = tfarima_ucarima_component([1.0_dp, -0.3_dp], [1.0_dp], &
      [1.0_dp], 1.5_dp)
   components(2) = tfarima_ucarima_component([1.0_dp, 0.4_dp], [1.0_dp], &
      [1.0_dp], 0.7_dp)
   ucarima = tfarima_build_ucarima(components, tolerance=1.0e-12_dp)
   matrix = reshape([1.0_dp, -0.3_dp, 1.0_dp, 0.4_dp], [2, 2])
   conversion = tfarima_arima_to_ucarima([1.0_dp, 0.1_dp, -0.12_dp], &
      [1.0_dp], ucarima%ma_polynomial, ucarima%innovation_variance, matrix, &
      [2, 2], [.true., .true.], canonical=.true., tolerance=1.0e-9_dp)
   call check(conversion%info == 0 .and. conversion%admissible .and. &
      conversion%reconstruction_error < 1.0e-7_dp, &
      'two-factor ARIMA to UCARIMA conversion')
   components(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp, -1.0_dp], &
      [1.0_dp], 1.0_dp)
   components(2) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 0.5_dp)
   ucarima = tfarima_build_ucarima(components, tolerance=1.0e-12_dp)
   values = [(0.1_dp*real(i, dp) + sin(0.3_dp*real(i, dp)), i=1,30)]
   ucarima_ssm = tfarima_ucarima_state_space(values, ucarima, &
      tolerance=1.0e-10_dp)
   call check(ucarima_ssm%info == 0 .and. ucarima_ssm%diffuse, &
      'integrated UCARIMA diffuse state-space construction')
   smoothing = tfarima_ucarima_smooth(values, ucarima, forecast_horizon=2, &
      tolerance=1.0e-10_dp)
   call check(smoothing%info == 0, 'integrated UCARIMA smoothing status')
   call check_close(sum(smoothing%component, dim=2), values, 1.0e-9_dp, &
      'integrated smoothed component reconstruction')
   ucarima_fit = tfarima_ucarima_fit(values, ucarima, &
      estimate_variance=[.true., .false.], max_iterations=150, &
      tolerance=1.0e-5_dp)
   call check(ucarima_fit%info == 0 .and. ucarima_fit%converged .and. &
      ucarima_fit%parameters(1) > 0.0_dp, &
      'exact diffuse UCARIMA variance fit')
   deallocate(components)

   allocate(components(1))
   components(1) = tfarima_ucarima_component([1.0_dp, -0.2_dp], [1.0_dp], &
      [1.0_dp], 0.4_dp)
   ucarima = tfarima_build_ucarima(components, tolerance=1.0e-12_dp)
   deallocate(values)
   allocate(values(120))
   values = 0.0_dp
   do i = 2, size(values)
      noise = modulo(sin(12.9898_dp*real(i, dp))*43758.5453_dp, 1.0_dp) - &
         0.5_dp
      values(i) = 0.6_dp*values(i - 1) + noise
   end do
   smoothing = tfarima_ucarima_smooth(values, ucarima, tolerance=1.0e-8_dp)
   ucarima_fit = tfarima_ucarima_fit(values, ucarima, &
      estimate_ar=reshape([.true.], [1, 1]), &
      estimate_variance=[.true.], max_iterations=200, tolerance=1.0e-6_dp)
   call check(ucarima_fit%info == 0 .and. ucarima_fit%converged, &
      'exact UCARIMA parameter fit convergence')
   call check(abs(ucarima_fit%parameters(1) + 0.6_dp) < 0.15_dp .and. &
      ucarima_fit%parameters(2) > 0.0_dp, &
      'exact UCARIMA AR and variance estimates')
   call check(ucarima_fit%log_likelihood > smoothing%log_likelihood .and. &
      size(ucarima_fit%residuals) == size(values) .and. &
      all(shape(ucarima_fit%covariance) == [2, 2]) .and. &
      all(ucarima_fit%standard_error >= 0.0_dp), &
      'exact UCARIMA fit likelihood and inference')
   components(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp, 0.1_dp], 0.3_dp)
   ucarima = tfarima_build_ucarima(components, tolerance=1.0e-12_dp)
   allocate(input(size(values)))
   do i = 1, size(values)
      input(i) = modulo(sin(7.123_dp*real(i, dp))*24634.6345_dp, 1.0_dp) - &
         0.5_dp
   end do
   values(1) = input(1)
   do i = 2, size(values)
      values(i) = input(i) + 0.5_dp*input(i - 1)
   end do
   ucarima_fit = tfarima_ucarima_fit(values, ucarima, &
      estimate_ma=reshape([.true.], [1, 1]), &
      estimate_variance=[.true.], max_iterations=200, tolerance=1.0e-6_dp)
   call check(ucarima_fit%info == 0 .and. &
      abs(ucarima_fit%parameters(1) - 0.5_dp) < 0.15_dp, &
      'exact UCARIMA MA parameter estimate')
   components(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 0.2_dp)
   ucarima = tfarima_build_ucarima(components, tolerance=1.0e-12_dp)
   values = 2.0_dp + input
   matrix = reshape([(1.0_dp, i=1,size(values))], [size(values), 1])
   ucarima_fit = tfarima_ucarima_fit(values, ucarima, &
      regressors=matrix, initial_regression=[0.0_dp], &
      estimate_variance=[.true.], max_iterations=200, tolerance=1.0e-6_dp)
   call check(ucarima_fit%info == 0 .and. &
      abs(ucarima_fit%regression_coefficients(1) - 2.0_dp) < 0.1_dp .and. &
      ucarima_fit%parameters(2) > 0.0_dp, &
      'exact UCARIMA regression and variance estimates')
   deallocate(components)

   lag_model = tfarima_lag_polynomial([1, 4], [0.0_dp, 0.0_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [2, 2]), &
      [0.5_dp, 0.2_dp], exponent=2)
   call check(lag_model%info == 0, 'restricted lag-polynomial status')
   call check_close(lag_model%base, [1.0_dp, -0.5_dp, 0.0_dp, 0.0_dp, &
      -0.4_dp], 1.0e-14_dp, 'restricted sparse lag polynomial')
   updated_lag_model = tfarima_update_lag_polynomial(lag_model, [0.4_dp, 0.1_dp])
   call check_close(updated_lag_model%base, [1.0_dp, -0.4_dp, 0.0_dp, 0.0_dp, &
      -0.2_dp], 1.0e-14_dp, 'restricted lag-polynomial update')
   call check(size(lag_model%expanded) == 9, 'raised lag-polynomial degree')
   call check(tfarima_polynomial_admissible([1.0_dp, -0.5_dp]) .and. &
      .not. tfarima_polynomial_admissible([1.0_dp, -1.1_dp]), &
      'lag-polynomial admissibility')

   if (.not. allocated(input)) allocate(input(120))
   allocate(output(120))
   do i = 1, size(input)
      input(i) = sin(0.17_dp*real(i, dp)) + 0.6_dp*cos(0.71_dp*real(i, dp)) + &
         0.4_dp*sin(1.13_dp*real(i, dp)) + &
         0.08_dp*real(modulo(37*i, 17) - 8, dp)
   end do
   transfer = tfarima_transfer([1.2_dp, -0.4_dp], [1.0_dp, -0.3_dp], delay=2)
   signal = tfarima_filter(input, transfer)
   do i = 1, size(output)
      output(i) = 0.5_dp + signal(i) + 0.03_dp*sin(1.37_dp*real(i, dp))
   end do
   transfer_fit = tfarima_transfer_fit(output, input, 2, 1, 1, [1.0_dp], &
      [1.0_dp], initial_parameters=[1.0_dp, -0.2_dp, -0.2_dp, 0.4_dp], &
      include_mean=.true., max_iterations=200)
   call check(transfer_fit%info == 0 .and. transfer_fit%converged, &
      'conditional transfer-function fit convergence')
   call check_close(transfer_fit%transfer%numerator, [1.2_dp, -0.4_dp], &
      3.0e-2_dp, 'conditional transfer numerator estimates')
   call check_close(transfer_fit%transfer%denominator, [1.0_dp, -0.3_dp], &
      3.0e-2_dp, 'conditional transfer denominator estimates')
   call check(abs(transfer_fit%mean - 0.5_dp) < 1.0e-2_dp .and. &
      transfer_fit%innovation_variance > 0.0_dp .and. &
      all(shape(transfer_fit%covariance) == [4, 4]), &
      'conditional transfer fit statistics')
   numerator_model = tfarima_lag_polynomial([1], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [0.2_dp])
   denominator_model = tfarima_lag_polynomial([1], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [0.2_dp])
   transfer_fit = tfarima_restricted_transfer_fit(output, input, 2, &
      numerator_model, denominator_model, [1.0_dp], [1.0_dp], 1.0_dp, &
      include_mean=.true., max_iterations=200)
   call check(transfer_fit%info == 0 .and. transfer_fit%converged, &
      'restricted transfer-function fit convergence')
   call check_close(transfer_fit%transfer%numerator, [1.2_dp, -0.4_dp], &
      3.0e-2_dp, 'restricted transfer numerator estimates')
   call check_close(transfer_fit%transfer%denominator, [1.0_dp, -0.3_dp], &
      3.0e-2_dp, 'restricted transfer denominator estimates')

   allocate(transfer_specs(2))
   transfer_specs(1) = tfarima_transfer_spec([1.0_dp], [1.0_dp], delay=0)
   transfer_specs(2) = tfarima_transfer_spec([-0.3_dp], &
      [1.0_dp, -0.1_dp], delay=2)
   matrix = reshape([(sin(0.17_dp*real(i, dp)) + &
      0.2_dp*cos(0.53_dp*real(i, dp)), i=1,120), &
      (cos(0.29_dp*real(i, dp)) + 0.3_dp*sin(0.71_dp*real(i, dp)), &
      i=1,120)], [120, 2])
   output = 0.8_dp + 1.3_dp*matrix(:, 1)
   transfer = tfarima_transfer([-0.6_dp], [1.0_dp, -0.2_dp], delay=2)
   output = output + tfarima_filter(matrix(:, 2), transfer)
   do i = 1, size(output)
      output(i) = output(i) + 0.04_dp*(modulo(real(37*i, dp), 17.0_dp) - 8.0_dp)
   end do
   transfer_noise = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 0.2_dp)
   exact_transfer_fit = tfarima_exact_transfer_fit(output, matrix, &
      transfer_specs, transfer_noise, &
      regressors=reshape([(1.0_dp, i=1,120)], [120, 1]), &
      initial_regression=[0.5_dp], max_iterations=200, tolerance=1.0e-6_dp)
   call check(exact_transfer_fit%info == 0 .and. &
      exact_transfer_fit%converged, 'exact multi-input transfer fit status')
   call check(abs(exact_transfer_fit%regression_coefficients(1) - 0.8_dp) < &
      0.05_dp .and. &
      abs(exact_transfer_fit%inputs(1)%transfer%numerator(1) - 1.3_dp) < &
      0.05_dp .and. &
      abs(exact_transfer_fit%inputs(2)%transfer%numerator(1) + 0.6_dp) < &
      0.05_dp .and. &
      abs(exact_transfer_fit%inputs(2)%transfer%denominator(2) + 0.2_dp) < &
      0.08_dp, 'exact multi-input transfer estimates')
   call check(size(exact_transfer_fit%residuals) == size(output) .and. &
      all(shape(exact_transfer_fit%signal) == [120, 2]) .and. &
      exact_transfer_fit%noise%innovation_variance > 0.0_dp .and. &
      all(shape(exact_transfer_fit%covariance) == [5, 5]), &
      'exact transfer fit statistics')
   future_matrix = reshape([0.2_dp, -0.1_dp, 0.4_dp, &
      0.3_dp, 0.6_dp, -0.2_dp], [3, 2])
   exact_transfer_forecast = tfarima_exact_transfer_forecast( &
      exact_transfer_fit, output, matrix, future_matrix, &
      future_regressors=reshape([1.0_dp, 1.0_dp, 1.0_dp], [3, 1]), &
      levels=[0.95_dp])
   call check(exact_transfer_forecast%info == 0 .and. &
      all(ieee_is_finite(exact_transfer_forecast%mean)) .and. &
      all(exact_transfer_forecast%variance > 0.0_dp) .and. &
      all(exact_transfer_forecast%lower(:, 1) < &
      exact_transfer_forecast%mean) .and. &
      all(exact_transfer_forecast%upper(:, 1) > &
      exact_transfer_forecast%mean), 'exact known-input transfer forecast')
   allocate(input_models(2))
   input_models(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 1.0_dp)
   input_models(2) = input_models(1)
   transfer_ccf = tfarima_exact_transfer_ccf(exact_transfer_fit, matrix, &
      input_models, 8)
   call check(size(transfer_ccf) == 2 .and. &
      all(transfer_ccf%info == 0), 'multi-input residual CCF status')
   call check(all(ieee_is_finite(transfer_ccf(1)%correlation)) .and. &
      size(transfer_ccf(2)%lag) == 17, 'multi-input residual CCF values')
   deallocate(input_models)
   deallocate(transfer_specs)

   allocate(transfer_specs(3))
   transfer_specs(1) = tfarima_transfer_spec([1.0_dp], [1.0_dp])
   transfer_specs(2) = tfarima_transfer_spec([-0.3_dp], &
      [1.0_dp, -0.1_dp], delay=2)
   transfer_specs(3) = tfarima_transfer_spec([0.05_dp], [1.0_dp])
   selection_matrix = reshape([(sin(0.17_dp*real(i, dp)) + &
      0.2_dp*cos(0.53_dp*real(i, dp)), i=1,120), &
      (cos(0.29_dp*real(i, dp)) + 0.3_dp*sin(0.71_dp*real(i, dp)), &
      i=1,120), (sin(0.91_dp*real(i, dp)) + &
      0.2_dp*cos(1.27_dp*real(i, dp)), i=1,120)], [120, 3])
   selection_regressors = reshape([(1.0_dp, i=1,120), &
      (real(modulo(13*i, 19) - 9, dp), i=1,120)], [120, 2])
   output = 0.8_dp + 1.3_dp*selection_matrix(:, 1)
   transfer = tfarima_transfer([-0.6_dp], [1.0_dp, -0.2_dp], delay=2)
   output = output + tfarima_filter(selection_matrix(:, 2), transfer)
   output = output + [(0.04_dp*real(modulo(37*i, 17) - 8, dp), &
      i=1,size(output))]
   transfer_selection = tfarima_select_transfer(output, selection_matrix, &
      transfer_specs, transfer_noise, regressors=selection_regressors, &
      initial_regression=[0.5_dp, 0.0_dp], significance_level=0.05_dp, &
      keep_regressors=[.true., .false.], max_iterations=200)
   call check(transfer_selection%info == 0 .and. &
      transfer_selection%converged .and. transfer_selection%steps == 2, &
      'backward transfer selection convergence')
   call check(all(transfer_selection%retained_regressor .eqv. &
      [.true., .false.]) .and. all(transfer_selection%retained_input .eqv. &
      [.true., .true., .false.]), 'backward transfer retained masks')
   call check(size(transfer_selection%fit%inputs) == 2 .and. &
      size(transfer_selection%fit%regression_coefficients) == 1 .and. &
      any(transfer_selection%removed_kind == tfarima_selection_regressor) .and. &
      any(transfer_selection%removed_kind == tfarima_selection_input), &
      'backward transfer reduced exact model')
   allocate(input_models(2))
   input_models(1) = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
      [1.0_dp], 1.0_dp)
   input_models(2) = input_models(1)
   transfer_diagnostics = tfarima_diagnose_transfer( &
      transfer_selection%fit, 20, selection_matrix(:, :2), input_models)
   call check(transfer_diagnostics%info == 0 .and. &
      transfer_diagnostics%ljung_box%info == 0 .and. &
      transfer_diagnostics%weighted_ljung_box%info == 0, &
      'transfer portmanteau diagnostics')
   call check(size(transfer_diagnostics%residual_acf) == 21 .and. &
      size(transfer_diagnostics%residual_pacf) == 20 .and. &
      size(transfer_diagnostics%input_ccf) == 2 .and. &
      abs(transfer_diagnostics%cumulative_periodogram( &
      size(transfer_diagnostics%cumulative_periodogram)) - 1.0_dp) < &
      1.0e-12_dp, 'transfer residual diagnostic values')
   call check(transfer_diagnostics%normality_p_value >= 0.0_dp .and. &
      transfer_diagnostics%normality_p_value <= 1.0_dp .and. &
      transfer_diagnostics%cumulative_periodogram_p_value >= 0.0_dp .and. &
      transfer_diagnostics%cumulative_periodogram_p_value <= 1.0_dp, &
      'transfer residual diagnostic probabilities')
   transfer_signal = tfarima_simulate_transfer(selection_matrix(:, :2), &
      transfer_selection%fit%inputs, selection_regressors(:, :1), &
      transfer_selection%fit%regression_coefficients)
   call check(transfer_signal%info == 0 .and. &
      maxval(abs(transfer_signal%total - transfer_signal%regression_signal - &
      sum(transfer_signal%input_signal, dim=2))) < 1.0e-12_dp, &
      'deterministic transfer simulation components')
   allocate(simulation_innovations(3, 1))
   simulation_innovations(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
   allocate(simulation_fit%inputs(0), simulation_fit%regression_coefficients(0))
   simulation_fit%noise = tfarima_ucarima_component([1.0_dp], &
      [1.0_dp, -1.0_dp], [1.0_dp], 1.0_dp)
   transfer_simulation = tfarima_simulate_exact_model(simulation_fit, 3, 1, &
      noise_history=[10.0_dp], noise_innovations=simulation_innovations)
   call check(transfer_simulation%info == 0, &
      'conditional integrated transfer simulation status')
   call check_close(transfer_simulation%noise(:, 1), &
      [11.0_dp, 11.0_dp, 11.0_dp], 1.0e-14_dp, &
      'conditional integrated transfer simulation')
   deallocate(simulation_innovations)
   allocate(simulation_innovations(35, 2))
   simulation_innovations = 0.0_dp
   transfer_simulation = tfarima_simulate_exact_model( &
      transfer_selection%fit, 30, 2, &
      regressors=selection_regressors(:35, :1), input_models=input_models, &
      noise_innovations=simulation_innovations, burn_in=5, seed=12345)
   repeated_simulation = tfarima_simulate_exact_model( &
      transfer_selection%fit, 30, 2, &
      regressors=selection_regressors(:35, :1), input_models=input_models, &
      noise_innovations=simulation_innovations, burn_in=5, seed=12345)
   call check(transfer_simulation%info == 0 .and. &
      all(shape(transfer_simulation%output) == [30, 2]) .and. &
      all(shape(transfer_simulation%input) == [30, 2, 2]), &
      'multi-path transfer simulation shape')
   call check_close(reshape(transfer_simulation%output, [60]), &
      reshape(repeated_simulation%output, [60]), 0.0_dp, &
      'seeded transfer simulation reproducibility')
   call check(maxval(abs(transfer_simulation%output - &
      transfer_simulation%regression_signal - &
      sum(transfer_simulation%input_signal, dim=2))) < 1.0e-12_dp, &
      'multi-path transfer simulation components')
   deallocate(input_models)
   deallocate(transfer_specs)

   deallocate(input, output)
   allocate(input(300), output(300))
   input = 0.0_dp
   do i = 2, size(input)
      noise = modulo(sin(19.123_dp*real(i, dp))*31415.9265_dp, 1.0_dp) - &
         0.5_dp
      input(i) = 0.65_dp*input(i - 1) + noise
   end do
   transfer = tfarima_transfer([1.1_dp, -0.3_dp], &
      [1.0_dp, -0.25_dp], delay=3)
   output = tfarima_filter(input, transfer)
   output = output + [(0.002_dp*sin(1.71_dp*real(i, dp)), i=1,size(output))]
   transfer_noise = tfarima_ucarima_component([1.0_dp, -0.65_dp], &
      [1.0_dp], [1.0_dp], 1.0_dp)
   prewhitened = tfarima_prewhitened_ccf(input, output, transfer_noise, 20)
   call check(prewhitened%info == 0 .and. &
      lbound(prewhitened%correlation, 1) == 1 .and. &
      ubound(prewhitened%correlation, 1) == 41 .and. &
      all(prewhitened%lag == [(i, i=-20,20)]), &
      'prewhitened CCF bounds and status')
   identified_transfer = tfarima_identify_transfer(output, input, &
      transfer_noise, 1, 1, 20)
   call check(identified_transfer%info == 0 .and. &
      identified_transfer%selected_delay == 3, &
      'automatic transfer delay identification')
   call check_close(identified_transfer%specification%transfer%numerator, &
      [1.1_dp, -0.3_dp], 0.12_dp, 'identified transfer numerator')
   call check_close(identified_transfer%specification%transfer%denominator, &
      [1.0_dp, -0.25_dp], 0.12_dp, 'identified transfer denominator')

   values = tfarima_outlier_response(tfarima_outlier_ao, 4, &
      [1.0_dp, -0.5_dp], [1.0_dp], [1.0_dp])
   call check_close(values, [1.0_dp, -0.5_dp, 0.0_dp, 0.0_dp], 1.0e-14_dp, &
      'additive-outlier residual response')
   values = tfarima_outlier_response(tfarima_outlier_ls, 4, &
      [1.0_dp, -0.5_dp], [1.0_dp], [1.0_dp])
   call check_close(values, [1.0_dp, 0.5_dp, 0.5_dp, 0.5_dp], 1.0e-14_dp, &
      'level-shift residual response')
   values = tfarima_outlier_response(tfarima_outlier_tc, 4, &
      [1.0_dp, -0.5_dp], [1.0_dp], [1.0_dp])
   call check_close(values, [1.0_dp, 0.2_dp, 0.14_dp, 0.098_dp], 1.0e-14_dp, &
      'temporary-change residual response')

   output = [(0.1_dp*sin(0.43_dp*real(i, dp)), i=1,80)]
   output(20) = output(20) + 6.0_dp
   outliers = tfarima_detect_outliers(output, [1.0_dp], [1.0_dp], [1.0_dp], &
      [tfarima_outlier_ao], candidate_positions=[20])
   call check(outliers%info == 0 .and. size(outliers%position) == 1 .and. &
      outliers%position(1) == 20 .and. &
      outliers%outlier_type(1) == tfarima_outlier_ao .and. &
      abs(outliers%effect(1) - 6.0_dp) < 0.1_dp .and. &
      outliers%t_ratio(1) > 10.0_dp, 'additive-outlier detection')
   outliers = tfarima_detect_outliers(output, [1.0_dp], [1.0_dp], [1.0_dp], &
      [tfarima_outlier_ao], cutoff=3.0_dp)
   call check(outliers%info == 0 .and. any(outliers%position == 20), &
      'automatic outlier timing scan')
   output = [(0.08_dp*cos(0.31_dp*real(i, dp)), i=1,80)]
   output(30:) = output(30:) + 3.0_dp
   outliers = tfarima_detect_outliers(output, [1.0_dp], [1.0_dp], [1.0_dp], &
      [tfarima_outlier_ls], candidate_positions=[30])
   call check(outliers%info == 0 .and. size(outliers%position) == 1 .and. &
      outliers%position(1) == 30 .and. &
      outliers%outlier_type(1) == tfarima_outlier_ls .and. &
      abs(outliers%effect(1) - 3.0_dp) < 0.1_dp, 'level-shift detection')

   values = tfarima_month_lengths(2024, 1, 3)
   call check_close(values, [31.0_dp, 29.0_dp, 31.0_dp], 1.0e-14_dp, &
      'monthly calendar lengths')
   matrix = tfarima_weekday_counts(2024, 1, 2)
   call check_close(matrix(1, :), [4.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, &
      4.0_dp, 4.0_dp, 4.0_dp], 1.0e-14_dp, 'January 2024 weekday counts')
   call check_close(matrix(2, :), [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, &
      5.0_dp, 4.0_dp, 4.0_dp], 1.0e-14_dp, 'February 2024 weekday counts')
   values = tfarima_leap_year_regressor(2024, 1, 14, &
      working_day_coding=.true.)
   call check(abs(values(2) - 0.75_dp) < 1.0e-14_dp .and. &
      abs(values(14) + 0.25_dp) < 1.0e-14_dp, 'leap-year contrasts')
   values = tfarima_easter_regressor(2024, 1, 4, 4)
   call check_close(values, [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], 1.0e-14_dp, &
      'March Easter aggregation')
   values = tfarima_easter_regressor(2025, 3, 2, 22)
   call check_close(values, [2.0_dp/22.0_dp, 20.0_dp/22.0_dp], 1.0e-14_dp, &
      'cross-month Easter aggregation')
   calendar = tfarima_calendar_regressors(2024, 1, 1, tfarima_calendar_dif, &
      reference_weekday=0, include_month_length=.false., &
      include_leap_year=.false.)
   call check(calendar%info == 0 .and. &
      all(shape(calendar%values) == [1, 6]), 'calendar contrast shape')
   call check_close(calendar%values(1, :), &
      [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 1.0e-14_dp, &
      'reference-weekday contrasts')
   calendar = tfarima_calendar_regressors(2024, 1, 1, tfarima_calendar_wd, &
      include_month_length=.false., include_leap_year=.false.)
   call check(abs(calendar%values(1, 1) - 3.0_dp) < 1.0e-14_dp, &
      'working-day versus weekend contrast')

contains

   subroutine check(condition, label)
      !! Stop the test when an assertion fails.
      logical, intent(in) :: condition !! Flag controlling condition.
      character(*), intent(in) :: label !! Label.

      if (.not. condition) then
         write (*, '(a)') 'FAILED: '//label
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      !! Compare two real vectors within an absolute tolerance.
      real(dp), intent(in) :: actual(:) !! Observed values used for evaluation.
      real(dp), intent(in) :: expected(:) !! Expected.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      character(*), intent(in) :: label !! Label.

      call check(size(actual) == size(expected), label//' size')
      if (size(actual) == size(expected)) then
         call check(maxval(abs(actual - expected)) <= tolerance, label)
      end if
   end subroutine check_close

end program test_tfarima
