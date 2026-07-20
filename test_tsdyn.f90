! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Numerical tests for algorithms translated from R tsDyn.
program test_tsdyn
   use kind_mod, only: dp
   use tsdyn_mod
   use random_mod, only: set_random_seed
   use resampling_mod, only: resample, block_resample, additive_resample, &
      wild_resample
   use forecast_mod, only: accuracy_result_t, rolling_forecast_accuracy
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   type(tsdyn_lstar_simulation_t) :: simulation, repeated
   type(tsdyn_lstar_model_t) :: model
   type(tsdyn_lstar_forecast_t) :: forecast
   type(tsdyn_lstar_selection_t) :: selection
   type(tsdyn_regime_test_t) :: regime_test
   type(tsdyn_star_model_t) :: star_model, expanded_star
   type(tsdyn_star_forecast_t) :: star_forecast
   type(tsdyn_llar_model_t) :: llar_model
   type(tsdyn_llar_forecast_t) :: llar_forecast
   type(tsdyn_aar_model_t) :: aar_model
   type(tsdyn_aar_forecast_t) :: aar_forecast
   type(tsdyn_nnet_model_t) :: nnet_model
   type(tsdyn_nnet_forecast_t) :: nnet_forecast
   type(tsdyn_nnet_selection_t) :: nnet_selection
   type(tsdyn_tvecm_model_t) :: generating_tvecm, fitted_tvecm, partial_tvecm
   type(tsdyn_tvecm_model_t) :: multivariate_tvecm, fitted_multivariate_tvecm
   type(tsdyn_tvecm_model_t) :: partial_multivariate_tvecm
   type(tsdyn_tvecm_model_t) :: advanced_tvecm, bounded_tvecm
   type(tsdyn_tvecm_simulation_t) :: tvecm_simulation
   type(tsdyn_tvecm_simulation_t) :: multivariate_tvecm_simulation
   type(tsdyn_tvecm_simulation_t) :: post_sample_tvecm_simulation
   type(tsdyn_tvecm_forecast_t) :: tvecm_forecast
   type(tsdyn_tvecm_forecast_t) :: multivariate_tvecm_forecast
   type(nts_tar_model_t) :: generating_setar, fitted_setar
   type(nts_tar_model_t) :: common_setar, symmetric_setar
   type(nts_tar_model_t) :: general_setar, differenced_setar, adf_setar
   type(nts_tar_simulation_t) :: setar_simulation, general_setar_simulation
   type(nts_tar_forecast_t) :: setar_forecast
   type(tsdyn_setar_selection_t) :: setar_selection
   type(nts_mtar_model_t) :: generating_tvar, fitted_tvar
   type(nts_mtar_model_t) :: common_tvar
   type(nts_mtar_model_t) :: trend_tvar
   type(nts_mtar_simulation_t) :: tvar_simulation, trend_tvar_simulation
   type(nts_mtar_forecast_t) :: tvar_forecast
   type(tsdyn_tvar_selection_t) :: tvar_selection
   type(tsdyn_threshold_test_t) :: hansen_seo, seo
   type(tsdyn_regime_count_test_t) :: setar_linearity, tvar_lr
   type(tsdyn_bbc_test_t) :: bbc_lr, bbc_wald, bbc_lm
   type(tsdyn_kapshin_test_t) :: kapshin
   type(tsdyn_girf_t) :: setar_girf, tvar_girf, tvecm_girf
   type(tsdyn_irf_t) :: setar_irf, tvar_irf, tvecm_irf
   type(tsdyn_fevd_t) :: nonlinear_fevd
   type(tsdyn_delta_test_t) :: delta_test
   type(tsdyn_delta_linear_test_t) :: delta_linear_test
   type(tsdyn_regime_diagnostics_t) :: ar_diagnostics, setar_diagnostics
   type(tsdyn_regime_diagnostics_t) :: lstar_diagnostics, star_diagnostics
   type(tsdyn_rank_test_t) :: rank_test
   type(tsdyn_rank_selection_t) :: rank_selection
   type(tsdyn_setar_variance_t) :: setar_variance, setar_ml_variance
   type(tsdyn_setar_inference_t) :: setar_inference
   type(tsdyn_tvar_inference_t) :: tvar_inference
   type(tsdyn_tvecm_inference_t) :: tvecm_inference
   type(tsdyn_regime_path_t) :: regime_path
   type(tsdyn_setar_bootstrap_t) :: setar_bootstrap_result
   type(tsdyn_tvar_bootstrap_t) :: tvar_bootstrap_result
   type(tsdyn_tvecm_bootstrap_t) :: tvecm_bootstrap_result
   type(tsdyn_tvecm_bootstrap_t) :: multivariate_tvecm_bootstrap
   type(tsdyn_forecast_distribution_t) :: forecast_distribution
   type(rolling_forecast_result_t) :: rolling_result
   type(accuracy_result_t), allocatable :: rolling_accuracy(:)
   real(dp) :: tvecm_start(2, 2), tvecm_innovations(600, 2)
   real(dp) :: setar_innovations(1000), setar_draws(5, 8)
   real(dp) :: tvar_draws(1000, 2), tvar_forecast_draws(5, 2, 8)
   real(dp) :: girf_innovations(4, 2, 10), shocks(2, 2)
   real(dp) :: setar_girf_innovations(4, 5), setar_shocks(2)
   real(dp), allocatable :: wild_multipliers(:, :)
   real(dp), allocatable :: setar_bootstrap(:, :), tvar_bootstrap(:, :, :)
   integer, allocatable :: residual_indices(:, :)
   integer, allocatable :: permutation_indices(:, :)
   real(dp), allocatable :: gradient(:, :)
   real(dp), allocatable :: delta_series(:), delta_draws(:, :)
   real(dp), allocatable :: bootstrap_values(:), bootstrap_matrix(:, :)
   real(dp), allocatable :: resampled_setar(:), resampled_tvar(:, :)
   real(dp), allocatable :: resampled_tvecm(:, :)
   real(dp), allocatable :: probability_innovations(:, :)
   real(dp), allocatable :: rolling_actual(:)
   real(dp), allocatable :: multivariate_start(:, :)
   real(dp), allocatable :: multivariate_innovations(:, :)
   real(dp), allocatable :: multivariate_girf_innovations(:, :, :)
   real(dp), allocatable :: multivariate_shocks(:, :)
   real(dp), allocatable :: advanced_cointegration_offset(:)
   real(dp) :: rank_system(100, 3), common, stationary_u, stationary_v
   integer, allocatable :: bootstrap_order(:)
   real(dp) :: delta_value, linear_delta_value
   integer :: time, particle

   call set_random_seed(4107)
   common = 0.0_dp
   stationary_u = 0.0_dp
   stationary_v = 0.0_dp
   do time = 1, size(rank_system, 1)
      common = common + 0.15_dp + sin(0.17_dp*real(time*time, dp))
      stationary_u = 0.45_dp*stationary_u + &
         cos(0.31_dp*real(time*time, dp))
      stationary_v = 0.30_dp*stationary_v + &
         sin(0.23_dp*real(time*time, dp))
      rank_system(time, :) = [common, common + stationary_u, &
         0.5_dp*common + stationary_v]
   end do
   rank_test = tsdyn_rank_test(rank_system, 2, 'H_lc', 'eigen')
   call check(rank_test%info == 0 .and. rank_test%sample_size == 97 .and. &
      maxval(abs(rank_test%eigenvalues - [0.2792715040519849_dp, &
      0.1247053499623679_dp, 0.007042824150148524_dp])) < 1.0e-10_dp, &
      'Doornik rank-test Johansen eigenvalues')
   call check(maxval(abs(rank_test%trace_statistic - &
      [45.372257007127544_dp, 13.605457485956569_dp, &
      0.68557096423149089_dp])) < 1.0e-9_dp .and. &
      maxval(abs(rank_test%eigen_statistic - &
      [31.766799521170977_dp, 12.919886521725079_dp, &
      0.68557096423149089_dp])) < 1.0e-9_dp, &
      'Doornik trace and maximum-eigenvalue statistics')
   call check(maxval(abs(rank_test%trace_p_value - &
      [0.00028309218001831837_dp, 0.09393084166381216_dp, &
      0.40767495796777553_dp])) < 1.0e-10_dp .and. &
      maxval(abs(rank_test%eigen_p_value - &
      [0.0006560182365625833_dp, 0.07960676378030307_dp, &
      0.40768186095866676_dp])) < 1.0e-10_dp, &
      'Doornik asymptotic gamma p-values')
   call check(maxval(abs(rank_test%adjusted_trace_p_value - &
      [0.0003647268147420890_dp, 0.09933759757861038_dp, &
      0.41471875628581101_dp])) < 1.0e-10_dp .and. &
      rank_test%selected_rank == 1, 'Doornik adjusted p-values and rank selection')
   rank_test = tsdyn_rank_test(rank_system, 2, 'H_lc', 'trace', null_rank=1)
   call check(rank_test%info == 0 .and. rank_test%tested_rank == 1 .and. &
      abs(rank_test%selected_p_value - rank_test%trace_p_value(1)) < &
      1.0e-14_dp, 'Doornik specified-rank trace test')
   rank_test = tsdyn_rank_test(rank_system, 2, 'H_z')
   call check(rank_test%info == 0 .and. all(ieee_is_finite( &
      rank_test%eigen_p_value)), 'Doornik no-deterministic case')
   rank_test = tsdyn_rank_test(rank_system, 2, 'H_c')
   call check(rank_test%info == 0 .and. all(ieee_is_finite( &
      rank_test%trace_p_value)), 'Doornik restricted-constant case')
   rank_test = tsdyn_rank_test(rank_system, 2, 'H_l')
   call check(rank_test%info == 0 .and. all(ieee_is_finite( &
      rank_test%trace_p_value)), 'Doornik restricted-trend case')
   rank_test = tsdyn_rank_test(rank_system, 2, 'H_ql')
   call check(rank_test%info == 0 .and. all(ieee_is_finite( &
      rank_test%trace_p_value)), 'Doornik unrestricted-trend case')
   rank_selection = tsdyn_rank_select(rank_system, 4, 3, 'const', 'LL')
   call check(rank_selection%info == 0 .and. &
      all(shape(rank_selection%aic) == [4, 4]) .and. &
      all(ieee_is_finite(rank_selection%aic)) .and. &
      all(ieee_is_finite(rank_selection%bic)) .and. &
      all(ieee_is_finite(rank_selection%hq)), &
      'cointegration rank-selection criterion grids')
   do time = 1, rank_selection%lag_max
      call check(all(rank_selection%log_likelihood(1:, time) >= &
         rank_selection%log_likelihood(:2, time)), &
         'cointegration likelihood monotonicity across ranks')
   end do
   call check(rank_selection%selected_rank(2) == 2 .and. &
      all(rank_selection%selected_lag >= 1) .and. &
      all(rank_selection%selected_lag <= 4), &
      'cointegration BIC rank and lag selection')
   rank_selection = tsdyn_rank_select(rank_system, 3, deterministic='none', &
      same_sample=.false.)
   call check(rank_selection%info == 0 .and. &
      trim(rank_selection%fit_measure) == 'SSR' .and. &
      .not. rank_selection%same_sample .and. &
      all(shape(rank_selection%hq) == [4, 3]), &
      'cointegration SSR selection with lag-specific samples')

   simulation = tsdyn_lstar_simulate(0.10_dp, [0.40_dp], -0.15_dp, &
      [0.25_dp], 8.0_dp, 0.0_dp, 1, 0.15_dp, 1200, 500)
   call check(simulation%info == 0 .and. size(simulation%series) == 1200 .and. &
      all(ieee_is_finite(simulation%series)), 'LSTAR simulation')
   call set_random_seed(4107)
   repeated = tsdyn_lstar_simulate(0.10_dp, [0.40_dp], -0.15_dp, &
      [0.25_dp], 8.0_dp, 0.0_dp, 1, 0.15_dp, 1200, 500)
   call check(maxval(abs(repeated%series - simulation%series)) < 1.0e-14_dp, &
      'LSTAR shared RNG reproducibility')

   model = tsdyn_lstar_fit(simulation%series, 1, 1, 1)
   call check(model%info == 0 .and. model%gamma > 0.0_dp .and. &
      ieee_is_finite(model%rss) .and. model%rss > 0.0_dp, 'LSTAR fit status')
   call check(abs(model%threshold) < 0.20_dp, 'LSTAR threshold recovery')
   call check(all(ieee_is_finite(model%low_coefficients)) .and. &
      all(ieee_is_finite(model%transition_coefficients)), &
      'LSTAR finite coefficient estimates')
   call check(model%rss < sum((model%data(model%first_fitted:) - &
      sum(model%data)/real(size(model%data), dp))**2), &
      'LSTAR improvement over constant mean')
   call check(all(model%transition_weight > 0.0_dp .and. &
      model%transition_weight < 1.0_dp), 'LSTAR transition weights')
   regime_path = tsdyn_lstar_regimes(model)
   call check(regime_path%info == 0 .and. regime_path%regime_count == 2 .and. &
      regime_path%first_valid == model%first_fitted .and. &
      all(.not. regime_path%valid(:model%first_fitted - 1)) .and. &
      all(regime_path%valid(model%first_fitted:)) .and. &
      maxval(abs(regime_path%transition_weight(model%first_fitted:, 1) - &
      model%transition_weight)) < 1.0e-12_dp, &
      'LSTAR full-series regime probabilities and labels')
   call check(all(shape(model%covariance) == [6, 6]) .and. &
      size(model%standard_error) == 6 .and. &
      all(model%confidence_lower <= model%confidence_upper), &
      'LSTAR Hessian inference')

   ar_diagnostics = tsdyn_ar_regime_diagnostics(0.2_dp, [0.5_dp])
   call check(ar_diagnostics%info == 0 .and. &
      abs(ar_diagnostics%equilibrium_mean(1) - 0.4_dp) < 1.0e-12_dp .and. &
      abs(ar_diagnostics%root_modulus(1, 1) - 2.0_dp) < 1.0e-10_dp .and. &
      ar_diagnostics%stable(1) .and. ar_diagnostics%mean_defined(1), &
      'linear AR equilibrium and stability diagnostics')
   ar_diagnostics = tsdyn_ar_regime_diagnostics(0.25_dp, &
      [0.75_dp, -0.125_dp])
   call check(abs(minval(ar_diagnostics%root_modulus(:, 1)) - 2.0_dp) < &
      1.0e-10_dp .and. abs(maxval(ar_diagnostics%root_modulus(:, 1)) - &
      4.0_dp) < 1.0e-10_dp, 'linear AR second-order characteristic roots')
   ar_diagnostics = tsdyn_ar_regime_diagnostics(4.0_dp, [0.3_dp], .false.)
   call check(ar_diagnostics%mean_defined(1) .and. &
      ar_diagnostics%equilibrium_mean(1) == 0.0_dp, &
      'linear AR no-constant equilibrium')
   ar_diagnostics = tsdyn_ar_regime_diagnostics(0.2_dp, [0.0_dp, 0.0_dp])
   call check(ar_diagnostics%info == 0 .and. ar_diagnostics%stable(1) .and. &
      ar_diagnostics%ar_order(1) == 0 .and. &
      abs(ar_diagnostics%equilibrium_mean(1) - 0.2_dp) < 1.0e-12_dp, &
      'linear AR zero-order effective regime')
   ar_diagnostics = tsdyn_ar_regime_diagnostics(0.2_dp, [1.0_dp])
   call check(ar_diagnostics%unit_root(1) .and. &
      .not. ar_diagnostics%stable(1) .and. &
      .not. ar_diagnostics%mean_defined(1), 'linear AR unit-root diagnostics')
   ar_diagnostics = tsdyn_ar_regime_diagnostics(0.2_dp, [0.5_dp], &
      include_trend=.true.)
   call check(.not. ar_diagnostics%mean_defined(1), &
      'linear AR trend has no equilibrium mean')

   lstar_diagnostics = tsdyn_lstar_regime_diagnostics(model)
   call check(lstar_diagnostics%info == 0 .and. &
      size(lstar_diagnostics%equilibrium_mean) == 2 .and. &
      all(lstar_diagnostics%mean_defined) .and. &
      trim(lstar_diagnostics%regime_label(1)) == 'low' .and. &
      trim(lstar_diagnostics%regime_label(2)) == 'high', &
      'LSTAR regime equilibrium and stability diagnostics')

   gradient = tsdyn_lstar_gradient(model)
   call check(all(shape(gradient) == [1199, 6]) .and. &
      all(ieee_is_finite(gradient)), 'LSTAR parameter gradient')
   regime_test = tsdyn_lstar_regime_test(model)
   call check(regime_test%info == 0 .and. regime_test%statistic >= 0.0_dp .and. &
      regime_test%p_value >= 0.0_dp .and. regime_test%p_value <= 1.0_dp, &
      'LSTAR additional-regime test')

   forecast = tsdyn_lstar_forecast(model, 6)
   call check(forecast%info == 0 .and. size(forecast%mean) == 6 .and. &
      all(ieee_is_finite(forecast%mean)) .and. &
      all(forecast%transition_weight > 0.0_dp .and. &
      forecast%transition_weight < 1.0_dp), 'LSTAR recursive forecast')

   selection = tsdyn_lstar_select(simulation%series, 2, 'bic')
   call check(selection%info == 0 .and. selection%selected >= 1 .and. &
      size(selection%aic) == 8 .and. selection%model%info == 0, &
      'LSTAR order and delay selection')

   star_model = tsdyn_star_fit(simulation%series, 1, 2, 1)
   call check(star_model%info == 0 .and. &
      all(shape(star_model%coefficients) == [2, 2]) .and. &
      size(star_model%gamma) == 1 .and. star_model%gamma(1) > 0.0_dp, &
      'two-regime general STAR fit')
   regime_path = tsdyn_star_regimes(star_model)
   call check(regime_path%info == 0 .and. &
      regime_path%first_valid == star_model%first_fitted .and. &
      maxval(abs(regime_path%transition_weight(star_model%first_fitted:, :) - &
      star_model%transition_weight)) < 1.0e-12_dp .and. &
      all(regime_path%regime(star_model%first_fitted:) >= 1 .and. &
      regime_path%regime(star_model%first_fitted:) <= star_model%regimes), &
      'STAR transition weights and dominant regimes')
   gradient = tsdyn_star_gradient(star_model)
   call check(all(shape(gradient) == [1199, 6]) .and. &
      all(ieee_is_finite(gradient)), 'general STAR parameter gradient')
   regime_test = tsdyn_star_regime_test(star_model)
   call check(regime_test%info == 0 .and. regime_test%p_value >= 0.0_dp .and. &
      regime_test%p_value <= 1.0_dp, 'general STAR additional-regime test')
   star_forecast = tsdyn_star_forecast(star_model, 6)
   call check(star_forecast%info == 0 .and. &
      all(shape(star_forecast%transition_weight) == [6, 1]) .and. &
      all(ieee_is_finite(star_forecast%mean)), 'general STAR recursive forecast')
   expanded_star = tsdyn_star_add_regime(star_model)
   call check(expanded_star%info == 0 .and. expanded_star%regimes == 3 .and. &
      all(shape(expanded_star%coefficients) == [3, 2]), &
      'general STAR regime expansion')
   star_diagnostics = tsdyn_star_regime_diagnostics(expanded_star)
   call check(star_diagnostics%info == 0 .and. &
      size(star_diagnostics%equilibrium_mean) == 3 .and. &
      all(star_diagnostics%ar_order == expanded_star%order), &
      'general STAR cumulative-regime diagnostics')

   llar_model = tsdyn_llar_fit(simulation%series(:350), 2, 1, 1, &
      [0.04_dp, 0.07_dp, 0.10_dp, 0.15_dp, 0.22_dp, 0.32_dp, 0.48_dp])
   call check(llar_model%info == 0 .and. llar_model%selected >= 1 .and. &
      llar_model%selected_epsilon > 0.0_dp .and. &
      count(ieee_is_finite(llar_model%normalized_rmse)) >= 2, &
      'local-linear AR epsilon selection')
   call check(size(llar_model%fitted) == 348 .and. &
      count(ieee_is_finite(llar_model%fitted)) > 100 .and. &
      all(llar_model%neighbor_count >= 0), 'local-linear AR fitted values')
   llar_forecast = tsdyn_llar_forecast(llar_model, 5, .true., 20.0_dp)
   call check(llar_forecast%info == 0 .and. size(llar_forecast%mean) == 5 .and. &
      all(ieee_is_finite(llar_forecast%mean)) .and. &
      all(llar_forecast%neighbor_count > 2*(llar_model%order + 1)), &
      'local-linear AR recursive forecasts')

   aar_model = tsdyn_aar_fit(simulation%series(:500), 2, 1, 1, 7, &
      [1.0e-3_dp, 1.0e-2_dp, 0.1_dp, 1.0_dp, 10.0_dp, 100.0_dp])
   call check(aar_model%info == 0 .and. &
      all(shape(aar_model%coefficients) == [7, 2]) .and. &
      all(aar_model%lambda >= 0.0_dp) .and. &
      all(aar_model%component_df > 0.0_dp), 'additive AR penalized fit')
   call check(size(tsdyn_aar_component(aar_model, 1, &
      [-0.2_dp, 0.0_dp, 0.2_dp])) == 3 .and. &
      all(ieee_is_finite(aar_model%fitted)), 'additive AR component effects')
   aar_forecast = tsdyn_aar_forecast(aar_model, 5)
   call check(aar_forecast%info == 0 .and. size(aar_forecast%mean) == 5 .and. &
      all(ieee_is_finite(aar_forecast%mean)), 'additive AR recursive forecast')

   nnet_model = tsdyn_nnet_fit(simulation%series(:350), 2, 3, 1, 1, 600, &
      1.0e-6_dp, 1.0e-5_dp)
   call check(nnet_model%info == 0 .and. nnet_model%rss > 0.0_dp .and. &
      all(ieee_is_finite(nnet_model%fitted)) .and. &
      nnet_model%rss < sum((simulation%series(3:350) - &
      sum(simulation%series(3:350))/348.0_dp)**2), &
      'neural-network autoregression fit')
   call check(size(tsdyn_nnet_predict(nnet_model, reshape( &
      simulation%series(349:350), [1, 2]))) == 1, &
      'neural-network one-step prediction')
   nnet_forecast = tsdyn_nnet_forecast(nnet_model, 5)
   call check(nnet_forecast%info == 0 .and. size(nnet_forecast%mean) == 5 .and. &
      all(ieee_is_finite(nnet_forecast%mean)), &
      'neural-network recursive forecast')
   nnet_selection = tsdyn_nnet_select(simulation%series(:250), 2, [1, 2], &
      'bic', 1, 1, 400, 1.0e-6_dp, 1.0e-5_dp)
   call check(nnet_selection%info == 0 .and. nnet_selection%selected >= 1 .and. &
      size(nnet_selection%bic) == 2 .and. &
      all(ieee_is_finite(nnet_selection%bic)), &
      'neural-network hidden-size selection')

   generating_tvecm%beta = 1.0_dp
   generating_tvecm%thresholds = [-0.05_dp, 0.05_dp]
   generating_tvecm%lag = 1
   generating_tvecm%threshold_count = 2
   generating_tvecm%include_constant = .false.
   generating_tvecm%include_trend = .false.
   generating_tvecm%only_error_correction = .false.
   allocate(generating_tvecm%coefficients(2, 9))
   generating_tvecm%coefficients = 0.0_dp
   generating_tvecm%coefficients(:, 1:3) = reshape( &
      [-0.40_dp, 0.30_dp, 0.20_dp, 0.02_dp, 0.02_dp, 0.20_dp], [2, 3])
   generating_tvecm%coefficients(:, 4:6) = reshape( &
      [-0.12_dp, 0.10_dp, 0.20_dp, 0.02_dp, 0.02_dp, 0.20_dp], [2, 3])
   generating_tvecm%coefficients(:, 7:9) = reshape( &
      [-0.35_dp, 0.35_dp, 0.20_dp, 0.02_dp, 0.02_dp, 0.20_dp], [2, 3])
   tvecm_start = reshape([0.0_dp, 0.0_dp, 0.02_dp, -0.01_dp], [2, 2])
   do time = 1, size(tvecm_innovations, 1)
      tvecm_innovations(time, 1) = 0.045_dp*sin(0.71_dp*real(time, dp))
      tvecm_innovations(time, 2) = 0.040_dp*cos(0.53_dp*real(time, dp))
   end do
   tvecm_simulation = tsdyn_tvecm_simulate_from_innovations(generating_tvecm, &
      tvecm_start, tvecm_innovations)
   call check(tvecm_simulation%info == 0 .and. &
      all(shape(tvecm_simulation%series) == [600, 2]) .and. &
      all([(count(tvecm_simulation%regime == particle) > 20, particle=1, 3)]), &
      'three-regime TVECM supplied-innovation simulation')
   fitted_tvecm = tsdyn_tvecm_fit(tvecm_simulation%series, 1, 2, 0.04_dp, &
      [1.0_dp], 17, .false., .false., .false.)
   call check(fitted_tvecm%info == 0 .and. size(fitted_tvecm%thresholds) == 2 .and. &
      all(shape(fitted_tvecm%coefficients) == [2, 9]) .and. &
      all([(count(fitted_tvecm%regime == particle) > 20, particle=1, 3)]), &
      'three-regime TVECM estimation')
   tvecm_inference = tsdyn_tvecm_inference(fitted_tvecm)
   call check(tvecm_inference%info == 0 .and. &
      trim(tvecm_inference%variance_mode) == 'pooled' .and. &
      trim(tvecm_inference%adjustment) == 'OLS' .and. &
      all(shape(tvecm_inference%standard_error) == [2, 9]) .and. &
      all(tvecm_inference%standard_error > 0.0_dp) .and. &
      all(ieee_is_finite(tvecm_inference%t_statistic)) .and. &
      maxval(abs(tvecm_inference%t_statistic - &
      fitted_tvecm%coefficients/tvecm_inference%standard_error)) < 1.0e-12_dp .and. &
      all(tvecm_inference%p_value >= 0.0_dp) .and. &
      all(tvecm_inference%p_value <= 1.0_dp) .and. &
      maxval(abs(tvecm_inference%coefficient_covariance - &
      transpose(tvecm_inference%coefficient_covariance))) < 1.0e-12_dp, &
      'TVECM pooled coefficient covariance and tests')
   tvecm_inference = tsdyn_tvecm_inference(fitted_tvecm, 'regime', 'ML')
   call check(tvecm_inference%info == 0 .and. &
      trim(tvecm_inference%variance_mode) == 'regime' .and. &
      trim(tvecm_inference%adjustment) == 'ML' .and. &
      all(tvecm_inference%regime_degrees_of_freedom == &
      tvecm_inference%regime_observations) .and. &
      maxval(abs(tvecm_inference%pooled_residual_covariance - &
      fitted_tvecm%covariance)) < 1.0e-12_dp, &
      'TVECM regime-specific ML covariance sandwich')
   regime_path = tsdyn_tvecm_regimes(fitted_tvecm)
   call check(regime_path%info == 0 .and. &
      regime_path%first_valid == fitted_tvecm%first_fitted .and. &
      all(regime_path%regime(fitted_tvecm%first_fitted:) == &
      fitted_tvecm%regime), 'TVECM error-correction regime path')
   tvecm_forecast = tsdyn_tvecm_forecast(fitted_tvecm, 6)
   call check(tvecm_forecast%info == 0 .and. &
      all(shape(tvecm_forecast%mean) == [6, 2]) .and. &
      all(ieee_is_finite(tvecm_forecast%mean)), 'TVECM recursive level forecast')
   partial_tvecm = tsdyn_tvecm_fit(tvecm_simulation%series, 1, 1, 0.05_dp, &
      [1.0_dp], 15, .true., .false., .false.)
   call check(partial_tvecm%info == 0 .and. &
      all(shape(partial_tvecm%coefficients) == [2, 4]), &
      'ECT-only switching TVECM estimation')
   tvecm_inference = tsdyn_tvecm_inference(partial_tvecm, 'regime')
   call check(tvecm_inference%info == 0 .and. &
      all(shape(tvecm_inference%standard_error) == [2, 4]) .and. &
      all(tvecm_inference%standard_error > 0.0_dp), &
      'partial TVECM shared-coefficient inference')

   multivariate_tvecm%cointegration = [1.0_dp, -1.0_dp, -0.5_dp]
   multivariate_tvecm%thresholds = [0.0_dp]
   multivariate_tvecm%lag = 1
   multivariate_tvecm%threshold_count = 1
   multivariate_tvecm%include_constant = .false.
   multivariate_tvecm%include_trend = .false.
   multivariate_tvecm%only_error_correction = .false.
   allocate(multivariate_tvecm%coefficients(3, 8))
   multivariate_tvecm%coefficients = 0.0_dp
   multivariate_tvecm%coefficients(:, 1) = [-0.25_dp, 0.10_dp, 0.05_dp]
   multivariate_tvecm%coefficients(1, 2) = 0.15_dp
   multivariate_tvecm%coefficients(2, 3) = 0.10_dp
   multivariate_tvecm%coefficients(3, 4) = 0.10_dp
   multivariate_tvecm%coefficients(:, 5) = [-0.20_dp, 0.08_dp, 0.04_dp]
   multivariate_tvecm%coefficients(1, 6) = 0.12_dp
   multivariate_tvecm%coefficients(2, 7) = 0.08_dp
   multivariate_tvecm%coefficients(3, 8) = 0.08_dp
   allocate(multivariate_start(2, 3), multivariate_innovations(500, 3))
   multivariate_start(1, :) = [0.0_dp, 0.0_dp, 0.0_dp]
   multivariate_start(2, :) = [0.02_dp, -0.01_dp, 0.01_dp]
   do time = 1, size(multivariate_innovations, 1)
      multivariate_innovations(time, 1) = 0.04_dp*sin(0.37_dp*real(time, dp))
      multivariate_innovations(time, 2) = 0.03_dp*cos(0.29_dp*real(time, dp))
      multivariate_innovations(time, 3) = 0.025_dp*sin(0.43_dp*real(time, dp))
   end do
   multivariate_tvecm_simulation = tsdyn_tvecm_simulate_from_innovations( &
      multivariate_tvecm, multivariate_start, multivariate_innovations)
   call check(multivariate_tvecm_simulation%info == 0 .and. &
      all(shape(multivariate_tvecm_simulation%series) == [500, 3]) .and. &
      count(multivariate_tvecm_simulation%regime == 1) > 20 .and. &
      count(multivariate_tvecm_simulation%regime == 2) > 20, &
      'three-variable TVECM supplied-innovation simulation')
   fitted_multivariate_tvecm = tsdyn_tvecm_fit( &
      multivariate_tvecm_simulation%series, 1, 1, 0.05_dp, &
      threshold_grid_size=15, only_error_correction=.false., &
      include_constant=.false., include_trend=.false., &
      cointegrating_coefficients=[1.0_dp, 0.5_dp])
   call check(fitted_multivariate_tvecm%info == 0 .and. &
      all(shape(fitted_multivariate_tvecm%coefficients) == [3, 8]) .and. &
      maxval(abs(fitted_multivariate_tvecm%cointegration - &
      [1.0_dp, -1.0_dp, -0.5_dp])) < 1.0e-12_dp, &
      'fixed-vector multivariate TVECM estimation')
   partial_multivariate_tvecm = tsdyn_tvecm_fit( &
      multivariate_tvecm_simulation%series, 1, 1, 0.05_dp, &
      threshold_grid_size=13, only_error_correction=.true., &
      include_constant=.false., include_trend=.false., &
      cointegrating_coefficients=[1.0_dp, 0.5_dp])
   call check(partial_multivariate_tvecm%info == 0 .and. &
      all(shape(partial_multivariate_tvecm%coefficients) == [3, 5]), &
      'multivariate ECT-only switching TVECM estimation')
   allocate(advanced_cointegration_offset( &
      size(multivariate_tvecm_simulation%series, 1)))
   advanced_cointegration_offset = 0.01_dp
   advanced_tvecm = tsdyn_tvecm_fit(multivariate_tvecm_simulation%series, &
      1, 1, 0.05_dp, threshold_grid_size=15, &
      only_error_correction=.false., include_constant=.true., &
      include_trend=.true., cointegrating_coefficients=[1.0_dp, 0.5_dp], &
      fixed_thresholds=[fitted_multivariate_tvecm%thresholds(1) - 0.01_dp], &
      cointegration_offset=advanced_cointegration_offset, &
      shared_deterministic=.true.)
   call check(advanced_tvecm%info == 0 .and. &
      advanced_tvecm%shared_deterministic .and. &
      all(shape(advanced_tvecm%coefficients) == [3, 10]) .and. &
      abs(advanced_tvecm%thresholds(1) - &
      (fitted_multivariate_tvecm%thresholds(1) - 0.01_dp)) < 1.0e-12_dp .and. &
      maxval(abs(advanced_tvecm%cointegration_offset - 0.01_dp)) < 1.0e-12_dp, &
      'fixed-threshold offset TVECM with shared deterministic terms')
   tvecm_inference = tsdyn_tvecm_inference(advanced_tvecm, 'regime')
   call check(tvecm_inference%info == 0 .and. &
      all(shape(tvecm_inference%standard_error) == [3, 10]), &
      'shared-deterministic offset TVECM inference')
   regime_path = tsdyn_tvecm_regimes(advanced_tvecm)
   multivariate_tvecm_forecast = tsdyn_tvecm_forecast(advanced_tvecm, 4)
   call check(regime_path%info == 0 .and. &
      all(regime_path%regime(advanced_tvecm%first_fitted:) == &
      advanced_tvecm%regime) .and. multivariate_tvecm_forecast%info == 0 .and. &
      all(ieee_is_finite(multivariate_tvecm_forecast%mean)), &
      'offset TVECM regime alignment and forecast continuation')
   post_sample_tvecm_simulation = tsdyn_tvecm_simulate_from_innovations( &
      advanced_tvecm, advanced_tvecm%data(size(advanced_tvecm%data, 1) - 1:, :), &
      multivariate_innovations(:4, :))
   call check(post_sample_tvecm_simulation%info == 0 .and. &
      all(shape(post_sample_tvecm_simulation%series) == [4, 3]) .and. &
      all(ieee_is_finite(post_sample_tvecm_simulation%series)), &
      'offset TVECM post-sample supplied-innovation simulation')
   bounded_tvecm = tsdyn_tvecm_fit(multivariate_tvecm_simulation%series, &
      1, 1, 0.05_dp, threshold_grid_size=15, &
      only_error_correction=.false., include_constant=.false., &
      include_trend=.false., cointegrating_coefficients=[1.0_dp, 0.5_dp], &
      threshold_bounds=reshape([fitted_multivariate_tvecm%thresholds(1) - &
      1.0e-10_dp, fitted_multivariate_tvecm%thresholds(1) + 1.0e-10_dp], &
      [1, 2]))
   call check(bounded_tvecm%info == 0 .and. &
      abs(bounded_tvecm%thresholds(1) - &
      fitted_multivariate_tvecm%thresholds(1)) < 1.0e-12_dp, &
      'bounded TVECM threshold search')
   regime_path = tsdyn_tvecm_regimes(fitted_multivariate_tvecm)
   call check(regime_path%info == 0 .and. &
      all(regime_path%regime(fitted_multivariate_tvecm%first_fitted:) == &
      fitted_multivariate_tvecm%regime), &
      'multivariate TVECM cointegration-regime path')
   multivariate_tvecm_forecast = tsdyn_tvecm_forecast( &
      fitted_multivariate_tvecm, 4)
   call check(multivariate_tvecm_forecast%info == 0 .and. &
      all(shape(multivariate_tvecm_forecast%mean) == [4, 3]), &
      'multivariate TVECM recursive forecast')
   tvecm_inference = tsdyn_tvecm_inference(fitted_multivariate_tvecm, 'regime')
   call check(tvecm_inference%info == 0 .and. &
      all(shape(tvecm_inference%standard_error) == [3, 8]) .and. &
      all(tvecm_inference%standard_error > 0.0_dp), &
      'multivariate TVECM regime covariance inference')
   multivariate_tvecm_bootstrap = tsdyn_tvecm_bootstrap( &
      fitted_multivariate_tvecm, fitted_multivariate_tvecm%residuals, .true., &
      threshold_grid_size=11, cointegrating_coefficients=[1.0_dp, 0.5_dp])
   call check(multivariate_tvecm_bootstrap%info == 0 .and. &
      multivariate_tvecm_bootstrap%refitted .and. &
      all(shape(multivariate_tvecm_bootstrap%series) == [500, 3]), &
      'multivariate TVECM residual bootstrap refit')
   multivariate_tvecm_bootstrap = tsdyn_tvecm_bootstrap(advanced_tvecm, &
      advanced_tvecm%residuals, .true., threshold_grid_size=11, &
      cointegrating_coefficients=[1.0_dp, 0.5_dp], &
      fixed_thresholds=advanced_tvecm%thresholds)
   call check(multivariate_tvecm_bootstrap%info == 0 .and. &
      multivariate_tvecm_bootstrap%refitted .and. &
      multivariate_tvecm_bootstrap%fitted_model%shared_deterministic .and. &
      maxval(abs(multivariate_tvecm_bootstrap%fitted_model%thresholds - &
      advanced_tvecm%thresholds)) < 1.0e-12_dp, &
      'fixed-threshold offset TVECM bootstrap refit')
   allocate(multivariate_girf_innovations(3, 3, 4), &
      multivariate_shocks(3, 3))
   multivariate_girf_innovations = 0.0_dp
   multivariate_shocks = 0.0_dp
   do particle = 1, 3
      multivariate_shocks(particle, particle) = 0.05_dp
   end do
   tvecm_girf = tsdyn_tvecm_girf_from_innovations( &
      fitted_multivariate_tvecm, &
      fitted_multivariate_tvecm%data(size(fitted_multivariate_tvecm%data, 1) - &
      1:, :), multivariate_girf_innovations, multivariate_shocks)
   call check(tvecm_girf%info == 0 .and. &
      all(shape(tvecm_girf%response) == [3, 3, 3]), &
      'multivariate TVECM generalized impulse responses')
   tvecm_girf = tsdyn_tvecm_girf_from_innovations(advanced_tvecm, &
      advanced_tvecm%data(size(advanced_tvecm%data, 1) - 1:, :), &
      multivariate_girf_innovations, multivariate_shocks)
   call check(tvecm_girf%info == 0 .and. &
      all(ieee_is_finite(tvecm_girf%response)), &
      'offset TVECM post-sample generalized impulse responses')

   generating_setar%coefficients = reshape([0.05_dp, 0.00_dp, -0.05_dp, &
      0.70_dp, 0.20_dp, 0.70_dp], [3, 2])
   generating_setar%thresholds = [-0.12_dp, 0.12_dp]
   generating_setar%innovation_sd = [0.08_dp, 0.08_dp, 0.08_dp]
   generating_setar%ar_order = [1, 1, 1]
   generating_setar%include_mean = [.true., .true., .true.]
   generating_setar%delay = 1
   setar_diagnostics = tsdyn_setar_regime_diagnostics(generating_setar)
   call check(setar_diagnostics%info == 0 .and. &
      all(setar_diagnostics%mean_defined) .and. &
      all(setar_diagnostics%stable) .and. &
      maxval(abs(setar_diagnostics%equilibrium_mean - &
      [1.0_dp/6.0_dp, 0.0_dp, -1.0_dp/6.0_dp])) < 1.0e-12_dp, &
      'SETAR regime equilibrium and stability diagnostics')
   do time = 1, size(setar_innovations)
      setar_innovations(time) = 0.08_dp*sin(0.73_dp*real(time, dp))
   end do
   setar_simulation = tsdyn_setar_simulate_from_innovations(generating_setar, &
      setar_innovations, 200)
   call check(setar_simulation%info == 0 .and. &
      size(setar_simulation%series) == 800, 'SETAR supplied-innovation simulation')
   fitted_setar = tsdyn_setar_fit(setar_simulation%series, [1, 1, 1], 2, 1, &
      0.08_dp, 15)
   call check(fitted_setar%info == 0 .and. size(fitted_setar%thresholds) == 2 .and. &
      all(fitted_setar%regime_observations > 30), 'three-regime SETAR estimation')
   regime_path = tsdyn_setar_regimes(fitted_setar)
   call check(regime_path%info == 0 .and. &
      all([(count(regime_path%regime == particle) == &
      fitted_setar%regime_observations(particle), particle=1, 3)]) .and. &
      all(regime_path%regime(:regime_path%first_valid - 1) == 0), &
      'SETAR aligned hard-regime classification')
   regime_path = tsdyn_setar_regimes(fitted_setar, fitted_setar%data, &
      fitted_setar%data)
   call check(regime_path%info == 0 .and. &
      all([(count(regime_path%regime == particle) == &
      fitted_setar%regime_observations(particle), particle=1, 3)]), &
      'SETAR replacement and external threshold classification')
   setar_variance = tsdyn_setar_residual_variance(fitted_setar)
   setar_ml_variance = tsdyn_setar_residual_variance(fitted_setar, 'ML')
   call check(setar_variance%info == 0 .and. &
      maxval(abs(setar_variance%regime - &
      fitted_setar%innovation_sd**2)) < 1.0e-12_dp .and. &
      setar_variance%pooled > 0.0_dp .and. &
      all(setar_variance%degrees_of_freedom > 0), &
      'SETAR OLS pooled and regime residual variances')
   setar_ml_variance = tsdyn_setar_residual_variance(fitted_setar, &
      threshold_variable=fitted_setar%data)
   call check(maxval(abs(setar_ml_variance%regime - &
      setar_variance%regime)) < 1.0e-12_dp, &
      'SETAR external threshold variance path')
   setar_ml_variance = tsdyn_setar_residual_variance(fitted_setar, 'ML')
   call check(setar_ml_variance%info == 0 .and. &
      all(setar_ml_variance%regime < setar_variance%regime) .and. &
      setar_ml_variance%pooled < setar_variance%pooled, &
      'SETAR maximum-likelihood residual variances')
   setar_inference = tsdyn_setar_inference(fitted_setar)
   call check(setar_inference%info == 0 .and. &
      all(setar_inference%active) .and. &
      all(setar_inference%standard_error > 0.0_dp) .and. &
      all([(maxval(abs(setar_inference%covariance(:, :, particle) - &
      transpose(setar_inference%covariance(:, :, particle)))) < &
      1.0e-12_dp, particle=1, 3)]), &
      'SETAR conditional coefficient covariance inference')
   setar_inference = tsdyn_setar_inference(fitted_setar, 'regime', 'ML')
   call check(setar_inference%info == 0 .and. &
      trim(setar_inference%variance_mode) == 'regime' .and. &
      trim(setar_inference%residual_variance%adjustment) == 'ML', &
      'SETAR regime-specific ML covariance inference')
   common_setar = tsdyn_setar_restricted_fit(setar_simulation%series, [1, 1, 1], &
      fitted_setar%thresholds, 1, 'both')
   call check(common_setar%info == 0 .and. &
      maxval(abs(common_setar%coefficients(1, :) - &
      common_setar%coefficients(2, :))) < 1.0e-12_dp .and. &
      maxval(abs(common_setar%coefficients(2, :) - &
      common_setar%coefficients(3, :))) < 1.0e-12_dp, &
      'SETAR shared deterministic and lag restrictions')
   symmetric_setar = tsdyn_setar_restricted_fit(setar_simulation%series, &
      [1, 1, 1], fitted_setar%thresholds, 1, 'none', .true.)
   call check(symmetric_setar%info == 0 .and. &
      maxval(abs(symmetric_setar%coefficients(1, :) - &
      symmetric_setar%coefficients(3, :))) < 1.0e-12_dp, &
      'SETAR outer-regime symmetry restriction')
   setar_draws = 0.0_dp
   setar_forecast = tsdyn_setar_forecast_draws(fitted_setar, &
      size(fitted_setar%data), setar_draws)
   call check(setar_forecast%info == 0 .and. &
      all(shape(setar_forecast%mean) == [5]), 'SETAR supplied-draw forecast')
   setar_selection = tsdyn_setar_select(setar_simulation%series, [1, 1, 1], &
      2, 2, 0.08_dp, 9)
   call check(setar_selection%info == 0 .and. setar_selection%selected >= 1 .and. &
      size(setar_selection%score) == 2, 'SETAR threshold-delay selection')
   setar_selection = tsdyn_setar_select(setar_simulation%series, [1, 1, 1], &
      2, 2, 0.08_dp, 7, criterion='BIC', order_candidates=[1, 2], &
      threshold_counts=[1, 2], same_lags=.true., same_sample=.true.)
   call check(setar_selection%info == 0 .and. &
      size(setar_selection%score) == 8 .and. &
      trim(setar_selection%criterion) == 'bic' .and. &
      setar_selection%same_sample .and. &
      all(setar_selection%threshold_count == [1, 1, 1, 1, 2, 2, 2, 2]), &
      'SETAR threshold-count, lag-order, delay, and BIC selection')
   setar_selection = tsdyn_setar_select(setar_simulation%series, [1, 1], &
      1, 1, 0.08_dp, 7, criterion='pooled-AIC', order_candidates=[1], &
      common='lags', threshold_variables=reshape([setar_simulation%series, &
      -setar_simulation%series], [size(setar_simulation%series), 2]))
   call check(setar_selection%info == 0 .and. &
      size(setar_selection%score) == 2 .and. &
      all(setar_selection%threshold_variable_index == [1, 2]) .and. &
      trim(setar_selection%common) == 'lags', &
      'SETAR external transition-variable and shared-lag selection')
   setar_selection = tsdyn_setar_select(setar_simulation%series, [1, 1, 1], &
      2, 1, 0.08_dp, 7, criterion='SSR', outer_symmetric=.true.)
   call check(setar_selection%info == 0 .and. &
      maxval(abs(setar_selection%model%coefficients(1, :) - &
      setar_selection%model%coefficients(3, :))) < 1.0e-12_dp, &
      'SETAR outer-symmetric model selection')
   general_setar = tsdyn_setar_fit(setar_simulation%series, [2, 1, 2], &
      2, 1, 0.08_dp, 7, include_trend=[.true., .true., .true.], &
      lag_spacing=2, forecast_step=1, lag_indices=reshape( &
      [1, 1, 1, 2, 0, 2], [3, 2]))
   call check(general_setar%info == 0 .and. &
      trim(general_setar%representation) == 'level' .and. &
      all(general_setar%include_trend) .and. &
      general_setar%lag_active(1, 1) .and. &
      .not. general_setar%lag_active(1, 2) .and. &
      general_setar%lag_active(1, 3), &
      'SETAR trends and nonconsecutive spaced lags')
   setar_forecast = tsdyn_setar_forecast_draws(general_setar, &
      size(general_setar%data), setar_draws)
   call check(setar_forecast%info == 0 .and. &
      all(ieee_is_finite(setar_forecast%mean)), &
      'generalized SETAR recursive forecast')
   general_setar_simulation = tsdyn_setar_simulate_from_innovations( &
      general_setar, setar_innovations, 200)
   call check(general_setar_simulation%info == 0 .and. &
      all(ieee_is_finite(general_setar_simulation%series)) .and. &
      trim(general_setar_simulation%representation) == 'level' .and. &
      general_setar_simulation%lag_active(1, 3), &
      'generalized SETAR supplied-innovation simulation')
   differenced_setar = tsdyn_setar_fit(setar_simulation%series, [1, 1], &
      1, 1, 0.08_dp, 7, representation='diff')
   adf_setar = tsdyn_setar_fit(setar_simulation%series, [1, 1], &
      1, 1, 0.08_dp, 7, representation='ADF')
   call check(differenced_setar%info == 0 .and. adf_setar%info == 0 .and. &
      trim(differenced_setar%representation) == 'diff' .and. &
      trim(adf_setar%representation) == 'adf' .and. &
      any(abs(adf_setar%adf_coefficient) > 0.0_dp), &
      'SETAR differenced and ADF representations')
   setar_selection = tsdyn_setar_select(setar_simulation%series, [1, 1], &
      1, 1, 0.08_dp, 6, criterion='BIC', &
      include_trend=[.true., .true.], representation='diff')
   call check(setar_selection%info == 0 .and. &
      trim(setar_selection%model%representation) == 'diff' .and. &
      all(setar_selection%model%include_trend), &
      'generalized SETAR model selection')

   allocate(generating_tvar%intercept(2, 3), generating_tvar%ar(2, 2, 1, 3))
   allocate(generating_tvar%covariance(2, 2, 3))
   generating_tvar%intercept = reshape([0.05_dp, 0.01_dp, 0.0_dp, 0.0_dp, &
      -0.05_dp, -0.01_dp], [2, 3])
   generating_tvar%ar = 0.0_dp
   generating_tvar%ar(:, :, 1, 1) = reshape([0.65_dp, 0.05_dp, &
      0.02_dp, 0.45_dp], [2, 2])
   generating_tvar%ar(:, :, 1, 2) = reshape([0.25_dp, 0.03_dp, &
      0.03_dp, 0.35_dp], [2, 2])
   generating_tvar%ar(:, :, 1, 3) = generating_tvar%ar(:, :, 1, 1)
   generating_tvar%covariance = 0.0_dp
   do particle = 1, 3
      generating_tvar%covariance(:, :, particle) = reshape( &
         [0.0064_dp, 0.001_dp, 0.001_dp, 0.0049_dp], [2, 2])
   end do
   generating_tvar%thresholds = [-0.10_dp, 0.10_dp]
   generating_tvar%ar_order = [1, 1, 1]
   generating_tvar%include_mean = [.true., .true., .true.]
   generating_tvar%threshold_component = 1
   generating_tvar%delay = 1
   do time = 1, size(tvar_draws, 1)
      tvar_draws(time, 1) = sin(0.61_dp*real(time, dp))
      tvar_draws(time, 2) = cos(0.47_dp*real(time, dp))
   end do
   tvar_simulation = tsdyn_tvar_simulate_from_standard(generating_tvar, &
      tvar_draws, 200)
   call check(tvar_simulation%info == 0 .and. &
      all(shape(tvar_simulation%series) == [800, 2]), &
      'TVAR supplied-normal simulation')
   fitted_tvar = tsdyn_tvar_fit(tvar_simulation%series, [1, 1, 1], 2, 1, 1, &
      0.06_dp, 13)
   call check(fitted_tvar%info == 0 .and. size(fitted_tvar%thresholds) == 2 .and. &
      all(fitted_tvar%regime_observations > 25), 'three-regime TVAR estimation')
   regime_path = tsdyn_tvar_regimes(fitted_tvar)
   call check(regime_path%info == 0 .and. &
      all([(count(regime_path%regime == particle) == &
      fitted_tvar%regime_observations(particle), particle=1, 3)]), &
      'TVAR aligned hard-regime classification')
   tvar_inference = tsdyn_tvar_inference(fitted_tvar, 'regime', 'ML')
   call check(tvar_inference%info == 0 .and. &
      maxval(abs(tvar_inference%regime_residual_covariance - &
      fitted_tvar%covariance)) < 1.0e-12_dp .and. &
      all(tvar_inference%active) .and. &
      all(tvar_inference%standard_error > 0.0_dp), &
      'TVAR regime ML covariance and standard errors')
   tvar_inference = tsdyn_tvar_inference(fitted_tvar)
   call check(tvar_inference%info == 0 .and. &
      trim(tvar_inference%variance_mode) == 'pooled' .and. &
      trim(tvar_inference%adjustment) == 'OLS' .and. &
      all([(maxval(abs(tvar_inference%coefficient_covariance(:, :, &
      particle) - transpose(tvar_inference%coefficient_covariance(:, :, &
      particle)))) < 1.0e-12_dp, particle=1, 3)]), &
      'TVAR pooled Kronecker coefficient covariance')
   tvar_inference = tsdyn_tvar_inference(fitted_tvar, &
      threshold_variable=fitted_tvar%data(:, fitted_tvar%threshold_component))
   call check(tvar_inference%info == 0 .and. &
      all(tvar_inference%parameter_count == 3), &
      'TVAR external threshold covariance path')
   common_tvar = tsdyn_tvar_restricted_fit(tvar_simulation%series, [1, 1, 1], &
      fitted_tvar%thresholds, 1, 1, .true.)
   call check(common_tvar%info == 0 .and. &
      maxval(abs(common_tvar%intercept(:, 1) - common_tvar%intercept(:, 2))) < &
      1.0e-12_dp .and. maxval(abs(common_tvar%intercept(:, 2) - &
      common_tvar%intercept(:, 3))) < 1.0e-12_dp, &
      'TVAR common-intercept restriction')
   tvar_forecast_draws = 0.0_dp
   tvar_forecast = tsdyn_tvar_forecast_draws(fitted_tvar, &
      size(fitted_tvar%data, 1), tvar_forecast_draws)
   call check(tvar_forecast%info == 0 .and. &
      all(shape(tvar_forecast%mean) == [5, 2]), 'TVAR supplied-draw forecast')
   tvar_selection = tsdyn_tvar_select(tvar_simulation%series, [1, 1, 1], 2, &
      2, [1, 2], 0.06_dp, 7)
   call check(tvar_selection%info == 0 .and. tvar_selection%selected >= 1 .and. &
      size(tvar_selection%score) == 4, 'TVAR variable and delay selection')
   tvar_selection = tsdyn_tvar_select(tvar_simulation%series, [1, 1, 1], 2, &
      2, [1, 2], 0.06_dp, 6, momentum=.true., criterion='BIC', &
      order_candidates=[1, 2], threshold_counts=[1, 2], same_lags=.true., &
      same_sample=.true.)
   call check(tvar_selection%info == 0 .and. &
      size(tvar_selection%score) == 16 .and. &
      trim(tvar_selection%criterion) == 'bic' .and. &
      tvar_selection%momentum .and. tvar_selection%same_sample .and. &
      count(tvar_selection%threshold_count == 1) == 8 .and. &
      count(tvar_selection%threshold_count == 2) == 8, &
      'TVAR threshold-count, lag-order, component, delay, and BIC selection')
   tvar_selection = tsdyn_tvar_select(tvar_simulation%series, [1, 1], 1, &
      1, [integer ::], 0.06_dp, 7, criterion='AIC', &
      common_intercept=.true., threshold_variables=reshape( &
      [tvar_simulation%series(:, 1), -tvar_simulation%series(:, 1)], &
      [size(tvar_simulation%series, 1), 2]))
   call check(tvar_selection%info == 0 .and. &
      size(tvar_selection%score) == 2 .and. &
      all(tvar_selection%threshold_component == 0) .and. &
      all(tvar_selection%threshold_variable_index == [1, 2]) .and. &
      tvar_selection%common_intercept .and. &
      maxval(abs(tvar_selection%model%intercept(:, 1) - &
      tvar_selection%model%intercept(:, 2))) < 1.0e-12_dp, &
      'TVAR external transition and common-intercept selection')
   trend_tvar = tsdyn_tvar_fit(tvar_simulation%series, [1, 1, 1], 2, &
      1, 1, 0.06_dp, 7, include_trend=[.true., .true., .true.])
   call check(trend_tvar%info == 0 .and. all(trend_tvar%include_trend) .and. &
      any(abs(trend_tvar%trend) > 0.0_dp), 'TVAR regime-specific trends')
   tvar_forecast = tsdyn_tvar_forecast_draws(trend_tvar, &
      size(trend_tvar%data, 1), tvar_forecast_draws)
   call check(tvar_forecast%info == 0 .and. &
      all(ieee_is_finite(tvar_forecast%mean)), 'TVAR trend forecast continuation')
   trend_tvar_simulation = tsdyn_tvar_simulate_from_standard(trend_tvar, &
      tvar_draws, 200)
   call check(trend_tvar_simulation%info == 0 .and. &
      all(ieee_is_finite(trend_tvar_simulation%series)) .and. &
      maxval(abs(trend_tvar_simulation%trend - trend_tvar%trend)) < &
      1.0e-12_dp, 'TVAR trend supplied-normal simulation')
   tvar_selection = tsdyn_tvar_select(tvar_simulation%series, [1, 1, 1], &
      2, 1, [1], 0.06_dp, 6, criterion='BIC', &
      include_trend=[.true., .true., .true.])
   call check(tvar_selection%info == 0 .and. &
      all(tvar_selection%model%include_trend), 'TVAR trend model selection')

   allocate(setar_bootstrap(size(setar_simulation%series), 2))
   setar_bootstrap(:, 1) = setar_simulation%series
   setar_bootstrap(:, 2) = -setar_simulation%series
   setar_linearity = tsdyn_setar_linearity_test(setar_simulation%series, 1, &
      [1, 2], 0.08_dp, 9, bootstrap_series=setar_bootstrap)
   call check(setar_linearity%info == 0 .and. &
      all(setar_linearity%statistic >= 0.0_dp) .and. &
      all(shape(setar_linearity%bootstrap_statistic) == [2, 3]), &
      'Hansen SETAR regime-count tests')
   bbc_lr = tsdyn_bbc_test(setar_simulation%series, 1, 'LR', 0.08_dp, &
      bootstrap_series=setar_bootstrap)
   bbc_wald = tsdyn_bbc_test(setar_simulation%series, 1, 'Wald', 0.08_dp)
   bbc_lm = tsdyn_bbc_test(setar_simulation%series, 1, 'LM', 0.08_dp)
   call check(bbc_lr%info == 0 .and. bbc_wald%info == 0 .and. &
      bbc_lm%info == 0 .and. bbc_lr%statistic >= 0.0_dp .and. &
      bbc_wald%statistic >= 0.0_dp .and. bbc_lm%statistic >= 0.0_dp .and. &
      size(bbc_lr%bootstrap_statistic) == 2, &
      'Bec-Ben Salem-Carrasco LR, Wald, and LM tests')
   kapshin = tsdyn_kapshin_test(setar_simulation%series, 1, 'const', &
      points=9, minimum_middle=20, bootstrap_series=setar_bootstrap)
   call check(kapshin%info == 0 .and. all(kapshin%statistic >= 0.0_dp) .and. &
      all(shape(kapshin%bootstrap_statistic) == [2, 3]) .and. &
      size(kapshin%thresholds, 2) == 2, 'Kapetanios-Shin unit-root tests')
   allocate(tvar_bootstrap(size(tvar_simulation%series, 1), 2, 2))
   tvar_bootstrap(:, :, 1) = tvar_simulation%series
   tvar_bootstrap(:, :, 2) = -tvar_simulation%series
   tvar_lr = tsdyn_tvar_lr_test(tvar_simulation%series, 1, 1, [1, 2], &
      0.06_dp, 7, bootstrap_series=tvar_bootstrap)
   call check(tvar_lr%info == 0 .and. all(tvar_lr%statistic >= 0.0_dp) .and. &
      all(shape(tvar_lr%bootstrap_statistic) == [2, 3]), &
      'TVAR likelihood-ratio regime-count tests')

   allocate(wild_multipliers(size(tvecm_simulation%series, 1) - 2, 5))
   do particle = 1, size(wild_multipliers, 2)
      do time = 1, size(wild_multipliers, 1)
         wild_multipliers(time, particle) = merge(1.0_dp, -1.0_dp, &
            mod(time + particle, 2) == 0)
      end do
   end do
   hansen_seo = tsdyn_hansen_seo_test(tvecm_simulation%series, 1, 1.0_dp, &
      0.05_dp, 11, wild_multipliers)
   call check(hansen_seo%info == 0 .and. hansen_seo%statistic >= 0.0_dp .and. &
      size(hansen_seo%bootstrap_statistic) == 5 .and. &
      hansen_seo%p_value >= 0.0_dp .and. hansen_seo%p_value <= 1.0_dp, &
      'Hansen-Seo threshold-cointegration test')
   allocate(residual_indices(size(wild_multipliers, 1), 3))
   do particle = 1, size(residual_indices, 2)
      do time = 1, size(residual_indices, 1)
         residual_indices(time, particle) = 1 + &
            mod(time + 17*particle, size(residual_indices, 1))
      end do
   end do
   seo = tsdyn_seo_test(tvecm_simulation%series, 1, 1.0_dp, 0.08_dp, 9, &
      residual_indices)
   call check(seo%info == 0 .and. seo%statistic >= 0.0_dp .and. &
      size(seo%thresholds, 2) == 2 .and. &
      size(seo%bootstrap_statistic) == 3, 'Seo threshold-cointegration test')

   girf_innovations = 0.0_dp
   shocks = 0.0_dp
   shocks(1, 1) = 0.10_dp
   shocks(2, 2) = 0.10_dp
   setar_irf = tsdyn_setar_regime_irf(fitted_setar, 1, 4)
   call check(setar_irf%info == 0 .and. &
      all(shape(setar_irf%response) == [5, 1, 1]) .and. &
      abs(setar_irf%response(1, 1, 1) - 1.0_dp) < 1.0e-12_dp .and. &
      abs(setar_irf%response(2, 1, 1) - &
      fitted_setar%coefficients(1, 2)) < 1.0e-12_dp, &
      'SETAR fixed-regime impulse response')
   setar_irf = tsdyn_setar_regime_irf(fitted_setar, 1, 4, .true.)
   call check(setar_irf%info == 0 .and. setar_irf%cumulative .and. &
      abs(setar_irf%response(2, 1, 1) - 1.0_dp - &
      fitted_setar%coefficients(1, 2)) < 1.0e-12_dp, &
      'SETAR cumulative impulse response')
   setar_girf_innovations = 0.0_dp
   setar_shocks = [0.10_dp, -0.10_dp]
   setar_girf = tsdyn_setar_girf_from_innovations(fitted_setar, &
      fitted_setar%data(size(fitted_setar%data) - 1:), &
      setar_girf_innovations, setar_shocks)
   call check(setar_girf%info == 0 .and. &
      all(shape(setar_girf%response) == [4, 1, 2]) .and. &
      abs(setar_girf%response(1, 1, 1) - 0.10_dp) < 1.0e-12_dp .and. &
      abs(setar_girf%response(1, 1, 2) + 0.10_dp) < 1.0e-12_dp, &
      'SETAR paired nonlinear GIRF')
   setar_irf = tsdyn_setar_irf_bootstrap(fitted_setar, &
      [fitted_setar, fitted_setar], 1, 3, 0.90_dp)
   call check(setar_irf%info == 0 .and. setar_irf%bootstrap_count == 2 .and. &
      maxval(abs(setar_irf%lower - setar_irf%response)) < 1.0e-12_dp .and. &
      maxval(abs(setar_irf%upper - setar_irf%response)) < 1.0e-12_dp, &
      'SETAR bootstrap impulse-response intervals')
   tvar_irf = tsdyn_tvar_regime_irf(fitted_tvar, 1, 3)
   call check(tvar_irf%info == 0 .and. &
      all(shape(tvar_irf%response) == [4, 2, 2]) .and. &
      abs(tvar_irf%response(1, 1, 1) - 1.0_dp) < 1.0e-12_dp .and. &
      abs(tvar_irf%response(1, 2, 2) - 1.0_dp) < 1.0e-12_dp .and. &
      abs(tvar_irf%response(1, 1, 2)) < 1.0e-12_dp, &
      'TVAR fixed-regime unit impulse responses')
   tvar_irf = tsdyn_tvar_irf_bootstrap(fitted_tvar, &
      [fitted_tvar, fitted_tvar], 1, 3, 0.90_dp, orthogonalized=.true.)
   call check(tvar_irf%info == 0 .and. tvar_irf%orthogonalized .and. &
      tvar_irf%bootstrap_count == 2 .and. &
      maxval(abs(tvar_irf%lower - tvar_irf%upper)) < 1.0e-12_dp, &
      'TVAR orthogonalized bootstrap impulse-response intervals')
   tvecm_irf = tsdyn_tvecm_regime_irf(fitted_tvecm, 1, 3)
   call check(tvecm_irf%info == 0 .and. &
      all(shape(tvecm_irf%response) == [4, 2, 2]) .and. &
      abs(tvecm_irf%response(1, 1, 1) - 1.0_dp) < 1.0e-12_dp .and. &
      abs(tvecm_irf%response(1, 2, 2) - 1.0_dp) < 1.0e-12_dp, &
      'TVECM fixed-regime level impulse responses')
   tvecm_irf = tsdyn_tvecm_irf_bootstrap(fitted_tvecm, &
      [fitted_tvecm, fitted_tvecm], 1, 3, 0.90_dp)
   call check(tvecm_irf%info == 0 .and. tvecm_irf%bootstrap_count == 2 .and. &
      maxval(abs(tvecm_irf%lower - tvecm_irf%upper)) < 1.0e-12_dp, &
      'TVECM bootstrap impulse-response intervals')
   tvar_girf = tsdyn_tvar_girf_from_innovations(fitted_tvar, &
      fitted_tvar%data(size(fitted_tvar%data, 1) - 1:, :), girf_innovations, shocks)
   call check(tvar_girf%info == 0 .and. &
      all(shape(tvar_girf%response) == [4, 2, 2]) .and. &
      abs(tvar_girf%response(1, 1, 1) - 0.10_dp) < 1.0e-12_dp, &
      'TVAR paired nonlinear GIRF')
   tvecm_girf = tsdyn_tvecm_girf_from_innovations(fitted_tvecm, &
      fitted_tvecm%data(size(fitted_tvecm%data, 1) - 1:, :), &
      girf_innovations, shocks)
   nonlinear_fevd = tsdyn_nonlinear_fevd(tvecm_girf)
   call check(tvecm_girf%info == 0 .and. nonlinear_fevd%info == 0 .and. &
      maxval(abs(sum(nonlinear_fevd%decomposition, dim=3) - 1.0_dp)) < &
      1.0e-12_dp, 'TVECM nonlinear GIRF and FEVD')

   allocate(probability_innovations(4, 6))
   probability_innovations = 0.0_dp
   setar_forecast = tsdyn_setar_forecast_draws(fitted_setar, &
      size(fitted_setar%data), probability_innovations)
   forecast_distribution = tsdyn_setar_forecast_distribution(fitted_setar, &
      probability_innovations)
   call check(forecast_distribution%info == 0 .and. &
      maxval(abs(forecast_distribution%mean - setar_forecast%mean)) < &
      1.0e-12_dp, 'SETAR zero-innovation forecast distribution')
   forecast = tsdyn_lstar_forecast(model, 4)
   forecast_distribution = tsdyn_lstar_forecast_distribution(model, &
      probability_innovations)
   call check(forecast_distribution%info == 0 .and. &
      maxval(abs(forecast_distribution%mean - forecast%mean)) < 1.0e-12_dp, &
      'LSTAR zero-innovation forecast distribution')
   star_forecast = tsdyn_star_forecast(star_model, 4)
   forecast_distribution = tsdyn_star_forecast_distribution(star_model, &
      probability_innovations)
   call check(forecast_distribution%info == 0 .and. &
      maxval(abs(forecast_distribution%mean - star_forecast%mean)) < &
      1.0e-12_dp, 'STAR zero-innovation forecast distribution')
   llar_forecast = tsdyn_llar_forecast(llar_model, 4)
   forecast_distribution = tsdyn_llar_forecast_distribution(llar_model, &
      probability_innovations)
   call check(forecast_distribution%info == 0 .and. &
      maxval(abs(forecast_distribution%mean - llar_forecast%mean)) < &
      1.0e-10_dp, 'local-linear AR zero-innovation forecast distribution')
   aar_forecast = tsdyn_aar_forecast(aar_model, 4)
   forecast_distribution = tsdyn_aar_forecast_distribution(aar_model, &
      probability_innovations)
   call check(forecast_distribution%info == 0 .and. &
      maxval(abs(forecast_distribution%mean - aar_forecast%mean)) < &
      1.0e-12_dp, 'additive AR zero-innovation forecast distribution')
   nnet_forecast = tsdyn_nnet_forecast(nnet_model, 4)
   forecast_distribution = tsdyn_nnet_forecast_distribution(nnet_model, &
      probability_innovations)
   call check(forecast_distribution%info == 0 .and. &
      maxval(abs(forecast_distribution%mean - nnet_forecast%mean)) < &
      1.0e-12_dp, 'neural AR zero-innovation forecast distribution')
   do particle = 1, size(probability_innovations, 2)
      do time = 1, size(probability_innovations, 1)
         probability_innovations(time, particle) = 0.03_dp* &
            sin(0.7_dp*real(time + 3*particle, dp))
      end do
   end do
   probability_innovations(1, :) = 1.0_dp
   forecast_distribution = tsdyn_lstar_forecast_distribution(model, &
      probability_innovations, 0.8_dp, .true.)
   call check(forecast_distribution%info == 0 .and. &
      forecast_distribution%first_innovation_zero .and. &
      abs(forecast_distribution%mean(1) - forecast%mean(1)) < 1.0e-12_dp .and. &
      forecast_distribution%standard_error(1) < 1.0e-12_dp .and. &
      any(forecast_distribution%standard_error(2:) > 0.0_dp) .and. &
      all(forecast_distribution%lower <= forecast_distribution%upper) .and. &
      all(ieee_is_finite(forecast_distribution%mean)), &
      'nonlinear empirical intervals and zero first innovation')

   rolling_actual = setar_simulation%series(795:800)
   rolling_result = tsdyn_setar_rolling_forecast(fitted_setar, rolling_actual, &
      [1, 2], 2)
   call check(rolling_result%info == 0 .and. all(rolling_result%valid) .and. &
      any(rolling_result%refitted) .and. &
      all(shape(rolling_result%forecast) == [6, 2]) .and. &
      maxval(abs(rolling_result%error - (rolling_result%actual - &
      rolling_result%forecast))) < 1.0e-14_dp, &
      'SETAR rolling forecasts and periodic refitting')
   rolling_accuracy = rolling_forecast_accuracy(rolling_result, &
      fitted_setar%data)
   call check(size(rolling_accuracy) == 2 .and. &
      ieee_is_finite(rolling_accuracy(1)%rmse) .and. &
      ieee_is_finite(rolling_accuracy(2)%rmse), &
      'rolling forecast accuracy by horizon')
   rolling_actual = simulation%series(size(simulation%series) - 5:)
   rolling_result = tsdyn_lstar_rolling_forecast(model, rolling_actual, [1, 2])
   call check(rolling_result%info == 0 .and. all(rolling_result%valid) .and. &
      all(rolling_result%origin(:, 1) == rolling_result%target - 1) .and. &
      all(rolling_result%origin(:, 2) == rolling_result%target - 2), &
      'LSTAR rolling forecast alignment')
   rolling_result = tsdyn_star_rolling_forecast(star_model, rolling_actual, &
      [1, 2])
   call check(rolling_result%info == 0 .and. all(rolling_result%valid), &
      'general STAR rolling forecasts')
   rolling_actual = simulation%series(351:356)
   rolling_result = tsdyn_llar_rolling_forecast(llar_model, rolling_actual, &
      [1, 2])
   call check(rolling_result%info == 0 .and. all(rolling_result%valid), &
      'local-linear AR rolling forecasts')
   rolling_actual = simulation%series(501:506)
   rolling_result = tsdyn_aar_rolling_forecast(aar_model, rolling_actual, &
      [1, 2])
   call check(rolling_result%info == 0 .and. all(rolling_result%valid), &
      'additive AR rolling forecasts')
   rolling_actual = simulation%series(351:356)
   rolling_result = tsdyn_nnet_rolling_forecast(nnet_model, rolling_actual, &
      [1, 2])
   call check(rolling_result%info == 0 .and. all(rolling_result%valid) .and. &
      all(ieee_is_finite(rolling_result%forecast)), &
      'neural AR rolling forecasts')

   bootstrap_values = resample([1.0_dp, 2.0_dp, 3.0_dp], [3, 1, 1])
   call check(all(bootstrap_values == [3.0_dp, 1.0_dp, 1.0_dp]), &
      'generic indexed vector resampling')
   bootstrap_values = block_resample([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [3, 1], 2)
   call check(all(bootstrap_values == [3.0_dp, 4.0_dp, 1.0_dp, 2.0_dp]), &
      'generic circular block resampling')
   bootstrap_matrix = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2])
   bootstrap_matrix = additive_resample(bootstrap_matrix, [1.0_dp, -1.0_dp])
   call check(all(bootstrap_matrix == reshape([2.0_dp, 1.0_dp, 4.0_dp, &
      3.0_dp], [2, 2])), 'tsDyn additive wild resampling')
   bootstrap_values = wild_resample([1.0_dp, 2.0_dp, 3.0_dp], &
      [1.0_dp, -1.0_dp, 1.0_dp])
   call check(all(bootstrap_values == [1.0_dp, -2.0_dp, 3.0_dp]), &
      'generic multiplicative wild resampling')

   bootstrap_order = [(time, time=1, size(fitted_setar%residuals))]
   resampled_setar = resample(fitted_setar%residuals, bootstrap_order)
   setar_bootstrap_result = tsdyn_setar_bootstrap(fitted_setar, &
      resampled_setar, .true., 0.08_dp, 15)
   call check(setar_bootstrap_result%info == 0 .and. &
      setar_bootstrap_result%refitted .and. &
      maxval(abs(setar_bootstrap_result%series - fitted_setar%data)) < &
      1.0e-10_dp, 'SETAR residual bootstrap and refit')

   bootstrap_order = [(time, time=1, size(fitted_tvar%residuals, 1))]
   resampled_tvar = resample(fitted_tvar%residuals, bootstrap_order)
   tvar_bootstrap_result = tsdyn_tvar_bootstrap(fitted_tvar, resampled_tvar, &
      .true., 0.06_dp, 13)
   call check(tvar_bootstrap_result%info == 0 .and. &
      tvar_bootstrap_result%refitted .and. &
      maxval(abs(tvar_bootstrap_result%series - fitted_tvar%data)) < &
      1.0e-10_dp, 'TVAR residual bootstrap and refit')

   bootstrap_order = [(time, time=1, size(fitted_tvecm%residuals, 1))]
   resampled_tvecm = resample(fitted_tvecm%residuals, bootstrap_order)
   tvecm_bootstrap_result = tsdyn_tvecm_bootstrap(fitted_tvecm, &
      resampled_tvecm, .true., [fitted_tvecm%beta], 11)
   call check(tvecm_bootstrap_result%info == 0 .and. &
      tvecm_bootstrap_result%refitted .and. &
      maxval(abs(tvecm_bootstrap_result%series - fitted_tvecm%data)) < &
      1.0e-8_dp, 'TVECM residual bootstrap and refit')

   allocate(delta_series(120), permutation_indices(120, 5))
   allocate(delta_draws(140, 5))
   do time = 1, size(delta_series)
      delta_series(time) = sin(0.17_dp*real(time, dp)) + &
         0.35_dp*cos(0.43_dp*real(time, dp))
   end do
   do particle = 1, size(permutation_indices, 2)
      do time = 1, size(permutation_indices, 1)
         permutation_indices(time, particle) = time
      end do
   end do
   permutation_indices(:, 1) = [(1 + mod(time + 6, size(delta_series)), &
      time=1, size(delta_series))]
   permutation_indices(:, 2) = permutation_indices(size(delta_series):1:-1, 1)
   permutation_indices(:, 3) = [(1 + mod(37*(time - 1), &
      size(delta_series)), time=1, size(delta_series))]
   permutation_indices(:, 4) = [(1 + mod(49*(time - 1), &
      size(delta_series)), time=1, size(delta_series))]
   permutation_indices(:, 5) = [(1 + mod(53*(time - 1), &
      size(delta_series)), time=1, size(delta_series))]
   do particle = 1, size(delta_draws, 2)
      do time = 1, size(delta_draws, 1)
         delta_draws(time, particle) = sin(0.71_dp*real(time + &
            13*particle, dp)) + cos(1.17_dp*real(time + 5*particle, dp))
      end do
   end do
   delta_value = tsdyn_delta(delta_series, 2, 1, 0.75_dp)
   linear_delta_value = tsdyn_delta_linear(delta_series, 2, 1)
   call check(ieee_is_finite(delta_value) .and. &
      ieee_is_finite(linear_delta_value), 'Manzan delta scalar statistics')
   delta_test = tsdyn_delta_test(delta_series, [2, 3], 1, &
      [0.65_dp, 0.90_dp], permutation_indices)
   call check(delta_test%info == 0 .and. &
      all(shape(delta_test%statistic) == [2, 2]) .and. &
      all(shape(delta_test%permutation_statistic) == [5, 2, 2]) .and. &
      all(delta_test%p_value >= 1.0_dp/6.0_dp .and. &
      delta_test%p_value <= 1.0_dp), 'Manzan delta permutation test')
   delta_linear_test = tsdyn_delta_linear_test(delta_series, [2, 3], 1, &
      [0.65_dp, 0.90_dp], delta_draws, 8)
   call check(delta_linear_test%info == 0 .and. &
      all(shape(delta_linear_test%statistic) == [2, 2]) .and. &
      all(shape(delta_linear_test%bootstrap_statistic) == [5, 2, 2]) .and. &
      delta_linear_test%ar_order >= 0 .and. &
      delta_linear_test%ar_order <= 8 .and. &
      delta_linear_test%innovation_variance > 0.0_dp .and. &
      all(delta_linear_test%p_value >= 1.0_dp/6.0_dp .and. &
      delta_linear_test%p_value <= 1.0_dp), &
      'Manzan delta AR-bootstrap linearity test')

   print '(a)', 'tsDyn tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_tsdyn
