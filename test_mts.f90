! SPDX-License-Identifier: Artistic-2.0
! SPDX-FileComment: Regression tests for the MTS translation.
program test_mts
   use kind_mod, only: dp
   use mts_mod
   use time_series_random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(mts_var_fit_t) :: fitted, selected
   type(mts_var_forecast_t) :: forecast
   type(mts_var_order_t) :: order_result
   type(mts_var_irf_t) :: irf
   type(mts_fevd_t) :: decomposition, generalized_decomposition
   type(mts_varx_fit_t) :: varx_fit
   type(mts_varx_irf_t) :: varx_response
   type(mts_varx_order_t) :: varx_order
   type(mts_varma_model_t) :: varma_model
   type(mts_varma_model_t) :: vma_model
   type(mts_varma_covariance_t) :: varma_covariance
   type(mts_varma_simulation_t) :: varma_simulation
   type(mts_varma_simulation_t) :: seasonal_simulation
   type(mts_vma_fit_t) :: vma_fit, selected_vma
   type(mts_varma_fit_t) :: varma_fit, selected_varma
   type(mts_varma_refinement_t) :: varma_refinement
   type(mts_svarma_model_t) :: seasonal_model
   type(mts_svarma_fit_t) :: seasonal_fit
   type(mts_svarma_refinement_t) :: seasonal_refinement
   type(mts_vecm_fit_t) :: vecm_fit, johansen_vecm
   type(mts_vecm_forecast_t) :: vecm_prediction
   type(mts_factor_model_t) :: factor_model
   type(mts_factor_forecast_t) :: factor_prediction
   type(mts_constrained_factor_t) :: constrained_factor
   type(mts_factor_forecast_t) :: constrained_prediction
   type(mts_bvar_prior_t) :: bvar_prior, strong_prior
   type(mts_bvar_fit_t) :: bvar_fit, shrunk_bvar
   type(mts_var_forecast_t) :: bvar_prediction
   type(mts_common_volatility_t) :: common_volatility, prewhitened_volatility
   type(mts_mch_diagnostic_t) :: mch_diagnostic
   type(mts_bekk_fit_t) :: bekk_fit
   type(mts_dcc_fit_t) :: dcc_fit
   type(mts_adcc_fit_t) :: adcc_fit, adcc_t_fit
   type(mts_tse_tsui_fit_t) :: tse_tsui_fit, tse_tsui_t_fit
   type(mts_ewma_fit_t) :: ewma_fit, estimated_ewma
   type(mts_mchol_fit_t) :: mchol_fit, prewhitened_mchol
   type(mts_sccor_t) :: grouped_correlation
   type(mts_arch_test_t) :: univariate_arch
   type(mts_march_test_t) :: multivariate_arch
   type(mts_copula_fit_t) :: copula_fit, estimated_copula, grouped_copula
   type(mts_missing_result_t) :: missing_result, partial_missing_result
   type(mts_granger_test_t) :: granger_result
   type(mts_mq_t) :: mq_result
   type(mts_diagnostic_t) :: diagnostic_result
   type(mts_var_backtest_t) :: backtest_result
   type(mts_scm_identification_t) :: scm_identification
   type(mts_scm_structure_t) :: scm_structure
   type(mts_scm_spec_t) :: scm_specification
   type(mts_scm_fit_t) :: scm_fit
   type(mts_scm_refinement_t) :: scm_refinement
   type(mts_transfer_fit_t) :: transfer_fit
   type(mts_transfer2_fit_t) :: transfer2_fit
   type(mts_transfer_forecast_t) :: transfer_prediction
   type(mts_transfer_backtest_t) :: transfer_backtest
   type(mts_regts_fit_t) :: regts_fit
   type(mts_regts_refinement_t) :: regts_refinement
   type(mts_ecm_known_fit_t) :: ecm_known, ecm_known_refined
   type(mts_ecm_normalized_fit_t) :: ecm_normalized, ecm_normalized_refined
   type(mts_reverse_mq_t) :: reverse_mq
   type(mts_eccm_t) :: extended_ccm
   type(mts_corner_t) :: corner_table
   type(mts_vma_exact_fit_t) :: exact_vma, exact_vma_free
   type(mts_vma_exact_refinement_t) :: exact_vma_refinement
   type(mts_apca_t) :: asymptotic_pca
   type(mts_diffusion_forecast_t) :: diffusion_forecast
   type(mts_multivariate_regression_t) :: multivariate_regression
   type(mts_var_chi_t) :: var_chi
   type(mts_kronecker_identification_t) :: kronecker_identification
   type(mts_kronecker_spec_t) :: kronecker_specification
   type(mts_kronecker_fit_t) :: kronecker_fit
   type(mts_kronecker_refinement_t) :: kronecker_refinement
   type(mts_vma_order_t) :: vma_order
   real(dp), allocatable :: psi(:, :, :), pi_weights(:, :, :)
   real(dp) :: series(120, 2), expected(2), exogenous(120, 1), future_exogenous(4, 1)
   real(dp) :: innovations(80, 2)
   real(dp) :: seasonal_innovations(100, 1)
   real(dp) :: component(1, 1, 1), seasonal_component(1, 1, 1)
   real(dp) :: zero_component(1, 1, 0), scalar_sigma(1, 1)
   real(dp) :: vecm_series(140, 2), cointegration(2, 1), difference_now(2)
   real(dp) :: walk, spread_value
   real(dp) :: factor_series(160, 6), factor_one, factor_two
   real(dp) :: identity_two(2, 2)
   real(dp) :: factor_constraint(6, 3), score_covariance(2, 2)
   real(dp) :: volatility_returns(240, 3), whitened_covariance(3, 3)
   real(dp) :: conditional_variance, common_shock
   real(dp), allocatable :: correlation_angles(:), recovered_correlation(:, :)
   real(dp) :: correlation_fixture(3, 3)
   real(dp) :: missing_series(50, 2), missing_pi(2, 2, 2), missing_sigma(2, 2)
   real(dp) :: missing_constant(2)
   real(dp) :: granger_series(180, 3)
   real(dp) :: mch_residuals(160, 2), mch_covariance(2, 2, 160), mch_standardized(160, 2)
   real(dp) :: first_variance, second_variance
   real(dp) :: transfer_response(100), transfer_input(100), transfer_filtered(100)
   real(dp) :: transfer_noise(100), transfer_innovation
   real(dp) :: transfer_input2(100), transfer_trend(100), transfer_equilibrium(100)
   real(dp) :: regts_response(100, 2), regts_input(100, 1), regts_error(100, 2)
   real(dp) :: regts_future(2, 1), regts_shock(2)
   real(dp) :: ecm_process(140, 1), filter_data(4, 2), filter_weights(2, 2, 1)
   real(dp), allocatable :: filtered_data(:, :)
   real(dp) :: apca_data(12, 30)
   integer :: j
   real(dp) :: bekk_returns(120, 2), bekk_covariance(2, 2, 120), bekk_residual(2)
   real(dp) :: bekk_constant(2, 2), bekk_arch(2, 2), bekk_garch(2, 2)
   real(dp) :: bekk_initial(13), bekk_shock(2), cholesky11, cholesky21, cholesky22
   real(dp), allocatable :: bekk_prediction(:, :, :)
   real(dp), allocatable :: dcc_prediction(:, :, :)
   real(dp), allocatable :: adcc_prediction(:, :, :)
   real(dp), allocatable :: tse_tsui_prediction(:, :, :)
   real(dp), allocatable :: ewma_prediction(:, :, :), mchol_prediction(:, :, :)
   real(dp) :: dcc_variance(160, 2), future_variance(3, 2), dcc_expected(2, 2)
   logical :: bekk_mask(13)
   logical :: vecm_mask(4, 2)
   integer :: t

   series(1, :) = [0.2_dp, -0.1_dp]
   do t = 2, 120
      series(t, 1) = 0.55_dp*series(t - 1, 1) + 0.15_dp*series(t - 1, 2) + &
         sin(1.3_dp*real(t, dp)) + 0.15_dp*sin(0.017_dp*real(t*t, dp))
      series(t, 2) = -0.2_dp*series(t - 1, 1) + 0.45_dp*series(t - 1, 2) + &
         0.7_dp*cos(0.71_dp*real(t, dp)) + 0.12_dp*cos(0.013_dp*real(t*t, dp))
   end do
   fitted = mts_var(series, 1)
   call check(fitted%info == 0 .and. all(shape(fitted%ar) == [2, 2, 1]) .and. &
      all(shape(fitted%residuals) == [119, 2]), 'VAR estimation dimensions')
   do t = 2, 120
      expected = fitted%intercept + matmul(fitted%ar(:, :, 1), series(t - 1, :))
      call check(maxval(abs(series(t, :) - expected - fitted%residuals(t - 1, :))) < 1.0e-11_dp, &
         'VAR fitted residual identity')
   end do
   call check(all([(fitted%sigma(t, t) > 0.0_dp, t=1, 2)]) .and. &
      fitted%aic < huge(1.0_dp), 'VAR covariance and criteria')

   psi = mts_var_psi(fitted%ar, 3)
   call check(maxval(abs(psi(:, :, 3) - matmul(fitted%ar(:, :, 1), fitted%ar(:, :, 1)))) < 1.0e-11_dp, &
      'VAR MA-representation recursion')
   forecast = mts_var_forecast(fitted, series, 4)
   expected = fitted%intercept + matmul(fitted%ar(:, :, 1), series(120, :))
   call check(forecast%info == 0 .and. maxval(abs(forecast%mean(1, :) - expected)) < 1.0e-11_dp, &
      'VAR forecast recursion')
   call check(all(forecast%standard_error > 0.0_dp) .and. &
      all(forecast%standard_error(4, :) >= forecast%standard_error(1, :)), &
      'VAR forecast uncertainty')

   irf = mts_var_irf(fitted, 3)
   call check(irf%info == 0 .and. maxval(abs(irf%psi(:, :, 1) - &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))) < 1.0e-12_dp, &
      'VAR impulse response initialization')
   call check(maxval(abs(irf%orthogonal(:, :, 1) - irf%shock_factor)) < 1.0e-12_dp .and. &
      maxval(abs(irf%cumulative_orthogonal(:, :, 2) - &
      irf%orthogonal(:, :, 1) - irf%orthogonal(:, :, 2))) < 1.0e-12_dp, &
      'orthogonalized cumulative impulse responses')
   decomposition = mts_fevd(fitted, 4)
   call check(decomposition%info == 0 .and. &
      maxval(abs(sum(decomposition%contribution, 2) - 1.0_dp)) < 1.0e-12_dp .and. &
      maxval(abs(decomposition%variance(:, 1) - &
      [fitted%sigma(1, 1), fitted%sigma(2, 2)])) < 1.0e-12_dp, &
      'orthogonal forecast-error variance decomposition')
   generalized_decomposition = mts_fevd(fitted, 4, .true.)
   call check(generalized_decomposition%info == 0 .and. &
      maxval(abs(sum(generalized_decomposition%contribution, 2) - 1.0_dp)) < 1.0e-12_dp, &
      'generalized forecast-error variance decomposition')

   selected = mts_vars(series, [1, 3])
   call check(selected%info == 0 .and. selected%max_lag == 3 .and. &
      maxval(abs(selected%ar(:, :, 2))) == 0.0_dp, 'selected-lag VAR')
   order_result = mts_var_order(series, 3)
   call check(order_result%info == 0 .and. order_result%aic_order >= 0 .and. &
      order_result%aic_order <= 3 .and. all(order_result%p_value >= 0.0_dp) .and. &
      all(order_result%p_value <= 1.0_dp), 'VAR order selection')

   do t = 1, 120
      exogenous(t, 1) = sin(0.37_dp*real(t, dp)) + 0.2_dp*cos(0.019_dp*real(t*t, dp))
   end do
   series(1, :) = [0.1_dp, -0.2_dp]
   do t = 2, 120
      series(t, 1) = 0.45_dp*series(t - 1, 1) + 0.1_dp*series(t - 1, 2) + &
         0.8_dp*exogenous(t, 1) - 0.2_dp*exogenous(t - 1, 1) + &
         0.08_dp*sin(0.023_dp*real(t*t, dp))
      series(t, 2) = -0.15_dp*series(t - 1, 1) + 0.35_dp*series(t - 1, 2) - &
         0.4_dp*exogenous(t, 1) + 0.3_dp*exogenous(t - 1, 1) + &
         0.07_dp*cos(0.031_dp*real(t*t, dp))
   end do
   varx_fit = mts_varx(series, exogenous, 1, 1)
   call check(varx_fit%info == 0 .and. all(shape(varx_fit%ar) == [2, 2, 1]) .and. &
      size(varx_fit%exogenous, 3) == 2, 'VARX estimation dimensions')
   do t = 2, 120
      expected = varx_fit%intercept + matmul(varx_fit%ar(:, :, 1), series(t - 1, :)) + &
         matmul(varx_fit%exogenous(:, :, 0), exogenous(t, :)) + &
         matmul(varx_fit%exogenous(:, :, 1), exogenous(t - 1, :))
      call check(maxval(abs(series(t, :) - expected - varx_fit%residuals(t - 1, :))) < 1.0e-11_dp, &
         'VARX fitted residual identity')
   end do
   future_exogenous(:, 1) = [0.3_dp, -0.1_dp, 0.2_dp, 0.0_dp]
   forecast = mts_varx_forecast(varx_fit, series, exogenous, future_exogenous, 4)
   expected = varx_fit%intercept + matmul(varx_fit%ar(:, :, 1), series(120, :)) + &
      matmul(varx_fit%exogenous(:, :, 0), future_exogenous(1, :)) + &
      matmul(varx_fit%exogenous(:, :, 1), exogenous(120, :))
   call check(forecast%info == 0 .and. maxval(abs(forecast%mean(1, :) - expected)) < 1.0e-11_dp, &
      'VARX forecast recursion')
   varx_response = mts_varx_irf(varx_fit, 3)
   call check(varx_response%info == 0 .and. maxval(abs(varx_response%exogenous(:, :, 1) - &
      varx_fit%exogenous(:, :, 0))) < 1.0e-12_dp .and. &
      maxval(abs(varx_response%exogenous(:, :, 2) - varx_fit%exogenous(:, :, 1) - &
      matmul(varx_fit%ar(:, :, 1), varx_fit%exogenous(:, :, 0)))) < 1.0e-12_dp, &
      'VARX exogenous impulse responses')
   varx_order = mts_varx_order(series, exogenous, 1, 1)
   call check(varx_order%info == 0 .and. all(varx_order%aic_order >= 0) .and. &
      all(varx_order%aic_order <= 1), 'VARX order selection')

   allocate(varma_model%ar(2, 2, 1), varma_model%ma(2, 2, 1))
   allocate(varma_model%intercept(2), varma_model%sigma(2, 2))
   varma_model%ar(:, :, 1) = reshape([0.45_dp, -0.1_dp, 0.15_dp, 0.35_dp], [2, 2])
   varma_model%ma(:, :, 1) = reshape([0.2_dp, 0.05_dp, -0.1_dp, 0.15_dp], [2, 2])
   varma_model%intercept = [0.1_dp, -0.05_dp]
   varma_model%sigma = reshape([1.0_dp, 0.25_dp, 0.25_dp, 0.8_dp], [2, 2])
   psi = mts_varma_psi(varma_model, 3)
   call check(maxval(abs(psi(:, :, 2) - varma_model%ar(:, :, 1) + &
      varma_model%ma(:, :, 1))) < 1.0e-12_dp .and. &
      maxval(abs(psi(:, :, 3) - matmul(varma_model%ar(:, :, 1), psi(:, :, 2)))) < 1.0e-12_dp, &
      'VARMA PSI weights')
   pi_weights = mts_varma_pi(varma_model, 3)
   call check(maxval(abs(pi_weights(:, :, 2) - varma_model%ar(:, :, 1) - &
      varma_model%ma(:, :, 1))) < 1.0e-12_dp .and. &
      maxval(abs(pi_weights(:, :, 3) - &
      matmul(varma_model%ma(:, :, 1), pi_weights(:, :, 2)))) < 1.0e-12_dp, &
      'VARMA PI weights')
   do t = 1, 80
      innovations(t, :) = [sin(1.11_dp*real(t, dp)), 0.7_dp*cos(0.83_dp*real(t, dp))]
   end do
   varma_simulation = mts_varma_simulate_from_innovations(varma_model, innovations, 10)
   call check(varma_simulation%info == 0 .and. all(shape(varma_simulation%series) == [70, 2]), &
      'VARMA supplied-innovation simulation')
   forecast = mts_varma_forecast(varma_model, varma_simulation%series, &
      varma_simulation%innovations, 4)
   expected = varma_model%intercept + matmul(varma_model%ar(:, :, 1), &
      varma_simulation%series(70, :)) - matmul(varma_model%ma(:, :, 1), &
      varma_simulation%innovations(70, :))
   call check(forecast%info == 0 .and. maxval(abs(forecast%mean(1, :) - expected)) < 1.0e-12_dp, &
      'VARMA forecast recursion')
   varma_covariance = mts_varma_covariance(varma_model, 3, 80)
   call check(varma_covariance%info == 0 .and. &
      maxval(abs([varma_covariance%correlation(1, 1, 1), &
      varma_covariance%correlation(2, 2, 1)] - 1.0_dp)) < 1.0e-12_dp, &
      'VARMA theoretical covariance')
   irf = mts_varma_irf(varma_model, 3)
   call check(irf%info == 0 .and. maxval(abs(irf%psi - psi)) < 1.0e-12_dp, &
      'VARMA impulse responses')
   varma_fit = mts_varma_fit(varma_simulation%series, 1, 1, &
      max_iterations=100, tolerance=1.0e-5_dp)
   call check(varma_fit%model%info == 0 .and. ieee_is_finite(varma_fit%log_likelihood) .and. &
      all(shape(varma_fit%model%ar) == [2, 2, 1]) .and. &
      all(shape(varma_fit%model%ma) == [2, 2, 1]), 'conditional VARMA estimation')
   do t = 2, size(varma_simulation%series, 1)
      expected = varma_simulation%series(t, :) - varma_fit%model%intercept - &
         matmul(varma_fit%model%ar(:, :, 1), varma_simulation%series(t - 1, :)) + &
         matmul(varma_fit%model%ma(:, :, 1), varma_fit%residuals(t - 1, :))
      call check(maxval(abs(varma_fit%residuals(t, :) - expected)) < 1.0e-11_dp, &
         'VARMA fitted residual identity')
   end do
   call check(varma_fit%stationarity_radius >= 0.0_dp .and. &
      varma_fit%invertibility_radius >= 0.0_dp .and. ieee_is_finite(varma_fit%aic) .and. &
      ieee_is_finite(varma_fit%bic), 'VARMA inference and diagnostics')
   varma_refinement = mts_refine_varma(varma_simulation%series, 1, 1, &
      max_steps=0, max_iterations=80, tolerance=1.0e-5_dp)
   call check(varma_refinement%info == 0 .and. varma_refinement%steps == 0 .and. &
      size(varma_refinement%active_count) == 1 .and. &
      maxval(varma_refinement%active_count) == size(varma_refinement%fit%coefficients), &
      'VARMA refinement initialization')
   selected_varma = mts_varmas_fit(varma_simulation%series, [2], [2], &
      max_iterations=80, tolerance=1.0e-5_dp)
   call check(selected_varma%model%info == 0 .and. &
      maxval(abs(selected_varma%model%ar(:, :, 1))) == 0.0_dp .and. &
      maxval(abs(selected_varma%model%ma(:, :, 1))) == 0.0_dp .and. &
      all(selected_varma%ar_lags == [2]) .and. all(selected_varma%ma_lags == [2]), &
      'selected-lag VARMA estimation')

   component(1, 1, 1) = 0.3_dp
   seasonal_component(1, 1, 1) = 0.2_dp
   scalar_sigma(1, 1) = 1.0_dp
   seasonal_model = mts_svarma_expand(component, seasonal_component, &
      zero_component, zero_component, [0.1_dp], scalar_sigma, 4)
   call check(seasonal_model%info == 0 .and. &
      abs(seasonal_model%expanded%ar(1, 1, 1) - 0.3_dp) < 1.0e-12_dp .and. &
      abs(seasonal_model%expanded%ar(1, 1, 4) - 0.2_dp) < 1.0e-12_dp .and. &
      abs(seasonal_model%expanded%ar(1, 1, 5) + 0.06_dp) < 1.0e-12_dp, &
      'multiplicative seasonal VARMA expansion')
   do t = 1, 100
      seasonal_innovations(t, 1) = sin(1.07_dp*real(t, dp)) + &
         0.13_dp*cos(0.019_dp*real(t*t, dp))
   end do
   seasonal_simulation = mts_varma_simulate_from_innovations(&
      seasonal_model%expanded, seasonal_innovations, 10)
   seasonal_fit = mts_svarma_fit(seasonal_simulation%series, 1, 0, 1, 0, 4, &
      max_iterations=100, tolerance=1.0e-6_dp)
   call check(seasonal_fit%model%info == 0 .and. ieee_is_finite(seasonal_fit%log_likelihood) .and. &
      all(shape(seasonal_fit%residuals) == [90, 1]) .and. &
      abs(seasonal_fit%model%expanded%ar(1, 1, 5) + &
      seasonal_fit%model%regular_ar(1, 1, 1)* &
      seasonal_fit%model%seasonal_ar(1, 1, 1)) < 1.0e-12_dp, &
      'multiplicative seasonal VARMA estimation')
   call check(seasonal_fit%stationarity_radius >= 0.0_dp .and. &
      seasonal_fit%invertibility_radius >= 0.0_dp .and. ieee_is_finite(seasonal_fit%aic), &
      'seasonal VARMA inference and diagnostics')
   seasonal_refinement = mts_refine_svarma(seasonal_simulation%series, 1, 0, 1, 0, 4, &
      threshold=1.0e6_dp, protected=[.true., .true., .false.], max_steps=1, &
      max_iterations=80, tolerance=1.0e-6_dp)
   call check(seasonal_refinement%info == 0 .and. seasonal_refinement%steps == 1 .and. &
      all(seasonal_refinement%active_count == [3, 2]) .and. &
      all(seasonal_refinement%removed_index == [3]) .and. &
      seasonal_refinement%fit%coefficients(3) == 0.0_dp, &
      'seasonal VARMA backward refinement')

   walk = 0.0_dp
   spread_value = 0.2_dp
   do t = 1, 140
      walk = walk + 0.15_dp*sin(0.47_dp*real(t, dp)) + &
         0.04_dp*cos(0.013_dp*real(t*t, dp))
      spread_value = 0.62_dp*spread_value + 0.12_dp*cos(0.83_dp*real(t, dp))
      vecm_series(t, :) = [walk + 0.5_dp*spread_value + 0.025_dp*sin(1.31_dp*real(t, dp)), &
         walk - 0.5_dp*spread_value + 0.02_dp*cos(1.17_dp*real(t, dp))]
   end do
   cointegration(:, 1) = [1.0_dp, -1.0_dp]
   vecm_fit = mts_vecm_fit(vecm_series, 2, 1, .true., cointegration)
   call check(vecm_fit%info == 0 .and. vecm_fit%rank == 1 .and. &
      all(shape(vecm_fit%gamma) == [2, 2, 1]) .and. &
      all(shape(vecm_fit%residuals) == [138, 2]), 'VECM supplied-vector estimation')
   do t = 3, 140
      difference_now = vecm_series(t, :) - vecm_series(t - 1, :) - vecm_fit%intercept - &
         matmul(vecm_fit%loading, matmul(transpose(vecm_fit%cointegration), &
         vecm_series(t - 1, :))) - matmul(vecm_fit%gamma(:, :, 1), &
         vecm_series(t - 1, :) - vecm_series(t - 2, :))
      call check(maxval(abs(difference_now - vecm_fit%residuals(t - 2, :))) < 1.0e-11_dp, &
         'VECM residual identity')
   end do
   expected = vecm_fit%level_var%intercept + &
      matmul(vecm_fit%level_var%ar(:, :, 1), vecm_series(140, :)) + &
      matmul(vecm_fit%level_var%ar(:, :, 2), vecm_series(139, :))
   vecm_prediction = mts_vecm_forecast(vecm_fit, vecm_series, 4)
   call check(vecm_prediction%info == 0 .and. &
      maxval(abs(vecm_prediction%level%mean(1, :) - expected)) < 1.0e-11_dp .and. &
      maxval(abs(vecm_prediction%difference(1, :) - expected + vecm_series(140, :))) < 1.0e-11_dp, &
      'VECM level and difference forecasts')
   johansen_vecm = mts_vecm_fit(vecm_series, 2, 1)
   call check(johansen_vecm%info == 0 .and. all(shape(johansen_vecm%cointegration) == [2, 1]), &
      'VECM Johansen-vector estimation')
   vecm_mask = .true.
   vecm_mask(4, 2) = .false.
   johansen_vecm = mts_vecm_fit(vecm_series, 2, 1, .true., cointegration, vecm_mask)
   call check(johansen_vecm%info == 0 .and. johansen_vecm%coefficients(4, 2) == 0.0_dp, &
      'VECM equation-specific fixed mask')
   filter_data = reshape([1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, &
      0.5_dp, -1.0_dp, 1.5_dp, 0.5_dp], [4, 2])
   filter_weights(:, :, 1) = reshape([0.2_dp, 0.0_dp, 0.1_dp, 0.3_dp], [2, 2])
   filtered_data = mts_matrix_filter(filter_data, filter_weights)
   call check(maxval(abs(filtered_data(1, :) - filter_data(1, :))) < 1.0e-12_dp .and. &
      maxval(abs(filtered_data(2, :) - filter_data(2, :) - &
      matmul(filter_weights(:, :, 1), filtered_data(1, :)))) < 1.0e-12_dp, &
      'multivariate matrix-polynomial filtering')
   ecm_process = matmul(vecm_series, cointegration)
   ecm_known = mts_ecm_known_fit(vecm_series, 2, ecm_process, .true.)
   call check(ecm_known%info == 0 .and. all(shape(ecm_known%loading) == [2, 1]) .and. &
      all(shape(ecm_known%gamma) == [2, 2, 1]) .and. &
      all(shape(ecm_known%residuals) == [139, 2]), 'known-process ECM estimation')
   ecm_known_refined = mts_ecm_known_refine(vecm_series, ecm_process, ecm_known, 0.0_dp)
   call check(ecm_known_refined%info == 0 .and. &
      maxval(abs(ecm_known_refined%coefficients - ecm_known%coefficients)) < 1.0e-10_dp, &
      'known-process ECM refinement')
   ecm_normalized = mts_ecm_normalized_fit(vecm_series, 2, cointegration, .true., &
      max_iterations=100, tolerance=1.0e-5_dp)
   call check((ecm_normalized%info == 0 .or. ecm_normalized%info >= 100) .and. &
      ecm_normalized%cointegration(1, 1) == 1.0_dp .and. &
      all(shape(ecm_normalized%loading) == [2, 1]) .and. &
      all(shape(ecm_normalized%gamma) == [2, 2, 1]) .and. &
      all(shape(ecm_normalized%residuals) == [139, 2]), 'normalized joint ECM estimation')
   ecm_normalized_refined = mts_ecm_normalized_refine(vecm_series, ecm_normalized, 0.0_dp, &
      max_iterations=40, tolerance=1.0e-5_dp)
   call check(ecm_normalized_refined%info == 0 .or. ecm_normalized_refined%info >= 100, &
      'normalized ECM refinement')

   factor_one = 0.0_dp
   factor_two = 0.0_dp
   do t = 1, 160
      factor_one = 0.7_dp*factor_one + sin(0.43_dp*real(t, dp))
      factor_two = -0.35_dp*factor_two + 0.8_dp*cos(0.67_dp*real(t, dp))
      factor_series(t, :) = [factor_one + 0.04_dp*sin(1.1_dp*real(t, dp)), &
         0.8_dp*factor_one + 0.3_dp*factor_two + 0.05_dp*cos(1.3_dp*real(t, dp)), &
         -0.6_dp*factor_one + factor_two + 0.04_dp*sin(1.7_dp*real(t, dp)), &
         0.9_dp*factor_two + 0.03_dp*cos(1.9_dp*real(t, dp)), &
         0.4_dp*factor_one - 0.7_dp*factor_two + 0.04_dp*sin(2.1_dp*real(t, dp)), &
         -0.3_dp*factor_one - 0.5_dp*factor_two + 0.03_dp*cos(2.3_dp*real(t, dp))]
   end do
   factor_model = mts_factor_fit(factor_series, 2, 4)
   identity_two = matmul(transpose(factor_model%loadings), factor_model%loadings)
   call check(factor_model%info == 0 .and. factor_model%factors == 2 .and. &
      maxval(abs(identity_two - reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))) < 1.0e-11_dp .and. &
      maxval(abs(factor_model%common + factor_model%residuals - factor_series)) < 1.0e-11_dp, &
      'principal-component factor extraction')
   call check(abs(sum(factor_model%explained) - 1.0_dp) < 1.0e-11_dp .and. &
      factor_model%ic1_factors >= 0 .and. factor_model%ic1_factors <= 4 .and. &
      factor_model%ic2_factors >= 0 .and. factor_model%ic2_factors <= 4 .and. &
      factor_model%ic3_factors >= 0 .and. factor_model%ic3_factors <= 4, &
      'factor explained variance and selection criteria')
   factor_prediction = mts_factor_forecast(factor_model, 1, 4)
   call check(factor_prediction%info == 0 .and. all(shape(factor_prediction%mean) == [4, 6]) .and. &
      maxval(abs(factor_prediction%mean(1, :) - factor_model%mean - &
      matmul(factor_prediction%factor_forecast%mean(1, :), &
      transpose(factor_model%loadings*spread(factor_model%scale, 2, 2))))) < 1.0e-11_dp .and. &
      maxval(abs(factor_prediction%covariance(:, :, 4) - &
      transpose(factor_prediction%covariance(:, :, 4)))) < 1.0e-11_dp .and. &
      all(factor_prediction%standard_error > 0.0_dp), 'factor VAR forecasting')
   factor_constraint = reshape([1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp], [6, 3])
   constrained_factor = mts_constrained_factor_fit(factor_series, factor_constraint, 2)
   score_covariance = matmul(transpose(constrained_factor%factor_model%scores), &
      constrained_factor%factor_model%scores)/159.0_dp
   call check(constrained_factor%info == 0 .and. &
      maxval(abs(constrained_factor%factor_model%loadings - &
      matmul(factor_constraint, constrained_factor%omega))) < 1.0e-11_dp .and. &
      maxval(abs(score_covariance - reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))) < 1.0e-11_dp .and. &
      maxval(abs(constrained_factor%factor_model%common + &
      constrained_factor%factor_model%residuals - factor_series)) < 1.0e-11_dp, &
      'constrained factor estimation identities')
   call check(constrained_factor%explained >= 0.0_dp .and. &
      maxval(abs(constrained_factor%psi - transpose(constrained_factor%psi))) < 1.0e-11_dp, &
      'constrained factor diagnostics')
   constrained_prediction = mts_constrained_factor_forecast(constrained_factor, 1, 3)
   call check(constrained_prediction%info == 0 .and. &
      all(shape(constrained_prediction%mean) == [3, 6]) .and. &
      all(constrained_prediction%standard_error > 0.0_dp), 'constrained factor forecasting')

   bvar_prior = mts_minnesota_prior(factor_model%scores, 1)
   call check(bvar_prior%info == 0 .and. &
      bvar_prior%mean(2, 1) == 1.0_dp .and. bvar_prior%mean(3, 2) == 1.0_dp .and. &
      all([(bvar_prior%precision(t, t) > 0.0_dp, t=1, 3)]), &
      'Minnesota-style BVAR prior')
   bvar_prior%mean = 0.0_dp
   bvar_prior%precision = 0.0_dp
   bvar_prior%scale = 0.0_dp
   do t = 1, 3
      bvar_prior%precision(t, t) = 1.0e-10_dp
   end do
   do t = 1, 2
      bvar_prior%scale(t, t) = 1.0e-10_dp
   end do
   selected = mts_var(factor_model%scores, 1)
   bvar_fit = mts_bvar_fit(factor_model%scores, 1, bvar_prior)
   call check(bvar_fit%info == 0 .and. &
      maxval(abs(bvar_fit%model%ar - selected%ar)) < 1.0e-9_dp .and. &
      maxval(abs(bvar_fit%model%intercept - selected%intercept)) < 1.0e-9_dp .and. &
      all(bvar_fit%standard_errors > 0.0_dp), 'diffuse-prior BVAR convergence to OLS')
   strong_prior = bvar_prior
   do t = 1, 3
      strong_prior%precision(t, t) = 1.0e6_dp
   end do
   shrunk_bvar = mts_bvar_fit(factor_model%scores, 1, strong_prior)
   call check(shrunk_bvar%info == 0 .and. &
      norm2(shrunk_bvar%coefficient_mean) < norm2(bvar_fit%coefficient_mean), &
      'strong-prior BVAR shrinkage')
   bvar_prediction = mts_var_forecast(bvar_fit%model, factor_model%scores, 3)
   call check(bvar_prediction%info == 0 .and. all(shape(bvar_prediction%mean) == [3, 2]) .and. &
      all(bvar_prediction%standard_error > 0.0_dp), 'BVAR shared VAR forecasting')

   conditional_variance = 1.0_dp
   do t = 1, 240
      common_shock = sin(0.71_dp*real(t, dp)) + 0.35_dp*cos(0.037_dp*real(t*t, dp))
      conditional_variance = 0.08_dp + 0.84_dp*conditional_variance + 0.11_dp*common_shock**2
      volatility_returns(t, :) = sqrt(conditional_variance)* &
         [common_shock + 0.2_dp*sin(1.31_dp*real(t, dp)), &
         0.7_dp*common_shock + 0.45_dp*cos(0.93_dp*real(t, dp)), &
         -0.5_dp*common_shock + 0.55_dp*sin(1.17_dp*real(t, dp))]
   end do
   common_volatility = mts_common_volatility(volatility_returns, 4, 0, .true., [2, 4])
   whitened_covariance = matmul(transpose(common_volatility%whitened), &
      common_volatility%whitened)/239.0_dp
   call check(common_volatility%info == 0 .and. &
      maxval(abs(whitened_covariance - reshape([1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3]))) < 1.0e-10_dp .and. &
      abs(sum(common_volatility%proportions) - 1.0_dp) < 1.0e-11_dp, &
      'common-volatility whitening and eigensystem')
   call check(maxval(abs([(norm2(common_volatility%directions(:, t)), t=1, 3)] - 1.0_dp)) < &
      1.0e-11_dp .and. all(shape(common_volatility%components) == [240, 3]) .and. &
      all(common_volatility%arch_statistic >= 0.0_dp) .and. &
      all(common_volatility%arch_p_value >= 0.0_dp) .and. &
      all(common_volatility%arch_p_value <= 1.0_dp), 'common-volatility directions and ARCH tests')
   prewhitened_volatility = mts_common_volatility(volatility_returns, 2, 1, arch_lags=[2])
   call check(prewhitened_volatility%info == 0 .and. &
      all(shape(prewhitened_volatility%residuals) == [239, 3]) .and. &
      prewhitened_volatility%prewhiten_order == 1, 'common-volatility VAR prewhitening')

   do t = 1, 160
      first_variance = 0.6_dp + 0.25_dp*sin(0.09_dp*real(t, dp))**2
      second_variance = 0.8_dp + 0.3_dp*cos(0.07_dp*real(t, dp))**2
      mch_standardized(t, :) = [sin(0.73_dp*real(t, dp)) + 0.2_dp*cos(0.031_dp*real(t*t, dp)), &
         0.8_dp*cos(0.61_dp*real(t, dp)) + 0.15_dp*sin(0.027_dp*real(t*t, dp))]
      mch_covariance(:, :, t) = 0.0_dp
      mch_covariance(1, 1, t) = first_variance
      mch_covariance(2, 2, t) = second_variance
      mch_residuals(t, :) = [sqrt(first_variance)*mch_standardized(t, 1), &
         sqrt(second_variance)*mch_standardized(t, 2)]
   end do
   mch_diagnostic = mts_mch_diagnostic(mch_residuals, mch_covariance, 5)
   call check(mch_diagnostic%info == 0 .and. &
      maxval(abs(mch_diagnostic%standardized_residuals - mch_standardized)) < 1.0e-11_dp .and. &
      maxval(abs(mch_diagnostic%radial_residual - sum(mch_standardized**2, 2) + 2.0_dp)) < 1.0e-11_dp, &
      'MCH covariance-path standardization')
   call check(mch_diagnostic%radial_q >= 0.0_dp .and. mch_diagnostic%rank_q >= 0.0_dp .and. &
      mch_diagnostic%multivariate_q >= 0.0_dp .and. mch_diagnostic%robust_q >= 0.0_dp .and. &
      all(mch_diagnostic%p_value >= 0.0_dp) .and. all(mch_diagnostic%p_value <= 1.0_dp) .and. &
      all(mch_diagnostic%degrees_of_freedom == [5, 5, 20, 20]) .and. &
      mch_diagnostic%robust_observations < 160, 'MCHdiag portmanteau diagnostics')

   bekk_constant = reshape([0.22_dp, 0.03_dp, 0.0_dp, 0.25_dp], [2, 2])
   bekk_arch = reshape([0.16_dp, 0.0_dp, 0.0_dp, 0.12_dp], [2, 2])
   bekk_garch = reshape([0.74_dp, 0.0_dp, 0.0_dp, 0.7_dp], [2, 2])
   bekk_covariance(:, :, 1) = matmul(bekk_constant, transpose(bekk_constant))
   do t = 1, 120
      cholesky11 = sqrt(bekk_covariance(1, 1, t))
      cholesky21 = bekk_covariance(2, 1, t)/cholesky11
      cholesky22 = sqrt(bekk_covariance(2, 2, t) - cholesky21**2)
      bekk_shock = [sin(0.79_dp*real(t, dp)) + 0.15_dp*cos(0.031_dp*real(t*t, dp)), &
         0.8_dp*cos(0.67_dp*real(t, dp)) + 0.12_dp*sin(0.027_dp*real(t*t, dp))]
      bekk_residual = [cholesky11*bekk_shock(1), &
         cholesky21*bekk_shock(1) + cholesky22*bekk_shock(2)]
      bekk_returns(t, :) = [0.02_dp, -0.01_dp] + bekk_residual
      if (t < 120) then
         bekk_shock = matmul(bekk_arch, bekk_residual)
         bekk_covariance(:, :, t + 1) = matmul(bekk_constant, transpose(bekk_constant)) + &
            spread(bekk_shock, 2, 2)*spread(bekk_shock, 1, 2) + &
            matmul(bekk_garch, matmul(bekk_covariance(:, :, t), transpose(bekk_garch)))
      end if
   end do
   bekk_initial = [0.02_dp, -0.01_dp, 0.22_dp, 0.03_dp, 0.25_dp, &
      0.16_dp, 0.0_dp, 0.0_dp, 0.12_dp, 0.74_dp, 0.0_dp, 0.0_dp, 0.7_dp]
   bekk_mask = .false.
   bekk_mask([1, 2, 3, 5]) = .true.
   bekk_fit = mts_bekk_fit(bekk_returns, .true., bekk_initial, bekk_mask, 80, 1.0e-6_dp)
   call check(allocated(bekk_fit%covariance) .and. ieee_is_finite(bekk_fit%log_likelihood) .and. &
      all(shape(bekk_fit%covariance) == [2, 2, 120]) .and. bekk_fit%persistence < 1.0_dp, &
      'BEKK Gaussian estimation and diagnostics')
   do t = 2, 120
      bekk_shock = matmul(bekk_fit%arch, bekk_fit%residuals(t - 1, :))
      call check(maxval(abs(bekk_fit%covariance(:, :, t) - &
         matmul(bekk_fit%constant, transpose(bekk_fit%constant)) - &
         spread(bekk_shock, 2, 2)*spread(bekk_shock, 1, 2) - &
         matmul(bekk_fit%garch, matmul(bekk_fit%covariance(:, :, t - 1), &
         transpose(bekk_fit%garch))))) < 1.0e-11_dp, 'BEKK covariance recursion')
   end do
   bekk_prediction = mts_bekk_forecast(bekk_fit, 3)
   call check(all(shape(bekk_prediction) == [2, 2, 3]) .and. &
      maxval(abs(bekk_prediction(:, :, 1) - transpose(bekk_prediction(:, :, 1)))) < 1.0e-12_dp .and. &
      all([(bekk_prediction(t, t, 3) > 0.0_dp, t=1, 2)]), 'BEKK covariance forecasting')
   mch_diagnostic = mts_mch_diagnostic(bekk_fit%residuals, bekk_fit%covariance, 4)
   call check(mch_diagnostic%info == 0, 'BEKK covariance-path diagnostics')

   do t = 1, 160
      dcc_variance(t, :) = [0.5_dp + 0.2_dp*sin(0.08_dp*real(t, dp))**2, &
         0.7_dp + 0.25_dp*cos(0.06_dp*real(t, dp))**2]
   end do
   dcc_fit = mts_dcc_fit(mch_standardized, dcc_variance, [0.04_dp, 0.9_dp], 100, 1.0e-6_dp)
   call check(allocated(dcc_fit%q) .and. ieee_is_finite(dcc_fit%log_likelihood) .and. &
      dcc_fit%arch >= 0.0_dp .and. dcc_fit%garch >= 0.0_dp .and. &
      dcc_fit%arch + dcc_fit%garch < 1.0_dp .and. &
      maxval(abs([(dcc_fit%correlation(t, t, 160), t=1, 2)] - 1.0_dp)) < 1.0e-12_dp .and. &
      maxval(abs([(dcc_fit%covariance(t, t, 160), t=1, 2)] - dcc_variance(160, :))) < 1.0e-11_dp, &
      'Gaussian DCC estimation and covariance path')
   do t = 2, 160
      dcc_expected = (1.0_dp - dcc_fit%arch - dcc_fit%garch)*dcc_fit%unconditional + &
         dcc_fit%arch*spread(mch_standardized(t - 1, :), 2, 2)* &
         spread(mch_standardized(t - 1, :), 1, 2) + dcc_fit%garch*dcc_fit%q(:, :, t - 1)
      call check(maxval(abs(dcc_fit%q(:, :, t) - dcc_expected)) < 1.0e-11_dp, &
         'DCC Q recursion')
   end do
   future_variance = reshape([0.8_dp, 0.85_dp, 0.9_dp, 0.7_dp, 0.75_dp, 0.8_dp], [3, 2])
   dcc_prediction = mts_dcc_forecast(dcc_fit, future_variance)
   call check(all(shape(dcc_prediction) == [2, 2, 3]) .and. &
      maxval(abs([(dcc_prediction(t, t, 3), t=1, 2)] - future_variance(3, :))) < 1.0e-11_dp .and. &
      maxval(abs(dcc_prediction(:, :, 3) - transpose(dcc_prediction(:, :, 3)))) < 1.0e-12_dp, &
      'DCC covariance forecasting')
   mch_diagnostic = mts_mch_diagnostic(dcc_fit%residuals, dcc_fit%covariance, 4)
   call check(mch_diagnostic%info == 0, 'DCC covariance-path diagnostics')

   adcc_fit = mts_adcc_fit(mch_standardized, dcc_variance, [0.04_dp, 0.85_dp, 0.04_dp], &
      100, 1.0e-6_dp)
   call check(allocated(adcc_fit%q) .and. ieee_is_finite(adcc_fit%log_likelihood) .and. &
      adcc_fit%arch >= 0.0_dp .and. adcc_fit%garch >= 0.0_dp .and. &
      adcc_fit%asymmetry >= 0.0_dp .and. &
      adcc_fit%arch + adcc_fit%garch + adcc_fit%asymmetry < 1.0_dp .and. &
      .not. adcc_fit%student_t .and. &
      maxval(abs([(adcc_fit%covariance(t, t, 160), t=1, 2)] - dcc_variance(160, :))) < 1.0e-11_dp, &
      'Gaussian ADCC estimation and covariance path')
   bekk_shock = min(mch_standardized(159, :), 0.0_dp)
   dcc_expected = (1.0_dp - adcc_fit%arch - adcc_fit%garch)*adcc_fit%unconditional - &
      adcc_fit%asymmetry*adcc_fit%negative_unconditional + adcc_fit%arch* &
      spread(mch_standardized(159, :), 2, 2)*spread(mch_standardized(159, :), 1, 2) + &
      adcc_fit%garch*adcc_fit%q(:, :, 159) + adcc_fit%asymmetry* &
      spread(bekk_shock, 2, 2)*spread(bekk_shock, 1, 2)
   call check(maxval(abs(adcc_fit%q(:, :, 160) - dcc_expected)) < 1.0e-11_dp, &
      'Gaussian ADCC negative-shock recursion')
   adcc_prediction = mts_adcc_forecast(adcc_fit, future_variance)
   call check(all(shape(adcc_prediction) == [2, 2, 3]) .and. &
      maxval(abs([(adcc_prediction(t, t, 3), t=1, 2)] - future_variance(3, :))) < 1.0e-11_dp, &
      'Gaussian ADCC covariance forecasting')

   adcc_t_fit = mts_adcc_t_fit(mch_standardized, dcc_variance, &
      [0.04_dp, 0.85_dp, 0.04_dp, 8.0_dp], 100, 1.0e-6_dp)
   call check(allocated(adcc_t_fit%q) .and. ieee_is_finite(adcc_t_fit%log_likelihood) .and. &
      adcc_t_fit%student_t .and. adcc_t_fit%degrees_of_freedom > 2.0_dp .and. &
      adcc_t_fit%arch + adcc_t_fit%garch + adcc_t_fit%asymmetry < 1.0_dp .and. &
      size(adcc_t_fit%standard_errors) == 4 .and. &
      maxval(abs([(adcc_t_fit%correlation(t, t, 160), t=1, 2)] - 1.0_dp)) < 1.0e-12_dp, &
      'Student-t ADCC estimation and inference')
   adcc_prediction = mts_adcc_forecast(adcc_t_fit, future_variance)
   call check(all(shape(adcc_prediction) == [2, 2, 3]) .and. &
      maxval(abs(adcc_prediction(:, :, 3) - transpose(adcc_prediction(:, :, 3)))) < 1.0e-12_dp, &
      'Student-t ADCC covariance forecasting')
   mch_diagnostic = mts_mch_diagnostic(adcc_t_fit%residuals, adcc_t_fit%covariance, 4)
   call check(mch_diagnostic%info == 0, 'Student-t ADCC covariance-path diagnostics')

   tse_tsui_fit = mts_tse_tsui_fit(mch_standardized, 10, dcc_variance, &
      [0.85_dp, 0.05_dp], 100, 1.0e-6_dp)
   call check(allocated(tse_tsui_fit%correlation) .and. &
      ieee_is_finite(tse_tsui_fit%log_likelihood) .and. &
      tse_tsui_fit%previous_weight >= 0.0_dp .and. tse_tsui_fit%rolling_weight >= 0.0_dp .and. &
      tse_tsui_fit%previous_weight + tse_tsui_fit%rolling_weight < 1.0_dp .and. &
      maxval(abs(tse_tsui_fit%correlation(:, :, 1) - tse_tsui_fit%unconditional)) < 1.0e-12_dp .and. &
      maxval(abs([(tse_tsui_fit%covariance(t, t, 160), t=1, 2)] - &
      dcc_variance(160, :))) < 1.0e-11_dp, 'Gaussian Tse-Tsui DCC estimation')
   dcc_expected = (1.0_dp - tse_tsui_fit%previous_weight - tse_tsui_fit%rolling_weight)* &
      tse_tsui_fit%unconditional + tse_tsui_fit%previous_weight* &
      tse_tsui_fit%correlation(:, :, 10) + tse_tsui_fit%rolling_weight* &
      tse_tsui_fit%local_correlation(:, :, 11)
   call check(maxval(abs(tse_tsui_fit%correlation(:, :, 11) - dcc_expected)) < 1.0e-11_dp, &
      'Tse-Tsui rolling-correlation recursion')
   tse_tsui_prediction = mts_tse_tsui_forecast(tse_tsui_fit, future_variance)
   call check(all(shape(tse_tsui_prediction) == [2, 2, 3]) .and. &
      maxval(abs([(tse_tsui_prediction(t, t, 3), t=1, 2)] - &
      future_variance(3, :))) < 1.0e-11_dp, 'Gaussian Tse-Tsui covariance forecasting')

   tse_tsui_t_fit = mts_tse_tsui_t_fit(mch_standardized, 10, dcc_variance, &
      [0.85_dp, 0.05_dp, 8.0_dp], 100, 1.0e-6_dp)
   call check(allocated(tse_tsui_t_fit%correlation) .and. &
      ieee_is_finite(tse_tsui_t_fit%log_likelihood) .and. tse_tsui_t_fit%student_t .and. &
      tse_tsui_t_fit%degrees_of_freedom > 2.0_dp .and. &
      size(tse_tsui_t_fit%standard_errors) == 3 .and. &
      maxval(abs([(tse_tsui_t_fit%correlation(t, t, 160), t=1, 2)] - 1.0_dp)) < 1.0e-12_dp, &
      'Student-t Tse-Tsui DCC estimation and inference')
   tse_tsui_prediction = mts_tse_tsui_forecast(tse_tsui_t_fit, future_variance)
   call check(all(shape(tse_tsui_prediction) == [2, 2, 3]) .and. &
      maxval(abs(tse_tsui_prediction(:, :, 3) - transpose(tse_tsui_prediction(:, :, 3)))) < 1.0e-12_dp, &
      'Student-t Tse-Tsui covariance forecasting')
   mch_diagnostic = mts_mch_diagnostic(tse_tsui_t_fit%residuals, tse_tsui_t_fit%covariance, 4)
   call check(mch_diagnostic%info == 0, 'Tse-Tsui covariance-path diagnostics')

   ewma_fit = mts_ewma_fit(volatility_returns, 0.94_dp)
   call check(ewma_fit%info == 0 .and. ewma_fit%decay == 0.94_dp .and. &
      all(shape(ewma_fit%covariance) == [3, 3, 240]) .and. ieee_is_finite(ewma_fit%log_likelihood), &
      'fixed-decay EWMA covariance estimation')
   call check(maxval(abs(ewma_fit%covariance(:, :, 2) - 0.06_dp* &
      spread(ewma_fit%residuals(1, :), 2, 3)*spread(ewma_fit%residuals(1, :), 1, 3) - &
      0.94_dp*ewma_fit%covariance(:, :, 1))) < 1.0e-11_dp, 'EWMA covariance recursion')
   estimated_ewma = mts_ewma_fit(volatility_returns, 0.0_dp, 100, 1.0e-7_dp)
   call check(allocated(estimated_ewma%covariance) .and. estimated_ewma%estimated .and. &
      estimated_ewma%decay > 0.0_dp .and. estimated_ewma%decay < 1.0_dp .and. &
      estimated_ewma%standard_error >= 0.0_dp, 'estimated-decay EWMA covariance')
   ewma_prediction = mts_ewma_forecast(ewma_fit, 3)
   call check(all(shape(ewma_prediction) == [3, 3, 3]) .and. &
      maxval(abs(ewma_prediction(:, :, 2) - ewma_prediction(:, :, 1))) < 1.0e-12_dp .and. &
      maxval(abs(ewma_prediction(:, :, 3) - transpose(ewma_prediction(:, :, 3)))) < 1.0e-12_dp, &
      'EWMA covariance forecasting')

   mchol_fit = mts_mchol_fit(volatility_returns, 30, 0.95_dp, 0, 80, 1.0e-5_dp)
   call check(mchol_fit%info == 0 .and. all(shape(mchol_fit%covariance) == [3, 3, 210]) .and. &
      all(shape(mchol_fit%smoothed_coefficients) == [210, 3]) .and. &
      all(mchol_fit%component_variance > 0.0_dp) .and. &
      all(mchol_fit%garch_parameters > 0.0_dp) .and. &
      all(sum(mchol_fit%garch_parameters(:, 2:3), 2) < 1.0_dp), &
      'moving-Cholesky GARCH estimation')
   call check(maxval(abs(mchol_fit%covariance(:, :, 210) - &
      transpose(mchol_fit%covariance(:, :, 210)))) < 1.0e-11_dp .and. &
      all([(mchol_fit%covariance(t, t, 210) > 0.0_dp, t=1, 3)]), &
      'moving-Cholesky covariance reconstruction')
   mchol_prediction = mts_mchol_forecast(mchol_fit, 3)
   call check(all(shape(mchol_prediction) == [3, 3, 3]) .and. &
      all([(mchol_prediction(t, t, 3) > 0.0_dp, t=1, 3)]) .and. &
      maxval(abs(mchol_prediction(:, :, 3) - transpose(mchol_prediction(:, :, 3)))) < 1.0e-11_dp, &
      'moving-Cholesky covariance forecasting')
   mch_diagnostic = mts_mch_diagnostic(mchol_fit%prewhitened(31:, :), mchol_fit%covariance, 4)
   call check(mch_diagnostic%info == 0, 'moving-Cholesky covariance-path diagnostics')
   prewhitened_mchol = mts_mchol_fit(volatility_returns, 30, 0.95_dp, 1, 60, 1.0e-5_dp)
   call check(prewhitened_mchol%info == 0 .and. prewhitened_mchol%var_order == 1 .and. &
      size(prewhitened_mchol%covariance, 3) == 209, 'moving-Cholesky VAR prewhitening')

   grouped_correlation = mts_sccor(factor_series, 160, 80, [3, 3])
   call check(grouped_correlation%info == 0 .and. grouped_correlation%start == 81 .and. &
      grouped_correlation%end == 160 .and. &
      abs(grouped_correlation%constrained(1, 2) - grouped_correlation%constrained(1, 3)) < &
      1.0e-12_dp .and. abs(grouped_correlation%constrained(1, 2) - &
      grouped_correlation%constrained(2, 3)) < 1.0e-12_dp .and. &
      maxval(abs(grouped_correlation%constrained(:3, 4:) - &
      grouped_correlation%constrained(1, 4))) < 1.0e-12_dp .and. &
      maxval(abs([(grouped_correlation%constrained(t, t), t=1, 6)] - 1.0_dp)) < 1.0e-12_dp, &
      'group-constrained sample correlation')
   univariate_arch = mts_arch_test(volatility_returns(:, 1), 5)
   call check(univariate_arch%info == 0 .and. all(univariate_arch%statistic >= 0.0_dp) .and. &
      all(univariate_arch%p_value >= 0.0_dp) .and. all(univariate_arch%p_value <= 1.0_dp), &
      'univariate ARCH serial-dependence tests')
   multivariate_arch = mts_march_test(volatility_returns, 4)
   call check(multivariate_arch%info == 0 .and. all(multivariate_arch%statistic >= 0.0_dp) .and. &
      all(multivariate_arch%p_value >= 0.0_dp) .and. all(multivariate_arch%p_value <= 1.0_dp) .and. &
      all(multivariate_arch%degrees_of_freedom == [4, 4, 36, 36]) .and. &
      multivariate_arch%robust_observations < 240, 'multivariate ARCH diagnostics')

   correlation_fixture = reshape([1.0_dp, 0.3_dp, -0.2_dp, &
      0.3_dp, 1.0_dp, 0.25_dp, -0.2_dp, 0.25_dp, 1.0_dp], [3, 3])
   correlation_angles = mts_correlation_to_angles(correlation_fixture)
   recovered_correlation = mts_angles_to_correlation(correlation_angles, 3)
   call check(size(correlation_angles) == 3 .and. &
      maxval(abs(recovered_correlation - correlation_fixture)) < 1.0e-11_dp, &
      'hyperspherical correlation angle round trip')
   copula_fit = mts_copula_fit(mch_standardized, [1, 1], 10, estimate_baseline=.false., &
      max_iterations=100, tolerance=1.0e-6_dp)
   call check(allocated(copula_fit%correlation) .and. ieee_is_finite(copula_fit%log_likelihood) .and. &
      copula_fit%degrees_of_freedom > 2.0_dp .and. copula_fit%previous_weight >= 0.0_dp .and. &
      copula_fit%local_weight >= 0.0_dp .and. &
      copula_fit%previous_weight + copula_fit%local_weight < 1.0_dp .and. &
      maxval(abs([(copula_fit%correlation(t, t, 160), t=1, 2)] - 1.0_dp)) < 1.0e-12_dp .and. &
      maxval(abs(copula_fit%angles(:10, :) - spread(copula_fit%baseline_angles, 1, 10))) < 1.0e-12_dp, &
      'fixed-baseline dynamic Student-t copula')
   call check(abs(copula_fit%angles(11, 1) - &
      ((1.0_dp - copula_fit%previous_weight - copula_fit%local_weight)* &
      copula_fit%baseline_angles(1) + copula_fit%previous_weight*copula_fit%angles(10, 1) + &
      copula_fit%local_weight*copula_fit%local_angles(11, 1))) < 1.0e-11_dp, &
      'dynamic copula angle recursion')
   estimated_copula = mts_copula_fit(mch_standardized, [1, 1], 10, &
      max_iterations=80, tolerance=1.0e-6_dp)
   call check(allocated(estimated_copula%correlation) .and. &
      estimated_copula%estimates_baseline .and. size(estimated_copula%coefficients) == 4 .and. &
      all(estimated_copula%baseline_angles > 0.0_dp) .and. &
      all(estimated_copula%baseline_angles < acos(-1.0_dp)), &
      'estimated-baseline grouped Student-t copula')
   grouped_copula = mts_copula_fit(common_volatility%whitened, [2, 1], 12, &
      estimate_baseline=.false., max_iterations=60, tolerance=1.0e-6_dp)
   call check(allocated(grouped_copula%correlation) .and. &
      size(grouped_copula%baseline_angles) == 2 .and. &
      maxval(abs(grouped_copula%correlation(1, 3, :) - &
      grouped_copula%correlation(2, 3, :))) < 1.0e-11_dp, &
      'multi-series grouped copula correlations')

   missing_pi(:, :, 1) = reshape([0.45_dp, -0.12_dp, 0.18_dp, 0.38_dp], [2, 2])
   missing_pi(:, :, 2) = reshape([0.08_dp, 0.03_dp, -0.04_dp, 0.06_dp], [2, 2])
   missing_sigma = reshape([1.0_dp, 0.2_dp, 0.2_dp, 0.8_dp], [2, 2])
   missing_constant = [0.1_dp, -0.05_dp]
   missing_series(1, :) = [0.2_dp, -0.1_dp]
   missing_series(2, :) = [0.15_dp, 0.05_dp]
   do t = 3, 50
      missing_series(t, :) = missing_constant + &
         matmul(missing_pi(:, :, 1), missing_series(t - 1, :)) + &
         matmul(missing_pi(:, :, 2), missing_series(t - 2, :))
   end do
   missing_result = mts_var_missing(missing_series, missing_pi, missing_sigma, 25, missing_constant)
   call check(missing_result%info == 0 .and. &
      maxval(abs(missing_result%estimate - missing_series(25, :))) < 1.0e-11_dp .and. &
      maxval(abs(missing_result%completed(25, :) - missing_series(25, :))) < 1.0e-11_dp .and. &
      maxval(abs(missing_result%covariance - transpose(missing_result%covariance))) < 1.0e-12_dp .and. &
      all([(missing_result%covariance(t, t) > 0.0_dp, t=1, 2)]), &
      'fully missing VAR observation smoothing')
   partial_missing_result = mts_var_partial_missing(missing_series, missing_pi, &
      missing_sigma, 25, [.true., .false.], missing_constant)
   call check(partial_missing_result%info == 0 .and. size(partial_missing_result%estimate) == 1 .and. &
      abs(partial_missing_result%estimate(1) - missing_series(25, 1)) < 1.0e-11_dp .and. &
      partial_missing_result%completed(25, 2) == missing_series(25, 2), &
      'partially missing VAR observation smoothing')

   granger_series(1, :) = [0.1_dp, -0.2_dp, 0.05_dp]
   do t = 2, 180
      granger_series(t, 1) = 0.55_dp*granger_series(t - 1, 1) + &
         sin(0.71_dp*real(t, dp)) + 0.08_dp*cos(0.017_dp*real(t*t, dp))
      granger_series(t, 2) = 0.35_dp*granger_series(t - 1, 2) + &
         0.8_dp*granger_series(t - 1, 1) + 0.15_dp*cos(1.13_dp*real(t, dp))
      granger_series(t, 3) = 0.4_dp*granger_series(t - 1, 3) + &
         0.7_dp*sin(0.93_dp*real(t, dp))
   end do
   granger_result = mts_granger_test(granger_series, 1, [2])
   call check(granger_result%info == 0 .and. granger_result%degrees_of_freedom == 2 .and. &
      all(granger_result%predictors == [1, 3]) .and. granger_result%statistic > 0.0_dp .and. &
      granger_result%p_value < 0.05_dp, 'VAR Granger-causality Wald test')
   mq_result = mts_mq(granger_result%unrestricted%residuals, 5)
   call check(mq_result%info == 0 .and. all(mq_result%statistic(2:) >= mq_result%statistic(:4)) .and. &
      all(mq_result%degrees_of_freedom == [9, 18, 27, 36, 45]) .and. &
      all(mq_result%p_value >= 0.0_dp) .and. all(mq_result%p_value <= 1.0_dp), &
      'multivariate Ljung-Box statistics')
   diagnostic_result = mts_diagnostic(granger_result%unrestricted%residuals, 5)
   call check(diagnostic_result%info == 0 .and. &
      maxval(abs(diagnostic_result%mq%statistic - mq_result%statistic)) < 1.0e-11_dp .and. &
      all(shape(diagnostic_result%cross_correlation) == [3, 3, 6]) .and. &
      all(diagnostic_result%p_value >= 0.0_dp) .and. all(diagnostic_result%p_value <= 1.0_dp), &
      'MTS residual model diagnostics')
   backtest_result = mts_var_backtest(series, 1, 80, 3, 2)
   call check(backtest_result%info == 0 .and. &
      all(shape(backtest_result%forecast) == [40, 3, 2]) .and. &
      abs(backtest_result%rmse(2, 1) - &
      sqrt(sum(backtest_result%error(:39, 2, 1)**2)/39.0_dp)) < 1.0e-12_dp .and. &
      abs(backtest_result%mean_absolute_error(3, 2) - &
      sum(abs(backtest_result%error(:38, 3, 2)))/38.0_dp) < 1.0e-12_dp, &
      'rolling-origin VAR backtest losses')
   scm_identification = mts_scm_identify(granger_series, 1, 2, 0, 0.05_dp)
   call check(scm_identification%info == 0 .and. &
      all(shape(scm_identification%zero_count) == [2, 3]) .and. &
      all(scm_identification%zero_count >= 0) .and. &
      all(scm_identification%zero_count <= 6) .and. &
      scm_identification%extra_lags == 0, 'scalar-component first-stage identification')
   do t = 1, 1
      call check(all(scm_identification%diagonal_difference(t, 1:2) == &
         min(3, scm_identification%zero_count(t, 1:2) - &
         scm_identification%zero_count(t - 1, 0:1))), &
         'scalar-component diagonal differences')
   end do
   scm_structure = mts_scm_identify_details(granger_series, 1, 1, 0, 0.05_dp)
   call check((scm_structure%info == 0 .or. scm_structure%info == 2) .and. &
      scm_structure%components >= 0 .and. scm_structure%components <= 3 .and. &
      all(shape(scm_structure%order) == [3, 2]) .and. &
      all(shape(scm_structure%transformation) == [3, 3]), &
      'scalar-component detailed identification')
   scm_specification = mts_scm_specification(reshape([0, 1, 1, 0, 0, 1], [3, 2]), [1, 2, 3])
   call check(scm_specification%info == 0 .and. scm_specification%ar_order == 1 .and. &
      scm_specification%ma_order == 1 .and. &
      all([(scm_specification%transformation(t, t) == 1, t=1, 3)]) .and. &
      all(scm_specification%ar(1, :, 1) == 0) .and. &
      all(scm_specification%ma(3, :, 1) == [0, 2, 2]), 'scalar-component model restrictions')
   scm_fit = mts_scm_fit(granger_series(:, :2), reshape([0, 0, 0, 0], [2, 2]), [1, 2], &
      initial=[0.1_dp, -0.1_dp, 0.0_dp], estimated=[(.false., t=1, 3)])
   call check(scm_fit%info == 0 .and. scm_fit%fit%converged .and. &
      all(shape(scm_fit%fit%model%ar) == [2, 2, 0]) .and. &
      maxval(abs(scm_fit%fit%structural_ar(:, :, 0) - &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))) < 1.0e-12_dp, &
      'conditional scalar-component model estimation')
   forecast = mts_scm_forecast(scm_fit, granger_series(:, :2), 2)
   call check(forecast%info == 0 .and. all(shape(forecast%mean) == [2, 2]), &
      'scalar-component model forecasting')
   scm_refinement = mts_scm_refine(granger_series(:, :2), scm_fit, threshold=0.0_dp, max_steps=1)
   call check(scm_refinement%info == 0 .and. scm_refinement%steps == 0, &
      'scalar-component model refinement')
   kronecker_specification = mts_kronecker_specification([0, 1])
   call check(kronecker_specification%info == 0 .and. kronecker_specification%order == 1 .and. &
      lbound(kronecker_specification%ar, 3) == 0 .and. &
      kronecker_specification%ar(1, 1, 0) == 1 .and. &
      kronecker_specification%ar(1, 2, 0) == 0 .and. &
      kronecker_specification%ar(2, 1, 1) == 0 .and. &
      all(kronecker_specification%ma(1, :, 1) == 0), &
      'Kronecker-index echelon restrictions')
   kronecker_identification = mts_kronecker_identify(granger_series, 2, 0.05_dp)
   call check(kronecker_identification%info == 0 .and. &
      all(kronecker_identification%index >= 0) .and. &
      all(kronecker_identification%index <= 2) .and. &
      all(kronecker_identification%p_value >= 0.0_dp) .and. &
      all(kronecker_identification%p_value <= 1.0_dp), &
      'Kronecker-index canonical-correlation identification')
   kronecker_fit = mts_kronecker_fit(granger_series(:, :2), [1, 1], &
      initial=[0.1_dp, -0.1_dp, 0.3_dp, 0.05_dp, -0.1_dp, 0.2_dp, &
      0.1_dp, 0.0_dp, 0.02_dp, 0.08_dp], estimated=[(.false., t=1, 10)])
   call check(kronecker_fit%info == 0 .and. kronecker_fit%converged .and. &
      all(shape(kronecker_fit%model%ar) == [2, 2, 1]) .and. &
      all(shape(kronecker_fit%residuals) == [180, 2]) .and. &
      maxval(abs(kronecker_fit%model%ar(:, :, 1) - &
      reshape([0.3_dp, 0.05_dp, -0.1_dp, 0.2_dp], [2, 2]))) < 1.0e-12_dp .and. &
      maxval(abs(kronecker_fit%model%ma(:, :, 1) + &
      reshape([0.1_dp, 0.0_dp, 0.02_dp, 0.08_dp], [2, 2]))) < 1.0e-12_dp .and. &
      abs(kronecker_fit%log_likelihood) < huge(1.0_dp), &
      'conditional Kronecker-index estimation')
   forecast = mts_kronecker_forecast(kronecker_fit, granger_series(:, :2), 3)
   call check(forecast%info == 0 .and. all(shape(forecast%mean) == [3, 2]) .and. &
      maxval(abs(forecast%mean(1, :) - kronecker_fit%model%intercept - &
      matmul(kronecker_fit%model%ar(:, :, 1), granger_series(180, :2)) + &
      matmul(kronecker_fit%model%ma(:, :, 1), kronecker_fit%residuals(180, :)))) < 1.0e-12_dp, &
      'Kronecker-index forecasting')
   kronecker_refinement = mts_kronecker_refine(granger_series(:, :2), kronecker_fit, &
      threshold=0.0_dp, max_steps=1, max_iterations=20)
   call check(kronecker_refinement%info == 0 .and. kronecker_refinement%steps == 0 .and. &
      size(kronecker_refinement%active_count) == 1, 'Kronecker-index refinement')
   call set_random_seed(11223)
   varma_simulation = mts_varma_simulate(varma_model, 20, 30)
   call check(varma_simulation%info == 0 .and. all(shape(varma_simulation%series) == [20, 2]), &
      'Gaussian VARMA simulation')

   allocate(vma_model%ar(2, 2, 0), vma_model%ma(2, 2, 1))
   allocate(vma_model%intercept(2), vma_model%sigma(2, 2))
   vma_model%ma(:, :, 1) = reshape([0.25_dp, 0.05_dp, -0.08_dp, 0.2_dp], [2, 2])
   vma_model%intercept = [0.2_dp, -0.1_dp]
   vma_model%sigma = varma_model%sigma
   do t = 1, 80
      innovations(t, :) = [sin(1.17_dp*real(t, dp)) + 0.1_dp*sin(0.021_dp*real(t*t, dp)), &
         0.7_dp*cos(0.79_dp*real(t, dp)) + 0.12_dp*cos(0.017_dp*real(t*t, dp))]
   end do
   varma_simulation = mts_varma_simulate_from_innovations(vma_model, innovations)
   vma_fit = mts_vma_fit(varma_simulation%series, 1, max_iterations=100, tolerance=1.0e-5_dp)
   call check(vma_fit%model%info == 0 .and. abs(vma_fit%log_likelihood) < huge(1.0_dp) .and. &
      all(shape(vma_fit%model%ma) == [2, 2, 1]) .and. all(shape(vma_fit%residuals) == [80, 2]), &
      'conditional VMA estimation')
   call check(vma_fit%aic < huge(1.0_dp) .and. vma_fit%bic < huge(1.0_dp) .and. &
      vma_fit%invertibility_radius >= 0.0_dp, 'VMA inference and diagnostics')
   selected_vma = mts_vmas_fit(varma_simulation%series, [2], max_iterations=80, tolerance=1.0e-5_dp)
   call check(selected_vma%model%info == 0 .and. all(selected_vma%lags == [2]) .and. &
      maxval(abs(selected_vma%model%ma(:, :, 1))) == 0.0_dp, 'selected-lag VMA estimation')
   vma_order = mts_vma_order(varma_simulation%series, 2, max_iterations=80, tolerance=1.0e-5_dp)
   call check(vma_order%info == 0 .and. vma_order%aic_order >= 1 .and. &
      vma_order%aic_order <= 2, 'VMA order selection')
   transfer_filtered = 0.0_dp
   transfer_noise = 0.0_dp
   do t = 1, 100
      transfer_input(t) = sin(0.19_dp*real(t, dp)) + 0.3_dp*cos(0.07_dp*real(t, dp))
      transfer_filtered(t) = transfer_input(t)
      if (t > 1) transfer_filtered(t) = transfer_filtered(t) + 0.35_dp*transfer_filtered(t - 1)
      transfer_innovation = 0.2_dp*sin(1.31_dp*real(t, dp))
      transfer_noise(t) = transfer_innovation
      if (t > 1) transfer_noise(t) = transfer_noise(t) + 0.4_dp*transfer_noise(t - 1)
      transfer_response(t) = 0.15_dp + transfer_noise(t)
      if (t > 1) transfer_response(t) = transfer_response(t) + 0.8_dp*transfer_filtered(t - 1)
      if (t > 2) transfer_response(t) = transfer_response(t) - 0.25_dp*transfer_filtered(t - 2)
   end do
   transfer_fit = mts_transfer_fit(transfer_response, transfer_input, [1, 0, 0], [1, 1, 1], &
      initial=[0.15_dp, 0.8_dp, -0.25_dp, 0.35_dp, 0.4_dp], &
      estimated=[(.false., t=1, 5)])
   call check(transfer_fit%info == 0 .and. transfer_fit%converged .and. &
      size(transfer_fit%residuals) == 98 .and. transfer_fit%delay == 1 .and. &
      maxval(abs(transfer_fit%filtered_input - transfer_filtered)) < 1.0e-11_dp .and. &
      transfer_fit%sigma2 > 0.0_dp, 'rational transfer-function likelihood')
   do t = 1, 100
      transfer_input2(t) = cos(0.11_dp*real(t, dp))
      transfer_trend(t) = real(t, dp)/100.0_dp
      transfer_equilibrium(t) = transfer_input(t) - transfer_input2(t)
      transfer_noise(t) = 0.12_dp*sin(1.07_dp*real(t, dp))
      if (t > 1) transfer_noise(t) = transfer_noise(t) + 0.3_dp*transfer_noise(t - 1)
      if (t > 12) transfer_noise(t) = transfer_noise(t) + 0.2_dp*transfer_noise(t - 12)
      if (t > 13) transfer_noise(t) = transfer_noise(t) - 0.06_dp*transfer_noise(t - 13)
      transfer_response(t) = 0.2_dp + 0.7_dp*transfer_input(t) - &
         0.4_dp*transfer_input2(t) + 0.15_dp*transfer_trend(t) + &
         0.25_dp*transfer_equilibrium(t) + transfer_noise(t)
   end do
   transfer2_fit = mts_transfer2_fit(transfer_response, transfer_input, [1, 0, 0], &
      [1, 0, 0], 12, [0, 0, 0], transfer_input2, [0, 0, 0], transfer_trend, &
      transfer_equilibrium, [0.2_dp, 0.7_dp, 0.15_dp, 0.25_dp, -0.4_dp, 0.3_dp, 0.2_dp], &
      [(.false., t=1, 7)])
   call check(transfer2_fit%info == 0 .and. transfer2_fit%converged .and. &
      transfer2_fit%has_second_input .and. transfer2_fit%has_deterministic .and. &
      transfer2_fit%has_equilibrium .and. size(transfer2_fit%residuals) == 88 .and. &
      transfer2_fit%sigma2 > 0.0_dp, 'two-input seasonal transfer-function likelihood')
   transfer_prediction = mts_transfer2_forecast(transfer2_fit, transfer_response, transfer_input, &
      [0.1_dp, 0.2_dp], [1, 0, 0], [1, 0, 0], [0, 0, 0], transfer_input2, &
      [0.3_dp, 0.4_dp], [0, 0, 0], transfer_trend, [1.01_dp, 1.02_dp], &
      transfer_equilibrium, [-0.2_dp, -0.2_dp])
   call check(transfer_prediction%info == 0 .and. size(transfer_prediction%mean) == 2 .and. &
      all(ieee_is_finite(transfer_prediction%mean)), 'two-input transfer-function forecasting')
   transfer_backtest = mts_transfer2_backtest(transfer_response, transfer_input, 90, &
      [1, 0, 0], [1, 0, 0], 12, [0, 0, 0], &
      [0.2_dp, 0.7_dp, 0.15_dp, 0.25_dp, -0.4_dp, 0.3_dp, 0.2_dp], &
      [(.false., t=1, 7)], transfer_input2, [0, 0, 0], transfer_trend, transfer_equilibrium)
   call check(transfer_backtest%info == 0 .and. size(transfer_backtest%error) == 10 .and. &
      transfer_backtest%rmse >= 0.0_dp .and. transfer_backtest%mean_absolute_error >= 0.0_dp, &
      'rolling transfer-function backtest')
   regts_error(1, :) = [0.1_dp, -0.05_dp]
   do t = 1, 100
      regts_input(t, 1) = sin(0.13_dp*real(t, dp))
      if (t > 1) then
         regts_shock = [0.15_dp*sin(1.17_dp*real(t, dp)), 0.12_dp*cos(0.83_dp*real(t, dp))]
         regts_error(t, :) = matmul(reshape([0.4_dp, -0.05_dp, 0.1_dp, 0.3_dp], [2, 2]), &
            regts_error(t - 1, :)) + regts_shock
      end if
      regts_response(t, 1) = 0.5_dp + 1.2_dp*regts_input(t, 1) + regts_error(t, 1)
      regts_response(t, 2) = -0.2_dp - 0.7_dp*regts_input(t, 1) + regts_error(t, 2)
   end do
   regts_fit = mts_regts_fit(regts_response, regts_input, 1, .true., &
      [0.5_dp, 1.2_dp, 0.4_dp, 0.1_dp, -0.2_dp, -0.7_dp, -0.05_dp, 0.3_dp], &
      reshape([(.false., t=1, 8)], [4, 2]))
   call check(regts_fit%info == 0 .and. regts_fit%converged .and. &
      all(shape(regts_fit%beta) == [2, 2]) .and. all(shape(regts_fit%ar) == [2, 2, 1]) .and. &
      all(shape(regts_fit%residuals) == [99, 2]), 'regression with VAR errors likelihood')
   regts_future(:, 1) = [0.2_dp, -0.1_dp]
   forecast = mts_regts_forecast(regts_fit, regts_response, regts_input, regts_future)
   expected = [0.5_dp + 1.2_dp*regts_future(1, 1), &
      -0.2_dp - 0.7_dp*regts_future(1, 1)] + &
      matmul(regts_fit%ar(:, :, 1), regts_error(100, :))
   call check(forecast%info == 0 .and. all(shape(forecast%mean) == [2, 2]) .and. &
      maxval(abs(forecast%mean(1, :) - expected)) < 1.0e-11_dp, &
      'regression with VAR errors forecasting')
   regts_refinement = mts_regts_refine(regts_response, regts_input, regts_fit, &
      threshold=0.0_dp, max_steps=1)
   call check(regts_refinement%info == 0 .and. regts_refinement%steps == 0, &
      'regression with VAR errors refinement')
   reverse_mq = mts_reverse_mq(granger_result%unrestricted%residuals, 5)
   call check(reverse_mq%info == 0 .and. &
      all(shape(reverse_mq%cross_correlation) == [3, 3, 5]) .and. &
      abs(reverse_mq%statistic(3) - reverse_mq%cumulative(5) + &
      reverse_mq%cumulative(2)) < 1.0e-12_dp .and. &
      all(reverse_mq%degrees_of_freedom == [45, 36, 27, 18, 9]) .and. &
      all(reverse_mq%p_value >= 0.0_dp) .and. all(reverse_mq%p_value <= 1.0_dp), &
      'reversed multivariate Ljung-Box tests')
   extended_ccm = mts_eccm(granger_series, 1, 2, .true., .true.)
   call check(extended_ccm%info == 0 .and. &
      all(shape(extended_ccm%cross_correlation) == [3, 3, 3, 2]) .and. &
      all(shape(extended_ccm%p_value) == [2, 3]) .and. &
      all(extended_ccm%p_value >= 0.0_dp) .and. all(extended_ccm%p_value <= 1.0_dp), &
      'extended cross-correlation order table')
   corner_table = mts_corner(transfer_response, transfer_input, 6, 4)
   call check(corner_table%info == 0 .and. all(shape(corner_table%value) == [6, 4]) .and. &
      maxval(abs(corner_table%value(:, 1) - &
      corner_table%normalized_cross_correlation(:6))) < 1.0e-12_dp .and. &
      abs(corner_table%value(1, 2) - corner_table%value(1, 1)**2) < 1.0e-12_dp .and. &
      all(corner_table%significant .eqv. abs(corner_table%value) > corner_table%threshold), &
      'transfer-function corner table')
   exact_vma = mts_vma_exact_fit(varma_simulation%series(:20, :), 1, .true., &
      [vma_model%intercept, reshape(vma_model%ma(:, :, 1), [4])], &
      [(.false., t=1, 6)])
   call check(exact_vma%info == 0 .and. exact_vma%converged .and. &
      all(shape(exact_vma%model%ma) == [2, 2, 1]) .and. &
      all(shape(exact_vma%residuals) == [20, 2]) .and. &
      ieee_is_finite(exact_vma%log_likelihood), 'exact finite-sample VMA likelihood')
   exact_vma_refinement = mts_vma_exact_refine(varma_simulation%series(:20, :), exact_vma, &
      threshold=0.0_dp, max_steps=1)
   call check(exact_vma_refinement%info == 0 .and. exact_vma_refinement%steps == 0, &
      'exact VMA refinement')
   exact_vma_free = mts_vma_exact_fit(varma_simulation%series(:20, 1:1), 1, .true., &
      max_iterations=40, tolerance=1.0e-5_dp)
   call check((exact_vma_free%info == 0 .or. exact_vma_free%info >= 100) .and. &
      size(exact_vma_free%coefficients) == 2 .and. &
      all(ieee_is_finite(exact_vma_free%coefficients)) .and. &
      ieee_is_finite(exact_vma_free%log_likelihood), 'exact VMA numerical estimation')
   do j = 1, 30
      do t = 1, 12
         apca_data(t, j) = sin(0.17_dp*real(t*j, dp)) + &
            0.3_dp*cos(0.11_dp*real(t + 2*j, dp))
      end do
   end do
   asymptotic_pca = mts_apca(apca_data, 3)
   call check(asymptotic_pca%info == 0 .and. .not. asymptotic_pca%transposed .and. &
      all(shape(asymptotic_pca%factors) == [12, 3]) .and. &
      all(shape(asymptotic_pca%loadings) == [30, 3]) .and. &
      all(asymptotic_pca%standard_deviation >= 0.0_dp), 'asymptotic principal components')
   diffusion_forecast = mts_diffusion_forecast(factor_series(:, 1), factor_series(:, 2:), 120, 2)
   call check(diffusion_forecast%info == 0 .and. &
      all(shape(diffusion_forecast%loadings) == [5, 2]) .and. &
      all(shape(diffusion_forecast%index) == [160, 2]) .and. &
      size(diffusion_forecast%forecast) == 40 .and. diffusion_forecast%mse >= 0.0_dp, &
      'Stock-Watson diffusion-index forecasting')
   multivariate_regression = mts_multivariate_regression(regts_response, regts_input, .true.)
   call check(multivariate_regression%info == 0 .and. &
      all(shape(multivariate_regression%coefficients) == [2, 2]) .and. &
      all(shape(multivariate_regression%covariance) == [4, 4]) .and. &
      maxval(abs(multivariate_regression%residuals - regts_response + &
      spread(multivariate_regression%coefficients(1, :), 1, 100) + &
      matmul(regts_input, multivariate_regression%coefficients(2:2, :)))) < 1.0e-11_dp, &
      'multivariate multiple linear regression')
   var_chi = mts_var_chi(granger_series, 1, .true., 1.0e12_dp)
   call check(var_chi%info == 0 .and. var_chi%degrees_of_freedom == 12 .and. &
      all(shape(var_chi%covariance) == [12, 12]) .and. var_chi%statistic >= 0.0_dp .and. &
      var_chi%p_value >= 0.0_dp .and. var_chi%p_value <= 1.0_dp, &
      'joint chi-square test of weak VAR coefficients')
   forecast = mts_var_fore(fitted, series, 3)
   order_result = mts_var_order(series, 1)
   call check(forecast%info == 0 .and. all(shape(forecast%mean) == [3, 2]) .and. &
      order_result%info == 0, 'VARfore compatibility entry point')
   print '(a)', 'All mts_mod tests passed.'

contains

   subroutine check(ok, name)
      ! Stop the test program when a named assertion fails.
      logical, intent(in) :: ok
      character(len=*), intent(in) :: name
      if (.not. ok) then
         print '(a)', 'FAILED: '//name
         error stop 1
      end if
   end subroutine check
end program test_mts
