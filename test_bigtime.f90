! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the bigtime translation.
program test_bigtime
   use kind_mod, only: dp
   use bigtime_mod
   use random_mod, only: set_random_seed
   implicit none
   type(bigtime_var_fit_t) :: lasso_fit, hlag_fit, zero_fit
   type(bigtime_var_path_t) :: path
   type(bigtime_varx_fit_t) :: varx_lasso, varx_hlag, varx_ridge, varx_zero
   type(bigtime_varx_path_t) :: varx_path
   type(bigtime_varx_grid_t) :: varx_grid
   type(bigtime_varma_fit_t) :: varma_supplied, varma_estimated
   type(bigtime_varma_path_t) :: varma_path
   type(bigtime_cv_result_t) :: var_cv, varx_cv
   type(bigtime_varma_cv_t) :: varma_cv
   type(bigtime_ic_result_t) :: var_ic, varx_ic, varma_ic
   type(bigtime_var_fit_t) :: selected_var
   type(bigtime_varx_fit_t) :: selected_varx
   type(bigtime_forecast_t) :: var_forecast, varx_forecast, varma_forecast
   type(bigtime_path_forecast_t) :: var_path_forecast, varx_path_forecast
   type(bigtime_stability_t) :: stable_result, unstable_result
   type(bigtime_lag_order_t) :: lag_orders
   type(bigtime_coefficient_t) :: dense_coefficients, l1_coefficients
   type(bigtime_coefficient_t) :: hlag_coefficients, random_first, random_second
   type(bigtime_simulation_t) :: supplied_simulation, random_simulation
   type(bigtime_simulation_t) :: repeated_simulation
   real(dp) :: series(90, 2), prox_input(6), expected(6)
   real(dp) :: exogenous(90, 2), varx_series(90, 2)
   real(dp) :: exogenous_forecast(93, 2), expected_forecast(3, 2)
   real(dp) :: varma_series(120, 2), innovations(120, 2)
   real(dp) :: proxy_expected(120, 2)
   real(dp) :: lag_primary(2, 6), lag_secondary(2, 4)
   real(dp) :: coefficient_draws(2, 4), simulation_innovations(4, 1)
   integer :: trailing_zeros(2, 2)
   real(dp), allocatable :: lambdas(:)
   real(dp) :: paired_phi(4), paired_beta(4)
   real(dp) :: prediction(2)
   integer :: time

   prox_input = [3.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   expected = [2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   call check(maxval(abs(bigtime_hlag_prox(prox_input, 1.0_dp, 2, 3) - &
      expected)) < 1.0e-12_dp, 'hierarchical lag proximal operator')
   call check(abs(bigtime_soft_threshold(-2.5_dp, 0.75_dp) + 1.75_dp) < &
      1.0e-12_dp, 'lasso soft threshold')

   series(1, :) = [0.3_dp, -0.2_dp]
   series(2, :) = [0.1_dp, 0.4_dp]
   do time = 3, 90
      series(time, 1) = 0.55_dp*series(time - 1, 1) - &
         0.20_dp*series(time - 1, 2) + &
         0.12_dp*series(time - 2, 1) + &
         0.03_dp*sin(0.7_dp*real(time, dp))
      series(time, 2) = 0.35_dp*series(time - 1, 2) + &
         0.18_dp*series(time - 1, 1) + &
         0.03_dp*cos(0.9_dp*real(time, dp))
   end do

   lasso_fit = bigtime_sparse_var(series, 2, 0.01_dp, &
      bigtime_penalty_l1, 1.0e-9_dp)
   call check(lasso_fit%info == 0 .and. lasso_fit%converged, &
      'lasso sparse VAR convergence')
   call check(all(shape(lasso_fit%phi) == [2, 4]) .and. &
      all(shape(lasso_fit%residuals) == [88, 2]), &
      'sparse VAR dimensions')
   do time = 3, 90
      prediction = lasso_fit%intercept + &
         matmul(lasso_fit%phi(:, 1:2), series(time - 1, :)) + &
         matmul(lasso_fit%phi(:, 3:4), series(time - 2, :))
      call check(maxval(abs(series(time, :) - prediction - &
         lasso_fit%residuals(time - 2, :))) < 1.0e-12_dp, &
         'sparse VAR residual identity')
   end do

   hlag_fit = bigtime_sparse_var(series, 2, 0.02_dp, &
      bigtime_penalty_hlag, 1.0e-9_dp)
   call check(hlag_fit%info == 0 .and. hlag_fit%converged .and. &
      hlag_fit%objective < huge(1.0_dp), &
      'hierarchical sparse VAR convergence')

   lambdas = bigtime_var_lambda_grid(series, 2, bigtime_penalty_l1, &
      100.0_dp, 5)
   call check(size(lambdas) == 5 .and. all(lambdas > 0.0_dp) .and. &
      all(lambdas(1:4) > lambdas(2:5)), 'geometric lambda grid')
   zero_fit = bigtime_sparse_var(series, 2, lambdas(1), &
      bigtime_penalty_l1, 1.0e-9_dp)
   call check(zero_fit%info == 0 .and. zero_fit%nonzero == 0, &
      'lambda upper bound gives zero coefficients')

   path = bigtime_sparse_var_path(series, 2, lambdas, &
      bigtime_penalty_l1, 1.0e-8_dp)
   call check(path%info == 0 .and. all(path%converged) .and. &
      all(shape(path%phi) == [2, 4, 5]), 'warm-started sparse VAR path')
   call check(path%nonzero(1) == 0 .and. &
      path%nonzero(5) >= path%nonzero(1), 'regularization path sparsity')
   var_cv = bigtime_var_cv(series, 2, lambdas(1:3), &
      bigtime_penalty_l1, training_fraction=0.85_dp, tolerance=1.0e-6_dp)
   call check(var_cv%info == 0 .and. &
      all(shape(var_cv%squared_error) == [14, 3]) .and. &
      var_cv%best_index >= 1 .and. var_cv%best_index <= 3 .and. &
      var_cv%one_se_index >= 1 .and. var_cv%one_se_index <= 3, &
      'sparse VAR expanding-window cross-validation')
   call check(maxval(abs(var_cv%mean_squared_error - &
      sum(var_cv%squared_error, dim=1)/14.0_dp)) < 1.0e-14_dp, &
      'VAR cross-validation score summaries')
   var_ic = bigtime_var_path_ic(series, path)
   call check(var_ic%info == 0 .and. all(shape(var_ic%value) == [5, 3]) .and. &
      all(var_ic%selected >= 1) .and. all(var_ic%selected <= 5), &
      'sparse VAR information criteria')
   selected_var = bigtime_select_var_path(series, path, &
      var_ic%selected(bigtime_ic_bic))
   call check(selected_var%info == 0 .and. maxval(abs(selected_var%phi - &
      path%phi(:, :, var_ic%selected(bigtime_ic_bic)))) == 0.0_dp, &
      'materialize selected sparse VAR')
   var_forecast = bigtime_var_forecast(lasso_fit, series, 3)
   expected_forecast(1, :) = lasso_fit%intercept + &
      matmul(lasso_fit%phi(:, 1:2), series(90, :)) + &
      matmul(lasso_fit%phi(:, 3:4), series(89, :))
   expected_forecast(2, :) = lasso_fit%intercept + &
      matmul(lasso_fit%phi(:, 1:2), expected_forecast(1, :)) + &
      matmul(lasso_fit%phi(:, 3:4), series(90, :))
   expected_forecast(3, :) = lasso_fit%intercept + &
      matmul(lasso_fit%phi(:, 1:2), expected_forecast(2, :)) + &
      matmul(lasso_fit%phi(:, 3:4), expected_forecast(1, :))
   call check(var_forecast%info == 0 .and. &
      maxval(abs(var_forecast%mean - expected_forecast)) < 1.0e-14_dp, &
      'recursive sparse VAR forecast')
   var_path_forecast = bigtime_var_path_forecast(path, series, 3)
   selected_var = bigtime_select_var_path(series, path, 3)
   var_forecast = bigtime_var_forecast(selected_var, series, 3)
   call check(var_path_forecast%info == 0 .and. &
      all(shape(var_path_forecast%mean) == [3, 2, 5]) .and. &
      maxval(abs(var_path_forecast%mean(:, :, 3) - var_forecast%mean)) < &
      1.0e-14_dp, 'sparse VAR path forecasts')

   do time = 1, 90
      exogenous(time, 1) = sin(0.21_dp*real(time, dp)) + &
         0.1_dp*cos(0.017_dp*real(time*time, dp))
      exogenous(time, 2) = cos(0.16_dp*real(time, dp)) - &
         0.08_dp*sin(0.013_dp*real(time*time, dp))
   end do
   varx_series(1, :) = [0.2_dp, -0.1_dp]
   varx_series(2, :) = [-0.05_dp, 0.25_dp]
   do time = 3, 90
      varx_series(time, 1) = 0.45_dp*varx_series(time - 1, 1) - &
         0.15_dp*varx_series(time - 1, 2) + &
         0.35_dp*exogenous(time - 1, 1) + &
         0.10_dp*exogenous(time - 2, 2) + &
         0.01_dp*sin(0.81_dp*real(time, dp))
      varx_series(time, 2) = 0.30_dp*varx_series(time - 1, 2) + &
         0.12_dp*varx_series(time - 1, 1) - &
         0.28_dp*exogenous(time - 1, 2) + &
         0.01_dp*cos(0.73_dp*real(time, dp))
   end do

   varx_lasso = bigtime_sparse_varx(varx_series, exogenous, 1, 2, &
      0.002_dp, 0.002_dp, bigtime_penalty_l1, tolerance=1.0e-6_dp, &
      max_iterations=1000)
   call check(varx_lasso%info == 0 .and. varx_lasso%converged .and. &
      all(shape(varx_lasso%phi) == [2, 2]) .and. &
      all(shape(varx_lasso%beta) == [2, 4]) .and. &
      all(shape(varx_lasso%residuals) == [88, 2]), &
      'lasso sparse VARX convergence and dimensions')
   do time = 3, 90
      prediction = varx_lasso%intercept + &
         matmul(varx_lasso%phi(:, 1:2), varx_series(time - 1, :)) + &
         matmul(varx_lasso%beta(:, 1:2), exogenous(time - 1, :)) + &
         matmul(varx_lasso%beta(:, 3:4), exogenous(time - 2, :))
      call check(maxval(abs(varx_series(time, :) - prediction - &
         varx_lasso%residuals(time - 2, :))) < 1.0e-12_dp, &
         'sparse VARX residual identity')
   end do

   varx_hlag = bigtime_sparse_varx(varx_series, exogenous, 2, 2, &
      0.01_dp, 0.01_dp, bigtime_penalty_hlag, tolerance=1.0e-6_dp, &
      max_iterations=1000)
   call check(varx_hlag%info == 0 .and. varx_hlag%converged .and. &
      all(shape(varx_hlag%phi) == [2, 4]) .and. &
      all(shape(varx_hlag%beta) == [2, 4]), &
      'hierarchical sparse VARX convergence')

   varx_ridge = bigtime_sparse_varx(varx_series, exogenous, 1, 2, &
      0.002_dp, 0.002_dp, bigtime_penalty_l1, alpha=0.2_dp, &
      tolerance=1.0e-6_dp, max_iterations=1000)
   call check(varx_ridge%info == 0 .and. &
      sum(varx_ridge%phi**2) + sum(varx_ridge%beta**2) < &
      sum(varx_lasso%phi**2) + sum(varx_lasso%beta**2), &
      'VARX ridge shrinkage')

   varx_grid = bigtime_varx_lambda_grid(varx_series, exogenous, 1, 2, &
      bigtime_penalty_l1, 100.0_dp, 4, 100.0_dp, 4)
   call check(varx_grid%info == 0 .and. &
      all(varx_grid%lambda_phi(1:3) > varx_grid%lambda_phi(2:4)) .and. &
      all(varx_grid%lambda_beta(1:3) > varx_grid%lambda_beta(2:4)), &
      'VARX geometric lambda grids')
   varx_zero = bigtime_sparse_varx(varx_series, exogenous, 1, 2, &
      varx_grid%lambda_phi(1), varx_grid%lambda_beta(1), &
      bigtime_penalty_l1, tolerance=1.0e-6_dp, max_iterations=1000)
   call check(varx_zero%info == 0 .and. varx_zero%nonzero == 0, &
      'VARX lambda upper bounds give zero coefficients')

   paired_phi = varx_grid%lambda_phi
   paired_beta = varx_grid%lambda_beta
   varx_path = bigtime_sparse_varx_path(varx_series, exogenous, 1, 2, &
      paired_phi, paired_beta, bigtime_penalty_l1, tolerance=1.0e-8_dp, &
      max_iterations=1000)
   call check(varx_path%info == 0 .and. all(varx_path%converged) .and. &
      all(shape(varx_path%phi) == [2, 2, 4]) .and. &
      all(shape(varx_path%beta) == [2, 4, 4]), &
      'warm-started sparse VARX path')
   call check(varx_path%nonzero(1) == 0 .and. &
      varx_path%nonzero(4) >= varx_path%nonzero(1), &
      'VARX regularization path sparsity')
   varx_cv = bigtime_varx_cv(varx_series, exogenous, 1, 2, paired_phi, &
      paired_beta, bigtime_penalty_l1, training_fraction=0.85_dp, &
      tolerance=1.0e-6_dp, max_iterations=1000)
   call check(varx_cv%info == 0 .and. &
      all(shape(varx_cv%squared_error) == [14, 4]) .and. &
      varx_cv%best_index >= 1 .and. varx_cv%best_index <= 4 .and. &
      paired_phi(varx_cv%one_se_index) >= paired_phi(varx_cv%best_index) .and. &
      paired_beta(varx_cv%one_se_index) >= paired_beta(varx_cv%best_index), &
      'sparse VARX expanding-window and one-SE selection')
   varx_ic = bigtime_varx_path_ic(varx_series, exogenous, varx_path)
   call check(varx_ic%info == 0 .and. &
      all(shape(varx_ic%value) == [4, 3]) .and. &
      all(varx_ic%selected >= 1) .and. all(varx_ic%selected <= 4), &
      'sparse VARX information criteria')
   selected_varx = bigtime_select_varx_path(varx_series, exogenous, varx_path, &
      varx_ic%selected(bigtime_ic_aic))
   call check(selected_varx%info == 0 .and. &
      maxval(abs(selected_varx%beta - varx_path%beta(:, :, &
      varx_ic%selected(bigtime_ic_aic)))) == 0.0_dp, &
      'materialize selected sparse VARX')
   exogenous_forecast(:90, :) = exogenous
   do time = 91, 93
      exogenous_forecast(time, 1) = sin(0.21_dp*real(time, dp)) + &
         0.1_dp*cos(0.017_dp*real(time*time, dp))
      exogenous_forecast(time, 2) = cos(0.16_dp*real(time, dp)) - &
         0.08_dp*sin(0.013_dp*real(time*time, dp))
   end do
   varx_forecast = bigtime_varx_forecast(varx_lasso, varx_series, &
      exogenous_forecast, 3)
   expected_forecast(1, :) = varx_lasso%intercept + &
      matmul(varx_lasso%phi, varx_series(90, :)) + &
      matmul(varx_lasso%beta(:, 1:2), exogenous_forecast(90, :)) + &
      matmul(varx_lasso%beta(:, 3:4), exogenous_forecast(89, :))
   expected_forecast(2, :) = varx_lasso%intercept + &
      matmul(varx_lasso%phi, expected_forecast(1, :)) + &
      matmul(varx_lasso%beta(:, 1:2), exogenous_forecast(91, :)) + &
      matmul(varx_lasso%beta(:, 3:4), exogenous_forecast(90, :))
   expected_forecast(3, :) = varx_lasso%intercept + &
      matmul(varx_lasso%phi, expected_forecast(2, :)) + &
      matmul(varx_lasso%beta(:, 1:2), exogenous_forecast(92, :)) + &
      matmul(varx_lasso%beta(:, 3:4), exogenous_forecast(91, :))
   call check(varx_forecast%info == 0 .and. &
      maxval(abs(varx_forecast%mean - expected_forecast)) < 1.0e-14_dp, &
      'recursive sparse VARX forecast')
   varx_path_forecast = bigtime_varx_path_forecast(varx_path, varx_series, &
      exogenous_forecast, 3)
   selected_varx = bigtime_select_varx_path(varx_series, exogenous, &
      varx_path, 2)
   varx_forecast = bigtime_varx_forecast(selected_varx, varx_series, &
      exogenous_forecast, 3)
   call check(varx_path_forecast%info == 0 .and. &
      all(shape(varx_path_forecast%mean) == [3, 2, 4]) .and. &
      maxval(abs(varx_path_forecast%mean(:, :, 2) - varx_forecast%mean)) < &
      1.0e-14_dp, 'sparse VARX path forecasts')

   innovations = 0.0_dp
   do time = 1, 120
      innovations(time, 1) = 0.04_dp*sin(0.83_dp*real(time, dp)) + &
         0.015_dp*cos(0.031_dp*real(time*time, dp))
      innovations(time, 2) = 0.035_dp*cos(0.67_dp*real(time, dp)) - &
         0.012_dp*sin(0.027_dp*real(time*time, dp))
   end do
   varma_series(1, :) = innovations(1, :)
   do time = 2, 120
      varma_series(time, 1) = 0.50_dp*varma_series(time - 1, 1) - &
         0.12_dp*varma_series(time - 1, 2) + &
         0.30_dp*innovations(time - 1, 1) + innovations(time, 1)
      varma_series(time, 2) = 0.35_dp*varma_series(time - 1, 2) + &
         0.10_dp*varma_series(time - 1, 1) - &
         0.25_dp*innovations(time - 1, 2) + innovations(time, 2)
   end do

   varma_supplied = bigtime_sparse_varma(varma_series, 3, 0.01_dp, &
      bigtime_penalty_l1, 1, 1, 0.001_dp, 0.001_dp, &
      bigtime_penalty_l1, innovations=innovations, tolerance=1.0e-6_dp, &
      phase2_max_iterations=1000)
   call check(varma_supplied%info == 0 .and. varma_supplied%converged .and. &
      .not. varma_supplied%estimated_innovations .and. &
      all(shape(varma_supplied%phi) == [2, 2]) .and. &
      all(shape(varma_supplied%theta) == [2, 2]) .and. &
      maxval(abs(varma_supplied%innovations - innovations)) == 0.0_dp, &
      'sparse VARMA with supplied innovations')
   do time = 2, 120
      prediction = varma_supplied%intercept + &
         matmul(varma_supplied%phi, varma_series(time - 1, :)) + &
         matmul(varma_supplied%theta, innovations(time - 1, :))
      call check(maxval(abs(varma_series(time, :) - prediction - &
         varma_supplied%residuals(time - 1, :))) < 1.0e-12_dp, &
         'sparse VARMA residual identity')
   end do

   varma_estimated = bigtime_sparse_varma(varma_series, 3, 0.002_dp, &
      bigtime_penalty_l1, 1, 1, 0.002_dp, 0.002_dp, &
      bigtime_penalty_l1, tolerance=1.0e-6_dp, &
      phase1_max_iterations=2000, phase2_max_iterations=1000)
   call check(varma_estimated%info == 0 .and. &
      varma_estimated%estimated_innovations .and. &
      varma_estimated%phase1%converged .and. varma_estimated%phase2%converged, &
      'two-stage sparse VARMA convergence')
   proxy_expected(4:, :) = varma_estimated%phase1%residuals
   proxy_expected(:3, :) = varma_estimated%phase1%residuals(:3, :)
   proxy_expected = proxy_expected - spread( &
      sum(proxy_expected, dim=1)/120.0_dp, 1, 120)
   call check(maxval(abs(varma_estimated%innovations - proxy_expected)) < &
      1.0e-14_dp .and. maxval(abs(sum(varma_estimated%innovations, dim=1))) < &
      1.0e-14_dp, 'Phase I innovation proxy extension and centering')

   paired_phi = [0.05_dp, 0.02_dp, 0.01_dp, 0.005_dp]
   paired_beta = [0.05_dp, 0.02_dp, 0.01_dp, 0.005_dp]
   varma_path = bigtime_sparse_varma_path(varma_series, 3, 0.002_dp, &
      bigtime_penalty_l1, 1, 1, paired_phi, paired_beta, &
      bigtime_penalty_l1, tolerance=1.0e-6_dp, &
      phase1_max_iterations=2000, phase2_max_iterations=1000)
   call check(varma_path%info == 0 .and. varma_path%estimated_innovations .and. &
      varma_path%phase1%converged .and. all(varma_path%phase2%converged) .and. &
      all(shape(varma_path%phase2%phi) == [2, 2, 4]) .and. &
      all(shape(varma_path%phase2%beta) == [2, 2, 4]), &
      'two-stage sparse VARMA regularization path')
   call check(maxval(abs(varma_path%innovations - &
      varma_estimated%innovations)) < 1.0e-14_dp, &
      'VARMA path shares one Phase I innovation proxy')
   varma_cv = bigtime_varma_cv(varma_series, 3, 0.002_dp, &
      bigtime_penalty_l1, 1, 1, paired_phi, paired_beta, &
      bigtime_penalty_l1, training_fraction=0.9_dp, tolerance=1.0e-6_dp, &
      phase1_max_iterations=2000, phase2_max_iterations=1000)
   call check(varma_cv%info == 0 .and. varma_cv%estimated_innovations .and. &
      varma_cv%phase1%converged .and. &
      all(shape(varma_cv%phase2%squared_error) == [12, 4]) .and. &
      maxval(abs(varma_cv%innovations - varma_path%innovations)) < 1.0e-14_dp, &
      'two-stage sparse VARMA cross-validation')
   varma_ic = bigtime_varma_path_ic(varma_series, varma_path)
   call check(varma_ic%info == 0 .and. &
      all(shape(varma_ic%value) == [4, 3]) .and. &
      all(varma_ic%selected >= 1) .and. all(varma_ic%selected <= 4), &
      'sparse VARMA information criteria')
   varma_forecast = bigtime_varma_forecast(varma_supplied, varma_series, 3)
   expected_forecast(1, :) = varma_supplied%intercept + &
      matmul(varma_supplied%phi, varma_series(120, :)) + &
      matmul(varma_supplied%theta, innovations(120, :))
   expected_forecast(2, :) = varma_supplied%intercept + &
      matmul(varma_supplied%phi, expected_forecast(1, :))
   expected_forecast(3, :) = varma_supplied%intercept + &
      matmul(varma_supplied%phi, expected_forecast(2, :))
   call check(varma_forecast%info == 0 .and. &
      maxval(abs(varma_forecast%mean - expected_forecast)) < 1.0e-14_dp, &
      'recursive sparse VARMA forecast')

   stable_result = bigtime_var_stability(reshape([0.8_dp], [1, 1]), 1)
   unstable_result = bigtime_var_stability(reshape([1.05_dp], [1, 1]), 1)
   call check(stable_result%info == 0 .and. stable_result%stable .and. &
      abs(stable_result%maximum_modulus - 0.8_dp) < 1.0e-12_dp .and. &
      unstable_result%info == 0 .and. .not. unstable_result%stable .and. &
      abs(unstable_result%maximum_modulus - 1.05_dp) < 1.0e-12_dp, &
      'sparse VAR companion stability')

   lag_primary = reshape([ &
      0.5_dp, 0.0_dp, 0.0_dp, 0.3_dp, 0.2_dp, 0.0_dp, &
      0.0_dp, 0.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.1_dp], [2, 6])
   lag_secondary = reshape([ &
      0.1_dp, 0.0_dp, 0.0_dp, 0.2_dp, &
      0.0_dp, 0.3_dp, 0.4_dp, 0.0_dp], [2, 4])
   lag_orders = bigtime_active_lags(lag_primary, 2, 3, lag_secondary, 2, 2)
   call check(lag_orders%info == 0 .and. &
      all(lag_orders%primary == reshape([2, 0, 0, 3], [2, 2])) .and. &
      all(lag_orders%secondary == reshape([1, 2, 2, 1], [2, 2])), &
      'active lag-order matrices')

   coefficient_draws = reshape([ &
      0.4_dp, 0.2_dp, -0.3_dp, 0.1_dp, &
      0.5_dp, -0.2_dp, 0.25_dp, 0.15_dp], [2, 4])
   dense_coefficients = bigtime_var_coefficients_from_draws( &
      coefficient_draws, 2, 0.99_dp, bigtime_sparsity_dense, 0.5_dp)
   call check(dense_coefficients%info == 0 .and. &
      maxval(abs(dense_coefficients%phi(:, :2) - &
      coefficient_draws(:, :2))) < 1.0e-14_dp .and. &
      maxval(abs(dense_coefficients%phi(:, 3:4) - &
      0.5_dp*coefficient_draws(:, 3:4))) < 1.0e-14_dp .and. &
      maxval(abs(dense_coefficients%companion(3:4, :2) - &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))) < 1.0e-14_dp, &
      'dense coefficient decay and companion layout')
   l1_coefficients = bigtime_var_coefficients_from_draws( &
      coefficient_draws, 2, 0.99_dp, bigtime_sparsity_l1, 1.0_dp, &
      zero_indices=[1, 6])
   call check(l1_coefficients%info == 0 .and. &
      l1_coefficients%phi(1, 1) == 0.0_dp .and. &
      l1_coefficients%phi(2, 3) == 0.0_dp, &
      'supplied L1 coefficient mask')
   trailing_zeros = reshape([0, 1, 2, 0], [2, 2])
   hlag_coefficients = bigtime_var_coefficients_from_draws( &
      coefficient_draws, 2, 0.99_dp, bigtime_sparsity_hlag, 1.0_dp, &
      trailing_zeros=trailing_zeros)
   call check(hlag_coefficients%info == 0 .and. &
      hlag_coefficients%phi(2, 1) /= 0.0_dp .and. &
      hlag_coefficients%phi(2, 3) == 0.0_dp .and. &
      hlag_coefficients%phi(1, 2) == 0.0_dp .and. &
      hlag_coefficients%phi(1, 4) == 0.0_dp, &
      'hierarchical trailing-lag mask')
   dense_coefficients = bigtime_var_coefficients_from_draws( &
      reshape([1.5_dp], [1, 1]), 1, 0.8_dp, &
      bigtime_sparsity_dense, 1.0_dp)
   call check(dense_coefficients%info == 0 .and. &
      dense_coefficients%scaling_iterations > 0 .and. &
      dense_coefficients%maximum_modulus <= 0.8_dp, &
      'coefficient stability contraction')

   simulation_innovations(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   supplied_simulation = bigtime_var_simulate_from_innovations( &
      reshape([0.5_dp], [1, 1]), [0.0_dp], simulation_innovations, &
      [0.0_dp], 1)
   call check(supplied_simulation%info == 0 .and. &
      all(shape(supplied_simulation%series) == [3, 1]) .and. &
      maxval(abs(supplied_simulation%series(:, 1) - &
      [0.5_dp, 0.25_dp, 0.125_dp])) < 1.0e-14_dp, &
      'supplied-innovation VAR simulation and burn-in')

   call set_random_seed(8137)
   random_first = bigtime_random_var_coefficients(2, 3, 0.8_dp, &
      bigtime_sparsity_hlag, 0.5_dp, zero_min=0, zero_max=2)
   call set_random_seed(8137)
   random_second = bigtime_random_var_coefficients(2, 3, 0.8_dp, &
      bigtime_sparsity_hlag, 0.5_dp, zero_min=0, zero_max=2)
   call check(random_first%info == 0 .and. random_second%info == 0, &
      'shared-RNG sparse coefficient generation')
   call check(maxval(abs(random_first%phi - random_second%phi)) == 0.0_dp, &
      'shared-RNG sparse coefficient reproducibility')
   call check(random_first%maximum_modulus <= 0.8_dp, &
      'random sparse coefficient stability')
   call set_random_seed(9173)
   random_simulation = bigtime_var_simulate(reshape([0.4_dp], [1, 1]), &
      [0.1_dp], 8, 4, innovation_scale=0.2_dp)
   call set_random_seed(9173)
   repeated_simulation = bigtime_var_simulate(reshape([0.4_dp], [1, 1]), &
      [0.1_dp], 8, 4, innovation_scale=0.2_dp)
   call check(random_simulation%info == 0 .and. &
      all(shape(random_simulation%series) == [8, 1]) .and. &
      maxval(abs(random_simulation%series - repeated_simulation%series)) == &
      0.0_dp, 'shared-RNG Gaussian VAR simulation')

   print '(a)', 'bigtime tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition is false.
      logical, intent(in) :: condition !! Flag controlling condition.
      character(len=*), intent(in) :: message !! Message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_bigtime
