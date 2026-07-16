! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the BigVAR translation.
program test_bigvar
   use kind_mod, only: dp
   use bigvar_mod
   use bigtime_mod, only: bigtime_forecast_t, bigtime_stability_t
   use bigtime_mod, only: bigtime_soft_threshold
   use time_series_random_mod, only: set_random_seed
   implicit none
   type(bigvar_fit_t) :: lag_fit, own_other_fit, sparse_lag_fit
   type(bigvar_fit_t) :: sparse_own_other_fit, zero_fit
   type(bigvar_fit_t) :: hlag_component_fit, hlag_own_other_fit
   type(bigvar_fit_t) :: hlag_element_fit
   type(bigvar_fit_t) :: basic_en_fit, tapered_fit
   type(bigvar_fit_t) :: basic_fit, separate_fit, relaxed_fit, direct_fit
   type(bigvar_fit_t) :: mcp_fit, scad_fit
   type(bigvar_fit_t) :: bgr_zero_fit, bgr_random_walk_fit
   type(bigvar_path_t) :: path
   type(bigvar_varx_fit_t) :: varx_lag, varx_own_other
   type(bigvar_varx_fit_t) :: varx_sparse_lag, varx_sparse_own_other
   type(bigvar_varx_fit_t) :: varx_basic_en, varx_mcp, varx_scad
   type(bigvar_varx_fit_t) :: varx_efx, varx_zero
   type(bigvar_varx_fit_t) :: relaxed_varx, contemporaneous_varx
   type(bigvar_varx_fit_t) :: current_only_varx, rejected_varx
   type(bigvar_varx_path_t) :: varx_path
   type(bigvar_dual_path_t) :: dual_path
   type(bigvar_varx_dual_path_t) :: varx_dual_path
   type(bigvar_separate_path_t) :: separate_path
   type(bigvar_validation_t) :: validation, evaluation, varx_validation
   type(bigvar_separate_validation_t) :: separate_validation
   type(bigvar_reselection_t) :: reselection
   type(bigvar_interval_forecast_t) :: interval, interval_90
   type(bigvar_ls_varx_fit_t) :: ls_varx, direct_ls_varx
   type(bigvar_ic_search_t) :: ic_search
   type(bigvar_ic_evaluation_t) :: ic_evaluation
   type(bigvar_benchmark_t) :: mean_benchmark, random_walk_benchmark
   type(bigvar_simulation_t) :: supplied_simulation
   type(bigvar_simulation_t) :: random_simulation, repeated_simulation
   type(bigtime_forecast_t) :: forecast
   type(bigtime_stability_t) :: stability
   real(dp) :: coefficients(2, 4), expected(2, 4), series(100, 2)
   real(dp) :: exogenous(100, 2), varx_series(100, 2)
   real(dp) :: transfer_series(100, 2)
   real(dp) :: contemporaneous_series(100, 2)
   real(dp) :: exogenous_forecast(103, 2)
   real(dp) :: sparse_norm, sparse_shrink
   real(dp) :: efx_input(2, 5), efx_expected(2, 5)
   real(dp) :: prediction(2)
   real(dp), allocatable :: direct_prediction(:)
   real(dp), allocatable :: innovation_covariance(:, :)
   real(dp), allocatable :: forecast_covariance(:, :, :)
   real(dp), allocatable :: companion(:, :)
   real(dp) :: simulation_innovations(4, 1)
   real(dp) :: manual_covariance(2, 2)
   real(dp) :: validation_lambdas(2), validation_alphas(2)
   real(dp) :: minnesota_target(2)
   real(dp), allocatable :: lambdas(:)
   real(dp), allocatable :: dual_lambdas(:, :)
   real(dp), allocatable :: linear_lambdas(:)
   real(dp), allocatable :: linear_lambda_matrix(:, :)
   real(dp), allocatable :: separate_lambdas(:, :)
   real(dp), allocatable :: bgr_grid(:)
   logical :: random_walk_prior(2)
   integer :: time

   call check(abs(bigvar_forecast_loss([1.0_dp, -3.0_dp], &
      bigvar_loss_l1) - 4.0_dp) < 1.0e-14_dp, 'L1 forecast loss')
   call check(abs(bigvar_forecast_loss([1.0_dp, -3.0_dp], &
      bigvar_loss_l2) - 10.0_dp) < 1.0e-14_dp, 'L2 forecast loss')
   call check(abs(bigvar_forecast_loss([1.0_dp, -3.0_dp], &
      bigvar_loss_huber, 2.0_dp) - 4.5_dp) < 1.0e-14_dp, &
      'Huber forecast loss')

   companion = bigvar_var_to_companion(reshape([0.5_dp, 0.2_dp], [1, 2]), 2)
   call check(all(shape(companion) == [2, 2]) .and. &
      maxval(abs(companion - reshape([0.5_dp, 1.0_dp, 0.2_dp, 0.0_dp], &
      [2, 2]))) < 1.0e-14_dp, 'BigVAR multiple companion form')
   simulation_innovations(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   supplied_simulation = bigvar_var_simulate_from_innovations( &
      reshape([0.5_dp, 0.2_dp], [1, 2]), simulation_innovations, &
      burnin=1)
   call check(supplied_simulation%info == 0 .and. &
      all(shape(supplied_simulation%series) == [3, 1]) .and. &
      maxval(abs(supplied_simulation%series(:, 1) - &
      [0.5_dp, 0.45_dp, 0.325_dp])) < 1.0e-14_dp .and. &
      maxval(abs(supplied_simulation%companion - companion)) < 1.0e-14_dp, &
      'supplied-innovation BigVAR simulation')
   call set_random_seed(7319)
   random_simulation = bigvar_var_simulate(reshape([0.4_dp, 0.0_dp, &
      0.0_dp, 0.3_dp], [2, 2]), reshape([0.2_dp, 0.08_dp, 0.08_dp, &
      0.3_dp], [2, 2]), 6)
   call set_random_seed(7319)
   repeated_simulation = bigvar_var_simulate(reshape([0.4_dp, 0.0_dp, &
      0.0_dp, 0.3_dp], [2, 2]), reshape([0.2_dp, 0.08_dp, 0.08_dp, &
      0.3_dp], [2, 2]), 6)
   call check(random_simulation%info == 0 .and. &
      all(shape(random_simulation%series) == [6, 2]) .and. &
      all(shape(random_simulation%innovations) == [506, 2]) .and. &
      random_simulation%burnin == 500 .and. &
      maxval(abs(random_simulation%series - &
      repeated_simulation%series)) == 0.0_dp .and. &
      maxval(abs(random_simulation%innovation_covariance - &
      reshape([0.2_dp, 0.08_dp, 0.08_dp, 0.3_dp], [2, 2]))) < 1.0e-14_dp, &
      'correlated Gaussian BigVAR simulation')
   random_simulation = bigvar_var_simulate(reshape([1.1_dp], [1, 1]), &
      reshape([1.0_dp], [1, 1]), 5, burnin=0)
   call check(random_simulation%info == 2, &
      'unstable BigVAR simulation rejection')

   coefficients = reshape([3.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 4])
   expected = coefficients
   expected(:, 1:2) = (1.0_dp - sqrt(2.0_dp)/5.0_dp)*expected(:, 1:2)
   call check(maxval(abs(bigvar_group_prox(coefficients, 1.0_dp, &
      bigvar_structure_lag, 2, 2) - expected)) < 1.0e-14_dp, &
      'Lag group proximal operator')

   coefficients = reshape([3.0_dp, 4.0_dp, 4.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 4])
   expected = 0.0_dp
   expected(1, 1) = 3.0_dp - sqrt(2.0_dp)
   expected(2, 2) = 0.0_dp
   expected(2, 1) = 3.0_dp
   expected(1, 2) = 3.0_dp
   call check(maxval(abs(bigvar_group_prox(coefficients, 1.0_dp, &
      bigvar_structure_own_other, 2, 2) - expected)) < 1.0e-14_dp, &
      'Own/Other group proximal operator')

   coefficients = reshape([3.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 4])
   expected = 0.0_dp
   sparse_norm = sqrt(2.75_dp**2 + 3.75_dp**2)
   sparse_shrink = 1.0_dp - 1.5_dp/sparse_norm
   expected(1, 1) = sparse_shrink*2.75_dp
   expected(2, 2) = sparse_shrink*3.75_dp
   call check(maxval(abs(bigvar_group_prox(coefficients, 1.0_dp, &
      bigvar_structure_sparse_lag, 2, 2, 0.25_dp) - expected)) < &
      1.0e-14_dp, 'SparseLag proximal operator')

   coefficients = 0.0_dp
   coefficients(1, 1:2) = [3.0_dp, 4.0_dp]
   expected = 0.0_dp
   expected(1, 1:2) = (1.0_dp - sqrt(2.0_dp)/5.0_dp)*[3.0_dp, 4.0_dp]
   call check(maxval(abs(bigvar_group_prox(coefficients, 1.0_dp, &
      bigvar_structure_hlag_component, 2, 2) - expected)) < 1.0e-14_dp, &
      'componentwise HLAG proximal operator')

   coefficients = 0.0_dp
   coefficients(1, 1:2) = [3.0_dp, 4.0_dp]
   expected = 0.0_dp
   expected(1, 1:2) = (1.0_dp - 1.0_dp/sqrt(18.0_dp))*[3.0_dp, 3.0_dp]
   call check(maxval(abs(bigvar_group_prox(coefficients, 1.0_dp, &
      bigvar_structure_hlag_own_other, 2, 2) - expected)) < 1.0e-14_dp, &
      'Own/Other HLAG proximal operator')

   coefficients = 0.0_dp
   coefficients(1, 1) = 3.0_dp
   coefficients(1, 3) = 4.0_dp
   expected = 0.0_dp
   expected(1, 1) = 3.0_dp*(1.0_dp - 1.0_dp/sqrt(18.0_dp))
   expected(1, 3) = expected(1, 1)
   call check(maxval(abs(bigvar_group_prox(coefficients, 1.0_dp, &
      bigvar_structure_hlag_element, 2, 2) - expected)) < 1.0e-14_dp, &
      'elementwise HLAG proximal operator')

   coefficients = reshape([3.0_dp, -2.0_dp, 0.5_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 4])
   expected = bigtime_soft_threshold(coefficients, 1.0_dp)
   call check(maxval(abs(bigvar_group_prox(coefficients, 1.0_dp, &
      bigvar_structure_basic, 2, 2) - expected)) < 1.0e-14_dp, &
      'Basic elementwise lasso proximal operator')
   expected = bigtime_soft_threshold(coefficients, 0.25_dp)/1.75_dp
   call check(maxval(abs(bigvar_group_prox(coefficients, 1.0_dp, &
      bigvar_structure_basic_en, 2, 2, 0.25_dp) - expected)) < &
      1.0e-14_dp, 'BasicEN proximal operator')

   coefficients = 2.0_dp
   expected(:, 1:2) = 1.0_dp
   expected(:, 3:4) = 0.0_dp
   call check(maxval(abs(bigvar_group_prox(coefficients, 1.0_dp, &
      bigvar_structure_tapered, 2, 2, 1.0_dp) - expected)) < &
      1.0e-14_dp, 'Tapered lag-weighted proximal operator')

   call check(bigvar_mcp_update(0.5_dp, 1.0_dp, 3.0_dp, 1.0_dp) == &
      0.0_dp .and. abs(bigvar_mcp_update(2.0_dp, 1.0_dp, 3.0_dp, &
      1.0_dp) - 0.6_dp) < 1.0e-14_dp .and. &
      abs(bigvar_mcp_update(7.0_dp, 1.0_dp, 3.0_dp, 1.0_dp) - &
      3.5_dp) < 1.0e-14_dp, 'MCP coordinate update regions')
   call check(bigvar_scad_update(0.5_dp, 1.0_dp, 3.0_dp, 1.0_dp) == &
      0.0_dp .and. abs(bigvar_scad_update(1.5_dp, 1.0_dp, 3.0_dp, &
      1.0_dp) - 0.25_dp) < 1.0e-14_dp .and. &
      abs(bigvar_scad_update(4.0_dp, 1.0_dp, 3.0_dp, 1.0_dp) - &
      5.0_dp/3.0_dp) < 1.0e-14_dp, 'SCAD coordinate update regions')

   efx_input = 0.0_dp
   efx_input(1, 1) = 3.0_dp
   efx_input(1, 5) = 4.0_dp
   efx_expected = 0.0_dp
   efx_expected(1, 1) = 3.0_dp*(1.0_dp - 1.0_dp/sqrt(18.0_dp))
   efx_expected(1, 5) = 3.0_dp*(1.0_dp - 1.0_dp/sqrt(18.0_dp))
   call check(maxval(abs(bigvar_efx_prox(efx_input, 1.0_dp, 2, 2, 1, 1) - &
      efx_expected)) < 1.0e-14_dp, 'EFX nested proximal operator')

   series(1, :) = [0.3_dp, -0.2_dp]
   series(2, :) = [0.1_dp, 0.4_dp]
   do time = 3, 100
      series(time, 1) = 0.60_dp*series(time - 1, 1) - &
         0.18_dp*series(time - 1, 2) + &
         0.10_dp*series(time - 2, 1) + &
         0.02_dp*sin(0.73_dp*real(time, dp))
      series(time, 2) = 0.42_dp*series(time - 1, 2) + &
         0.16_dp*series(time - 1, 1) + &
         0.02_dp*cos(0.61_dp*real(time, dp))
   end do
   minnesota_target = [1.0_dp, 0.8_dp]

   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_basic, &
      100.0_dp, 3, include_intercept=.false.)
   zero_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_basic, 1.0e-8_dp, include_intercept=.false.)
   call check(zero_fit%info == 0 .and. .not. zero_fit%include_intercept .and. &
      zero_fit%nonzero == 0 .and. &
      maxval(abs(zero_fit%intercept)) < 1.0e-14_dp, &
      'structured VAR without intercept')
   path = bigvar_structured_var_path(series, 2, lambdas, &
      bigvar_structure_basic, 1.0e-8_dp, include_intercept=.false.)
   call check(path%info == 0 .and. .not. path%include_intercept .and. &
      maxval(abs(path%intercept)) < 1.0e-14_dp, &
      'structured VAR path without intercept')
   relaxed_fit = bigvar_structured_var(series, 2, lambdas(2), &
      bigvar_structure_basic, 1.0e-8_dp, refit_fraction=0.5_dp, &
      include_intercept=.false.)
   call check(relaxed_fit%info == 0 .and. relaxed_fit%relaxed .and. &
      .not. relaxed_fit%include_intercept .and. &
      maxval(abs(relaxed_fit%intercept)) < 1.0e-14_dp, &
      'relaxed VAR refit without intercept')
   separate_lambdas = bigvar_separate_lambda_grid(series, 2, &
      bigvar_structure_basic, 100.0_dp, 3, include_intercept=.false.)
   separate_path = bigvar_structured_var_separate_path(series, 2, &
      separate_lambdas, bigvar_structure_basic, 1.0e-8_dp, &
      include_intercept=.false.)
   call check(separate_path%info == 0 .and. &
      .not. separate_path%include_intercept .and. &
      maxval(abs(separate_path%intercept)) < 1.0e-14_dp, &
      'response-specific VAR path without intercept')
   validation = bigvar_var_validate(series, 2, lambdas(1:2), &
      bigvar_structure_basic, 94, 96, 1, tolerance=1.0e-7_dp, &
      max_iterations=5000, include_intercept=.false.)
   call check((validation%info == 0 .or. validation%info == 4) .and. &
      all(validation%valid), 'rolling VAR validation without intercept')
   validation = bigvar_var_validate_loo(series, 2, lambdas(1:2), &
      bigvar_structure_basic, tolerance=1.0e-7_dp, max_iterations=5000, &
      first_observation=97, last_observation=98, &
      include_intercept=.false.)
   call check((validation%info == 0 .or. validation%info == 4) .and. &
      all(validation%valid), 'leave-one-out VAR validation without intercept')
   validation = bigvar_var_validate(series, 2, [0.1_dp], &
      bigvar_structure_bgr, 94, 96, 1, include_intercept=.false.)
   call check(validation%info == 1, &
      'BGR rejects intercept suppression')

   mean_benchmark = bigvar_mean_benchmark(series, 94, 96, 2)
   call check(mean_benchmark%info == 0 .and. &
      all(shape(mean_benchmark%forecasts) == [3, 2]) .and. &
      all(mean_benchmark%valid) .and. &
      maxval(abs(mean_benchmark%forecasts(1, :) - &
      sum(series(:94, :), dim=1)/94.0_dp)) < 1.0e-14_dp, &
      'expanding unconditional-mean benchmark')
   call check(abs(mean_benchmark%loss(1) - sum((series(96, :) - &
      mean_benchmark%forecasts(1, :))**2)) < 1.0e-14_dp .and. &
      abs(mean_benchmark%mean_loss - sum(mean_benchmark%loss)/3.0_dp) < &
      1.0e-14_dp, 'mean benchmark loss and target alignment')
   mean_benchmark = bigvar_mean_benchmark(series, 94, 96, 1, &
      bigvar_loss_huber, 0.05_dp, 20)
   call check(mean_benchmark%info == 0 .and. &
      mean_benchmark%window_size == 20 .and. &
      maxval(abs(mean_benchmark%forecasts(1, :) - &
      sum(series(75:94, :), dim=1)/20.0_dp)) < 1.0e-14_dp .and. &
      abs(mean_benchmark%loss(1) - bigvar_forecast_loss(series(95, :) - &
      mean_benchmark%forecasts(1, :), bigvar_loss_huber, 0.05_dp)) < &
      1.0e-14_dp, 'fixed-window Huber mean benchmark')
   random_walk_benchmark = bigvar_random_walk_benchmark(series, 94, 96, 3, &
      bigvar_loss_l1)
   call check(random_walk_benchmark%info == 0 .and. &
      maxval(abs(random_walk_benchmark%forecasts(1, :) - series(94, :))) < &
      1.0e-14_dp .and. abs(random_walk_benchmark%loss(1) - &
      sum(abs(series(97, :) - series(94, :)))) < 1.0e-14_dp .and. &
      random_walk_benchmark%standard_error >= 0.0_dp, &
      'multi-step random-walk benchmark')

   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_basic, &
      100.0_dp, 4, minnesota_target=minnesota_target)
   zero_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_basic, 1.0e-9_dp, &
      minnesota_target=minnesota_target)
   call check(zero_fit%info == 0 .and. zero_fit%minnesota .and. &
      zero_fit%nonzero == 0 .and. &
      maxval(abs(zero_fit%intercept)) < 1.0e-14_dp .and. &
      abs(zero_fit%phi(1, 1) - 1.0_dp) < 1.0e-12_dp .and. &
      abs(zero_fit%phi(2, 2) - 0.8_dp) < 1.0e-12_dp .and. &
      maxval(abs(zero_fit%phi(:, 3:4))) < 1.0e-12_dp, &
      'Minnesota zero-deviation VAR target')
   do time = 3, 100
      prediction = matmul(zero_fit%phi(:, 1:2), series(time - 1, :)) + &
         matmul(zero_fit%phi(:, 3:4), series(time - 2, :))
      call check(maxval(abs(series(time, :) - prediction - &
         zero_fit%residuals(time - 2, :))) < 1.0e-12_dp, &
         'Minnesota VAR residual identity')
   end do

   path = bigvar_structured_var_path(series, 2, lambdas, &
      bigvar_structure_basic, 1.0e-8_dp, &
      minnesota_target=minnesota_target)
   call check(path%info == 0 .and. path%minnesota .and. &
      path%nonzero(1) == 0 .and. &
      maxval(abs(path%shrinkage_target - minnesota_target)) < 1.0e-14_dp, &
      'Minnesota VAR warm-started path')
   relaxed_fit = bigvar_structured_var(series, 2, lambdas(2), &
      bigvar_structure_basic, 1.0e-8_dp, refit_fraction=0.5_dp, &
      minnesota_target=minnesota_target)
   call check(relaxed_fit%info == 0 .and. relaxed_fit%minnesota .and. &
      maxval(abs(relaxed_fit%intercept)) < 1.0e-14_dp, &
      'Minnesota relaxed VAR refit')
   validation = bigvar_var_validate(series, 2, lambdas(1:2), &
      bigvar_structure_basic, 94, 96, 1, tolerance=1.0e-7_dp, &
      max_iterations=5000, minnesota_target=minnesota_target)
   call check((validation%info == 0 .or. validation%info == 4) .and. &
      all(validation%valid), 'Minnesota rolling VAR validation')
   validation = bigvar_var_validate_loo(series, 2, lambdas(1:2), &
      bigvar_structure_basic, tolerance=1.0e-7_dp, max_iterations=5000, &
      first_observation=97, last_observation=98, &
      minnesota_target=minnesota_target)
   call check((validation%info == 0 .or. validation%info == 4) .and. &
      all(validation%valid), 'Minnesota leave-one-out VAR validation')
   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_basic, &
      100.0_dp, 3, direct_horizon=3, minnesota_target=minnesota_target)
   direct_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_basic, 1.0e-8_dp, direct_horizon=3, &
      minnesota_target=minnesota_target)
   call check(direct_fit%info == 0 .and. direct_fit%minnesota .and. &
      direct_fit%direct .and. direct_fit%forecast_horizon == 3 .and. &
      direct_fit%nonzero == 0, 'direct Minnesota VAR target')
   separate_lambdas = bigvar_separate_lambda_grid(series, 2, &
      bigvar_structure_basic, 100.0_dp, 3, &
      minnesota_target=minnesota_target)
   separate_fit = bigvar_structured_var_separate(series, 2, &
      separate_lambdas(1, :), bigvar_structure_basic, 1.0e-8_dp, &
      minnesota_target=minnesota_target)
   call check(separate_fit%info == 0 .and. separate_fit%minnesota .and. &
      separate_fit%nonzero == 0 .and. &
      maxval(abs(separate_fit%intercept)) < 1.0e-14_dp, &
      'response-specific Minnesota VAR target')
   separate_path = bigvar_structured_var_separate_path(series, 2, &
      separate_lambdas, bigvar_structure_basic, 1.0e-8_dp, &
      minnesota_target=minnesota_target)
   call check(separate_path%info == 0 .and. separate_path%minnesota .and. &
      separate_path%nonzero(1) == 0, &
      'response-specific Minnesota VAR path')
   separate_validation = bigvar_var_validate_separate(series, 2, &
      separate_lambdas(1:2, :), bigvar_structure_basic, 94, 96, 1, &
      tolerance=1.0e-7_dp, max_iterations=5000, &
      minnesota_target=minnesota_target)
   call check((separate_validation%info == 0 .or. &
      separate_validation%info == 4) .and. &
      all(separate_validation%valid), &
      'response-specific Minnesota rolling validation')
   reselection = bigvar_var_reselect_separate(series, 2, &
      separate_lambdas(1:2, :), bigvar_structure_basic, 94, 97, 98, 1, &
      selection_window=3, tolerance=1.0e-7_dp, max_iterations=5000)
   call check((reselection%info == 0 .or. reselection%info == 4) .and. &
      all(reselection%valid) .and. &
      all(shape(reselection%selected_lambda) == [2, 2]) .and. &
      all(reselection%selected_index >= 1), &
      'response-specific VAR rolling reselection')

   lag_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_lag, 1.0e-9_dp)
   call check(lag_fit%info == 0 .and. lag_fit%converged .and. &
      all(shape(lag_fit%phi) == [2, 4]) .and. &
      all(shape(lag_fit%residuals) == [98, 2]), &
      'Lag structured VAR convergence and dimensions')
   do time = 3, 100
      prediction = lag_fit%intercept + &
         matmul(lag_fit%phi(:, 1:2), series(time - 1, :)) + &
         matmul(lag_fit%phi(:, 3:4), series(time - 2, :))
      call check(maxval(abs(series(time, :) - prediction - &
         lag_fit%residuals(time - 2, :))) < 1.0e-12_dp, &
         'structured VAR residual identity')
   end do

   validation_lambdas = [0.05_dp, 0.005_dp]
   validation_alphas = [0.2_dp, 0.8_dp]
   dual_lambdas = bigvar_lambda_alpha_grid(series, 2, &
      bigvar_structure_sparse_lag, 100.0_dp, 3, validation_alphas)
   linear_lambda_matrix = bigvar_lambda_alpha_grid(series, 2, &
      bigvar_structure_sparse_lag, 100.0_dp, 3, validation_alphas, &
      linear=.true.)
   call check(all(shape(dual_lambdas) == [3, 2]) .and. &
      maxval(abs(dual_lambdas(:, 1) - dual_lambdas(:, 2))) > 1.0e-10_dp, &
      'alpha-specific VAR lambda grids')
   call check(maxval(abs(linear_lambda_matrix(1, :) - &
      dual_lambdas(1, :))) < 1.0e-14_dp .and. &
      maxval(abs(linear_lambda_matrix(3, :) - &
      dual_lambdas(3, :))) < 1.0e-14_dp .and. &
      maxval(abs(linear_lambda_matrix(2, :) - 0.5_dp*( &
      linear_lambda_matrix(1, :) + linear_lambda_matrix(3, :)))) < &
      1.0e-14_dp, 'linearly spaced dual VAR grids')
   dual_path = bigvar_structured_var_dual_path(series, 2, dual_lambdas, &
      validation_alphas, bigvar_structure_sparse_lag, 1.0e-8_dp)
   call check(dual_path%info == 0 .and. &
      all(shape(dual_path%phi) == [2, 4, 3, 2]) .and. &
      all(dual_path%nonzero(1, :) == 0), 'dual VAR path layout')
   validation = bigvar_var_validate_dual(series, 2, dual_lambdas, &
      validation_alphas, bigvar_structure_sparse_lag, 94, 96, 1, &
      one_standard_error=.true., tolerance=1.0e-7_dp, max_iterations=5000)
   call check((validation%info == 0 .or. validation%info == 4) .and. &
      all(shape(validation%mean_loss_surface) == [3, 2]) .and. &
      validation%selected_lambda == validation%lambda_grid( &
      validation%selected_lambda_index, validation%selected_alpha_index) .and. &
      validation%selected_alpha == &
      validation_alphas(validation%selected_alpha_index), &
      'dual rolling VAR validation surface')
   validation = bigvar_var_validate_loo_dual(series, 2, dual_lambdas, &
      validation_alphas, bigvar_structure_sparse_lag, tolerance=1.0e-7_dp, &
      max_iterations=5000, first_observation=97, last_observation=98)
   call check((validation%info == 0 .or. validation%info == 4) .and. &
      all(validation%valid), 'dual leave-one-out VAR validation')
   reselection = bigvar_var_reselect_dual(series, 2, dual_lambdas, &
      validation_alphas, bigvar_structure_sparse_lag, 94, 97, 98, 2, &
      one_standard_error=.true., selection_window=3, tolerance=1.0e-7_dp, &
      max_iterations=5000, recursive=.false.)
   call check((reselection%info == 0 .or. reselection%info == 4) .and. &
      all(reselection%valid) .and. &
      all(shape(reselection%forecasts) == [2, 2]) .and. &
      all(reselection%selected_lambda_index >= 1) .and. &
      all(reselection%selected_alpha_index >= 1), &
      'dual direct VAR rolling reselection')
   call check(abs(reselection%loss(1, 1) - sum((series(99, :) - &
      reselection%forecasts(1, :))**2)) < 1.0e-14_dp, &
      'rolling reselection forecast alignment')
   reselection = bigvar_var_reselect(series, 2, validation_lambdas, &
      bigvar_structure_sparse_lag, 94, 97, 98, 1, &
      alphas=validation_alphas, selection_window=3, tolerance=1.0e-7_dp, &
      max_iterations=5000)
   call check((reselection%info == 0 .or. reselection%info == 4) .and. &
      all(reselection%valid) .and. all(reselection%selected_lambda > 0.0_dp), &
      'ordinary VAR rolling reselection')
   validation = bigvar_var_validate(series, 2, validation_lambdas, &
      bigvar_structure_sparse_lag, 94, 96, 2, bigvar_loss_huber, 0.05_dp, &
      .true., 40, validation_alphas, 1.0e-7_dp, 5000)
   call check((validation%info == 0 .or. validation%info == 4) .and. &
      all(shape(validation%loss) == [3, 4]) .and. &
      all(shape(validation%forecasts) == [3, 4, 2]) .and. &
      all(validation%valid), 'rolling VAR validation dimensions')
   call check(maxval(abs(validation%lambda - &
      [0.05_dp, 0.005_dp, 0.05_dp, 0.005_dp])) < 1.0e-14_dp .and. &
      maxval(abs(validation%alpha - &
      [0.2_dp, 0.2_dp, 0.8_dp, 0.8_dp])) < 1.0e-14_dp, &
      'joint lambda-alpha candidate ordering')
   call check(abs(validation%loss(1, 1) - bigvar_forecast_loss( &
      series(96, :) - validation%forecasts(1, 1, :), &
      bigvar_loss_huber, 0.05_dp)) < 1.0e-14_dp, &
      'multi-step VAR validation origin alignment')
   call check(validation%best_index >= 1 .and. validation%best_index <= 4 .and. &
      validation%one_se_index >= 1 .and. validation%one_se_index <= 4 .and. &
      validation%selected_index == validation%one_se_index .and. &
      validation%mean_loss(validation%one_se_index) < &
      minval(validation%mean_loss) + validation%selection_standard_error, &
      'one-standard-error penalty selection')
   evaluation = bigvar_var_evaluate(series, 2, &
      validation%selected_lambda, bigvar_structure_sparse_lag, 97, 98, 1, &
      bigvar_loss_l1, alpha=validation%selected_alpha, tolerance=1.0e-7_dp, &
      max_iterations=5000)
   call check((evaluation%info == 0 .or. evaluation%info == 4) .and. &
      all(shape(evaluation%forecasts) == [2, 1, 2]) .and. &
      evaluation%selected_index == 1 .and. all(evaluation%valid), &
      'rolling evaluation with selected VAR penalty')
   evaluation = bigvar_var_validate(series, 2, [0.01_dp], &
      bigvar_structure_hlag_element, 95, 96, 1, tolerance=1.0e-7_dp, &
      max_iterations=5000)
   call check((evaluation%info == 0 .or. evaluation%info == 4) .and. &
      all(evaluation%valid) .and. evaluation%selected_index == 1, &
      'elementwise HLAG rolling validation')
   separate_lambdas = bigvar_separate_lambda_grid(series, 2, &
      bigvar_structure_basic, 100.0_dp, 5)
   separate_validation = bigvar_var_validate_separate(series, 2, &
      separate_lambdas([1, 3, 5], :), bigvar_structure_basic, 94, 96, 1, &
      bigvar_loss_l2, one_standard_error=.true., tolerance=1.0e-7_dp, &
      max_iterations=5000)
   call check((separate_validation%info == 0 .or. &
      separate_validation%info == 4) .and. &
      all(shape(separate_validation%loss) == [3, 3, 2]) .and. &
      all(separate_validation%valid) .and. &
      all(separate_validation%selected_index >= 1) .and. &
      all(separate_validation%selected_index <= 3), &
      'response-specific rolling validation')
   call check(abs(separate_validation%loss(1, 2, 1) - &
      (series(95, 1) - separate_validation%forecasts(1, 2, 1))**2) < &
      1.0e-14_dp .and. all(separate_validation%selected_lambda > 0.0_dp), &
      'per-response validation loss and selected penalties')
   evaluation = bigvar_var_validate(series, 2, [0.01_dp, 0.005_dp], &
      bigvar_structure_basic, 94, 96, 1, tolerance=1.0e-7_dp, &
      max_iterations=5000, refit_fraction=1.0_dp)
   call check((evaluation%info == 0 .or. evaluation%info == 4) .and. &
      all(evaluation%valid), 'relaxed rolling VAR validation')
   evaluation = bigvar_var_validate(series, 2, [0.01_dp], &
      bigvar_structure_basic, 94, 96, 3, tolerance=1.0e-7_dp, &
      max_iterations=5000, recursive=.false., window_size=40)
   call check((evaluation%info == 0 .or. evaluation%info == 4) .and. &
      all(evaluation%valid) .and. &
      all(shape(evaluation%forecasts) == [3, 1, 2]), &
      'direct rolling VAR validation')
   direct_fit = bigvar_structured_var(series(52:94, :), 2, 0.01_dp, &
      bigvar_structure_basic, 1.0e-7_dp, 5000, direct_horizon=3)
   direct_prediction = bigvar_direct_forecast(direct_fit, series(52:94, :))
   call check(maxval(abs(evaluation%forecasts(1, 1, :) - &
      direct_prediction)) < 1.0e-12_dp .and. &
      abs(evaluation%loss(1, 1) - sum((series(97, :) - &
      direct_prediction)**2)) < 1.0e-12_dp, &
      'direct rolling origin and window alignment')
   separate_lambdas = spread([0.01_dp, 0.005_dp], 2, 2)
   separate_validation = bigvar_var_validate_separate(series, 2, &
      separate_lambdas, bigvar_structure_basic, 94, 96, 3, &
      recursive=.false., tolerance=1.0e-7_dp, max_iterations=5000)
   call check((separate_validation%info == 0 .or. &
      separate_validation%info == 4) .and. &
      all(separate_validation%valid), &
      'direct response-specific rolling validation')
   evaluation = bigvar_var_evaluate(series, 2, 0.01_dp, &
      bigvar_structure_basic, 94, 96, 3, recursive=.false., &
      tolerance=1.0e-7_dp, max_iterations=5000)
   call check((evaluation%info == 0 .or. evaluation%info == 4) .and. &
      all(evaluation%valid), 'direct rolling VAR evaluation')
   validation = bigvar_var_validate_loo(series, 2, &
      [0.01_dp, 0.005_dp], bigvar_structure_basic, &
      one_standard_error=.true., tolerance=1.0e-7_dp, &
      max_iterations=5000, refit_fraction=0.5_dp, &
      first_observation=96, last_observation=98)
   call check((validation%info == 0 .or. validation%info == 4) .and. &
      all(shape(validation%loss) == [3, 2]) .and. all(validation%valid) .and. &
      validation%selected_index >= 1 .and. validation%selected_index <= 2, &
      'bounded leave-one-out VAR validation')
   call check(abs(validation%loss(1, 1) - bigvar_forecast_loss( &
      series(96, :) - validation%forecasts(1, 1, :), bigvar_loss_l2)) < &
      1.0e-14_dp, 'leave-one-out VAR loss identity')
   validation = bigvar_var_validate_loo(series, 2, [0.01_dp], &
      bigvar_structure_basic, tolerance=1.0e-7_dp, max_iterations=5000, &
      horizon=2, recursive=.false., first_observation=96, &
      last_observation=97)
   call check((validation%info == 0 .or. validation%info == 4) .and. &
      all(validation%valid) .and. validation%horizon == 2, &
      'direct multi-step leave-one-out VAR validation')
   separate_lambdas = spread([0.01_dp, 0.005_dp], 2, 2)
   separate_validation = bigvar_var_validate_separate_loo(series, 2, &
      separate_lambdas, bigvar_structure_basic, tolerance=1.0e-7_dp, &
      max_iterations=5000, first_observation=96, last_observation=97)
   call check((separate_validation%info == 0 .or. &
      separate_validation%info == 4) .and. &
      all(shape(separate_validation%loss) == [2, 2, 2]) .and. &
      all(separate_validation%valid), &
      'response-specific leave-one-out VAR validation')

   own_other_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_own_other, 1.0e-9_dp)
   call check(own_other_fit%info == 0 .and. own_other_fit%converged .and. &
      own_other_fit%active_groups >= 1 .and. &
      own_other_fit%active_groups <= 4, &
      'Own/Other structured VAR convergence')

   sparse_lag_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_sparse_lag, 1.0e-9_dp, alpha=0.25_dp)
   call check(sparse_lag_fit%info == 0 .and. &
      sparse_lag_fit%converged .and. &
      abs(sparse_lag_fit%alpha - 0.25_dp) < 1.0e-14_dp .and. &
      sparse_lag_fit%nonzero >= sparse_lag_fit%active_groups, &
      'SparseLag structured VAR convergence')
   sparse_own_other_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_sparse_own_other, 1.0e-9_dp)
   call check(sparse_own_other_fit%info == 0 .and. &
      sparse_own_other_fit%converged .and. &
      abs(sparse_own_other_fit%alpha - 1.0_dp/3.0_dp) < 1.0e-14_dp, &
      'SparseOO default alpha and convergence')
   hlag_component_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_hlag_component, 1.0e-9_dp)
   call check(hlag_component_fit%info == 0 .and. &
      hlag_component_fit%converged .and. &
      hlag_component_fit%active_groups >= 1 .and. &
      hlag_component_fit%active_groups <= 4, &
      'componentwise HLAG convergence')
   hlag_own_other_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_hlag_own_other, 1.0e-9_dp)
   call check(hlag_own_other_fit%info == 0 .and. &
      hlag_own_other_fit%converged .and. &
      hlag_own_other_fit%active_groups >= 1 .and. &
      hlag_own_other_fit%active_groups <= 8, &
      'Own/Other HLAG convergence')
   hlag_element_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_hlag_element, 1.0e-9_dp)
   call check(hlag_element_fit%info == 0 .and. &
      hlag_element_fit%converged .and. &
      hlag_element_fit%active_groups >= 1 .and. &
      hlag_element_fit%active_groups <= 8 .and. &
      all(shape(hlag_element_fit%phi) == [2, 4]), &
      'elementwise HLAG convergence and diagnostics')
   basic_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_basic, 1.0e-9_dp)
   separate_fit = bigvar_structured_var_separate(series, 2, &
      [0.01_dp, 0.01_dp], bigvar_structure_basic, 1.0e-9_dp)
   call check(basic_fit%info == 0 .and. separate_fit%info == 0 .and. &
      maxval(abs(basic_fit%phi - separate_fit%phi)) < 1.0e-12_dp .and. &
      maxval(abs(separate_fit%lambda_by_response - 0.01_dp)) < 1.0e-14_dp, &
      'equal separate penalties reproduce the Basic fit')
   relaxed_fit = bigvar_relaxed_var(series, basic_fit, 0.0_dp)
   call check(relaxed_fit%info == 0 .and. relaxed_fit%relaxed .and. &
      relaxed_fit%refit_fraction == 0.0_dp .and. &
      maxval(abs(relaxed_fit%phi - basic_fit%phi)) < 1.0e-14_dp, &
      'zero-fraction relaxed VAR preserves penalized coefficients')
   relaxed_fit = bigvar_relaxed_var(series, basic_fit, 1.0_dp)
   call check(relaxed_fit%info == 0 .and. relaxed_fit%relaxed .and. &
      relaxed_fit%refit_fraction == 1.0_dp .and. &
      all((abs(relaxed_fit%phi) > 0.0_dp) .eqv. &
      (abs(basic_fit%phi) > 1.0e-8_dp)), &
      'relaxed VAR preserves the selected support')
   direct_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_basic, 1.0e-9_dp, direct_horizon=3)
   call check(direct_fit%info == 0 .and. direct_fit%direct .and. &
      direct_fit%forecast_horizon == 3 .and. &
      all(shape(direct_fit%residuals) == [96, 2]), &
      'direct horizon-specific VAR dimensions')
   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_basic, &
      100.0_dp, 5, direct_horizon=3)
   zero_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_basic, 1.0e-9_dp, direct_horizon=3)
   call check(zero_fit%info == 0 .and. zero_fit%nonzero == 0, &
      'direct VAR lambda bound gives the zero model')
   do time = 5, 100
      prediction = direct_fit%intercept + &
         matmul(direct_fit%phi(:, 1:2), series(time - 3, :)) + &
         matmul(direct_fit%phi(:, 3:4), series(time - 4, :))
      call check(maxval(abs(series(time, :) - prediction - &
         direct_fit%residuals(time - 4, :))) < 1.0e-12_dp, &
         'direct VAR shifted residual identity')
   end do
   direct_prediction = bigvar_direct_forecast(direct_fit, series)
   prediction = direct_fit%intercept + &
      matmul(direct_fit%phi(:, 1:2), series(100, :)) + &
      matmul(direct_fit%phi(:, 3:4), series(99, :))
   call check(size(direct_prediction) == 2 .and. &
      maxval(abs(direct_prediction - prediction)) < 1.0e-14_dp, &
      'direct VAR latest-origin forecast')
   forecast = bigvar_forecast(direct_fit, series, 3)
   call check(forecast%info == 0 .and. &
      maxval(abs(forecast%mean(3, :) - prediction)) < 1.0e-14_dp, &
      'direct VAR forecast adapter')
   relaxed_fit = bigvar_relaxed_var(series, direct_fit, 0.5_dp)
   call check(relaxed_fit%info == 0 .and. relaxed_fit%direct .and. &
      relaxed_fit%forecast_horizon == 3 .and. relaxed_fit%relaxed, &
      'relaxed direct VAR retains horizon metadata')
   path = bigvar_structured_var_path(series, 2, [0.02_dp, 0.01_dp], &
      bigvar_structure_basic, 1.0e-8_dp, refit_fraction=0.5_dp)
   call check(path%info == 0 .and. all(path%converged), &
      'relaxed scalar-lambda path')
   path = bigvar_structured_var_path(series, 2, [0.02_dp, 0.01_dp], &
      bigvar_structure_basic, 1.0e-8_dp, direct_horizon=3)
   call check(path%info == 0 .and. path%direct .and. &
      path%forecast_horizon == 3 .and. all(path%converged), &
      'direct scalar-lambda path')
   separate_fit = bigvar_structured_var_separate(series, 2, &
      [0.01_dp, 0.02_dp], bigvar_structure_basic, 1.0e-8_dp, &
      refit_fraction=0.5_dp)
   call check(separate_fit%info == 0 .and. separate_fit%relaxed .and. &
      abs(separate_fit%refit_fraction - 0.5_dp) < 1.0e-14_dp, &
      'integrated relaxed response-specific fit')
   separate_fit = bigvar_structured_var_separate(series, 2, &
      [0.01_dp, 0.02_dp], bigvar_structure_basic_en, 1.0e-8_dp, &
      alpha=0.25_dp)
   call check(separate_fit%info == 0 .and. separate_fit%converged, &
      'separate BasicEN penalties')
   separate_fit = bigvar_structured_var_separate(series, 2, &
      [0.01_dp, 0.02_dp], bigvar_structure_hlag_component, 1.0e-8_dp)
   call check(separate_fit%info == 0 .and. separate_fit%converged, &
      'separate componentwise HLAG penalties')
   separate_fit = bigvar_structured_var_separate(series, 2, &
      [0.01_dp, 0.02_dp], bigvar_structure_hlag_own_other, 1.0e-8_dp)
   call check(separate_fit%info == 0 .and. separate_fit%converged, &
      'separate Own/Other HLAG penalties')
   separate_fit = bigvar_structured_var_separate(series, 2, &
      [0.01_dp, 0.02_dp], bigvar_structure_hlag_element, 1.0e-8_dp)
   call check(separate_fit%info == 0 .and. separate_fit%converged, &
      'separate elementwise HLAG penalties')
   separate_fit = bigvar_structured_var_separate(series, 2, &
      [0.005_dp, 0.01_dp], bigvar_structure_mcp, 1.0e-8_dp, gamma=3.0_dp)
   call check(separate_fit%info == 0 .and. separate_fit%converged, &
      'separate MCP penalties')
   separate_fit = bigvar_structured_var_separate(series, 2, &
      [0.005_dp, 0.01_dp], bigvar_structure_scad, 1.0e-8_dp, gamma=3.5_dp)
   call check(separate_fit%info == 0 .and. separate_fit%converged, &
      'separate SCAD penalties')
   basic_en_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_basic_en, 1.0e-9_dp, alpha=0.25_dp)
   call check(basic_en_fit%info == 0 .and. basic_en_fit%converged .and. &
      abs(basic_en_fit%alpha - 0.25_dp) < 1.0e-14_dp .and. &
      basic_en_fit%objective < huge(1.0_dp) .and. &
      all(shape(basic_en_fit%residuals) == [98, 2]), &
      'BasicEN convergence and diagnostics')
   tapered_fit = bigvar_structured_var(series, 2, 0.01_dp, &
      bigvar_structure_tapered, 1.0e-9_dp, alpha=1.0_dp)
   call check(tapered_fit%info == 0 .and. tapered_fit%converged .and. &
      abs(tapered_fit%alpha - 1.0_dp) < 1.0e-14_dp, &
      'Tapered structured VAR convergence')
   mcp_fit = bigvar_structured_var(series, 2, 0.005_dp, &
      bigvar_structure_mcp, 1.0e-8_dp, gamma=3.0_dp)
   call check(mcp_fit%info == 0 .and. mcp_fit%converged .and. &
      abs(mcp_fit%gamma - 3.0_dp) < 1.0e-14_dp .and. &
      mcp_fit%objective < huge(1.0_dp), 'MCP coordinate-descent fit')
   scad_fit = bigvar_structured_var(series, 2, 0.005_dp, &
      bigvar_structure_scad, 1.0e-8_dp, gamma=3.5_dp)
   call check(scad_fit%info == 0 .and. scad_fit%converged .and. &
      abs(scad_fit%gamma - 3.5_dp) < 1.0e-14_dp .and. &
      all(shape(scad_fit%residuals) == [98, 2]), &
      'SCAD coordinate-descent fit')

   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_lag, 100.0_dp, 5)
   linear_lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_lag, &
      100.0_dp, 5, linear=.true.)
   call check(abs(linear_lambdas(1) - lambdas(1)) < 1.0e-14_dp .and. &
      abs(linear_lambdas(5) - lambdas(5)) < 1.0e-14_dp .and. &
      maxval(abs((linear_lambdas(2:5) - linear_lambdas(1:4)) - &
      (linear_lambdas(2) - linear_lambdas(1)))) < 1.0e-14_dp, &
      'linearly spaced VAR lambda grid')
   call check(size(lambdas) == 5 .and. all(lambdas > 0.0_dp) .and. &
      all(lambdas(1:4) > lambdas(2:5)), 'geometric lambda grid')
   zero_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_lag, 1.0e-9_dp)
   call check(zero_fit%info == 0 .and. zero_fit%nonzero == 0 .and. &
      zero_fit%active_groups == 0, 'lambda bound gives the zero model')

   path = bigvar_structured_var_path(series, 2, lambdas, &
      bigvar_structure_lag, 1.0e-8_dp)
   call check(path%info == 0 .and. all(path%converged) .and. &
      all(shape(path%phi) == [2, 4, 5]) .and. &
      path%active_groups(1) == 0 .and. &
      path%active_groups(5) >= path%active_groups(1), &
      'warm-started structured VAR path')
   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_sparse_lag, &
      100.0_dp, 5, 0.25_dp)
   zero_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_sparse_lag, 1.0e-9_dp, alpha=0.25_dp)
   call check(zero_fit%info == 0 .and. zero_fit%nonzero == 0, &
      'SparseLag lambda bound gives the zero model')
   path = bigvar_structured_var_path(series, 2, lambdas, &
      bigvar_structure_sparse_lag, 1.0e-8_dp, alpha=0.25_dp)
   call check(path%info == 0 .and. all(path%converged) .and. &
      abs(path%alpha - 0.25_dp) < 1.0e-14_dp .and. &
      path%nonzero(1) == 0, 'SparseLag warm-started path')
   lambdas = bigvar_lambda_grid(series, 2, &
      bigvar_structure_hlag_component, 100.0_dp, 5)
   zero_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_hlag_component, 1.0e-9_dp)
   call check(zero_fit%info == 0 .and. zero_fit%nonzero == 0, &
      'componentwise HLAG lambda bound gives the zero model')
   path = bigvar_structured_var_path(series, 2, lambdas, &
      bigvar_structure_hlag_own_other, 1.0e-8_dp)
   call check(path%info == 0 .and. all(path%converged) .and. &
      all(shape(path%phi) == [2, 4, 5]), &
      'Own/Other HLAG warm-started path')
   lambdas = bigvar_lambda_grid(series, 2, &
      bigvar_structure_hlag_element, 100.0_dp, 5)
   zero_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_hlag_element, 1.0e-9_dp)
   call check(zero_fit%info == 0 .and. zero_fit%nonzero == 0 .and. &
      zero_fit%active_groups == 0, &
      'elementwise HLAG lambda bound gives the zero model')
   path = bigvar_structured_var_path(series, 2, lambdas, &
      bigvar_structure_hlag_element, 1.0e-8_dp)
   call check(path%info == 0 .and. all(path%converged) .and. &
      all(shape(path%phi) == [2, 4, 5]) .and. &
      path%active_groups(1) == 0, &
      'elementwise HLAG warm-started path')
   separate_lambdas = bigvar_separate_lambda_grid(series, 2, &
      bigvar_structure_basic, 100.0_dp, 5)
   linear_lambda_matrix = bigvar_separate_lambda_grid(series, 2, &
      bigvar_structure_basic, 100.0_dp, 5, linear=.true.)
   call check(maxval(abs(linear_lambda_matrix(1, :) - &
      separate_lambdas(1, :))) < 1.0e-14_dp .and. &
      maxval(abs(linear_lambda_matrix(5, :) - &
      separate_lambdas(5, :))) < 1.0e-14_dp .and. &
      maxval(abs((linear_lambda_matrix(2:5, :) - &
      linear_lambda_matrix(1:4, :)) - spread( &
      linear_lambda_matrix(2, :) - linear_lambda_matrix(1, :), 1, 4))) < &
      1.0e-14_dp, 'linearly spaced response-specific grids')
   call check(all(shape(separate_lambdas) == [5, 2]) .and. &
      all(separate_lambdas > 0.0_dp) .and. &
      all(separate_lambdas(1:4, :) > separate_lambdas(2:5, :)), &
      'response-specific lambda grids')
   separate_fit = bigvar_structured_var_separate(series, 2, &
      separate_lambdas(1, :), bigvar_structure_basic, 1.0e-9_dp)
   call check(separate_fit%info == 0 .and. separate_fit%nonzero == 0, &
      'separate lambda bounds give the zero model')
   separate_path = bigvar_structured_var_separate_path(series, 2, &
      separate_lambdas, bigvar_structure_basic, 1.0e-8_dp)
   call check(separate_path%info == 0 .and. &
      all(separate_path%converged) .and. &
      all(shape(separate_path%phi) == [2, 4, 5]) .and. &
      all(shape(separate_path%lambda) == [5, 2]) .and. &
      separate_path%nonzero(1) == 0, &
      'response-specific warm-started path')
   separate_path = bigvar_structured_var_separate_path(series, 2, &
      separate_lambdas, bigvar_structure_basic, 1.0e-8_dp, &
      refit_fraction=0.5_dp)
   call check(separate_path%info == 0 .and. all(separate_path%converged), &
      'relaxed response-specific warm-started path')
   separate_lambdas = bigvar_separate_lambda_grid(series, 2, &
      bigvar_structure_basic, 100.0_dp, 5, direct_horizon=3)
   separate_fit = bigvar_structured_var_separate(series, 2, &
      separate_lambdas(1, :), bigvar_structure_basic, 1.0e-9_dp, &
      direct_horizon=3)
   call check(separate_fit%info == 0 .and. separate_fit%direct .and. &
      separate_fit%nonzero == 0, 'direct separate lambda zero model')
   separate_path = bigvar_structured_var_separate_path(series, 2, &
      separate_lambdas, bigvar_structure_basic, 1.0e-8_dp, &
      direct_horizon=3)
   call check(separate_path%info == 0 .and. separate_path%direct .and. &
      separate_path%forecast_horizon == 3 .and. &
      all(separate_path%converged), 'direct response-specific path')
   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_basic_en, &
      100.0_dp, 5, 0.25_dp)
   zero_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_basic_en, 1.0e-9_dp, alpha=0.25_dp)
   call check(zero_fit%info == 0 .and. zero_fit%nonzero == 0, &
      'BasicEN lambda bound gives the zero model')
   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_tapered, &
      100.0_dp, 5, 1.0_dp)
   path = bigvar_structured_var_path(series, 2, lambdas, &
      bigvar_structure_tapered, 1.0e-8_dp, alpha=1.0_dp)
   call check(path%info == 0 .and. all(path%converged) .and. &
      path%nonzero(1) == 0 .and. abs(path%alpha - 1.0_dp) < 1.0e-14_dp, &
      'Tapered warm-started path')
   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_mcp, &
      100.0_dp, 5)
   zero_fit = bigvar_structured_var(series, 2, lambdas(1), &
      bigvar_structure_mcp, 1.0e-8_dp, gamma=3.0_dp)
   call check(zero_fit%info == 0 .and. zero_fit%nonzero == 0, &
      'MCP lambda bound gives the zero model')
   path = bigvar_structured_var_path(series, 2, lambdas, &
      bigvar_structure_mcp, 1.0e-8_dp, gamma=3.0_dp)
   call check(path%info == 0 .and. all(path%converged) .and. &
      path%nonzero(1) == 0 .and. abs(path%gamma - 3.0_dp) < 1.0e-14_dp, &
      'MCP warm-started non-convex path')
   lambdas = bigvar_lambda_grid(series, 2, bigvar_structure_scad, &
      100.0_dp, 5)
   path = bigvar_structured_var_path(series, 2, lambdas, &
      bigvar_structure_scad, 1.0e-8_dp, gamma=3.5_dp)
   call check(path%info == 0 .and. all(path%converged) .and. &
      path%nonzero(1) == 0 .and. abs(path%gamma - 3.5_dp) < 1.0e-14_dp, &
      'SCAD warm-started non-convex path')

   bgr_zero_fit = bigvar_bgr(series, 2, 3.0_dp)
   random_walk_prior = [.true., .true.]
   bgr_random_walk_fit = bigvar_bgr(series, 2, 3.0_dp, random_walk_prior)
   call check(bgr_zero_fit%info == 0 .and. bgr_zero_fit%converged .and. &
      all(shape(bgr_zero_fit%phi) == [2, 4]) .and. &
      all(shape(bgr_zero_fit%residuals) == [98, 2]) .and. &
      bgr_zero_fit%structure == bigvar_structure_bgr, &
      'BGR fixed fit and dimensions')
   do time = 3, 100
      prediction = bgr_zero_fit%intercept + &
         matmul(bgr_zero_fit%phi(:, 1:2), series(time - 1, :)) + &
         matmul(bgr_zero_fit%phi(:, 3:4), series(time - 2, :))
      call check(maxval(abs(series(time, :) - prediction - &
         bgr_zero_fit%residuals(time - 2, :))) < 1.0e-12_dp, &
         'BGR residual identity')
   end do
   call check(bgr_random_walk_fit%info == 0 .and. &
      maxval(abs(bgr_random_walk_fit%phi - bgr_zero_fit%phi)) > &
      1.0e-8_dp, 'BGR random-walk dummy changes posterior coefficients')
   bgr_grid = bigvar_bgr_default_grid(2, 2)
   call check(size(bgr_grid) == 161 .and. &
      abs(bgr_grid(1) - 2.0_dp) < 1.0e-14_dp .and. &
      abs(bgr_grid(161) - 10.0_dp) < 1.0e-14_dp, &
      'BGR package-default tightness grid')
   path = bigvar_bgr_path(series, 2, bgr_grid([1, 41, 81]), &
      random_walk_prior)
   call check(path%info == 0 .and. all(path%converged) .and. &
      all(shape(path%phi) == [2, 4, 3]) .and. &
      path%structure == bigvar_structure_bgr, 'BGR tightness path')
   forecast = bigvar_forecast(bgr_random_walk_fit, series, 2)
   prediction = bgr_random_walk_fit%intercept + &
      matmul(bgr_random_walk_fit%phi(:, 1:2), series(100, :)) + &
      matmul(bgr_random_walk_fit%phi(:, 3:4), series(99, :))
   call check(forecast%info == 0 .and. &
      maxval(abs(forecast%mean(1, :) - prediction)) < 1.0e-14_dp, &
      'BGR forecast adapter')

   do time = 1, 100
      exogenous(time, 1) = sin(0.17_dp*real(time, dp)) + &
         0.08_dp*cos(0.013_dp*real(time*time, dp))
      exogenous(time, 2) = cos(0.11_dp*real(time, dp)) - &
         0.06_dp*sin(0.019_dp*real(time*time, dp))
   end do
   exogenous_forecast(:100, :) = exogenous
   do time = 101, 103
      exogenous_forecast(time, 1) = sin(0.17_dp*real(time, dp)) + &
         0.08_dp*cos(0.013_dp*real(time*time, dp))
      exogenous_forecast(time, 2) = cos(0.11_dp*real(time, dp)) - &
         0.06_dp*sin(0.019_dp*real(time*time, dp))
   end do
   transfer_series(1, :) = [0.2_dp, -0.1_dp]
   transfer_series(2, :) = [-0.05_dp, 0.25_dp]
   do time = 3, 100
      transfer_series(time, 1) = 0.35_dp + &
         0.62_dp*exogenous(time - 1, 1) - &
         0.18_dp*exogenous(time - 2, 2) + &
         0.004_dp*sin(0.41_dp*real(time, dp))
      transfer_series(time, 2) = -0.22_dp - &
         0.48_dp*exogenous(time - 1, 2) + &
         0.16_dp*exogenous(time - 2, 1) + &
         0.004_dp*cos(0.37_dp*real(time, dp))
   end do
   varx_series(1, :) = [0.2_dp, -0.1_dp]
   varx_series(2, :) = [-0.05_dp, 0.25_dp]
   do time = 3, 100
      varx_series(time, 1) = 0.48_dp*varx_series(time - 1, 1) - &
         0.14_dp*varx_series(time - 1, 2) + &
         0.30_dp*exogenous(time - 1, 1) + &
         0.12_dp*exogenous(time - 2, 2) + &
         0.01_dp*sin(0.71_dp*real(time, dp))
      varx_series(time, 2) = 0.34_dp*varx_series(time - 1, 2) + &
         0.11_dp*varx_series(time - 1, 1) - &
         0.25_dp*exogenous(time - 1, 2) + &
         0.01_dp*cos(0.67_dp*real(time, dp))
   end do
   contemporaneous_series(1, :) = [0.1_dp, -0.2_dp]
   do time = 2, 100
      contemporaneous_series(time, 1) = &
         0.35_dp*contemporaneous_series(time - 1, 1) + &
         0.42_dp*exogenous(time, 1) + &
         0.14_dp*exogenous(time - 1, 2) + &
         0.005_dp*sin(0.53_dp*real(time, dp))
      contemporaneous_series(time, 2) = &
         0.30_dp*contemporaneous_series(time - 1, 2) - &
         0.38_dp*exogenous(time, 2) + &
         0.10_dp*exogenous(time - 1, 1) + &
         0.005_dp*cos(0.47_dp*real(time, dp))
   end do

   lambdas = bigvar_varx_lambda_grid(transfer_series, exogenous, 0, 2, &
      bigvar_structure_basic, 100.0_dp, 3)
   call check(size(lambdas) == 3 .and. all(lambdas > 0.0_dp) .and. &
      all(lambdas(1:2) > lambdas(2:3)), &
      'transfer-function lambda grid')
   varx_zero = bigvar_structured_varx(transfer_series, exogenous, 0, 2, &
      lambdas(1), bigvar_structure_basic, 1.0e-8_dp)
   call check(varx_zero%info == 0 .and. varx_zero%transfer_function .and. &
      all(shape(varx_zero%phi) == [2, 0]) .and. &
      all(shape(varx_zero%beta) == [2, 4]) .and. varx_zero%nonzero == 0, &
      'Basic transfer-function zero model')
   varx_path = bigvar_structured_varx_path(transfer_series, exogenous, 0, 2, &
      lambdas, bigvar_structure_basic, 1.0e-8_dp)
   call check(varx_path%info == 0 .and. varx_path%transfer_function .and. &
      all(shape(varx_path%phi) == [2, 0, 3]) .and. &
      all(shape(varx_path%beta) == [2, 4, 3]), &
      'Basic transfer-function path')
   varx_lag = bigvar_structured_varx(transfer_series, exogenous, 0, 2, &
      lambdas(3), bigvar_structure_basic, 1.0e-8_dp)
   call check(varx_lag%info == 0 .and. varx_lag%transfer_function .and. &
      varx_lag%nonzero > 0, 'Basic transfer-function fit')
   do time = 3, 100
      prediction = varx_lag%intercept + &
         matmul(varx_lag%beta(:, 1:2), exogenous(time - 1, :)) + &
         matmul(varx_lag%beta(:, 3:4), exogenous(time - 2, :))
      call check(maxval(abs(transfer_series(time, :) - prediction - &
         varx_lag%residuals(time - 2, :))) < 1.0e-12_dp, &
         'transfer-function residual identity')
   end do
   relaxed_varx = bigvar_relaxed_varx(transfer_series, exogenous, varx_lag, &
      0.5_dp)
   call check(relaxed_varx%info == 0 .and. relaxed_varx%relaxed .and. &
      relaxed_varx%transfer_function .and. &
      all(shape(relaxed_varx%phi) == [2, 0]), &
      'relaxed transfer-function refit')
   forecast = bigvar_varx_forecast(varx_lag, transfer_series, exogenous, 1)
   prediction = varx_lag%intercept + &
      matmul(varx_lag%beta(:, 1:2), exogenous(100, :)) + &
      matmul(varx_lag%beta(:, 3:4), exogenous(99, :))
   call check(forecast%info == 0 .and. &
      maxval(abs(forecast%mean(1, :) - prediction)) < 1.0e-14_dp, &
      'lagged transfer-function forecast')
   interval = bigvar_varx_forecast_interval(varx_lag, transfer_series, &
      exogenous, 1)
   call check(interval%info == 0 .and. &
      maxval(abs(interval%covariance(:, :, 1) - &
      bigvar_innovation_covariance(varx_lag%residuals))) < 1.0e-14_dp, &
      'transfer-function forecast interval')
   dual_lambdas = bigvar_varx_lambda_alpha_grid(transfer_series, exogenous, &
      0, 2, bigvar_structure_basic_en, 100.0_dp, 3, validation_alphas)
   varx_dual_path = bigvar_structured_varx_dual_path(transfer_series, &
      exogenous, 0, 2, dual_lambdas, validation_alphas, &
      bigvar_structure_basic_en, 1.0e-8_dp)
   call check(varx_dual_path%info == 0 .and. &
      varx_dual_path%transfer_function .and. &
      all(shape(varx_dual_path%phi) == [2, 0, 3, 2]), &
      'BasicEN transfer-function dual path')
   varx_mcp = bigvar_structured_varx(transfer_series, exogenous, 0, 2, &
      0.002_dp, bigvar_structure_mcp, 1.0e-8_dp)
   varx_scad = bigvar_structured_varx(transfer_series, exogenous, 0, 2, &
      0.002_dp, bigvar_structure_scad, 1.0e-8_dp)
   call check(varx_mcp%info == 0 .and. varx_scad%info == 0 .and. &
      varx_mcp%transfer_function .and. varx_scad%transfer_function, &
      'MCP and SCAD transfer-function fits')
   rejected_varx = bigvar_structured_varx(transfer_series, exogenous, 0, 2, &
      0.01_dp, bigvar_structure_lag)
   call check(rejected_varx%info == 1, &
      'grouped penalty rejects transfer function')
   current_only_varx = bigvar_structured_varx(transfer_series, exogenous, &
      0, 0, 0.001_dp, bigvar_structure_basic, 1.0e-8_dp, &
      contemporaneous=.true.)
   forecast = bigvar_varx_forecast(current_only_varx, transfer_series, &
      exogenous_forecast, 1)
   prediction = current_only_varx%intercept + &
      matmul(current_only_varx%beta, exogenous_forecast(101, :))
   call check(current_only_varx%info == 0 .and. &
      current_only_varx%transfer_function .and. forecast%info == 0 .and. &
      maxval(abs(forecast%mean(1, :) - prediction)) < 1.0e-14_dp, &
      'contemporaneous transfer-function forecast')
   varx_validation = bigvar_varx_validate(transfer_series, exogenous, 0, 2, &
      lambdas(1:2), bigvar_structure_basic, 94, 96, 1, &
      tolerance=1.0e-7_dp, max_iterations=5000)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. all(varx_validation%valid), &
      'rolling transfer-function validation')
   varx_validation = bigvar_varx_validate_loo(transfer_series, exogenous, &
      0, 2, lambdas(1:2), bigvar_structure_basic, tolerance=1.0e-7_dp, &
      max_iterations=5000, first_observation=97, last_observation=98)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. all(varx_validation%valid), &
      'leave-one-out transfer-function validation')
   reselection = bigvar_varx_reselect(transfer_series, exogenous, 0, 2, &
      lambdas(1:2), bigvar_structure_basic, 94, 97, 98, 1, &
      tolerance=1.0e-7_dp, max_iterations=5000)
   call check((reselection%info == 0 .or. reselection%info == 4) .and. &
      all(reselection%valid), 'rolling transfer-function reselection')

   lambdas = bigvar_varx_lambda_grid(varx_series, exogenous, 1, 2, &
      bigvar_structure_basic, 100.0_dp, 3, include_intercept=.false.)
   varx_zero = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      lambdas(1), bigvar_structure_basic, 1.0e-8_dp, &
      include_intercept=.false.)
   call check(varx_zero%info == 0 .and. &
      .not. varx_zero%include_intercept .and. varx_zero%nonzero == 0 .and. &
      maxval(abs(varx_zero%intercept)) < 1.0e-14_dp, &
      'structured VARX without intercept')
   varx_path = bigvar_structured_varx_path(varx_series, exogenous, 1, 2, &
      lambdas, bigvar_structure_basic, 1.0e-8_dp, &
      include_intercept=.false.)
   call check(varx_path%info == 0 .and. &
      .not. varx_path%include_intercept .and. &
      maxval(abs(varx_path%intercept)) < 1.0e-14_dp, &
      'structured VARX path without intercept')
   relaxed_varx = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      lambdas(2), bigvar_structure_basic, 1.0e-8_dp, &
      refit_fraction=0.5_dp, include_intercept=.false.)
   call check(relaxed_varx%info == 0 .and. relaxed_varx%relaxed .and. &
      .not. relaxed_varx%include_intercept .and. &
      maxval(abs(relaxed_varx%intercept)) < 1.0e-14_dp, &
      'relaxed VARX refit without intercept')
   varx_validation = bigvar_varx_validate(varx_series, exogenous, 1, 2, &
      lambdas(1:2), bigvar_structure_basic, 94, 96, 1, &
      tolerance=1.0e-7_dp, max_iterations=5000, &
      include_intercept=.false.)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. all(varx_validation%valid), &
      'rolling VARX validation without intercept')
   varx_validation = bigvar_varx_validate_loo(varx_series, exogenous, 1, 2, &
      lambdas(1:2), bigvar_structure_basic, tolerance=1.0e-7_dp, &
      max_iterations=5000, first_observation=97, last_observation=98, &
      include_intercept=.false.)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. all(varx_validation%valid), &
      'leave-one-out VARX validation without intercept')
   ls_varx = bigvar_least_squares_varx(varx_series, exogenous, 1, 2, &
      information_criterion=bigvar_ic_bic, include_intercept=.false.)
   call check(ls_varx%info == 0 .and. .not. ls_varx%include_intercept .and. &
      maxval(abs(ls_varx%intercept)) < 1.0e-14_dp, &
      'least-squares VARX without intercept')
   ic_search = bigvar_varx_ic_select(varx_series, exogenous, 2, 2, &
      bigvar_ic_bic, include_intercept=.false.)
   call check(ic_search%info == 0 .and. &
      .not. ic_search%fit%include_intercept .and. &
      maxval(abs(ic_search%fit%intercept)) < 1.0e-14_dp, &
      'information-criterion VARX selection without intercept')

   dual_lambdas = bigvar_varx_lambda_alpha_grid(varx_series, exogenous, 1, 2, &
      bigvar_structure_sparse_lag, 100.0_dp, 3, validation_alphas)
   linear_lambda_matrix = bigvar_varx_lambda_alpha_grid(varx_series, &
      exogenous, 1, 2, bigvar_structure_sparse_lag, 100.0_dp, 3, &
      validation_alphas, linear=.true.)
   call check(all(shape(dual_lambdas) == [3, 2]) .and. &
      maxval(abs(dual_lambdas(:, 1) - dual_lambdas(:, 2))) > 1.0e-10_dp, &
      'alpha-specific VARX lambda grids')
   call check(maxval(abs(linear_lambda_matrix(1, :) - &
      dual_lambdas(1, :))) < 1.0e-14_dp .and. &
      maxval(abs(linear_lambda_matrix(3, :) - &
      dual_lambdas(3, :))) < 1.0e-14_dp .and. &
      maxval(abs(linear_lambda_matrix(2, :) - 0.5_dp*( &
      linear_lambda_matrix(1, :) + linear_lambda_matrix(3, :)))) < &
      1.0e-14_dp, 'linearly spaced dual VARX grids')
   varx_dual_path = bigvar_structured_varx_dual_path(varx_series, exogenous, &
      1, 2, dual_lambdas, validation_alphas, bigvar_structure_sparse_lag, &
      1.0e-8_dp)
   call check(varx_dual_path%info == 0 .and. &
      all(shape(varx_dual_path%phi) == [2, 2, 3, 2]) .and. &
      all(shape(varx_dual_path%beta) == [2, 4, 3, 2]) .and. &
      all(varx_dual_path%nonzero(1, :) == 0), 'dual VARX path layout')
   varx_validation = bigvar_varx_validate_dual(varx_series, exogenous, 1, 2, &
      dual_lambdas, validation_alphas, bigvar_structure_sparse_lag, 94, 96, 1, &
      tolerance=1.0e-7_dp, max_iterations=5000)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. &
      all(shape(varx_validation%mean_loss_surface) == [3, 2]), &
      'dual rolling VARX validation surface')
   varx_validation = bigvar_varx_validate_loo_dual(varx_series, exogenous, &
      1, 2, dual_lambdas, validation_alphas, bigvar_structure_sparse_lag, &
      tolerance=1.0e-7_dp, max_iterations=5000, first_observation=97, &
      last_observation=98)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. all(varx_validation%valid), &
      'dual leave-one-out VARX validation')
   reselection = bigvar_varx_reselect_dual(varx_series, exogenous, 1, 2, &
      dual_lambdas, validation_alphas, bigvar_structure_sparse_lag, 94, 97, &
      98, 1, selection_window=3, tolerance=1.0e-7_dp, max_iterations=5000)
   call check((reselection%info == 0 .or. reselection%info == 4) .and. &
      all(reselection%valid) .and. all(reselection%selected_lambda > 0.0_dp) &
      .and. all(reselection%selected_alpha_index >= 1), &
      'dual VARX rolling reselection')
   lambdas = bigvar_varx_lambda_grid(varx_series, exogenous, 1, 2, &
      bigvar_structure_basic_en, 100.0_dp, 3, alpha=0.4_dp)
   linear_lambdas = bigvar_varx_lambda_grid(varx_series, exogenous, 1, 2, &
      bigvar_structure_basic_en, 100.0_dp, 3, alpha=0.4_dp, linear=.true.)
   call check(abs(linear_lambdas(1) - lambdas(1)) < 1.0e-14_dp .and. &
      abs(linear_lambdas(3) - lambdas(3)) < 1.0e-14_dp .and. &
      abs(linear_lambdas(2) - 0.5_dp*(linear_lambdas(1) + &
      linear_lambdas(3))) < 1.0e-14_dp, 'linearly spaced VARX lambda grid')
   reselection = bigvar_varx_reselect(varx_series, exogenous, 1, 2, &
      lambdas(1:2), bigvar_structure_basic_en, 94, 97, 98, 1, &
      alphas=[0.4_dp], selection_window=3, tolerance=1.0e-7_dp, &
      max_iterations=5000)
   call check((reselection%info == 0 .or. reselection%info == 4) .and. &
      all(reselection%valid) .and. &
      maxval(abs(reselection%selected_alpha - 0.4_dp)) < 1.0e-14_dp, &
      'ordinary VARX rolling reselection')

   lambdas = bigvar_varx_lambda_grid(varx_series, exogenous, 1, 2, &
      bigvar_structure_basic, 100.0_dp, 3, &
      minnesota_target=minnesota_target)
   varx_zero = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      lambdas(1), bigvar_structure_basic, 1.0e-8_dp, &
      minnesota_target=minnesota_target)
   call check(varx_zero%info == 0 .and. varx_zero%minnesota .and. &
      varx_zero%nonzero == 0 .and. &
      maxval(abs(varx_zero%intercept)) < 1.0e-14_dp .and. &
      abs(varx_zero%phi(1, 1) - 1.0_dp) < 1.0e-12_dp .and. &
      abs(varx_zero%phi(2, 2) - 0.8_dp) < 1.0e-12_dp .and. &
      maxval(abs(varx_zero%beta)) < 1.0e-12_dp, &
      'Minnesota zero-deviation VARX target')
   varx_path = bigvar_structured_varx_path(varx_series, exogenous, 1, 2, &
      lambdas, bigvar_structure_basic, 1.0e-8_dp, &
      minnesota_target=minnesota_target)
   call check(varx_path%info == 0 .and. varx_path%minnesota .and. &
      varx_path%nonzero(1) == 0, 'Minnesota VARX path')
   relaxed_varx = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      lambdas(2), bigvar_structure_basic, 1.0e-8_dp, &
      refit_fraction=0.5_dp, minnesota_target=minnesota_target)
   call check(relaxed_varx%info == 0 .and. relaxed_varx%minnesota .and. &
      relaxed_varx%relaxed .and. &
      maxval(abs(relaxed_varx%intercept)) < 1.0e-14_dp, &
      'Minnesota relaxed VARX refit')
   varx_validation = bigvar_varx_validate(varx_series, exogenous, 1, 2, &
      lambdas(1:2), bigvar_structure_basic, 94, 96, 1, &
      tolerance=1.0e-7_dp, max_iterations=5000, &
      minnesota_target=minnesota_target)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. all(varx_validation%valid), &
      'Minnesota rolling VARX validation')
   varx_validation = bigvar_varx_validate_loo(varx_series, exogenous, 1, 2, &
      lambdas(1:2), bigvar_structure_basic, tolerance=1.0e-7_dp, &
      max_iterations=5000, first_observation=97, last_observation=98, &
      minnesota_target=minnesota_target)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. all(varx_validation%valid), &
      'Minnesota leave-one-out VARX validation')
   lambdas = bigvar_varx_lambda_grid(contemporaneous_series, exogenous, 1, 1, &
      bigvar_structure_basic, 100.0_dp, 3, contemporaneous=.true., &
      minnesota_target=minnesota_target)
   contemporaneous_varx = bigvar_structured_varx(contemporaneous_series, &
      exogenous, 1, 1, lambdas(1), bigvar_structure_basic, 1.0e-8_dp, &
      contemporaneous=.true., minnesota_target=minnesota_target)
   call check(contemporaneous_varx%info == 0 .and. &
      contemporaneous_varx%minnesota .and. &
      contemporaneous_varx%contemporaneous .and. &
      contemporaneous_varx%nonzero == 0 .and. &
      maxval(abs(contemporaneous_varx%beta)) < 1.0e-12_dp, &
      'contemporaneous Minnesota VARX target')

   ls_varx = bigvar_least_squares_varx(varx_series, exogenous, 1, 2, &
      information_criterion=bigvar_ic_bic)
   call check(ls_varx%info == 0 .and. &
      all(shape(ls_varx%phi) == [2, 2]) .and. &
      all(shape(ls_varx%beta) == [2, 4]) .and. &
      all(shape(ls_varx%residuals) == [98, 2]) .and. &
      ls_varx%criterion < huge(1.0_dp), &
      'least-squares VARX dimensions and BIC')
   do time = 3, 100
      prediction = ls_varx%intercept + &
         matmul(ls_varx%phi, varx_series(time - 1, :)) + &
         matmul(ls_varx%beta(:, 1:2), exogenous(time - 1, :)) + &
         matmul(ls_varx%beta(:, 3:4), exogenous(time - 2, :))
      call check(maxval(abs(varx_series(time, :) - prediction - &
         ls_varx%residuals(time - 2, :))) < 1.0e-12_dp, &
         'least-squares VARX residual identity')
   end do
   call check(maxval(abs(ls_varx%innovation_covariance - &
      matmul(transpose(ls_varx%residuals), ls_varx%residuals)/98.0_dp)) < &
      1.0e-14_dp, 'least-squares VARX innovation covariance')
   ic_search = bigvar_varx_ic_select(varx_series, exogenous, 3, 3, &
      bigvar_ic_bic)
   call check(ic_search%info == 0 .and. &
      all(shape(ic_search%criterion) == [4, 4]) .and. &
      abs(ic_search%fit%criterion - minval(ic_search%criterion)) < &
      1.0e-14_dp .and. ic_search%selected_ar_order >= 0 .and. &
      ic_search%selected_ar_order <= 3 .and. &
      ic_search%selected_exogenous_order >= 0 .and. &
      ic_search%selected_exogenous_order <= 3, &
      'joint BIC lag-order selection')
   direct_ls_varx = bigvar_least_squares_varx(varx_series, exogenous, &
      1, 2, 3, bigvar_ic_aic)
   call check(direct_ls_varx%info == 0 .and. direct_ls_varx%direct .and. &
      direct_ls_varx%forecast_horizon == 3 .and. &
      all(shape(direct_ls_varx%residuals) == [96, 2]), &
      'direct least-squares VARX alignment dimensions')
   prediction = direct_ls_varx%intercept + &
      matmul(direct_ls_varx%phi, varx_series(2, :)) + &
      matmul(direct_ls_varx%beta(:, 1:2), exogenous(2, :)) + &
      matmul(direct_ls_varx%beta(:, 3:4), exogenous(1, :))
   call check(maxval(abs(varx_series(5, :) - prediction - &
      direct_ls_varx%residuals(1, :))) < 1.0e-12_dp, &
      'direct least-squares VARX residual identity')
   forecast = bigvar_ls_varx_forecast(direct_ls_varx, varx_series, &
      exogenous, 3)
   prediction = direct_ls_varx%intercept + &
      matmul(direct_ls_varx%phi, varx_series(100, :)) + &
      matmul(direct_ls_varx%beta(:, 1:2), exogenous(100, :)) + &
      matmul(direct_ls_varx%beta(:, 3:4), exogenous(99, :))
   call check(forecast%info == 0 .and. &
      maxval(abs(forecast%mean(3, :) - prediction)) < 1.0e-14_dp, &
      'direct least-squares VARX forecast')
   ic_evaluation = bigvar_varx_ic_evaluate(varx_series, exogenous, 2, 2, &
      bigvar_ic_bic, 94, 96, 2, iterated=.false.)
   call check(ic_evaluation%info == 0 .and. &
      all(shape(ic_evaluation%forecasts) == [3, 2]) .and. &
      all(ic_evaluation%valid) .and. &
      all(ic_evaluation%ar_order >= 0) .and. &
      all(ic_evaluation%exogenous_order >= 0), &
      'rolling direct BIC VARX benchmark')
   call check(abs(ic_evaluation%loss(1) - sum((varx_series(96, :) - &
      ic_evaluation%forecasts(1, :))**2)) < 1.0e-14_dp, &
      'rolling BIC VARX loss alignment')
   ic_evaluation = bigvar_varx_ic_evaluate(varx_series, exogenous, 2, 0, &
      bigvar_ic_aic, 94, 96, 2, iterated=.true.)
   call check(ic_evaluation%info == 0 .and. all(ic_evaluation%valid) .and. &
      all(ic_evaluation%exogenous_order == 0), &
      'rolling iterated AIC VAR benchmark')

   varx_lag = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      0.002_dp, bigvar_structure_lag, 1.0e-8_dp)
   call check(varx_lag%info == 0 .and. varx_lag%converged .and. &
      all(shape(varx_lag%phi) == [2, 2]) .and. &
      all(shape(varx_lag%beta) == [2, 4]) .and. &
      all(shape(varx_lag%residuals) == [98, 2]), &
      'Lag structured VARX convergence and dimensions')
   relaxed_varx = bigvar_relaxed_varx(varx_series, exogenous, varx_lag, &
      0.5_dp)
   call check(relaxed_varx%info == 0 .and. relaxed_varx%relaxed .and. &
      abs(relaxed_varx%refit_fraction - 0.5_dp) < 1.0e-14_dp .and. &
      all((abs(relaxed_varx%phi) > 0.0_dp) .eqv. &
      (abs(varx_lag%phi) > 1.0e-8_dp)) .and. &
      all((abs(relaxed_varx%beta) > 0.0_dp) .eqv. &
      (abs(varx_lag%beta) > 1.0e-8_dp)), &
      'relaxed VARX support and blending')
   varx_path = bigvar_structured_varx_path(varx_series, exogenous, 1, 2, &
      [0.004_dp, 0.002_dp], bigvar_structure_lag, 1.0e-7_dp, &
      refit_fraction=0.5_dp)
   call check(varx_path%info == 0 .and. all(varx_path%converged), &
      'relaxed VARX path')
   do time = 3, 100
      prediction = varx_lag%intercept + &
         matmul(varx_lag%phi, varx_series(time - 1, :)) + &
         matmul(varx_lag%beta(:, 1:2), exogenous(time - 1, :)) + &
         matmul(varx_lag%beta(:, 3:4), exogenous(time - 2, :))
      call check(maxval(abs(varx_series(time, :) - prediction - &
         varx_lag%residuals(time - 2, :))) < 1.0e-12_dp, &
         'structured VARX residual identity')
   end do

   varx_validation = bigvar_varx_validate(varx_series, exogenous, 1, 2, &
      validation_lambdas, bigvar_structure_lag, 95, 97, 2, &
      bigvar_loss_l2, window_size=45, tolerance=1.0e-7_dp, &
      max_iterations=5000)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. &
      all(shape(varx_validation%loss) == [3, 2]) .and. &
      all(varx_validation%valid), 'rolling VARX validation dimensions')
   call check(abs(varx_validation%loss(1, 2) - sum((varx_series(97, :) - &
      varx_validation%forecasts(1, 2, :))**2)) < 1.0e-14_dp, &
      'VARX validation future-exogenous alignment')
   varx_validation = bigvar_varx_validate_loo(varx_series, exogenous, 1, 2, &
      [0.01_dp, 0.002_dp], bigvar_structure_lag, &
      one_standard_error=.true., tolerance=1.0e-7_dp, &
      max_iterations=5000, horizon=2, first_observation=96, &
      last_observation=98)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. &
      all(shape(varx_validation%loss) == [3, 2]) .and. &
      all(varx_validation%valid), 'leave-one-out VARX validation')
   call check(abs(varx_validation%loss(1, 2) - sum((varx_series(96, :) - &
      varx_validation%forecasts(1, 2, :))**2)) < 1.0e-14_dp, &
      'leave-one-out VARX loss identity')
   varx_own_other = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      0.002_dp, bigvar_structure_own_other, 1.0e-8_dp)
   varx_sparse_lag = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      0.002_dp, bigvar_structure_sparse_lag, 1.0e-8_dp, alpha=0.25_dp)
   varx_sparse_own_other = bigvar_structured_varx(varx_series, exogenous, &
      1, 2, 0.002_dp, bigvar_structure_sparse_own_other, 1.0e-8_dp)
   varx_basic_en = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      0.002_dp, bigvar_structure_basic_en, 1.0e-8_dp, alpha=0.25_dp)
   call check(varx_own_other%info == 0 .and. varx_own_other%converged .and. &
      varx_sparse_lag%info == 0 .and. varx_sparse_lag%converged .and. &
      varx_sparse_own_other%info == 0 .and. &
      varx_sparse_own_other%converged .and. varx_basic_en%info == 0 .and. &
      varx_basic_en%converged, 'convex structured VARX families')
   varx_mcp = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      0.001_dp, bigvar_structure_mcp, 1.0e-8_dp, gamma=3.0_dp)
   varx_scad = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      0.001_dp, bigvar_structure_scad, 1.0e-8_dp, gamma=3.5_dp)
   call check(varx_mcp%info == 0 .and. varx_mcp%converged .and. &
      varx_scad%info == 0 .and. varx_scad%converged, &
      'non-convex structured VARX families')
   varx_efx = bigvar_structured_varx(varx_series, exogenous, 2, 2, &
      0.002_dp, bigvar_structure_efx, 1.0e-4_dp)
   call check(varx_efx%info == 0 .and. varx_efx%converged .and. &
      all(shape(varx_efx%phi) == [2, 4]) .and. &
      all(shape(varx_efx%beta) == [2, 4]) .and. &
      varx_efx%active_groups >= 1, 'EFX structured VARX fit')

   lambdas = bigvar_varx_lambda_grid(varx_series, exogenous, 1, 2, &
      bigvar_structure_sparse_lag, 100.0_dp, 5, 0.25_dp)
   varx_zero = bigvar_structured_varx(varx_series, exogenous, 1, 2, &
      lambdas(1), bigvar_structure_sparse_lag, 1.0e-8_dp, alpha=0.25_dp)
   call check(varx_zero%info == 0 .and. varx_zero%nonzero == 0, &
      'structured VARX lambda bound gives the zero model')
   varx_path = bigvar_structured_varx_path(varx_series, exogenous, 1, 2, &
      lambdas, bigvar_structure_sparse_lag, 1.0e-8_dp, alpha=0.25_dp)
   call check(varx_path%info == 0 .and. all(varx_path%converged) .and. &
      all(shape(varx_path%phi) == [2, 2, 5]) .and. &
      all(shape(varx_path%beta) == [2, 4, 5]) .and. &
      varx_path%nonzero(1) == 0, 'structured VARX warm-started path')
   lambdas = bigvar_varx_lambda_grid(varx_series, exogenous, 2, 2, &
      bigvar_structure_efx, 100.0_dp, 5)
   varx_zero = bigvar_structured_varx(varx_series, exogenous, 2, 2, &
      lambdas(1), bigvar_structure_efx, 1.0e-4_dp)
   call check(varx_zero%info == 0 .and. varx_zero%nonzero == 0, &
      'EFX lambda bound gives the zero model')
   varx_path = bigvar_structured_varx_path(varx_series, exogenous, 2, 2, &
      lambdas, bigvar_structure_efx, 1.0e-4_dp)
   call check(varx_path%info == 0 .and. all(varx_path%converged) .and. &
      varx_path%nonzero(1) == 0 .and. &
      all(shape(varx_path%phi) == [2, 4, 5]), &
      'EFX warm-started path')

   contemporaneous_varx = bigvar_structured_varx( &
      contemporaneous_series, exogenous, 1, 1, 0.001_dp, &
      bigvar_structure_basic, 1.0e-8_dp, contemporaneous=.true.)
   call check(contemporaneous_varx%info == 0 .and. &
      contemporaneous_varx%converged .and. &
      contemporaneous_varx%contemporaneous .and. &
      all(shape(contemporaneous_varx%phi) == [2, 2]) .and. &
      all(shape(contemporaneous_varx%beta) == [2, 4]) .and. &
      all(shape(contemporaneous_varx%residuals) == [99, 2]), &
      'contemporaneous VARX dimensions and metadata')
   do time = 2, 100
      prediction = contemporaneous_varx%intercept + &
         matmul(contemporaneous_varx%phi, &
         contemporaneous_series(time - 1, :)) + &
         matmul(contemporaneous_varx%beta(:, 1:2), &
         exogenous(time, :)) + &
         matmul(contemporaneous_varx%beta(:, 3:4), &
         exogenous(time - 1, :))
      call check(maxval(abs(contemporaneous_series(time, :) - prediction - &
         contemporaneous_varx%residuals(time - 1, :))) < 1.0e-12_dp, &
         'contemporaneous VARX residual identity')
   end do
   current_only_varx = bigvar_structured_varx( &
      contemporaneous_series, exogenous, 1, 0, 0.001_dp, &
      bigvar_structure_basic, 1.0e-8_dp, contemporaneous=.true.)
   call check(current_only_varx%info == 0 .and. &
      current_only_varx%converged .and. &
      all(shape(current_only_varx%beta) == [2, 2]), &
      'current-only VARX permits zero exogenous lag order')
   relaxed_varx = bigvar_relaxed_varx(contemporaneous_series, exogenous, &
      contemporaneous_varx, 0.5_dp)
   call check(relaxed_varx%info == 0 .and. relaxed_varx%relaxed .and. &
      relaxed_varx%contemporaneous .and. &
      all(shape(relaxed_varx%beta) == [2, 4]), &
      'relaxed contemporaneous VARX preserves layout')
   lambdas = bigvar_varx_lambda_grid(contemporaneous_series, exogenous, &
      1, 1, bigvar_structure_basic, 100.0_dp, 4, &
      contemporaneous=.true.)
   varx_zero = bigvar_structured_varx(contemporaneous_series, exogenous, &
      1, 1, lambdas(1), bigvar_structure_basic, 1.0e-8_dp, &
      contemporaneous=.true.)
   call check(size(lambdas) == 4 .and. varx_zero%info == 0 .and. &
      varx_zero%nonzero == 0, &
      'contemporaneous VARX lambda bound gives the zero model')
   varx_path = bigvar_structured_varx_path(contemporaneous_series, &
      exogenous, 1, 1, lambdas, bigvar_structure_basic, 1.0e-8_dp, &
      contemporaneous=.true.)
   call check(varx_path%info == 0 .and. all(varx_path%converged) .and. &
      varx_path%contemporaneous .and. &
      all(shape(varx_path%beta) == [2, 4, 4]), &
      'contemporaneous VARX warm-started path')
   rejected_varx = bigvar_structured_varx(contemporaneous_series, &
      exogenous, 1, 1, 0.001_dp, bigvar_structure_efx, &
      contemporaneous=.true.)
   call check(rejected_varx%info /= 0, &
      'EFX rejects contemporaneous exogenous predictors')
   varx_validation = bigvar_varx_validate(contemporaneous_series, &
      exogenous, 1, 1, [0.01_dp, 0.001_dp], bigvar_structure_basic, &
      95, 97, 2, tolerance=1.0e-7_dp, max_iterations=5000, &
      contemporaneous=.true.)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. all(varx_validation%valid), &
      'rolling contemporaneous VARX validation')
   varx_validation = bigvar_varx_validate_loo(contemporaneous_series, &
      exogenous, 1, 1, [0.01_dp], bigvar_structure_basic, &
      tolerance=1.0e-7_dp, max_iterations=5000, contemporaneous=.true., &
      first_observation=98, last_observation=99)
   call check((varx_validation%info == 0 .or. &
      varx_validation%info == 4) .and. all(varx_validation%valid), &
      'contemporaneous leave-one-out VARX validation')
   evaluation = bigvar_varx_evaluate(contemporaneous_series, exogenous, &
      1, 1, 0.001_dp, bigvar_structure_basic, 96, 97, 2, &
      tolerance=1.0e-7_dp, max_iterations=5000, contemporaneous=.true.)
   call check((evaluation%info == 0 .or. evaluation%info == 4) .and. &
      all(evaluation%valid), 'contemporaneous VARX rolling evaluation')

   exogenous_forecast(:100, :) = exogenous
   do time = 101, 103
      exogenous_forecast(time, 1) = sin(0.17_dp*real(time, dp)) + &
         0.08_dp*cos(0.013_dp*real(time*time, dp))
      exogenous_forecast(time, 2) = cos(0.11_dp*real(time, dp)) - &
         0.06_dp*sin(0.019_dp*real(time*time, dp))
   end do
   forecast = bigvar_varx_forecast(varx_lag, varx_series, &
      exogenous_forecast, 3)
   prediction = varx_lag%intercept + &
      matmul(varx_lag%phi, varx_series(100, :)) + &
      matmul(varx_lag%beta(:, 1:2), exogenous_forecast(100, :)) + &
      matmul(varx_lag%beta(:, 3:4), exogenous_forecast(99, :))
   call check(forecast%info == 0 .and. &
      maxval(abs(forecast%mean(1, :) - prediction)) < 1.0e-14_dp, &
      'structured VARX forecast adapter')
   interval = bigvar_varx_forecast_interval(varx_lag, varx_series, &
      exogenous_forecast, 3)
   innovation_covariance = &
      bigvar_innovation_covariance(varx_lag%residuals)
   manual_covariance = innovation_covariance + matmul(matmul( &
      varx_lag%phi(:, 1:2), innovation_covariance), &
      transpose(varx_lag%phi(:, 1:2)))
   call check(interval%info == 0 .and. &
      maxval(abs(interval%mean - forecast%mean)) < 1.0e-14_dp .and. &
      maxval(abs(interval%covariance(:, :, 1) - &
      innovation_covariance)) < 1.0e-14_dp .and. &
      maxval(abs(interval%covariance(:, :, 2) - manual_covariance)) < &
      1.0e-12_dp, 'conditional VARX forecast covariance')
   forecast = bigvar_varx_forecast(contemporaneous_varx, &
      contemporaneous_series, exogenous_forecast, 3)
   prediction = contemporaneous_varx%intercept + &
      matmul(contemporaneous_varx%phi, contemporaneous_series(100, :)) + &
      matmul(contemporaneous_varx%beta(:, 1:2), &
      exogenous_forecast(101, :)) + &
      matmul(contemporaneous_varx%beta(:, 3:4), &
      exogenous_forecast(100, :))
   call check(forecast%info == 0 .and. &
      maxval(abs(forecast%mean(1, :) - prediction)) < 1.0e-14_dp, &
      'contemporaneous VARX forecast alignment')
   interval = bigvar_varx_forecast_interval(contemporaneous_varx, &
      contemporaneous_series, exogenous_forecast, 3)
   call check(interval%info == 0 .and. &
      maxval(abs(interval%mean - forecast%mean)) < 1.0e-14_dp, &
      'contemporaneous VARX forecast intervals')

   forecast = bigvar_forecast(lag_fit, series, 3)
   prediction = lag_fit%intercept + &
      matmul(lag_fit%phi(:, 1:2), series(100, :)) + &
      matmul(lag_fit%phi(:, 3:4), series(99, :))
   call check(forecast%info == 0 .and. &
      maxval(abs(forecast%mean(1, :) - prediction)) < 1.0e-14_dp, &
      'forecast adapter reuses bigtime')
   innovation_covariance = bigvar_innovation_covariance(lag_fit%residuals)
   forecast_covariance = bigvar_forecast_covariance(lag_fit%phi, 2, &
      innovation_covariance, 3)
   manual_covariance = innovation_covariance + matmul(matmul( &
      lag_fit%phi(:, 1:2), innovation_covariance), &
      transpose(lag_fit%phi(:, 1:2)))
   call check(maxval(abs(forecast_covariance(:, :, 1) - &
      innovation_covariance)) < 1.0e-14_dp .and. &
      maxval(abs(forecast_covariance(:, :, 2) - manual_covariance)) < &
      1.0e-12_dp, 'recursive VAR forecast covariance identities')
   interval = bigvar_var_forecast_interval(lag_fit, series, 3)
   interval_90 = bigvar_var_forecast_interval(lag_fit, series, 3, 0.90_dp)
   call check(interval%info == 0 .and. interval%level == 0.95_dp .and. &
      all(shape(interval%covariance) == [2, 2, 3]) .and. &
      maxval(abs(interval%mean - forecast%mean)) < 1.0e-14_dp .and. &
      maxval(abs(interval%lower + interval%upper - &
      2.0_dp*interval%mean)) < 1.0e-14_dp, &
      'recursive VAR marginal confidence intervals')
   call check(all(interval_90%upper - interval_90%lower < &
      interval%upper - interval%lower), &
      'lower confidence level gives narrower VAR intervals')
   interval = bigvar_var_forecast_interval(direct_fit, series, 3)
   call check(interval%info == 0 .and. &
      maxval(abs(interval%covariance(:, :, 3) - &
      bigvar_innovation_covariance(direct_fit%residuals))) < 1.0e-14_dp .and. &
      maxval(abs(interval%covariance(:, :, :2))) == 0.0_dp, &
      'direct VAR horizon-specific forecast covariance')
   stability = bigvar_stability(lag_fit)
   call check(stability%info == 0 .and. stability%maximum_modulus >= 0.0_dp, &
      'stability adapter reuses bigtime')

   print '(a)', 'BigVAR tests passed'

contains

   subroutine check(condition, message)
      ! Stop the test program when a condition is false.
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_bigvar
