! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Numerical tests for algorithms translated from R bvarsv.
program test_bvarsv
   use kind_mod, only: dp
   use bvarsv_mod, only: bvarsv_state_draw_t, bvarsv_carter_kohn, &
      bvarsv_data_t, bvarsv_prepare, bvarsv_coefficient_update, &
      bvarsv_contemporaneous_matrix, bvarsv_covariance_from_state, &
      bvarsv_volatility_draw_t, bvarsv_log_volatility_update, &
      bvarsv_contemporaneous_draw_t, bvarsv_contemporaneous_update, &
      bvarsv_covariance_draw_t, bvarsv_random_walk_covariance, &
      bvarsv_contemporaneous_covariance_update, bvarsv_prior_t, &
      bvarsv_ols_prior, bvarsv_draws_t, bvarsv_gibbs
   use bvarsv_mod, only: bvarsv_predictive_t, bvarsv_forecast
   use bvarsv_mod, only: bvarsv_irf_t, bvarsv_irf
   use bvarsv_mod, only: bvarsv_predictive_density
   use bvarsv_mod, only: bvarsv_simulation_t, bvarsv_simulate_var1
   use bvarsv_mod, only: bvarsv_predictive_draws_t, bvarsv_predictive_draws
   use bvarsv_mod, only: bvarsv_parameter_path_t, bvarsv_parameter_draws
   use random_mod, only: set_random_seed
   implicit none

   real(dp) :: response(1, 4), design(1, 1, 4), observation_covariance(1, 1, 4)
   real(dp) :: state_covariance(1, 1), initial_mean(1), initial_covariance(1, 1)
   real(dp), allocatable :: contemporaneous(:, :), covariance(:, :)
   type(bvarsv_state_draw_t) :: draw
   type(bvarsv_volatility_draw_t) :: volatility
   real(dp) :: residual(2, 4), current_log_variance(2, 4)
   real(dp) :: volatility_covariance(2, 2), volatility_initial_mean(2)
   real(dp) :: volatility_initial_covariance(2, 2)
   real(dp) :: structural_residual(3, 4), structural_log_variance(3, 4)
   real(dp) :: structural_state_covariance(3, 3), structural_initial_mean(3)
   real(dp) :: structural_initial_covariance(3, 3)
   type(bvarsv_contemporaneous_draw_t) :: structural_draw
   type(bvarsv_covariance_draw_t) :: covariance_draw, block_covariance_draw
   real(dp) :: covariance_prior(3, 3)
   real(dp) :: var_series(2, 7), coefficient_observation_covariance(2, 2, 6)
   real(dp) :: prior_series(2, 15)
   real(dp) :: coefficient_state_covariance(6, 6), coefficient_initial_mean(6)
   real(dp) :: coefficient_initial_covariance(6, 6)
   type(bvarsv_data_t) :: prepared
   type(bvarsv_state_draw_t) :: coefficient_draw
   integer :: time
   type(bvarsv_prior_t) :: prior
   type(bvarsv_draws_t) :: posterior
   type(bvarsv_predictive_t) :: prediction
   type(bvarsv_irf_t) :: irf
   real(dp), allocatable :: density(:), distribution(:)
   real(dp) :: simulation_coefficient(2, 3), simulation_q(6, 6)
   real(dp) :: simulation_s(1, 1), simulation_w(2, 2)
   type(bvarsv_simulation_t) :: simulation
   type(bvarsv_predictive_draws_t) :: predictive_selection
   type(bvarsv_parameter_path_t) :: parameter_selection

   response = reshape([0.1_dp, 0.2_dp, 0.15_dp, 0.25_dp], [1, 4])
   design = 1.0_dp
   observation_covariance = 0.1_dp
   state_covariance = 0.01_dp
   initial_mean = 0.0_dp
   initial_covariance = 1.0_dp
   call set_random_seed(910)
   draw = bvarsv_carter_kohn(response, design, observation_covariance, &
      state_covariance, initial_mean, initial_covariance)
   call assert_true(draw%info == 0 .and. all(shape(draw%state) == [1, 4]) .and. &
      abs(draw%log_likelihood) < huge(1.0_dp), 'Carter-Kohn state draw')
   contemporaneous = bvarsv_contemporaneous_matrix([0.5_dp], 2)
   call assert_true(all(abs(contemporaneous - reshape([1.0_dp, 0.5_dp, &
      0.0_dp, 1.0_dp], [2, 2])) < 1.0e-14_dp), &
      'triangular contemporaneous matrix')
   covariance = bvarsv_covariance_from_state([0.5_dp], [0.0_dp, 0.0_dp])
   call assert_true(abs(covariance(1, 1) - 1.0_dp) < 1.0e-14_dp .and. &
      abs(covariance(1, 2) + 0.5_dp) < 1.0e-14_dp .and. &
      abs(covariance(2, 2) - 1.25_dp) < 1.0e-14_dp, &
      'structural covariance transformation')
   residual = reshape([0.5_dp, -0.2_dp, 0.4_dp, 0.3_dp, 0.6_dp, &
      -0.1_dp, 0.2_dp, 0.5_dp], [2, 4])
   current_log_variance = 0.0_dp
   volatility_covariance = 0.0_dp
   volatility_covariance(1, 1) = 0.02_dp
   volatility_covariance(2, 2) = 0.02_dp
   volatility_initial_mean = 0.0_dp
   volatility_initial_covariance = 0.0_dp
   volatility_initial_covariance(1, 1) = 1.0_dp
   volatility_initial_covariance(2, 2) = 1.0_dp
   call set_random_seed(911)
   volatility = bvarsv_log_volatility_update(residual, current_log_variance, &
      volatility_covariance, volatility_initial_mean, &
      volatility_initial_covariance)
   call assert_true(volatility%info == 0 .and. &
      all(shape(volatility%log_variance) == [2, 4]) .and. &
      all(volatility%component >= 1) .and. all(volatility%component <= 7) .and. &
      all(volatility%standard_deviation > 0.0_dp), &
      'KSC log-volatility update')
   structural_residual = reshape([0.2_dp, -0.1_dp, 0.3_dp, &
      0.1_dp, 0.2_dp, -0.2_dp, 0.3_dp, 0.1_dp, 0.2_dp, &
      -0.1_dp, 0.2_dp, 0.4_dp], [3, 4])
   structural_log_variance = 0.0_dp
   structural_state_covariance = 0.0_dp
   structural_initial_covariance = 0.0_dp
   structural_state_covariance(1, 1) = 0.01_dp
   structural_state_covariance(2, 2) = 0.01_dp
   structural_state_covariance(3, 3) = 0.01_dp
   structural_initial_covariance(1, 1) = 1.0_dp
   structural_initial_covariance(2, 2) = 1.0_dp
   structural_initial_covariance(3, 3) = 1.0_dp
   structural_initial_mean = 0.0_dp
   call set_random_seed(912)
   structural_draw = bvarsv_contemporaneous_update(structural_residual, &
      structural_log_variance, structural_state_covariance, &
      structural_initial_mean, structural_initial_covariance)
   call assert_true(structural_draw%info == 0 .and. &
      all(shape(structural_draw%alpha) == [3, 4]), &
      'time-varying contemporaneous-state update')
   covariance_prior = 0.0_dp
   covariance_prior(1, 1) = 0.1_dp
   covariance_prior(2, 2) = 0.1_dp
   covariance_prior(3, 3) = 0.1_dp
   call set_random_seed(913)
   covariance_draw = bvarsv_random_walk_covariance(structural_draw%alpha, &
      covariance_prior, 5.0_dp, degrees_offset=1)
   call assert_true(covariance_draw%info == 0 .and. &
      all(shape(covariance_draw%covariance) == [3, 3]) .and. &
      all(diagonal_values(covariance_draw%covariance) > 0.0_dp), &
      'random-walk innovation covariance update')
   call set_random_seed(914)
   block_covariance_draw = bvarsv_contemporaneous_covariance_update( &
      structural_draw%alpha, covariance_prior, [3.0_dp, 4.0_dp])
   call assert_true(block_covariance_draw%info == 0 .and. &
      all(abs(block_covariance_draw%covariance(1, 2:3)) < 1.0e-14_dp), &
      'block contemporaneous covariance update')
   do time = 1, 7
      var_series(1, time) = 0.1_dp*real(time, dp)
      var_series(2, time) = 0.2_dp + 0.5_dp*var_series(1, time) + &
         0.03_dp*sin(real(time, dp))
   end do
   prepared = bvarsv_prepare(var_series, 1)
   call assert_true(prepared%info == 0 .and. prepared%observations == 6 .and. &
      prepared%states == 6, 'TVP-VAR design preparation')
   coefficient_observation_covariance = 0.0_dp
   coefficient_observation_covariance(1, 1, :) = 0.1_dp
   coefficient_observation_covariance(2, 2, :) = 0.1_dp
   coefficient_state_covariance = 0.0_dp
   coefficient_initial_covariance = 0.0_dp
   do time = 1, 6
      coefficient_state_covariance(time, time) = 0.001_dp
      coefficient_initial_covariance(time, time) = 1.0_dp
   end do
   coefficient_initial_mean = 0.0_dp
   call set_random_seed(915)
   coefficient_draw = bvarsv_coefficient_update(prepared, &
      coefficient_observation_covariance(:, :, :prepared%observations), &
      coefficient_state_covariance, coefficient_initial_mean, &
      coefficient_initial_covariance)
   call assert_true(coefficient_draw%info == 0 .and. &
      all(shape(coefficient_draw%state) == [6, 6]), &
      'time-varying VAR coefficient update')
   call set_random_seed(916)
   do time = 1, 15
      prior_series(1, time) = 0.1_dp*real(time, dp) + &
         0.05_dp*sin(real(time, dp))
      prior_series(2, time) = 0.3_dp + 0.4_dp*prior_series(1, time) + &
         0.04_dp*cos(0.7_dp*real(time, dp))
   end do
   prior = bvarsv_ols_prior(prior_series, 1, covariance_draws=30)
   call assert_true(prior%info == 0 .and. size(prior%coefficient_mean) == 6 .and. &
      all(shape(prior%coefficient_covariance) == [6, 6]) .and. &
      size(prior%contemporaneous_mean) == 1 .and. &
      all(shape(prior%contemporaneous_covariance) == [1, 1]) .and. &
      size(prior%log_variance_mean) == 2, 'OLS prior construction')
   call set_random_seed(917)
   posterior = bvarsv_gibbs(prepared, prior, 1)
   call assert_true(posterior%info == 0 .and. posterior%retained_draws == 1 .and. &
      all(shape(posterior%coefficient) == [6, 6, 1]) .and. &
      all(shape(posterior%covariance) == [2, 2, 6, 1]), &
      'complete Primiceri Gibbs sampler')
   call set_random_seed(918)
   prediction = bvarsv_forecast(posterior, prepared, 2, parameter_drift=.true.)
   call assert_true(prediction%info == 0 .and. &
      all(shape(prediction%path) == [2, 2, 1]) .and. &
      all(prediction%covariance(1, 1, :, :) > 0.0_dp), &
      'parameter-drift posterior forecast')
   call set_random_seed(919)
   prediction = bvarsv_forecast(posterior, prepared, 2, parameter_drift=.false.)
   call assert_true(prediction%info == 0 .and. &
      all(shape(prediction%mean) == [2, 2, 1]), &
      'fixed-parameter posterior forecast')
   irf = bvarsv_irf(posterior, prepared, 2, scenario=1)
   call assert_true(irf%info == 0 .and. &
      all(shape(irf%response) == [2, 2, 3, 1]) .and. &
      all(abs(irf%response(:, :, 1, 1) - &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])) < 1.0e-14_dp), &
      'unorthogonalized time-varying impulse responses')
   irf = bvarsv_irf(posterior, prepared, 1, scenario=2)
   call assert_true(irf%info == 0 .and. &
      irf%response(1, 1, 1, 1) > 0.0_dp .and. &
      irf%response(2, 2, 1, 1) > 0.0_dp, &
      'Cholesky time-varying impulse responses')
   density = bvarsv_predictive_density(prediction, [0.0_dp, 1.0_dp], 1, 1)
   distribution = bvarsv_predictive_density(prediction, &
      [0.0_dp, 1.0_dp], 1, 1, cumulative=.true.)
   call assert_true(size(density) == 2 .and. all(density > 0.0_dp) .and. &
      all(distribution >= 0.0_dp) .and. all(distribution <= 1.0_dp) .and. &
      distribution(2) >= distribution(1), 'Gaussian-mixture predictive density')
   predictive_selection = bvarsv_predictive_draws(prediction, 2, 1)
   call assert_true(predictive_selection%info == 0 .and. &
      size(predictive_selection%observation) == 1 .and. &
      abs(predictive_selection%observation(1) - prediction%path(2, 1, 1)) < &
      1.0e-14_dp .and. predictive_selection%variance(1) > 0.0_dp, &
      'predictive draw extraction')
   parameter_selection = bvarsv_parameter_draws(posterior, prepared, &
      'intercept', row=2)
   call assert_true(parameter_selection%info == 0 .and. &
      all(shape(parameter_selection%value) == [1, 6]) .and. &
      all(abs(parameter_selection%value(1, :) - &
      posterior%coefficient(2, :, 1)) < 1.0e-14_dp), &
      'intercept draw extraction')
   parameter_selection = bvarsv_parameter_draws(posterior, prepared, 'lag', &
      row=1, column=2, lag=1)
   call assert_true(parameter_selection%info == 0 .and. &
      all(abs(parameter_selection%value(1, :) - &
      posterior%coefficient(4, :, 1)) < 1.0e-14_dp), &
      'lag-coefficient draw extraction')
   parameter_selection = bvarsv_parameter_draws(posterior, prepared, 'vcv', &
      row=1, column=2)
   call assert_true(parameter_selection%info == 0 .and. &
      all(abs(parameter_selection%value(1, :) - &
      posterior%covariance(1, 2, :, 1)) < 1.0e-14_dp), &
      'covariance draw extraction')
   simulation_coefficient = reshape([0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
      0.0_dp, 0.5_dp], [2, 3])
   simulation_q = 0.0_dp
   simulation_s = 0.0_dp
   simulation_w = 0.0_dp
   call set_random_seed(917)
   simulation = bvarsv_simulate_var1(simulation_coefficient, [0.0_dp], &
      [0.0_dp, 0.0_dp], simulation_q, simulation_s, simulation_w, 5, burnin=2)
   call assert_true(simulation%info == 0 .and. &
      all(shape(simulation%series) == [2, 5]) .and. &
      all(shape(simulation%coefficient) == [2, 3, 5]) .and. &
      all(abs(simulation%coefficient(:, :, 5) - simulation_coefficient) < &
      1.0e-14_dp) .and. &
      all(abs(simulation%covariance(:, :, 5) - &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])) < 1.0e-14_dp), &
      'TVP-VAR(1) stochastic-volatility simulation')
   print '(a)', 'bvarsv tests passed'
contains
   subroutine assert_true(condition, message)
      !! Stop the test program when a logical assertion fails.
      logical, intent(in) :: condition !! Assertion condition.
      character(len=*), intent(in) :: message !! Failure description.
      if (.not. condition) error stop 'FAIL: '//trim(message)
   end subroutine assert_true

   pure function diagonal_values(matrix) result(values)
      !! Return the diagonal of a square test matrix.
      real(dp), intent(in) :: matrix(:, :) !! Square input matrix.
      real(dp) :: values(size(matrix, 1))
      integer :: item
      do item = 1, size(values)
         values(item) = matrix(item, item)
      end do
   end function diagonal_values
end program test_bvarsv
