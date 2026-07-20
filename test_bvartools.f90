! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Numerical tests for algorithms translated from R bvartools.
program test_bvartools
   use kind_mod, only: dp
   use bvartools_mod
   use random_mod, only: set_random_seed
   implicit none

   real(dp) :: covariance(2, 2), y(2, 5), x(2, 5), precision(2, 2)
   real(dp) :: prior_mean(4), prior_precision(4, 4), normals(4), uniforms(4)
   real(dp) :: errors(2, 5), states(2, 5), initial_state(2)
   real(dp) :: z(10, 4), time_precision(2, 2, 5), draw(4), gamma_draw(2)
   type(bvartools_minnesota_prior_t) :: minnesota
   type(bvartools_loglik_normal_t) :: normal_loglik, tvp_normal_loglik
   type(bvartools_kalman_dk_t) :: kalman_draw
   type(bvartools_inclusion_prior_t) :: inclusion
   type(bvartools_ssvs_result_t) :: ssvs
   type(bvartools_ssvs_prior_t) :: fixed_ssvs_prior, automatic_ssvs_prior
   type(bvartools_normal_posterior_t) :: posterior, sur_posterior
   type(bvartools_gamma_posterior_t) :: measurement, state
   type(bvartools_covariance_data_t) :: covariance_data
   type(bvartools_normal_posterior_t) :: covariance_posterior, tvp_posterior
   type(bvartools_bvs_result_t) :: bvs
   type(bvartools_bvs_result_t) :: tvp_bvs, covariance_bvs
   type(bvartools_bvar_draws_t) :: bvar_draws, diagonal_bvar_draws
   type(bvartools_bvar_draws_t) :: structural_bvar_draws
   type(bvartools_bvar_prior_bundle_t) :: bvar_prior_bundle
   type(bvartools_bvar_draws_t) :: bundled_bvar_draws
   type(bvartools_bvar_fit_t) :: bvar_fit
   type(bvartools_cointegration_draw_t) :: cointegration, cointegration_x
   type(bvartools_cointegration_draw_t) :: cointegration_sur, cointegration_sur_x
   type(bvartools_bvec_draws_t) :: bvec_draws, diagonal_bvec_draws
   type(bvartools_bvec_draws_t) :: structural_bvec_draws
   type(bvartools_bvec_prior_bundle_t) :: bvec_prior_bundle
   type(bvartools_bvec_draws_t) :: bundled_bvec_draws
   type(bvartools_bvec_prior_bundle_t) :: rank_zero_prior_bundle
   type(bvartools_bvec_fit_t) :: bvec_fit, rank_zero_fit
   type(bvartools_tvp_bvar_prior_bundle_t) :: tvp_bvar_prior, tvp_bvar_sv_prior
   type(bvartools_tvp_bvar_fit_t) :: tvp_bvar_fit, tvp_bvar_sv_fit
   type(bvartools_structural_tvp_bvar_prior_t) :: structural_tvp_prior
   type(bvartools_structural_tvp_bvar_prior_t) :: structural_constant_prior
   type(bvartools_structural_tvp_bvar_fit_t) :: structural_tvp_fit
   type(bvartools_structural_tvp_bvar_fit_t) :: structural_constant_fit
   type(bvartools_tvp_bvec_prior_bundle_t) :: tvp_bvec_prior, tvp_bvec_sv_prior
   type(bvartools_tvp_bvec_fit_t) :: tvp_bvec_fit, tvp_bvec_sv_fit
   type(bvartools_tvp_bvec_fit_t) :: structural_tvp_bvec_fit
   type(bvartools_tvp_bvar_draws_t) :: tvp_bvar_draws, diagonal_tvp_bvar_draws
   type(bvartools_tvp_covariance_draws_t) :: tvp_covariance_draws
   type(bvartools_joint_tvp_bvar_draws_t) :: joint_tvp_draws
   type(bvartools_joint_tvp_bvar_draws_t) :: joint_constant_draws
   type(bvartools_stochastic_volatility_t) :: volatility_draw, ocsn_volatility_draw
   type(bvartools_tvp_bvar_draws_t) :: volatility_tvp_bvar_draws
   type(bvartools_tvp_bvec_draws_t) :: tvp_bvec_draws, tvp_bvec_x_draws
   type(bvartools_tvp_bvec_draws_t) :: tvp_bvec_covariance_draws
   type(bvartools_tvp_bvec_draws_t) :: tvp_bvec_sv_draws
   type(bvartools_normal_posterior_t) :: dfm_factor_posterior
   type(bvartools_dfm_draws_t) :: dfm_draws
   type(bvartools_dfm_prior_t) :: dfm_prior
   type(bvartools_dfm_grid_draws_t) :: dfm_grid_draws
   type(bvartools_model_comparison_t) :: model_comparison, tvp_model_comparison
   type(bvartools_model_likelihood_data_t) :: comparison_data(2)
   type(bvartools_model_comparison_set_t) :: comparison_set
   type(bvartools_normal_posterior_t) :: random_walk_posterior
   type(bvartools_bvar_draws_t) :: predictive_draws, structural_draws
   type(bvartools_tvp_bvar_draws_t) :: predictive_tvp_draws
   type(bvartools_predictive_t) :: prediction, structural_prediction, tvp_prediction
   type(bvartools_bvar_draws_t) :: irf_draws
   type(bvartools_tvp_bvar_draws_t) :: irf_tvp_draws
   type(bvartools_irf_t) :: feir, oir, gir, sir, sgir, cumulative_feir, tvp_irf
   type(bvartools_fevd_t) :: oir_fevd, gir_fevd, normalized_gir_fevd
   type(bvartools_fevd_t) :: sir_fevd, sgir_fevd, tvp_fevd
   type(bvartools_bvec_draws_t) :: conversion_draws
   type(bvartools_tvp_bvec_draws_t) :: conversion_tvp_draws
   type(bvartools_level_var_draws_t) :: level_draws
   type(bvartools_tvp_level_var_draws_t) :: tvp_level_draws
   type(bvartools_var_data_t) :: prepared_var
   type(bvartools_vecm_data_t) :: prepared_vecm
   type(bvartools_dfm_data_t) :: prepared_dfm
   type(bvartools_dfm_data_t) :: dfm_grid_data
   type(bvartools_var_data_t) :: prior_var_data
   type(bvartools_var_data_t) :: structural_var_data
   type(bvartools_vecm_data_t) :: prior_vecm_data
   type(bvartools_vecm_data_t) :: structural_vecm_data
   real(dp) :: covariance_prior_mean(1), covariance_prior_precision(1, 1)
   real(dp) :: covariance_state_precision(1, 1, 5), psi(5)
   real(dp) :: state_precision4(4, 4, 5)
   real(dp), allocatable :: covariance_matrix(:, :), expanded(:, :)
   logical :: initial_inclusion(4)
   real(dp) :: initial_beta(2, 1), cointegration_precision(2, 2)
   real(dp) :: loading_normals(2), loading_gamma_normals(6), beta_normals(2)
   real(dp) :: initial_alpha(2, 1), initial_gamma(2, 2)
   real(dp) :: state_prior6(6), state_prior_precision6(6, 6)
   integer :: t
   real(dp) :: predictive_normals(1, 2, 2), predictive_history(1, 1)
   real(dp) :: structural_normals(1, 1, 2), structural_impact(1, 1, 2)
   real(dp) :: structural_fit_normals(2, 1, 2), structural_fit_history(2, 1)
   real(dp) :: structural_fit_future(1, 1)
   real(dp) :: future_constant(1, 1)
   real(dp) :: irf_impact(2, 2, 2)
   real(dp) :: conversion_pi(2, 3), conversion_gamma(2, 3)
   real(dp) :: exogenous_pi(2, 1, 2), exogenous_difference(2, 1, 2, 2)
   real(dp), allocatable :: level_ar(:, :, :), level_exogenous(:, :, :)
   real(dp), allocatable :: reconstructed(:, :)
   real(dp), allocatable :: zero_gamma(:, :, :)
   real(dp) :: raw_series(2, 8), raw_exogenous(1, 8)
   real(dp) :: comparison_residual(1, 2, 2)
   real(dp) :: comparison_covariance(1, 1, 2)
   real(dp) :: comparison_tvp_covariance(1, 1, 2, 2)
   real(dp) :: expected_log_likelihood
   real(dp) :: ssvs_y(2, 6), ssvs_x(2, 6)
   real(dp) :: prior_series(2, 10)
   real(dp) :: likelihood_residual(1, 2), likelihood_covariance(2, 1)
   real(dp) :: kalman_y(1, 3), kalman_z(3, 1), kalman_covariance(3, 1)
   real(dp) :: kalman_state_covariance(3, 1), kalman_transition(3, 1)

   covariance = reshape([1.0_dp, 0.2_dp, 0.2_dp, 0.8_dp], [2, 2])
   likelihood_residual = reshape([0.0_dp, 1.0_dp], [1, 2])
   normal_loglik = bvartools_loglik_normal(likelihood_residual, &
      reshape([1.0_dp], [1, 1]))
   call assert_true(normal_loglik%info == 0, &
      'constant Gaussian log likelihood status')
   call assert_close(normal_loglik%value(1), &
      -0.5_dp*log(2.0_dp*acos(-1.0_dp)), 1.0e-14_dp, &
      'constant Gaussian log likelihood value')
   likelihood_covariance(:, 1) = [1.0_dp, 4.0_dp]
   tvp_normal_loglik = bvartools_loglik_normal(likelihood_residual, &
      likelihood_covariance)
   call assert_close(tvp_normal_loglik%value(2), &
      -0.5_dp*(log(2.0_dp*acos(-1.0_dp)) + log(4.0_dp) + 0.25_dp), &
      1.0e-14_dp, 'time-varying Gaussian log likelihood value')

   kalman_y = reshape([0.2_dp, -0.1_dp, 0.3_dp], [1, 3])
   kalman_z = 1.0_dp
   kalman_covariance(:, 1) = [0.2_dp, 0.3_dp, 0.4_dp]
   kalman_state_covariance(:, 1) = [0.05_dp, 0.04_dp, 0.03_dp]
   kalman_transition(:, 1) = [1.0_dp, 0.9_dp, 0.8_dp]
   call set_random_seed(700)
   kalman_draw = bvartools_kalman_dk(kalman_y, kalman_z, &
      kalman_covariance, kalman_state_covariance, kalman_transition, &
      [0.0_dp], reshape([1.0_dp], [1, 1]))
   call assert_true(kalman_draw%info == 0 .and. &
      all(shape(kalman_draw%state) == [1, 4]), &
      'Durbin-Koopman compatibility draw')
   minnesota = bvartools_minnesota_prior(covariance, 2, deterministic_count=1, &
      cointegrated_var=.true.)
   call assert_true(minnesota%info == 0, 'Minnesota prior status')
   call assert_true(all(shape(minnesota%variance) == [2, 5]), 'Minnesota prior shape')
   call assert_close(minnesota%mean(1, 1), 1.0_dp, 1.0e-14_dp, &
      'Minnesota random-walk mean')

   inclusion = bvartools_inclusion_prior(2, 2, minnesota_like=.true., &
      deterministic_count=1)
   call assert_true(inclusion%info == 0, 'inclusion prior status')
   call assert_true(size(inclusion%include) == 8, 'deterministic exclusion')
   call assert_close(inclusion%probability(1, 1), 0.8_dp, 1.0e-14_dp, &
      'own-lag inclusion probability')

   ssvs_x(1, :) = 1.0_dp
   ssvs_x(2, :) = [(real(t, dp), t=1, 6)]
   ssvs_y(1, :) = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp]
   ssvs_y(2, :) = [2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 5.0_dp, 4.0_dp]
   fixed_ssvs_prior = bvartools_ssvs_prior(ssvs_y, ssvs_x, &
      tau=[0.1_dp, 2.0_dp], covariance_count=1)
   call assert_true(fixed_ssvs_prior%info == 0 .and. &
      size(fixed_ssvs_prior%tau0) == 5 .and. &
      all(abs(fixed_ssvs_prior%tau0 - 0.1_dp) < 1.0e-14_dp) .and. &
      all(abs(fixed_ssvs_prior%tau1 - 2.0_dp) < 1.0e-14_dp), &
      'fixed SSVS prior scales')
   automatic_ssvs_prior = bvartools_ssvs_prior(ssvs_y, ssvs_x, &
      semiautomatic=[0.1_dp, 10.0_dp], covariance_count=1)
   call assert_true(automatic_ssvs_prior%info == 0 .and. &
      automatic_ssvs_prior%semiautomatic .and. &
      automatic_ssvs_prior%regression_parameters == 4 .and. &
      automatic_ssvs_prior%covariance_parameters == 1, &
      'semiautomatic SSVS prior status')
   call assert_true(all(shape(automatic_ssvs_prior%coefficients) == [2, 2]) .and. &
      all(shape(automatic_ssvs_prior%residual_covariance) == [2, 2]) .and. &
      size(automatic_ssvs_prior%regression_standard_error) == 4, &
      'semiautomatic SSVS OLS shapes')
   call assert_close(automatic_ssvs_prior%coefficients(1, 2), &
      17.0_dp/35.0_dp, 1.0e-13_dp, 'semiautomatic SSVS OLS slope')
   call assert_true(all(abs(automatic_ssvs_prior%tau1(:4) - &
      100.0_dp*automatic_ssvs_prior%tau0(:4)) < 1.0e-12_dp), &
      'semiautomatic SSVS scale ratio')
   call assert_close(automatic_ssvs_prior%tau0(5), 0.05_dp, &
      1.0e-14_dp, 'SSVS covariance fallback spike scale')
   call assert_close(automatic_ssvs_prior%tau1(5), 10.0_dp, &
      1.0e-14_dp, 'SSVS covariance fallback slab scale')

   prior_mean = 0.0_dp
   prior_precision = 0.0_dp
   do t = 1, 4
      prior_precision(t, t) = 0.2_dp
   end do
   x = reshape([(0.1_dp*real(t, dp), t=1, 10)], [2, 5])
   y(1, :) = 0.4_dp*x(1, :) - 0.2_dp*x(2, :)
   y(2, :) = 0.1_dp*x(1, :) + 0.5_dp*x(2, :)
   precision = 0.0_dp
   precision(1, 1) = 1.0_dp
   precision(2, 2) = 1.0_dp
   posterior = bvartools_normal_posterior(y, x, precision, prior_mean, prior_precision)
   call assert_true(posterior%info == 0, 'normal posterior status')
   normals = [0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp]
   draw = bvartools_normal_draw(posterior, normals)
   call assert_true(all(abs(draw) < huge(1.0_dp)), 'normal posterior draw')

   z = 0.0_dp
   do t = 1, 5
      z(2*t - 1, [1, 3]) = x(:, t)
      z(2*t, [2, 4]) = x(:, t)
      time_precision(:, :, t) = precision
   end do
   sur_posterior = bvartools_sur_normal_posterior(y, z, time_precision, &
      prior_mean, prior_precision)
   call assert_true(sur_posterior%info == 0, 'SUR posterior status')
   call assert_true(maxval(abs(sur_posterior%mean - posterior%mean)) < 1.0e-10_dp, &
      'SUR and Kronecker posterior agreement')

   uniforms = [0.1_dp, 0.8_dp, 0.2_dp, 0.9_dp]
   ssvs = bvartools_ssvs(draw, spread(0.05_dp, 1, 4), spread(10.0_dp, 1, 4), &
      spread(0.5_dp, 1, 4), uniforms)
   call assert_true(ssvs%info == 0, 'SSVS status')
   call assert_true(all(ssvs%posterior_probability >= 0.0_dp .and. &
      ssvs%posterior_probability <= 1.0_dp), 'SSVS probabilities')

   errors = reshape([(0.02_dp*real(t - 5, dp), t=1, 10)], [2, 5])
   states = reshape([(0.01_dp*real(t, dp), t=1, 10)], [2, 5])
   initial_state = 0.0_dp
   measurement = bvartools_measurement_variance_posterior(errors, &
      [1.0_dp, 1.0_dp], [0.1_dp, 0.1_dp])
   state = bvartools_state_variance_posterior(states, initial_state, &
      [1.0_dp, 1.0_dp], [0.1_dp, 0.1_dp])
   call assert_true(measurement%info == 0 .and. state%info == 0, &
      'gamma posterior status')
   call set_random_seed(625)
   gamma_draw = bvartools_gamma_precision_draw(measurement)
   call assert_true(all(gamma_draw > 0.0_dp), 'gamma precision draws')

   covariance_data = bvartools_covar_prepare_data(y, time_precision)
   call assert_true(covariance_data%info == 0, 'covariance data status')
   call assert_true(all(shape(covariance_data%design) == [5, 1]), &
      'constant covariance design shape')
   covariance_prior_mean = 0.0_dp
   covariance_prior_precision(1, 1) = 0.1_dp
   covariance_posterior = bvartools_covar_const_posterior(y, time_precision, &
      covariance_prior_mean, covariance_prior_precision)
   call assert_true(covariance_posterior%info == 0, 'constant covariance posterior')
   covariance_state_precision = 2.0_dp
   tvp_posterior = bvartools_covar_tvp_posterior(y, time_precision, &
      covariance_state_precision, covariance_prior_mean)
   call assert_true(tvp_posterior%info == 0, 'TVP covariance posterior')
   psi = [(0.1_dp*real(t, dp), t=1, 5)]
   covariance_matrix = bvartools_covar_vector_to_matrix(psi, 2, 5)
   call assert_true(all(shape(covariance_matrix) == [10, 10]), &
      'covariance block matrix shape')
   expanded = bvartools_sur_const_to_tvp(z, 2, 5)
   call assert_true(all(shape(expanded) == [10, 20]), 'TVP SUR expansion shape')

   initial_inclusion = .true.
   bvs = bvartools_bvs(y, z, draw, initial_inclusion, time_precision, &
      spread(0.5_dp, 1, 4), [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp], &
      [0.7_dp, 0.6_dp, 0.5_dp, 0.4_dp])
   call assert_true(bvs%info == 0, 'BVS status')
   call assert_true(size(bvs%included) == 4, 'BVS indicator count')

   state_precision4 = 0.0_dp
   do t = 1, 5
      state_precision4(1, 1, t) = 2.0_dp
      state_precision4(2, 2, t) = 2.0_dp
      state_precision4(3, 3, t) = 2.0_dp
      state_precision4(4, 4, t) = 2.0_dp
   end do
   random_walk_posterior = bvartools_random_walk_posterior(y, z, time_precision, &
      state_precision4, prior_mean)
   call assert_true(random_walk_posterior%info == 0, 'random-walk posterior status')
   covariance_state_precision = 2.0_dp
   call set_random_seed(734)
   bvar_draws = bvartools_bvar_gibbs(y, x, prior_mean, prior_precision, &
      covariance, 8, 2, covariance_prior_scale=covariance, &
      covariance_prior_df=4.0_dp)
   call assert_true(bvar_draws%info == 0, 'inverse-Wishart BVAR status')
   call assert_true(bvar_draws%retained_draws == 8, 'inverse-Wishart BVAR draws')
   call set_random_seed(735)
   diagonal_bvar_draws = bvartools_bvar_gibbs(y, x, prior_mean, prior_precision, &
      covariance, 8, 2, gamma_shape=[1.0_dp, 1.0_dp], &
      gamma_rate=[0.1_dp, 0.1_dp], tau0=spread(0.05_dp, 1, 4), &
      tau1=spread(10.0_dp, 1, 4), &
      inclusion_probability=spread(0.5_dp, 1, 4))
   call assert_true(diagonal_bvar_draws%info == 0, 'diagonal SSVS BVAR status')
   call assert_true(allocated(diagonal_bvar_draws%included), 'BVAR inclusion draws')

   initial_beta(:, 1) = [1.0_dp, 0.0_dp]
   cointegration_precision = 0.0_dp
   cointegration_precision(1, 1) = 1.0_dp
   cointegration_precision(2, 2) = 1.0_dp
   loading_normals = [0.1_dp, -0.2_dp]
   loading_gamma_normals = [0.1_dp, -0.2_dp, 0.05_dp, -0.1_dp, 0.2_dp, -0.05_dp]
   beta_normals = [0.15_dp, -0.1_dp]
   cointegration = bvartools_cointegration_draw(y, initial_beta, x, precision, &
      0.2_dp, cointegration_precision, precision, loading_normals, beta_normals)
   call assert_true(cointegration%info == 0, 'KLS cointegration draw status')
   call assert_close(sum(cointegration%beta(:, 1)**2), 1.0_dp, 1.0e-10_dp, &
      'KLS normalized beta')
   call assert_true(maxval(abs(cointegration%pi - matmul(cointegration%alpha, &
      transpose(cointegration%beta)))) < 1.0e-10_dp, 'KLS Pi identity')
   cointegration_sur = bvartools_cointegration_sur_draw(y, initial_beta, x, &
      time_precision, 0.2_dp, cointegration_precision, precision, &
      loading_normals, beta_normals)
   call assert_true(cointegration_sur%info == 0, 'KLS SUR cointegration draw status')
   call assert_true(maxval(abs(cointegration_sur%alpha - cointegration%alpha)) < &
      1.0e-10_dp .and. maxval(abs(cointegration_sur%beta - cointegration%beta)) < &
      1.0e-10_dp, 'KLS SUR constant-precision agreement')
   cointegration_x = bvartools_cointegration_draw(y, initial_beta, x, precision, &
      0.2_dp, cointegration_precision, precision, loading_gamma_normals, &
      beta_normals, x=x, gamma_prior_mean=prior_mean, &
      gamma_prior_precision=prior_precision)
   call assert_true(cointegration_x%info == 0, 'KLS regression draw status')
   call assert_true(all(shape(cointegration_x%gamma) == [2, 2]), &
      'KLS unrestricted coefficient shape')
   cointegration_sur_x = bvartools_cointegration_sur_draw(y, initial_beta, x, &
      time_precision, 0.2_dp, cointegration_precision, precision, &
      loading_gamma_normals, beta_normals, x=x, gamma_prior_mean=prior_mean, &
      gamma_prior_precision=prior_precision)
   call assert_true(cointegration_sur_x%info == 0, &
      'KLS SUR regression draw status')
   call assert_true(all(shape(cointegration_sur_x%gamma) == [2, 2]) .and. &
      maxval(abs(cointegration_sur_x%pi - matmul(cointegration_sur_x%alpha, &
      transpose(cointegration_sur_x%beta)))) < 1.0e-10_dp, &
      'KLS SUR regression structure')
   time_precision(:, :, 5) = 2.0_dp*precision
   cointegration_sur = bvartools_cointegration_sur_draw(y, initial_beta, x, &
      time_precision, 0.2_dp, cointegration_precision, precision, &
      loading_normals, beta_normals)
   call assert_true(cointegration_sur%info == 0 .and. &
      maxval(abs(cointegration_sur%pi - cointegration%pi)) > 1.0e-8_dp, &
      'KLS observation-specific precision draw')
   time_precision(:, :, 5) = precision

   call set_random_seed(736)
   bvec_draws = bvartools_bvec_gibbs(y, x, initial_beta, covariance, 6, 2, &
      0.2_dp, cointegration_precision, precision, covariance_prior_df=4.0_dp, &
      x=x, gamma_prior_mean=prior_mean, gamma_prior_precision=prior_precision, &
      tau0=spread(0.05_dp, 1, 4), tau1=spread(10.0_dp, 1, 4), &
      inclusion_probability=spread(0.5_dp, 1, 4))
   call assert_true(bvec_draws%info == 0, 'inverse-Wishart BVEC status')
   call assert_true(bvec_draws%retained_draws == 6, 'inverse-Wishart BVEC draws')
   call assert_true(allocated(bvec_draws%included), 'BVEC inclusion draws')
   call assert_true(all(shape(bvec_draws%pi) == [4, 6]), 'BVEC Pi draw shape')
   call assert_true(maxval(abs(reshape(bvec_draws%pi(:, 1), [2, 2]) - &
      matmul(reshape(bvec_draws%alpha(:, 1), [2, 1]), transpose( &
      reshape(bvec_draws%beta(:, 1), [2, 1]))))) < 1.0e-10_dp, &
      'BVEC retained Pi identity')
   call set_random_seed(737)
   diagonal_bvec_draws = bvartools_bvec_gibbs(y, x, initial_beta, covariance, &
      6, 2, 0.2_dp, cointegration_precision, precision, &
      gamma_shape=[1.0_dp, 1.0_dp], gamma_rate=[0.1_dp, 0.1_dp])
   call assert_true(diagonal_bvec_draws%info == 0, 'diagonal BVEC status')
   call assert_true(all(abs(diagonal_bvec_draws%covariance(1, 2, :)) < &
      1.0e-14_dp), 'diagonal BVEC covariance')

   tvp_bvs = bvartools_tvp_bvs(y, z, reshape([(0.01_dp*real(t, dp), &
      t=1, 20)], [4, 5]), [.true., .true., .true., .true.], &
      time_precision, spread(0.5_dp, 1, 4), spread(0.0_dp, 1, 4), &
      spread(0.5_dp, 1, 4))
   call assert_true(tvp_bvs%info == 0 .and. size(tvp_bvs%included) == 4, &
      'grouped TVP coefficient BVS')
   covariance_bvs = bvartools_tvp_covariance_bvs(errors, &
      reshape([(0.01_dp*real(t, dp), t=1, 5)], [1, 5]), [.true.], &
      time_precision, [0.5_dp], [0.0_dp], [0.5_dp])
   call assert_true(covariance_bvs%info == 0 .and. &
      size(covariance_bvs%included) == 1, 'grouped TVP covariance BVS')

   call set_random_seed(738)
   tvp_bvar_draws = bvartools_tvp_bvar_gibbs(y, z, prior_mean, covariance, &
      4, 2, prior_mean, prior_precision, spread(1.0_dp, 1, 4), &
      spread(0.1_dp, 1, 4), covariance_prior_scale=covariance, &
      covariance_prior_df=4.0_dp, &
      inclusion_probability=spread(0.5_dp, 1, 4))
   call assert_true(tvp_bvar_draws%info == 0, 'inverse-Wishart TVP-BVAR status')
   call assert_true(tvp_bvar_draws%retained_draws == 4, 'TVP-BVAR retained draws')
   call assert_true(all(shape(tvp_bvar_draws%states) == [20, 4]), &
      'TVP-BVAR state-path shape')
   call assert_true(all(tvp_bvar_draws%state_variance > 0.0_dp), &
      'TVP-BVAR state variances')
   call assert_true(allocated(tvp_bvar_draws%included), &
      'TVP-BVAR inclusion draws')
   call set_random_seed(739)
   diagonal_tvp_bvar_draws = bvartools_tvp_bvar_gibbs(y, z, prior_mean, &
      covariance, 4, 2, prior_mean, prior_precision, spread(1.0_dp, 1, 4), &
      spread(0.1_dp, 1, 4), measurement_shape=[1.0_dp, 1.0_dp], &
      measurement_rate=[0.1_dp, 0.1_dp])
   call assert_true(diagonal_tvp_bvar_draws%info == 0, 'diagonal TVP-BVAR status')
   call assert_true(all(abs(diagonal_tvp_bvar_draws%covariance(1, 2, :)) < &
      1.0e-14_dp), 'diagonal TVP-BVAR covariance')

   volatility_draw = bvartools_stochastic_volatility_draw(errors, &
      spread([0.0_dp, log(0.8_dp)], 2, 5), [0.05_dp, 0.05_dp], &
      [0.0_dp, log(0.8_dp)], [0.0001_dp, 0.0001_dp], &
      spread(spread(0.5_dp, 1, 2), 2, 5), &
      spread(spread(0.0_dp, 1, 2), 2, 5))
   call assert_true(volatility_draw%info == 0, 'KSC volatility draw status')
   call assert_true(all(volatility_draw%component >= 1 .and. &
      volatility_draw%component <= 7), 'KSC mixture components')
   ocsn_volatility_draw = bvartools_stochastic_volatility_ocsn_draw(errors, &
      spread([0.0_dp, log(0.8_dp)], 2, 5), [0.05_dp, 0.05_dp], &
      [0.0_dp, log(0.8_dp)], [0.0001_dp, 0.0001_dp], &
      spread(spread(0.5_dp, 1, 2), 2, 5), &
      spread(spread(0.0_dp, 1, 2), 2, 5))
   call assert_true(ocsn_volatility_draw%info == 0, &
      'OCSN volatility draw status')
   call assert_true(all(ocsn_volatility_draw%component >= 1 .and. &
      ocsn_volatility_draw%component <= 10), 'OCSN mixture components')
   call set_random_seed(741)
   volatility_tvp_bvar_draws = bvartools_tvp_bvar_gibbs(y, z, prior_mean, &
      covariance, 3, 1, prior_mean, prior_precision, spread(1.0_dp, 1, 4), &
      spread(0.1_dp, 1, 4), &
      initial_log_variance=spread([0.0_dp, log(0.8_dp)], 2, 5), &
      initial_log_variance_level=[0.0_dp, log(0.8_dp)], &
      initial_log_variance_state_variance=[0.05_dp, 0.05_dp], &
      log_variance_state_shape=[1.0_dp, 1.0_dp], &
      log_variance_state_rate=[0.1_dp, 0.1_dp], &
      log_variance_initial_prior_mean=[0.0_dp, log(0.8_dp)], &
      log_variance_initial_prior_precision=precision, &
      log_variance_offset=[0.0001_dp, 0.0001_dp], volatility_method='ocsn')
   call assert_true(volatility_tvp_bvar_draws%info == 0, &
      'OCSN stochastic-volatility TVP-BVAR status')
   call assert_true(all(shape(volatility_tvp_bvar_draws%log_variance) == [10, 3]), &
      'stochastic-volatility retained paths')
   call assert_true(all(volatility_tvp_bvar_draws%log_variance_state_variance > &
      0.0_dp), 'stochastic-volatility state variances')
   call assert_true(all([(volatility_tvp_bvar_draws%time_covariance(t, t, 1, 1) > &
      0.0_dp, t=1, 2)]), 'stochastic-volatility timestamp covariance')

   initial_alpha(:, 1) = [0.1_dp, -0.1_dp]
   call set_random_seed(742)
   tvp_bvec_draws = bvartools_tvp_bvec_gibbs(y, x, initial_alpha, initial_beta, &
      covariance, 4, 2, initial_alpha(:, 1), precision, [1.0_dp, 1.0_dp], &
      [0.1_dp, 0.1_dp], 0.9_dp, precision, initial_beta(:, 1), precision, &
      covariance_prior_scale=covariance, covariance_prior_df=4.0_dp)
   call assert_true(tvp_bvec_draws%info == 0, 'TVP-BVEC status')
   call assert_true(all(shape(tvp_bvec_draws%alpha) == [10, 4]) .and. &
      all(shape(tvp_bvec_draws%beta) == [10, 4]), 'TVP-BVEC path shapes')
   call assert_true(maxval(abs(reshape(tvp_bvec_draws%pi(1:4, 1), [2, 2]) - &
      matmul(reshape(tvp_bvec_draws%alpha(1:2, 1), [2, 1]), transpose( &
      reshape(tvp_bvec_draws%beta(1:2, 1), [2, 1]))))) < 1.0e-10_dp, &
      'TVP-BVEC retained Pi identity')
   initial_gamma = 0.0_dp
   state_prior6 = 0.0_dp
   state_prior6(1:2) = initial_alpha(:, 1)
   state_prior_precision6 = 0.0_dp
   do t = 1, 6
      state_prior_precision6(t, t) = 0.2_dp
   end do
   call set_random_seed(743)
   tvp_bvec_x_draws = bvartools_tvp_bvec_gibbs(y, x, initial_alpha, &
      initial_beta, covariance, 3, 1, state_prior6, state_prior_precision6, &
      spread(1.0_dp, 1, 6), spread(0.1_dp, 1, 6), 0.9_dp, precision, &
      initial_beta(:, 1), precision, measurement_shape=[1.0_dp, 1.0_dp], &
      measurement_rate=[0.1_dp, 0.1_dp], x=x, initial_gamma=initial_gamma)
   call assert_true(tvp_bvec_x_draws%info == 0, 'TVP-BVEC regression status')
   call assert_true(all(shape(tvp_bvec_x_draws%gamma) == [20, 3]), &
      'TVP-BVEC unrestricted path shape')
   call set_random_seed(744)
   tvp_bvec_x_draws = bvartools_tvp_bvec_gibbs(y, x, initial_alpha, &
      initial_beta, covariance, 3, 1, state_prior6, state_prior_precision6, &
      spread(1.0_dp, 1, 6), spread(0.1_dp, 1, 6), 0.9_dp, precision, &
      initial_beta(:, 1), precision, measurement_shape=[1.0_dp, 1.0_dp], &
      measurement_rate=[0.1_dp, 0.1_dp], x=x, initial_gamma=initial_gamma, &
      inclusion_probability=spread(0.5_dp, 1, 4))
   call assert_true(tvp_bvec_x_draws%info == 0 .and. &
      allocated(tvp_bvec_x_draws%included), 'TVP-BVEC BVS draws')
   call set_random_seed(745)
   tvp_bvec_covariance_draws = bvartools_tvp_bvec_gibbs(y, x, initial_alpha, &
      initial_beta, covariance, 3, 1, initial_alpha(:, 1), precision, &
      [1.0_dp, 1.0_dp], [0.1_dp, 0.1_dp], 0.9_dp, precision, &
      initial_beta(:, 1), precision, measurement_shape=[1.0_dp, 1.0_dp], &
      measurement_rate=[0.1_dp, 0.1_dp], initial_covariance_state=[0.0_dp], &
      covariance_state_initial_prior_mean=[0.0_dp], &
      covariance_state_initial_prior_precision=reshape([0.2_dp], [1, 1]), &
      covariance_state_shape=[1.0_dp], covariance_state_rate=[0.1_dp], &
      covariance_inclusion_probability=[0.5_dp])
   call assert_true(tvp_bvec_covariance_draws%info == 0, &
      'TVP-BVEC covariance-state status')
   call assert_true(allocated(tvp_bvec_covariance_draws%covariance_state) .and. &
      allocated(tvp_bvec_covariance_draws%covariance_included), &
      'TVP-BVEC covariance-state draws')
   call assert_true(all(shape(tvp_bvec_covariance_draws%covariance_state) == [5, 3]), &
      'TVP-BVEC covariance-state path shape')
   call assert_true(all(tvp_bvec_covariance_draws%covariance_state_variance > 0.0_dp), &
      'TVP-BVEC covariance-state variances')
   call assert_true(abs(tvp_bvec_covariance_draws%time_covariance(1, 2, 1, 1) - &
      tvp_bvec_covariance_draws%time_covariance(2, 1, 1, 1)) < 1.0e-12_dp, &
      'TVP-BVEC time covariance symmetry')
   call assert_true(all([(tvp_bvec_covariance_draws%time_covariance(t, t, 1, 1) > &
      0.0_dp, t=1, 2)]), 'TVP-BVEC time covariance positive diagonal')
   call set_random_seed(746)
   tvp_bvec_sv_draws = bvartools_tvp_bvec_gibbs(y, x, initial_alpha, &
      initial_beta, covariance, 3, 1, initial_alpha(:, 1), precision, &
      [1.0_dp, 1.0_dp], [0.1_dp, 0.1_dp], 0.9_dp, precision, &
      initial_beta(:, 1), precision, &
      initial_log_variance=spread([0.0_dp, log(0.8_dp)], 2, 5), &
      initial_log_variance_level=[0.0_dp, log(0.8_dp)], &
      initial_log_variance_state_variance=[0.05_dp, 0.05_dp], &
      log_variance_state_shape=[1.0_dp, 1.0_dp], &
      log_variance_state_rate=[0.1_dp, 0.1_dp], &
      log_variance_initial_prior_mean=[0.0_dp, log(0.8_dp)], &
      log_variance_initial_prior_precision=precision, &
      log_variance_offset=[0.0001_dp, 0.0001_dp])
   call assert_true(tvp_bvec_sv_draws%info == 0, &
      'stochastic-volatility TVP-BVEC status')
   call assert_true(all(shape(tvp_bvec_sv_draws%log_variance) == [10, 3]), &
      'stochastic-volatility TVP-BVEC retained paths')
   call assert_true(all(tvp_bvec_sv_draws%log_variance_state_variance > 0.0_dp), &
      'stochastic-volatility TVP-BVEC state variances')
   call assert_true(abs(tvp_bvec_sv_draws%time_covariance(1, 1, 1, 1) - &
      exp(tvp_bvec_sv_draws%log_variance(1, 1))) < 1.0e-12_dp, &
      'stochastic-volatility TVP-BVEC covariance path')

   call set_random_seed(740)
   tvp_covariance_draws = bvartools_tvp_covariance_gibbs(errors, &
      [1.0_dp, 1.0_dp], [0.0_dp], 4, 2, [0.0_dp], &
      reshape([0.2_dp], [1, 1]), [1.0_dp], [0.1_dp], &
      [1.0_dp, 1.0_dp], [0.1_dp, 0.1_dp], &
      inclusion_probability=[0.5_dp])
   call assert_true(tvp_covariance_draws%info == 0, &
      'TVP covariance Gibbs status')
   call assert_true(tvp_covariance_draws%retained_draws == 4, &
      'TVP covariance retained draws')
   call assert_true(all(shape(tvp_covariance_draws%states) == [5, 4]), &
      'TVP covariance state-path shape')
   call assert_true(all(tvp_covariance_draws%state_variance > 0.0_dp) .and. &
      all(tvp_covariance_draws%diagonal_variance > 0.0_dp), &
      'TVP covariance innovation variances')
   call assert_true(allocated(tvp_covariance_draws%included), &
      'TVP covariance inclusion draws')
   call assert_true(abs(tvp_covariance_draws%covariance(1, 2, 1, 1) - &
      tvp_covariance_draws%covariance(2, 1, 1, 1)) < 1.0e-12_dp, &
      'TVP covariance symmetry')
   call assert_true(all([(tvp_covariance_draws%covariance(t, t, 1, 1) > &
      0.0_dp, t=1, 2)]), 'TVP covariance positive diagonal')

   call set_random_seed(751)
   joint_tvp_draws = bvartools_joint_tvp_bvar_gibbs(y, z, prior_mean, &
      [1.0_dp, 1.0_dp], [0.0_dp], 3, 1, prior_mean, prior_precision, &
      spread(1.0_dp, 1, 4), spread(0.1_dp, 1, 4), [0.0_dp], &
      covariance_prior_precision, [1.0_dp], [0.1_dp], &
      [1.0_dp, 1.0_dp], [0.1_dp, 0.1_dp], &
      inclusion_probability=spread(0.5_dp, 1, 4), &
      covariance_inclusion_probability=[0.5_dp])
   call assert_true(joint_tvp_draws%info == 0 .and. &
      joint_tvp_draws%retained_draws == 3, 'joint TVP-BVAR sampler status')
   call assert_true(all(shape(joint_tvp_draws%coefficient_states) == [20, 3]) .and. &
      all(shape(joint_tvp_draws%covariance_states) == [5, 3]), &
      'joint TVP-BVAR retained state shapes')
   call assert_true(allocated(joint_tvp_draws%coefficient_included) .and. &
      allocated(joint_tvp_draws%covariance_included) .and. &
      allocated(joint_tvp_draws%covariance_state_variance), &
      'joint TVP-BVAR selection and state variance draws')
   call assert_true(abs(joint_tvp_draws%covariance(1, 2, 1, 1) - &
      joint_tvp_draws%covariance(2, 1, 1, 1)) < 1.0e-12_dp .and. &
      all([(joint_tvp_draws%covariance(t, t, 1, 1) > 0.0_dp, t=1, 2)]), &
      'joint TVP-BVAR covariance validity')

   call set_random_seed(752)
   joint_constant_draws = bvartools_joint_tvp_bvar_gibbs(y, z, prior_mean, &
      [1.0_dp, 1.0_dp], [0.0_dp], 3, 1, prior_mean, prior_precision, &
      spread(1.0_dp, 1, 4), spread(0.1_dp, 1, 4), [0.0_dp], &
      covariance_prior_precision, [1.0_dp], [0.1_dp], &
      [1.0_dp, 1.0_dp], [0.1_dp, 0.1_dp], &
      time_varying_covariance=.false., covariance_tau0=[0.05_dp], &
      covariance_tau1=[10.0_dp], covariance_ssvs_probability=[0.5_dp])
   call assert_true(joint_constant_draws%info == 0 .and. &
      .not. joint_constant_draws%time_varying_covariance .and. &
      allocated(joint_constant_draws%covariance_included), &
      'joint constant-covariance BVAR sampler status')
   call assert_true(.not. allocated(joint_constant_draws%covariance_state_variance) .and. &
      all(abs(joint_constant_draws%covariance_states(1, :) - &
      joint_constant_draws%covariance_states(5, :)) < 1.0e-12_dp), &
      'joint constant-covariance repeated state paths')

   dfm_factor_posterior = bvartools_dfm_factor_posterior(y, &
      reshape([1.0_dp, 0.4_dp], [2, 1]), [1.0_dp, 1.0_dp], &
      reshape([0.5_dp], [1, 1]), [1.0_dp])
   call assert_true(dfm_factor_posterior%info == 0, 'DFM factor posterior status')
   call assert_true(size(dfm_factor_posterior%mean) == 5 .and. &
      all(shape(dfm_factor_posterior%covariance) == [5, 5]), &
      'DFM factor posterior shape')
   call set_random_seed(747)
   dfm_draws = bvartools_dfm_gibbs(y, 1, 1, 3, 1, &
      reshape([0.1_dp], [1, 1]), [1.0_dp, 1.0_dp], &
      [0.1_dp, 0.1_dp], [0.0_dp], reshape([0.1_dp], [1, 1]), &
      [1.0_dp], [0.1_dp])
   call assert_true(dfm_draws%info == 0 .and. dfm_draws%retained_draws == 3, &
      'DFM Gibbs status')
   call assert_true(all(shape(dfm_draws%loadings) == [2, 3]) .and. &
      all(shape(dfm_draws%factors) == [5, 3]) .and. &
      all(shape(dfm_draws%transition) == [1, 3]), 'DFM retained shapes')
   call assert_true(all(abs(dfm_draws%loadings(1, :) - 1.0_dp) < 1.0e-14_dp), &
      'DFM unit-loading identification')
   call assert_true(all(dfm_draws%measurement_variance > 0.0_dp) .and. &
      all(dfm_draws%factor_variance > 0.0_dp), 'DFM positive variances')

   predictive_draws%retained_draws = 2
   allocate(predictive_draws%coefficients(1, 2), &
      predictive_draws%covariance(1, 1, 2))
   predictive_draws%coefficients(1, :) = [0.5_dp, 1.0_dp]
   predictive_draws%covariance = 0.0_dp
   predictive_history = 2.0_dp
   predictive_normals = 0.0_dp
   prediction = bvartools_bvar_predictive_from_random(predictive_draws, &
      predictive_history, 1, predictive_normals, probability=0.5_dp)
   call assert_true(prediction%info == 0, 'constant BVAR prediction status')
   call assert_true(all(shape(prediction%paths) == [1, 2, 2]), &
      'constant BVAR predictive shape')
   call assert_close(prediction%paths(1, 2, 1), 0.5_dp, 1.0e-14_dp, &
      'constant BVAR recursive path')
   call assert_close(prediction%mean(1, 1), 1.5_dp, 1.0e-14_dp, &
      'constant BVAR predictive mean')
   call assert_close(prediction%lower(1, 1), 1.25_dp, 1.0e-14_dp, &
      'constant BVAR lower quantile')
   call assert_close(prediction%upper(1, 1), 1.75_dp, 1.0e-14_dp, &
      'constant BVAR upper quantile')

   structural_draws%retained_draws = 2
   allocate(structural_draws%coefficients(2, 2), &
      structural_draws%covariance(1, 1, 2))
   structural_draws%coefficients = 0.0_dp
   structural_draws%coefficients(2, :) = 1.0_dp
   structural_draws%covariance = 0.0_dp
   structural_normals = 0.0_dp
   structural_impact = 2.0_dp
   future_constant = 1.0_dp
   structural_prediction = bvartools_bvar_predictive_from_random( &
      structural_draws, predictive_history, 1, structural_normals, &
      future_constant, structural_impact)
   call assert_true(structural_prediction%info == 0, &
      'structural BVAR prediction status')
   call assert_close(structural_prediction%paths(1, 1, 1), 2.0_dp, &
      1.0e-14_dp, 'structural transformation of conditional mean')

   predictive_tvp_draws%retained_draws = 2
   allocate(predictive_tvp_draws%states(2, 2), &
      predictive_tvp_draws%time_covariance(1, 1, 2, 2))
   predictive_tvp_draws%states(:, 1) = [0.1_dp, 0.5_dp]
   predictive_tvp_draws%states(:, 2) = [0.2_dp, 1.0_dp]
   predictive_tvp_draws%time_covariance = 0.0_dp
   tvp_prediction = bvartools_tvp_bvar_predictive_from_random( &
      predictive_tvp_draws, predictive_history, 1, 2, predictive_normals)
   call assert_true(tvp_prediction%info == 0, 'TVP-BVAR prediction status')
   call assert_close(tvp_prediction%paths(1, 2, 1), 0.5_dp, 1.0e-14_dp, &
      'TVP-BVAR terminal-state recursion')
   call assert_close(tvp_prediction%paths(1, 2, 2), 2.0_dp, 1.0e-14_dp, &
      'TVP-BVAR draw-specific terminal state')

   irf_draws%retained_draws = 2
   allocate(irf_draws%coefficients(4, 2), irf_draws%covariance(2, 2, 2))
   irf_draws%coefficients(:, 1) = [0.5_dp, 0.0_dp, 0.1_dp, 0.25_dp]
   irf_draws%coefficients(:, 2) = irf_draws%coefficients(:, 1)
   irf_draws%covariance(:, :, 1) = reshape( &
      [4.0_dp, 1.0_dp, 1.0_dp, 9.0_dp], [2, 2])
   irf_draws%covariance(:, :, 2) = irf_draws%covariance(:, :, 1)
   feir = bvartools_bvar_irf(irf_draws, 1, 1, 'feir')
   call assert_true(feir%info == 0 .and. &
      all(shape(feir%paths) == [2, 2, 2, 2]), 'FEIR posterior shape')
   call assert_close(feir%paths(1, 1, 2, 1), 0.5_dp, 1.0e-14_dp, &
      'FEIR first-lag response')
   call assert_close(feir%paths(1, 2, 2, 1), 0.1_dp, 1.0e-14_dp, &
      'FEIR cross response')

   gir = bvartools_bvar_irf(irf_draws, 1, 0, 'gir')
   call assert_true(gir%info == 0, 'GIR posterior status')
   call assert_close(gir%paths(2, 1, 1, 1), 0.25_dp, 1.0e-14_dp, &
      'GIR covariance normalization')
   oir = bvartools_bvar_irf(irf_draws, 1, 0, 'oir', &
      scale_by_standard_deviation=.true.)
   call assert_true(oir%info == 0, 'OIR posterior status')
   call assert_close(oir%paths(1, 1, 1, 1), 2.0_dp, 1.0e-14_dp, &
      'OIR standard-deviation shock')
   call assert_close(oir%paths(2, 1, 1, 1), 0.5_dp, 1.0e-14_dp, &
      'OIR Cholesky transmission')

   irf_impact = 0.0_dp
   irf_impact(1, 1, :) = 2.0_dp
   irf_impact(2, 2, :) = 3.0_dp
   sir = bvartools_bvar_irf(irf_draws, 1, 0, 'sir', &
      structural_impact=irf_impact)
   call assert_true(sir%info == 0, 'SIR posterior status')
   call assert_close(sir%paths(2, 2, 1, 1), 3.0_dp, 1.0e-14_dp, &
      'SIR structural impact')
   sgir = bvartools_bvar_irf(irf_draws, 1, 0, 'sgir', &
      structural_impact=irf_impact)
   call assert_true(sgir%info == 0, 'SGIR posterior status')
   call assert_close(sgir%paths(2, 1, 1, 1), 0.75_dp, 1.0e-14_dp, &
      'SGIR structural covariance normalization')
   cumulative_feir = bvartools_bvar_irf(irf_draws, 1, 1, 'feir', &
      cumulative=.true.)
   call assert_close(cumulative_feir%paths(1, 1, 2, 1), 1.5_dp, &
      1.0e-14_dp, 'cumulative FEIR response')

   irf_tvp_draws%retained_draws = 2
   allocate(irf_tvp_draws%states(8, 2), &
      irf_tvp_draws%time_covariance(2, 2, 2, 2))
   irf_tvp_draws%states = 0.0_dp
   irf_tvp_draws%states(5:8, 1) = irf_draws%coefficients(:, 1)
   irf_tvp_draws%states(5:8, 2) = irf_draws%coefficients(:, 2)
   irf_tvp_draws%time_covariance(:, :, 1, :) = irf_draws%covariance
   irf_tvp_draws%time_covariance(:, :, 2, :) = irf_draws%covariance
   tvp_irf = bvartools_tvp_bvar_irf(irf_tvp_draws, 2, 1, 2, 2, 1, 'feir')
   call assert_true(tvp_irf%info == 0, 'TVP-BVAR IRF status')
   call assert_close(tvp_irf%paths(1, 1, 2, 1), 0.5_dp, 1.0e-14_dp, &
      'TVP-BVAR period coefficient response')

   oir_fevd = bvartools_bvar_fevd(irf_draws, 1, 1, 'oir')
   call assert_true(oir_fevd%info == 0 .and. &
      all(shape(oir_fevd%paths) == [2, 2, 2, 2]), 'OIR FEVD posterior shape')
   call assert_close(oir_fevd%paths(1, 1, 1, 1), 1.0_dp, 1.0e-14_dp, &
      'OIR FEVD contemporaneous own share')
   call assert_close(sum(oir_fevd%paths(2, :, 2, 1)), 1.0_dp, 1.0e-14_dp, &
      'OIR FEVD shares sum to one')

   gir_fevd = bvartools_bvar_fevd(irf_draws, 1, 0, 'gir')
   call assert_true(gir_fevd%info == 0 .and. .not. gir_fevd%normalized, &
      'GIR FEVD posterior status')
   call assert_close(gir_fevd%paths(1, 2, 1, 1), 1.0_dp/36.0_dp, &
      1.0e-14_dp, 'GIR FEVD impulse-variance scaling')
   normalized_gir_fevd = bvartools_bvar_fevd(irf_draws, 1, 0, 'gir', &
      normalize_generalized=.true.)
   call assert_true(normalized_gir_fevd%info == 0 .and. &
      normalized_gir_fevd%normalized, 'normalized GIR FEVD status')
   call assert_close(sum(normalized_gir_fevd%paths(1, :, 1, 1)), 1.0_dp, &
      1.0e-14_dp, 'normalized GIR FEVD shares')

   sir_fevd = bvartools_bvar_fevd(irf_draws, 1, 0, 'sir', &
      structural_impact=irf_impact)
   call assert_true(sir_fevd%info == 0, 'SIR FEVD posterior status')
   call assert_close(sir_fevd%paths(2, 2, 1, 1), 1.0_dp, 1.0e-14_dp, &
      'SIR FEVD structural share')
   sgir_fevd = bvartools_bvar_fevd(irf_draws, 1, 0, 'sgir', &
      normalize_generalized=.true., structural_impact=irf_impact)
   call assert_true(sgir_fevd%info == 0, 'SGIR FEVD posterior status')
   call assert_close(sum(sgir_fevd%paths(1, :, 1, 1)), 1.0_dp, 1.0e-14_dp, &
      'normalized SGIR FEVD shares')

   tvp_fevd = bvartools_tvp_bvar_fevd(irf_tvp_draws, 2, 1, 2, 2, 1, 'oir')
   call assert_true(tvp_fevd%info == 0, 'TVP-BVAR FEVD status')
   call assert_close(sum(tvp_fevd%paths(1, :, 2, 1)), 1.0_dp, 1.0e-14_dp, &
      'TVP-BVAR FEVD selected-period shares')

   conversion_pi = reshape([-0.2_dp, 0.05_dp, 0.1_dp, -0.3_dp, &
      0.7_dp, -0.4_dp], [2, 3])
   conversion_gamma = reshape([0.4_dp, 0.2_dp, 0.1_dp, 0.3_dp, &
      1.5_dp, -0.5_dp], [2, 3])
   level_ar = bvartools_vecm_level_ar(conversion_pi(:, :2), &
      reshape(conversion_gamma(:, :2), [2, 2, 1]))
   call assert_true(all(shape(level_ar) == [2, 2, 2]), &
      'VECM-to-VAR coefficient shape')
   call assert_close(level_ar(1, 1, 1), 1.2_dp, 1.0e-14_dp, &
      'VECM-to-VAR first level lag')
   call assert_close(level_ar(2, 1, 2), -0.2_dp, 1.0e-14_dp, &
      'VECM-to-VAR final level lag')
   allocate(zero_gamma(2, 2, 0))
   level_ar = bvartools_vecm_level_ar(reshape([0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp], [2, 2]), zero_gamma)
   call assert_close(level_ar(1, 1, 1), 1.0_dp, 1.0e-14_dp, &
      'rank-zero one-lag level VAR')

   level_exogenous = bvartools_vecm_level_exogenous( &
      reshape([0.2_dp, -0.1_dp], [2, 1]), &
      reshape([0.5_dp, 0.2_dp, 0.1_dp, 0.05_dp], [2, 1, 2]))
   call assert_close(level_exogenous(1, 1, 1), 0.5_dp, 1.0e-14_dp, &
      'level exogenous contemporaneous coefficient')
   call assert_close(level_exogenous(1, 1, 2), -0.2_dp, 1.0e-14_dp, &
      'level exogenous middle coefficient')
   call assert_close(level_exogenous(2, 1, 3), -0.05_dp, 1.0e-14_dp, &
      'level exogenous final coefficient')

   conversion_draws%retained_draws = 2
   allocate(conversion_draws%pi(6, 2), conversion_draws%gamma(6, 2))
   conversion_draws%pi(:, 1) = reshape(conversion_pi, [6])
   conversion_draws%pi(:, 2) = conversion_draws%pi(:, 1)
   conversion_draws%gamma(:, 1) = reshape(conversion_gamma, [6])
   conversion_draws%gamma(:, 2) = conversion_draws%gamma(:, 1)
   exogenous_pi(:, 1, 1) = [0.2_dp, -0.1_dp]
   exogenous_pi(:, 1, 2) = exogenous_pi(:, 1, 1)
   exogenous_difference(:, 1, :, 1) = reshape( &
      [0.5_dp, 0.2_dp, 0.1_dp, 0.05_dp], [2, 2])
   exogenous_difference(:, 1, :, 2) = exogenous_difference(:, 1, :, 1)
   level_draws = bvartools_bvec_to_level_var(conversion_draws, 2, 2, &
      exogenous_pi, exogenous_difference)
   call assert_true(level_draws%info == 0 .and. &
      all(shape(level_draws%ar) == [2, 2, 2, 2]), &
      'constant BVEC-to-VAR status and shape')
   call assert_close(level_draws%restricted_deterministic(1, 1, 1), &
      0.7_dp, 1.0e-14_dp, 'restricted deterministic preservation')
   call assert_close(level_draws%unrestricted_deterministic(1, 1, 1), &
      1.5_dp, 1.0e-14_dp, 'unrestricted deterministic preservation')
   call assert_close(level_draws%exogenous(1, 1, 2, 1), -0.2_dp, &
      1.0e-14_dp, 'constant BVEC exogenous conversion')

   conversion_tvp_draws%retained_draws = 1
   allocate(conversion_tvp_draws%pi(12, 1), conversion_tvp_draws%gamma(12, 1))
   conversion_tvp_draws%pi(:, 1) = [reshape(conversion_pi, [6]), &
      2.0_dp*reshape(conversion_pi, [6])]
   conversion_tvp_draws%gamma(:, 1) = [reshape(conversion_gamma, [6]), &
      2.0_dp*reshape(conversion_gamma, [6])]
   tvp_level_draws = bvartools_tvp_bvec_to_level_var( &
      conversion_tvp_draws, 2, 2, 2)
   call assert_true(tvp_level_draws%info == 0 .and. &
      all(shape(tvp_level_draws%ar) == [2, 2, 2, 2, 1]), &
      'TVP-BVEC-to-VAR status and shape')
   call assert_close(tvp_level_draws%restricted_deterministic(1, 1, 2, 1), &
      1.4_dp, 1.0e-14_dp, 'TVP restricted deterministic period')

   reconstructed = bvartools_reconstruct_levels([10.0_dp, 20.0_dp], &
      reshape([1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp], [2, 2]))
   call assert_true(all(shape(reconstructed) == [2, 3]), &
      'reconstructed level shape')
   call assert_close(reconstructed(1, 3), 11.5_dp, 1.0e-14_dp, &
      'reconstructed cumulative level')

   raw_series(1, :) = [(real(t, dp), t=1, 8)]
   raw_series(2, :) = [(10.0_dp*real(t, dp), t=1, 8)]
   raw_exogenous(1, :) = [(2.0_dp*real(t, dp), t=1, 8)]
   prepared_var = bvartools_prepare_var(raw_series, 2, raw_exogenous, 1, &
      include_constant=.true., include_trend=.true., seasonal_period=4, &
      holdout=2, time_varying=.true., structural=.true.)
   call assert_true(prepared_var%info == 0 .and. &
      all(shape(prepared_var%y) == [2, 4]) .and. &
      all(shape(prepared_var%x) == [11, 4]), 'prepared BVAR matrix shapes')
   call assert_close(prepared_var%y(1, 1), 3.0_dp, 1.0e-14_dp, &
      'prepared BVAR aligned response')
   call assert_close(prepared_var%x(1, 1), 2.0_dp, 1.0e-14_dp, &
      'prepared BVAR first endogenous lag')
   call assert_close(prepared_var%x(5, 1), 6.0_dp, 1.0e-14_dp, &
      'prepared BVAR current exogenous value')
   call assert_close(prepared_var%x(7, 1), 1.0_dp, 1.0e-14_dp, &
      'prepared BVAR constant')
   call assert_close(prepared_var%x(8, 1), 1.0_dp, 1.0e-14_dp, &
      'prepared BVAR aligned trend')
   call assert_close(prepared_var%x(11, 1), 1.0_dp, 1.0e-14_dp, &
      'prepared BVAR seasonal indicator')
   call assert_true(all(shape(prepared_var%holdout_y) == [2, 2]) .and. &
      all(shape(prepared_var%sur) == [8, 22]) .and. &
      all(shape(prepared_var%tvp_sur) == [8, 88]), &
      'prepared BVAR holdout and SUR shapes')
   call assert_true(all(shape(prepared_var%structural) == [8, 1]), &
      'prepared structural BVAR design shape')
   call assert_close(prepared_var%structural(2, 1), -3.0_dp, 1.0e-14_dp, &
      'prepared structural BVAR off-diagonal regressor')

   prepared_vecm = bvartools_prepare_vecm(raw_series, 2, raw_exogenous, 2, &
      constant='restricted', trend='unrestricted', seasonal='unrestricted', &
      seasonal_period=4, holdout=1, time_varying=.true., structural=.true.)
   call assert_true(prepared_vecm%info == 0 .and. &
      all(shape(prepared_vecm%y) == [2, 5]) .and. &
      all(shape(prepared_vecm%w) == [4, 5]) .and. &
      all(shape(prepared_vecm%x) == [8, 5]), 'prepared BVEC matrix shapes')
   call assert_close(prepared_vecm%y(2, 1), 10.0_dp, 1.0e-14_dp, &
      'prepared BVEC differenced response')
   call assert_close(prepared_vecm%w(1, 1), 2.0_dp, 1.0e-14_dp, &
      'prepared BVEC lagged endogenous level')
   call assert_close(prepared_vecm%w(3, 1), 4.0_dp, 1.0e-14_dp, &
      'prepared BVEC lagged exogenous level')
   call assert_close(prepared_vecm%w(4, 1), 1.0_dp, 1.0e-14_dp, &
      'prepared BVEC restricted constant')
   call assert_close(prepared_vecm%x(3, 1), 2.0_dp, 1.0e-14_dp, &
      'prepared BVEC current exogenous difference')
   call assert_close(prepared_vecm%x(5, 1), 1.0_dp, 1.0e-14_dp, &
      'prepared BVEC unrestricted trend')
   call assert_close(prepared_vecm%x(8, 1), 1.0_dp, 1.0e-14_dp, &
      'prepared BVEC unrestricted seasonal indicator')
   call assert_true(all(shape(prepared_vecm%holdout_y) == [2, 1]) .and. &
      all(shape(prepared_vecm%sur) == [10, 16]) .and. &
      all(shape(prepared_vecm%tvp_sur) == [10, 80]), &
      'prepared BVEC holdout and SUR shapes')
   call assert_true(all(shape(prepared_vecm%structural) == [10, 1]), &
      'prepared structural BVEC design shape')

   prior_series(1, :) = [1.0_dp, 2.0_dp, 4.0_dp, 3.0_dp, 7.0_dp, &
      5.0_dp, 8.0_dp, 6.0_dp, 10.0_dp, 9.0_dp]
   prior_series(2, :) = [2.0_dp, 1.0_dp, 3.0_dp, 6.0_dp, 4.0_dp, &
      8.0_dp, 5.0_dp, 9.0_dp, 7.0_dp, 11.0_dp]
   prior_var_data = bvartools_prepare_var(prior_series, 1, &
      include_constant=.true.)
   bvar_prior_bundle = bvartools_prepare_bvar_prior(prior_var_data, &
      ssvs_tau=[0.1_dp, 2.0_dp], inclusion_probability=0.4_dp, &
      exclude_deterministics=.true.)
   call assert_true(bvar_prior_bundle%info == 0 .and. &
      bvar_prior_bundle%use_ssvs .and. &
      all(shape(bvar_prior_bundle%initial_coefficients) == [2, 3]) .and. &
      all(shape(bvar_prior_bundle%coefficient_precision) == [6, 6]), &
      'BVAR prior-bundle status and dimensions')
   call assert_true(size(bvar_prior_bundle%selectable) == 4 .and. &
      all(bvar_prior_bundle%selectable == [1, 2, 3, 4]), &
      'BVAR prior-bundle deterministic exclusion')
   call assert_close(bvar_prior_bundle%coefficient_precision(1, 1), &
      0.25_dp, 1.0e-14_dp, 'BVAR SSVS initial slab precision')
   call assert_true(all(shape(bvar_prior_bundle%initial_covariance) == [2, 2]) .and. &
      abs(bvar_prior_bundle%covariance_df - 2.0_dp) < 1.0e-14_dp, &
      'BVAR covariance prior and initialization')
   call set_random_seed(749)
   bundled_bvar_draws = bvartools_bvar_gibbs(prior_var_data%y, &
      prior_var_data%x, bvar_prior_bundle%coefficient_mean, &
      bvar_prior_bundle%coefficient_precision, &
      bvar_prior_bundle%initial_covariance, 2, 1, &
      bvar_prior_bundle%covariance_scale, bvar_prior_bundle%covariance_df, &
      tau0=bvar_prior_bundle%tau0, tau1=bvar_prior_bundle%tau1, &
      inclusion_probability=bvar_prior_bundle%inclusion_probability, &
      selectable=bvar_prior_bundle%selectable)
   call assert_true(bundled_bvar_draws%info == 0 .and. &
      bundled_bvar_draws%retained_draws == 2, &
      'BVAR prior-bundle Gibbs compatibility')
   call set_random_seed(751)
   bvar_fit = bvartools_fit_bvar(prior_var_data, bvar_prior_bundle, 2, 1)
   call assert_true(bvar_fit%info == 0 .and. &
      bvar_fit%draws%retained_draws == 2 .and. bvar_fit%variables == 2 .and. &
      bvar_fit%lag_order == 1, 'end-to-end BVAR fit wrapper')

   prior_vecm_data = bvartools_prepare_vecm(prior_series, 2, &
      constant='unrestricted')
   bvec_prior_bundle = bvartools_prepare_bvec_prior(prior_vecm_data, 1, &
      ssvs_tau=[0.1_dp, 2.0_dp], inclusion_probability=0.4_dp, &
      exclude_deterministics=.true.)
   call assert_true(bvec_prior_bundle%info == 0 .and. &
      bvec_prior_bundle%use_ssvs .and. bvec_prior_bundle%rank == 1, &
      'BVEC prior-bundle status')
   call assert_true(all(shape(bvec_prior_bundle%initial_beta) == [2, 1]) .and. &
      all(shape(bvec_prior_bundle%initial_gamma) == [2, 3]) .and. &
      all(shape(bvec_prior_bundle%gamma_precision) == [6, 6]), &
      'BVEC prior-bundle dimensions')
   call assert_close(bvec_prior_bundle%initial_beta(1, 1), 1.0_dp, &
      1.0e-14_dp, 'BVEC normalized initial beta')
   call assert_true(size(bvec_prior_bundle%selectable) == 4 .and. &
      all(bvec_prior_bundle%selectable == [1, 2, 3, 4]), &
      'BVEC prior-bundle deterministic exclusion')
   call assert_close(bvec_prior_bundle%covariance_df, 3.0_dp, &
      1.0e-14_dp, 'BVEC rank-adjusted covariance degrees of freedom')
   call set_random_seed(750)
   bundled_bvec_draws = bvartools_bvec_gibbs(prior_vecm_data%y, &
      prior_vecm_data%w, bvec_prior_bundle%initial_beta, &
      bvec_prior_bundle%initial_covariance, 2, 1, bvec_prior_bundle%shrinkage, &
      bvec_prior_bundle%cointegration_precision, &
      bvec_prior_bundle%loading_precision, &
      covariance_prior_df=bvec_prior_bundle%covariance_df, x=prior_vecm_data%x, &
      gamma_prior_mean=bvec_prior_bundle%gamma_mean, &
      gamma_prior_precision=bvec_prior_bundle%gamma_precision, &
      tau0=bvec_prior_bundle%tau0, tau1=bvec_prior_bundle%tau1, &
      inclusion_probability=bvec_prior_bundle%inclusion_probability, &
      selectable=bvec_prior_bundle%selectable)
   call assert_true(bundled_bvec_draws%info == 0 .and. &
      bundled_bvec_draws%retained_draws == 2, &
      'BVEC prior-bundle Gibbs compatibility')
   call set_random_seed(752)
   bvec_fit = bvartools_fit_bvec(prior_vecm_data, bvec_prior_bundle, 2, 1)
   call assert_true(bvec_fit%info == 0 .and. &
      bvec_fit%draws%retained_draws == 2 .and. bvec_fit%rank == 1 .and. &
      .not. bvec_fit%rank_zero, 'end-to-end positive-rank BVEC fit wrapper')

   rank_zero_prior_bundle = bvartools_prepare_bvec_prior(prior_vecm_data, 0, &
      covariance_prior='gamma', error_shape=2.0_dp, error_rate=0.5_dp, &
      ssvs_tau=[0.1_dp, 2.0_dp])
   call assert_true(rank_zero_prior_bundle%info == 0 .and. &
      all(shape(rank_zero_prior_bundle%initial_beta) == [2, 0]), &
      'rank-zero BVEC prior bundle')
   call set_random_seed(753)
   rank_zero_fit = bvartools_fit_bvec(prior_vecm_data, &
      rank_zero_prior_bundle, 2, 1)
   call assert_true(rank_zero_fit%info == 0 .and. rank_zero_fit%rank_zero .and. &
      rank_zero_fit%rank == 0 .and. &
      rank_zero_fit%rank_zero_draws%retained_draws == 2, &
      'rank-zero gamma BVEC fit wrapper')

   tvp_bvar_prior = bvartools_prepare_tvp_bvar_prior(prior_var_data, &
      bvs_probability=0.5_dp, exclude_deterministics=.true.)
   call assert_true(tvp_bvar_prior%info == 0 .and. tvp_bvar_prior%use_bvs .and. &
      size(tvp_bvar_prior%initial_state) == 6 .and. &
      size(tvp_bvar_prior%selectable) == 4, 'TVP-BVAR prior bundle')
   call set_random_seed(754)
   tvp_bvar_fit = bvartools_fit_tvp_bvar(prior_var_data, tvp_bvar_prior, 2, 1)
   call assert_true(tvp_bvar_fit%info == 0 .and. &
      tvp_bvar_fit%draws%retained_draws == 2 .and. &
      allocated(tvp_bvar_fit%draws%included), 'TVP-BVAR Wishart and BVS fit')

   tvp_bvar_sv_prior = bvartools_prepare_tvp_bvar_prior(prior_var_data, &
      observation_prior='sv', volatility_method='ocsn')
   call assert_true(tvp_bvar_sv_prior%info == 0 .and. &
      all(shape(tvp_bvar_sv_prior%initial_log_variance) == [2, 9]), &
      'TVP-BVAR stochastic-volatility prior bundle')
   call set_random_seed(755)
   tvp_bvar_sv_fit = bvartools_fit_tvp_bvar(prior_var_data, &
      tvp_bvar_sv_prior, 2, 1)
   call assert_true(tvp_bvar_sv_fit%info == 0 .and. &
      tvp_bvar_sv_fit%draws%retained_draws == 2 .and. &
      allocated(tvp_bvar_sv_fit%draws%log_variance), &
      'TVP-BVAR stochastic-volatility fit')

   structural_var_data = bvartools_prepare_var(prior_series, 1, &
      include_constant=.true., structural=.true.)
   bvar_prior_bundle = bvartools_prepare_bvar_prior(structural_var_data)
   call set_random_seed(761)
   structural_bvar_draws = bvartools_structural_bvar_gibbs( &
      structural_var_data%y, structural_var_data%sur, &
      structural_var_data%structural, bvar_prior_bundle%coefficient_mean, &
      bvar_prior_bundle%coefficient_precision, [0.0_dp], &
      reshape([1.0_dp], [1, 1]), bvar_prior_bundle%initial_covariance, 2, 1, &
      covariance_prior_scale=bvar_prior_bundle%covariance_scale, &
      covariance_prior_df=bvar_prior_bundle%covariance_df, &
      structural_tau0=[0.05_dp], structural_tau1=[10.0_dp], &
      structural_inclusion_probability=[0.5_dp])
   call assert_true(structural_bvar_draws%info == 0 .and. &
      all(shape(structural_bvar_draws%structural) == [1, 2]) .and. &
      all(shape(structural_bvar_draws%structural_impact) == [2, 2, 2]) .and. &
      allocated(structural_bvar_draws%structural_included), &
      'constant structural BVAR draws and impact matrices')
   call set_random_seed(764)
   bvar_fit = bvartools_fit_structural_bvar(structural_var_data, &
      bvar_prior_bundle, 2, 1, structural_inclusion_probability=0.5_dp)
   call assert_true(bvar_fit%info == 0 .and. &
      bvar_fit%draws%retained_draws == 2 .and. &
      allocated(bvar_fit%draws%structural_impact), &
      'end-to-end constant structural BVAR fit')
   structural_fit_normals = 0.0_dp
   structural_fit_history(:, 1) = prior_series(:, 10)
   structural_fit_future = 1.0_dp
   structural_prediction = bvartools_bvar_predictive_from_random( &
      bvar_fit%draws, structural_fit_history, 1, structural_fit_normals, &
      future_regressors=structural_fit_future, &
      structural_impact=bvar_fit%draws%structural_impact)
   call assert_true(structural_prediction%info == 0, &
      'constant structural BVAR impacts usable by prediction')
   sir = bvartools_bvar_irf(bvar_fit%draws, 1, 1, 'sir', &
      structural_impact=bvar_fit%draws%structural_impact)
   call assert_true(sir%info == 0, &
      'constant structural BVAR impacts usable by IRF')
   sir_fevd = bvartools_bvar_fevd(bvar_fit%draws, 1, 1, 'sir', &
      structural_impact=bvar_fit%draws%structural_impact)
   call assert_true(sir_fevd%info == 0, &
      'constant structural BVAR impacts usable by FEVD')
   structural_tvp_prior = bvartools_prepare_structural_tvp_bvar_prior( &
      structural_var_data, coefficient_bvs_probability=0.5_dp, &
      structural_bvs_probability=0.4_dp, covariance_bvs_probability=0.3_dp)
   call assert_true(structural_tvp_prior%info == 0 .and. &
      structural_tvp_prior%reduced_form_states == 6 .and. &
      structural_tvp_prior%structural_states == 1 .and. &
      size(structural_tvp_prior%selectable) == 7, &
      'structural TVP-BVAR prior bundle')
   call set_random_seed(758)
   structural_tvp_fit = bvartools_fit_structural_tvp_bvar( &
      structural_var_data, structural_tvp_prior, 2, 1)
   call assert_true(structural_tvp_fit%info == 0 .and. &
      structural_tvp_fit%draws%retained_draws == 2 .and. &
      all(shape(structural_tvp_fit%reduced_form_states) == [54, 2]) .and. &
      all(shape(structural_tvp_fit%structural_states) == [9, 2]), &
      'structural TVP-BVAR split state draws')
   call assert_true(all(shape(structural_tvp_fit%structural_impact) == &
      [2, 2, 9, 2]), 'structural TVP-BVAR impact matrices')
   irf_tvp_draws%states = structural_tvp_fit%reduced_form_states
   irf_tvp_draws%time_covariance = structural_tvp_fit%draws%covariance
   irf_tvp_draws%retained_draws = structural_tvp_fit%draws%retained_draws
   irf_tvp_draws%info = structural_tvp_fit%draws%info
   tvp_irf = bvartools_tvp_bvar_irf(irf_tvp_draws, 2, 1, 9, 9, 1, 'sir', &
      structural_impact=structural_tvp_fit%structural_impact(:, :, 9, :))
   call assert_true(tvp_irf%info == 0, &
      'period-specific structural TVP-BVAR impacts usable by IRF')
   tvp_fevd = bvartools_tvp_bvar_fevd(irf_tvp_draws, 2, 1, 9, 9, 1, 'sir', &
      structural_impact=structural_tvp_fit%structural_impact(:, :, 9, :))
   call assert_true(tvp_fevd%info == 0, &
      'period-specific structural TVP-BVAR impacts usable by FEVD')
   call assert_true(allocated(structural_tvp_fit%structural_included) .and. &
      allocated(structural_tvp_fit%draws%covariance_included), &
      'structural and covariance trajectory BVS draws')

   structural_constant_prior = bvartools_prepare_structural_tvp_bvar_prior( &
      structural_var_data, time_varying_covariance=.false., &
      structural_bvs_probability=0.5_dp, covariance_ssvs_probability=0.5_dp)
   call assert_true(structural_constant_prior%info == 0 .and. &
      structural_constant_prior%use_covariance_ssvs, &
      'structural constant-covariance prior bundle')
   call set_random_seed(759)
   structural_constant_fit = bvartools_fit_structural_tvp_bvar( &
      structural_var_data, structural_constant_prior, 2, 1)
   call assert_true(structural_constant_fit%info == 0 .and. &
      .not. structural_constant_fit%draws%time_varying_covariance .and. &
      allocated(structural_constant_fit%draws%covariance_included), &
      'structural constant-covariance SSVS fit')

   tvp_bvec_prior = bvartools_prepare_tvp_bvec_prior(prior_vecm_data, 1, &
      observation_prior='gamma', bvs_probability=0.5_dp, &
      exclude_deterministics=.true., time_varying_covariance=.true., &
      covariance_bvs_probability=0.5_dp)
   call assert_true(tvp_bvec_prior%info == 0 .and. tvp_bvec_prior%use_bvs .and. &
      tvp_bvec_prior%use_covariance_state .and. &
      tvp_bvec_prior%use_covariance_bvs, 'TVP-BVEC covariance-state prior bundle')
   call set_random_seed(756)
   tvp_bvec_fit = bvartools_fit_tvp_bvec(prior_vecm_data, tvp_bvec_prior, 2, 1)
   call assert_true(tvp_bvec_fit%info == 0 .and. &
      tvp_bvec_fit%draws%retained_draws == 2 .and. &
      allocated(tvp_bvec_fit%draws%covariance_state) .and. &
      allocated(tvp_bvec_fit%draws%covariance_included), &
      'TVP-BVEC covariance-state and BVS fit')

   structural_vecm_data = bvartools_prepare_vecm(prior_series, 2, &
      constant='unrestricted', structural=.true.)
   bvec_prior_bundle = bvartools_prepare_bvec_prior(structural_vecm_data, 1, &
      covariance_prior='gamma', error_shape=2.0_dp, error_rate=0.5_dp)
   call set_random_seed(762)
   structural_bvec_draws = bvartools_bvec_gibbs(structural_vecm_data%y, &
      structural_vecm_data%w, bvec_prior_bundle%initial_beta, &
      bvec_prior_bundle%initial_covariance, 2, 1, bvec_prior_bundle%shrinkage, &
      bvec_prior_bundle%cointegration_precision, &
      bvec_prior_bundle%loading_precision, x=structural_vecm_data%x, &
      gamma_prior_mean=bvec_prior_bundle%gamma_mean, &
      gamma_prior_precision=bvec_prior_bundle%gamma_precision, &
      gamma_shape=bvec_prior_bundle%error_shape, &
      gamma_rate=bvec_prior_bundle%error_rate, &
      structural_design=structural_vecm_data%structural, &
      structural_prior_mean=[0.0_dp], &
      structural_prior_precision=reshape([1.0_dp], [1, 1]), &
      structural_tau0=[0.05_dp], structural_tau1=[10.0_dp], &
      structural_inclusion_probability=[0.5_dp])
   call assert_true(structural_bvec_draws%info == 0 .and. &
      all(shape(structural_bvec_draws%structural) == [1, 2]) .and. &
      all(shape(structural_bvec_draws%structural_impact) == [2, 2, 2]) .and. &
      allocated(structural_bvec_draws%structural_included), &
      'constant structural BVEC draws and SSVS indicators')
   call set_random_seed(765)
   bvec_fit = bvartools_fit_structural_bvec(structural_vecm_data, &
      bvec_prior_bundle, 2, 1, structural_inclusion_probability=0.5_dp)
   call assert_true(bvec_fit%info == 0 .and. &
      bvec_fit%draws%retained_draws == 2 .and. &
      allocated(bvec_fit%draws%structural_impact), &
      'end-to-end constant structural BVEC fit')
   tvp_bvec_prior = bvartools_prepare_tvp_bvec_prior(structural_vecm_data, 1, &
      observation_prior='gamma', bvs_probability=0.5_dp, &
      time_varying_covariance=.true., covariance_bvs_probability=0.5_dp)
   call set_random_seed(760)
   structural_tvp_bvec_fit = bvartools_fit_structural_tvp_bvec( &
      structural_vecm_data, tvp_bvec_prior, 2, 1, &
      structural_bvs_probability=0.5_dp)
   call assert_true(structural_tvp_bvec_fit%info == 0 .and. &
      structural_tvp_bvec_fit%draws%retained_draws == 2 .and. &
      all(shape(structural_tvp_bvec_fit%draws%structural) == [8, 2]) .and. &
      allocated(structural_tvp_bvec_fit%draws%structural_included), &
      'structural TVP-BVEC state and BVS draws')
   call assert_true(all(shape(structural_tvp_bvec_fit%draws%structural_impact) == &
      [2, 2, 8, 2]), 'structural TVP-BVEC impact matrices')
   tvp_bvec_prior = bvartools_prepare_tvp_bvec_prior(structural_vecm_data, 1, &
      observation_prior='gamma', time_varying_covariance=.true.)
   call set_random_seed(763)
   structural_tvp_bvec_fit = bvartools_fit_structural_tvp_bvec( &
      structural_vecm_data, tvp_bvec_prior, 2, 1)
   call assert_true(structural_tvp_bvec_fit%info == 0 .and. &
      structural_tvp_bvec_fit%draws%retained_draws == 2, &
      'unselected structural TVP-BVEC fit')

   tvp_bvec_prior = bvartools_prepare_tvp_bvec_prior(structural_vecm_data, 1, &
      observation_prior='gamma')
   call set_random_seed(766)
   structural_tvp_bvec_fit = bvartools_fit_structural_tvp_bvec( &
      structural_vecm_data, tvp_bvec_prior, 2, 1)
   call assert_true(structural_tvp_bvec_fit%info == 0 .and. &
      allocated(structural_tvp_bvec_fit%draws%structural_impact) .and. &
      .not. allocated(structural_tvp_bvec_fit%draws%covariance_state), &
      'structural TVP-BVEC diagonal gamma fit')

   tvp_bvec_prior = bvartools_prepare_tvp_bvec_prior(structural_vecm_data, 1, &
      observation_prior='wishart')
   call set_random_seed(767)
   structural_tvp_bvec_fit = bvartools_fit_structural_tvp_bvec( &
      structural_vecm_data, tvp_bvec_prior, 2, 1)
   call assert_true(structural_tvp_bvec_fit%info == 0 .and. &
      allocated(structural_tvp_bvec_fit%draws%structural_impact), &
      'structural TVP-BVEC inverse-Wishart fit')

   tvp_bvec_prior = bvartools_prepare_tvp_bvec_prior(structural_vecm_data, 1, &
      observation_prior='sv', volatility_method='ksc')
   call set_random_seed(768)
   structural_tvp_bvec_fit = bvartools_fit_structural_tvp_bvec( &
      structural_vecm_data, tvp_bvec_prior, 2, 1)
   call assert_true(structural_tvp_bvec_fit%info == 0 .and. &
      allocated(structural_tvp_bvec_fit%draws%structural_impact) .and. &
      allocated(structural_tvp_bvec_fit%draws%log_variance), &
      'structural TVP-BVEC stochastic-volatility fit')

   tvp_bvec_sv_prior = bvartools_prepare_tvp_bvec_prior(prior_vecm_data, 1, &
      observation_prior='sv', volatility_method='ksc')
   call assert_true(tvp_bvec_sv_prior%info == 0 .and. &
      all(shape(tvp_bvec_sv_prior%initial_log_variance) == [2, 8]), &
      'TVP-BVEC stochastic-volatility prior bundle')
   call set_random_seed(757)
   tvp_bvec_sv_fit = bvartools_fit_tvp_bvec(prior_vecm_data, &
      tvp_bvec_sv_prior, 2, 1)
   call assert_true(tvp_bvec_sv_fit%info == 0 .and. &
      tvp_bvec_sv_fit%draws%retained_draws == 2 .and. &
      allocated(tvp_bvec_sv_fit%draws%log_variance), &
      'TVP-BVEC stochastic-volatility fit')

   prepared_dfm = bvartools_prepare_dfm(raw_series, [0, 2], [1, 2], &
      iterations=100, burnin=20)
   call assert_true(prepared_dfm%info == 0 .and. &
      all(shape(prepared_dfm%x) == [2, 8]), 'prepared DFM matrix shape')
   call assert_close(prepared_dfm%mean(1), 4.5_dp, 1.0e-14_dp, &
      'prepared DFM series mean')
   call assert_close(prepared_dfm%standard_deviation(1), sqrt(6.0_dp), &
      1.0e-14_dp, 'prepared DFM sample standard deviation')
   call assert_close(prepared_dfm%x(1, 1), -3.5_dp/sqrt(6.0_dp), &
      1.0e-14_dp, 'prepared DFM standardized value')
   call assert_true(all(prepared_dfm%factor_count == [1, 1, 2, 2]) .and. &
      all(prepared_dfm%lag_order == [0, 2, 0, 2]), &
      'prepared DFM model-grid order')
   call assert_true(prepared_dfm%iterations == 100 .and. &
      prepared_dfm%burnin == 20 .and. prepared_dfm%observations == 8 .and. &
      prepared_dfm%variables == 2, 'prepared DFM metadata')

   dfm_prior = bvartools_dfm_prior(2, 1, 2)
   call assert_true(dfm_prior%info == 0 .and. &
      all(shape(dfm_prior%loading_precision) == [1, 1]) .and. &
      all(shape(dfm_prior%transition_precision) == [2, 2]), &
      'DFM prior dimensions')
   call assert_close(dfm_prior%loading_precision(1, 1), 0.01_dp, &
      1.0e-14_dp, 'DFM default loading precision')
   call assert_close(dfm_prior%measurement_shape(1), 5.0_dp, &
      1.0e-14_dp, 'DFM default measurement shape')
   call assert_close(dfm_prior%measurement_rate(1), 4.0_dp, &
      1.0e-14_dp, 'DFM default measurement rate')
   call assert_close(dfm_prior%transition_precision(1, 1), 0.01_dp, &
      1.0e-14_dp, 'DFM default transition precision')

   dfm_grid_data = bvartools_prepare_dfm(raw_series, [0, 1], [1], &
      iterations=2, burnin=1)
   call set_random_seed(748)
   dfm_grid_draws = bvartools_dfm_grid_gibbs(dfm_grid_data)
   call assert_true(dfm_grid_draws%info == 0 .and. &
      dfm_grid_draws%failed_model == 0 .and. &
      size(dfm_grid_draws%draws) == 2, 'DFM model-grid Gibbs status')
   call assert_true(all(dfm_grid_draws%factor_count == [1, 1]) .and. &
      all(dfm_grid_draws%lag_order == [0, 1]), 'DFM fitted model-grid metadata')
   call assert_true(dfm_grid_draws%draws(1)%retained_draws == 2 .and. &
      dfm_grid_draws%draws(2)%retained_draws == 2, &
      'DFM model-grid retained draws')
   call assert_true(all(shape(dfm_grid_draws%draws(1)%transition) == [0, 2]) .and. &
      all(shape(dfm_grid_draws%draws(2)%transition) == [1, 2]), &
      'DFM model-grid lag-dependent transition shapes')

   comparison_residual(:, :, 1) = 0.0_dp
   comparison_residual(:, :, 2) = 1.0_dp
   comparison_covariance = 1.0_dp
   comparison_tvp_covariance = 1.0_dp
   expected_log_likelihood = -log(2.0_dp*acos(-1.0_dp)) - 0.5_dp
   model_comparison = bvartools_model_comparison(comparison_residual, &
      comparison_covariance, 1)
   call assert_true(model_comparison%info == 0 .and. &
      model_comparison%observations == 2 .and. model_comparison%draws == 2, &
      'constant-covariance model-comparison status')
   call assert_close(model_comparison%log_likelihood, expected_log_likelihood, &
      1.0e-13_dp, 'posterior-mean multivariate log likelihood')
   call assert_close(model_comparison%aic, 2.0_dp - &
      2.0_dp*expected_log_likelihood, 1.0e-13_dp, 'posterior model AIC')
   call assert_close(model_comparison%bic, log(2.0_dp) - &
      2.0_dp*expected_log_likelihood, 1.0e-13_dp, 'posterior model BIC')
   tvp_model_comparison = bvartools_model_comparison(comparison_residual, &
      comparison_tvp_covariance, 1)
   call assert_close(tvp_model_comparison%log_likelihood, &
      model_comparison%log_likelihood, 1.0e-13_dp, &
      'time-varying covariance model likelihood')

   comparison_data(1)%residual = comparison_residual
   comparison_data(1)%covariance = reshape(comparison_covariance, [1, 1, 1, 2])
   comparison_data(1)%parameter_count = 1
   comparison_data(2)%residual = comparison_residual
   comparison_data(2)%covariance = comparison_tvp_covariance
   comparison_data(2)%parameter_count = 2
   comparison_set = bvartools_compare_models(comparison_data)
   call assert_true(comparison_set%info == 0 .and. &
      size(comparison_set%model) == 2, 'posterior model-collection status')
   call assert_true(comparison_set%best_aic == 1 .and. &
      comparison_set%best_bic == 1 .and. comparison_set%best_hq == 2, &
      'posterior information-criterion minimizing models')

   print '(a)', 'All bvartools tests passed.'

contains

   subroutine assert_true(condition, label)
      !! Stop when a logical test condition is false.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: label !! Failure label.

      if (.not. condition) error stop 'FAIL: '//label
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, label)
      !! Stop when two scalar values differ beyond tolerance.
      real(dp), intent(in) :: actual !! Computed value.
      real(dp), intent(in) :: expected !! Reference value.
      real(dp), intent(in) :: tolerance !! Absolute tolerance.
      character(len=*), intent(in) :: label !! Failure label.

      if (abs(actual - expected) > tolerance) error stop 'FAIL: '//label
   end subroutine assert_close

end program test_bvartools
