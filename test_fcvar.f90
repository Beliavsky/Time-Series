! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the FCVAR translation.
program test_fcvar
   use kind_mod, only: dp
   use fcvar_mod
   use time_series_random_mod, only: set_random_seed
   use time_series_diagnostics_mod, only: multivariate_white_noise_test_t, &
      multivariate_white_noise_test
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   type(fcvar_options_t) :: options, level_options
   type(fcvar_transform_t) :: transformed
   type(fcvar_fit_t) :: fit, zero_rank, full_rank, restricted_fit, unrestricted_fit
   type(fcvar_fit_t) :: level_fit, shifted_level_fit, fixed_level_fit
   type(fcvar_grid_t) :: grid, local_grid
   type(fcvar_estimation_t) :: estimation, equal_estimation
   type(fcvar_path_t) :: forecast, simulation, bootstrap, random_path
   type(fcvar_lr_test_t) :: lr_test
   type(fcvar_rank_tests_t) :: rank_tests
   type(fcvar_lag_selection_t) :: lag_selection
   type(fcvar_bootstrap_rank_t) :: bootstrap_rank, random_bootstrap_rank
   type(fcvar_standard_errors_t) :: standard_errors
   type(fcvar_roots_t) :: roots
   type(fcvar_restrictions_t) :: restrictions
   type(fcvar_restricted_estimation_t) :: restricted_estimation
   type(fcvar_switching_estimation_t) :: switching_estimation
   type(fcvar_order_estimation_t) :: constrained_orders
   type(fcvar_bootstrap_hypothesis_t) :: bootstrap_hypothesis
   type(fcvar_bootstrap_hypothesis_t) :: random_bootstrap_hypothesis
   type(fcvar_fit_t) :: unpacked_fit
   type(multivariate_white_noise_test_t) :: white_noise
   real(dp) :: x(12, 2), shifted_x(12, 2), expected_fd(12, 2)
   real(dp) :: expected_alpha(2)
   real(dp) :: expected_beta(2), expected_covariance(2, 2)
   real(dp) :: expected_gamma(2, 2)
   real(dp), allocatable :: differenced(:, :), lagged(:, :)
   real(dp), allocatable :: packed_parameters(:)
   real(dp) :: supplied_innovations(3, 2), expected_path(3, 2)
   real(dp) :: rank_signs(10, 3)
   real(dp) :: synthetic_surface(4, 3)
   logical :: synthetic_feasible(4, 3)
   integer, allocatable :: local_rows(:), local_columns(:)

   x = reshape([ &
      1.0_dp, 1.3_dp, 1.1_dp, 0.8_dp, 0.6_dp, 0.5_dp, &
      0.2_dp, 0.1_dp, -0.1_dp, -0.2_dp, -0.15_dp, -0.3_dp, &
      2.0_dp, 1.7_dp, 1.4_dp, 1.2_dp, 1.0_dp, 0.7_dp, &
      0.6_dp, 0.3_dp, 0.2_dp, 0.0_dp, -0.1_dp, -0.05_dp], [12, 2])
   expected_fd = reshape([ &
      1.0_dp, 0.7_dp, 0.2_dp, -0.072_dp, -0.1184_dp, -0.084128_dp, &
      -0.3002176_dp, -0.20032_dp, -0.296012544_dp, -0.2432754688_dp, &
      -0.091706283008_dp, -0.2374742153216_dp, &
      2.0_dp, 0.5_dp, 0.14_dp, 0.044_dp, -0.0504_dp, -0.225216_dp, &
      -0.126592_dp, -0.32664192_dp, -0.218328832_dp, -0.3087525888_dp, &
      -0.252432807936_dp, -0.09830102528_dp], [12, 2])

   differenced = fcvar_fractional_difference(x, 0.6_dp)
   call check(maxval(abs(differenced - expected_fd)) < 2.0e-13_dp, &
      'level-preserving fractional difference')
   lagged = fcvar_fractional_lags(x, 0.4_dp, 2)
   call check(all(shape(lagged) == [12, 4]), 'fractional lag shape')
   call check(maxval(abs(lagged(1, :))) < 2.0e-13_dp, 'fractional lag initial row')
   call check(abs(lagged(12, 4) - 0.2530548764672_dp) < 2.0e-13_dp, &
      'second fractional lag power')

   options%initial_values = 2
   transformed = fcvar_transform(x, 1, 0.6_dp, 0.4_dp, options)
   call check(transformed%info == 0, 'transformation status')
   call check(all(shape(transformed%z0) == [10, 2]), 'transformation shape')
   call check(maxval(abs(transformed%z0 - expected_fd(3:, :))) < 2.0e-13_dp, &
      'transformed response')
   call check(abs(transformed%z1(1, 1) - 0.56_dp) < 2.0e-13_dp, &
      'transformed levels')
   call check(abs(transformed%z2(10, 2) + 0.14830102528_dp) < 2.0e-13_dp, &
      'transformed lag')

   fit = fcvar_estimate_fixed(x, 1, 1, 0.6_dp, 0.4_dp, options)
   expected_alpha = [1.14464886240144_dp, 0.92552452198536_dp]
   expected_beta = [1.0_dp, -1.15514226977979_dp]
   expected_covariance = reshape([ &
      0.00496844105003716_dp, -0.00100102895686122_dp, &
      -0.00100102895686122_dp, 0.00114282276243961_dp], [2, 2])
   expected_gamma = reshape([ &
      -1.20590412163072_dp, 1.09098288393801_dp, &
      2.20575792390292_dp, 0.0747900879978058_dp], [2, 2])
   call check(fit%info == 0, 'fixed FCVAR estimation status')
   call check(maxval(abs(fit%alpha(:, 1) - expected_alpha)) < 2.0e-10_dp, &
      'adjustment coefficients')
   call check(maxval(abs(fit%beta(:, 1) - expected_beta)) < 2.0e-10_dp, &
      'identified cointegrating vector')
   call check(maxval(abs(fit%covariance - expected_covariance)) < 2.0e-10_dp, &
      'innovation covariance')
   call check(maxval(abs(fit%gamma - expected_gamma)) < 2.0e-10_dp, &
      'short-run coefficients')
   call check(abs(fit%log_likelihood - 32.9865791474851_dp) < 2.0e-9_dp, &
      'concentrated log likelihood')
   call check(maxval(abs(fit%residuals(:, 1))) > 0.0_dp, 'estimated residuals')

   level_options = options
   level_options%level_parameter = .true.
   shifted_x = x + spread([3.0_dp, -2.0_dp], 1, size(x, 1))
   level_fit = fcvar_estimate_fixed(x, 1, 2, 0.6_dp, 0.4_dp, level_options)
   shifted_level_fit = fcvar_estimate_fixed(shifted_x, 1, 2, 0.6_dp, 0.4_dp, &
      level_options)
   call check(level_fit%info == 0 .and. shifted_level_fit%info == 0 .and. &
      size(level_fit%level) == 2, 'profiled FCVAR level status')
   call check(maxval(abs(shifted_level_fit%level - level_fit%level - &
      [3.0_dp, -2.0_dp])) < 2.0e-4_dp, 'profiled level translation')
   call check(abs(shifted_level_fit%log_likelihood - &
      level_fit%log_likelihood) < 2.0e-6_dp .and. &
      maxval(abs(shifted_level_fit%alpha - level_fit%alpha)) < 2.0e-8_dp .and. &
      maxval(abs(shifted_level_fit%beta - level_fit%beta)) < 2.0e-4_dp, &
      'level-profile likelihood translation invariance')
   fixed_level_fit = fcvar_estimate_fixed(shifted_x, 1, 1, 0.6_dp, 0.4_dp, &
      level_options, [3.0_dp, -2.0_dp])
   call check(fixed_level_fit%info == 0 .and. &
      abs(fixed_level_fit%log_likelihood - fit%log_likelihood) < 2.0e-11_dp, &
      'supplied FCVAR level equivalence')
   packed_parameters = fcvar_pack_parameters(level_fit, .false.)
   unpacked_fit = fcvar_unpack_parameters(packed_parameters, level_fit, .false.)
   call check(size(packed_parameters) == 12 .and. unpacked_fit%info == 0 .and. &
      maxval(abs(unpacked_fit%level - level_fit%level)) < 1.0e-14_dp, &
      'level parameter packing')
   call check(fcvar_free_parameter_count(2, 1, 2, level_options, .false.) == &
      12, 'level parameter count')
   call check(abs(fcvar_full_log_likelihood(x, level_fit, level_options) - &
      level_fit%log_likelihood) < 2.0e-10_dp, 'level-aware full likelihood')
   forecast = fcvar_forecast(shifted_x, shifted_level_fit, 2, level_options)
   simulation = fcvar_forecast(x, level_fit, 2, level_options)
   call check(forecast%info == 0 .and. simulation%info == 0 .and. &
      maxval(abs(forecast%series - simulation%series - &
      spread([3.0_dp, -2.0_dp], 1, 2))) < 2.0e-4_dp, &
      'automatic forecast level restoration')

   zero_rank = fcvar_estimate_fixed(x, 1, 0, 0.6_dp, 0.4_dp, options)
   full_rank = fcvar_estimate_fixed(x, 1, 2, 0.6_dp, 0.4_dp, options)
   call check(zero_rank%info == 0, 'zero-rank estimation')
   call check(full_rank%info == 0, 'full-rank estimation')
   call check(full_rank%log_likelihood >= fit%log_likelihood - 1.0e-9_dp, &
      'nested rank likelihood')
   call check(abs(fcvar_log_likelihood_fixed(x, 1, 1, 0.6_dp, 0.4_dp, &
      options) - fit%log_likelihood) < 1.0e-13_dp, 'likelihood convenience function')

   options%restricted_constant = .true.
   restricted_fit = fcvar_estimate_fixed(x, 1, 1, 0.6_dp, 0.4_dp, options)
   call check(restricted_fit%info == 0, 'restricted-constant estimation')
   call check(abs(restricted_fit%rho(1) - 0.113411162484186_dp) < 2.0e-10_dp, &
      'restricted constant')
   call check(abs(restricted_fit%log_likelihood - 33.5587692132949_dp) < &
      2.0e-9_dp, 'restricted-constant likelihood')

   options%restricted_constant = .false.
   options%unrestricted_constant = .true.
   unrestricted_fit = fcvar_estimate_fixed(x, 1, 1, 0.6_dp, 0.4_dp, options)
   call check(unrestricted_fit%info == 0, 'unrestricted-constant estimation')
   call check(maxval(abs(unrestricted_fit%unrestricted_constant - &
      [-0.478556750046639_dp, 0.00450523485517007_dp])) < 2.0e-10_dp, &
      'unrestricted constant')
   call check(abs(unrestricted_fit%log_likelihood - 35.6798115350941_dp) < &
      2.0e-9_dp, 'unrestricted-constant likelihood')

   options%unrestricted_constant = .false.
   grid = fcvar_likelihood_grid(x, 1, 1, &
      [0.4_dp, 0.5_dp, 0.6_dp, 0.7_dp, 0.8_dp], &
      [0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp], options, .false.)
   call check(grid%info == 0, 'likelihood grid status')
   call check(abs(grid%best_d - 0.4_dp) < 1.0e-14_dp .and. &
      abs(grid%best_b - 0.2_dp) < 1.0e-14_dp, 'likelihood grid maximum')
   call check(abs(grid%best_log_likelihood - 38.4424656478953_dp) < &
      2.0e-9_dp, 'likelihood grid reference')
   call check(abs(grid%global_log_likelihood - grid%best_log_likelihood) < &
      1.0e-14_dp .and. size(grid%local_d) > 0 .and. &
      .not. grid%used_local_maximum, 'global grid diagnostics')

   synthetic_surface = reshape([10.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 9.0_dp], [4, 3])
   synthetic_feasible = .true.
   call fcvar_find_local_maxima(synthetic_surface, synthetic_feasible, &
      local_rows, local_columns)
   call check(size(local_rows) == 2 .and. all(local_rows == [1, 4]) .and. &
      all(local_columns == [1, 3]), 'strict grid local maxima')
   options%prefer_high_b_local_max = .true.
   local_grid = fcvar_likelihood_grid(x, 1, 1, &
      [0.4_dp, 0.5_dp, 0.6_dp, 0.7_dp, 0.8_dp], &
      [0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp], options, .false.)
   call check(local_grid%info == 0 .and. local_grid%used_local_maximum .and. &
      local_grid%best_b == maxval(local_grid%local_b), &
      'highest-b local grid selection')
   call check(abs(local_grid%global_log_likelihood - &
      grid%global_log_likelihood) < 1.0e-13_dp, &
      'local selection retains global diagnostics')
   options%prefer_high_b_local_max = .false.

   level_options = options
   level_options%level_parameter = .true.
   local_grid = fcvar_likelihood_grid(x, 1, 2, [0.5_dp, 0.6_dp], &
      [0.3_dp, 0.4_dp], level_options, .true.)
   call check(local_grid%info == 0 .and. &
      all(shape(local_grid%level) == [2, 2, 2]) .and. &
      all(ieee_is_finite(local_grid%level)), 'stored grid level profiles')

   estimation = fcvar_estimate(x, 1, 1, [0.4_dp, 0.8_dp], &
      [0.2_dp, 0.6_dp], options, .false., 5, 100)
   call check(estimation%info == 0, 'fractional-order estimation status')
   call check(estimation%fit%log_likelihood >= grid%best_log_likelihood - &
      1.0e-10_dp, 'refinement preserves grid maximum')
   call check(estimation%fit%log_likelihood > 38.46_dp, &
      'refinement improves the boundary grid point')
   call check(abs(estimation%fit%b - 0.224476582546426_dp) < 2.0e-3_dp, &
      'refined fractional cointegration order')
   call check(estimation%fit%d >= 0.4_dp .and. estimation%fit%d <= 0.8_dp, &
      'estimated d bounds')
   call check(estimation%fit%b >= 0.2_dp .and. estimation%fit%b <= 0.6_dp, &
      'estimated b bounds')

   equal_estimation = fcvar_estimate_equal_orders(x, 1, 1, &
      [0.2_dp, 0.8_dp], options, 7)
   call check(equal_estimation%info == 0, 'equal-order estimation status')
   call check(abs(equal_estimation%fit%d - equal_estimation%fit%b) < &
      1.0e-14_dp, 'equal fractional orders')
   call check(abs(equal_estimation%fit%d - 0.339604099274659_dp) < &
      2.0e-7_dp, 'equal-order estimate')
   call check(abs(equal_estimation%fit%log_likelihood - &
      38.6624396730209_dp) < 2.0e-8_dp, 'equal-order likelihood')

   restrictions = fcvar_restrictions_t()
   allocate(restrictions%order_matrix(1, 2), restrictions%order_value(1))
   restrictions%order_matrix = reshape([1.0_dp, -1.0_dp], [1, 2])
   restrictions%order_value = 0.0_dp
   constrained_orders = fcvar_estimate_constrained_orders(x, 1, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), restrictions, &
      options, .true., 7, 200)
   call check(constrained_orders%info == 0 .and. &
      constrained_orders%equality_rank == 1 .and. &
      constrained_orders%free_dimensions == 1, &
      'one-dimensional order restriction status')
   call check(abs(constrained_orders%fit%d - constrained_orders%fit%b) < &
      2.0e-12_dp .and. abs(constrained_orders%fit%d - &
      equal_estimation%fit%d) < 3.0e-7_dp, 'general d-equals-b restriction')

   restrictions%order_matrix = reshape([1.0_dp, -2.0_dp], [1, 2])
   constrained_orders = fcvar_estimate_constrained_orders(x, 1, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), restrictions, &
      options, .true., 7, 200)
   call check(constrained_orders%info == 0 .and. &
      abs(constrained_orders%fit%d - 2.0_dp*constrained_orders%fit%b) < &
      2.0e-12_dp, 'general affine order equality')

   restrictions = fcvar_restrictions_t()
   allocate(restrictions%order_matrix(2, 2), restrictions%order_value(2))
   restrictions%order_matrix = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
      [2, 2])
   restrictions%order_value = [0.6_dp, 0.4_dp]
   constrained_orders = fcvar_estimate_constrained_orders(x, 1, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), restrictions, &
      options, .true., 7, 200)
   call check(constrained_orders%info == 0 .and. &
      constrained_orders%free_dimensions == 0 .and. &
      abs(constrained_orders%fit%log_likelihood - fit%log_likelihood) < &
      2.0e-11_dp, 'exactly identified fractional orders')

   restrictions = fcvar_restrictions_t()
   allocate(restrictions%order_inequality_matrix(1, 2))
   allocate(restrictions%order_inequality_value(1))
   restrictions%order_inequality_matrix = reshape([0.0_dp, 1.0_dp], [1, 2])
   restrictions%order_inequality_value = 0.4_dp
   constrained_orders = fcvar_estimate_constrained_orders(x, 1, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), restrictions, &
      options, .true., 9, 200)
   call check(constrained_orders%info == 0 .and. &
      constrained_orders%free_dimensions == 2 .and. &
      constrained_orders%fit%b >= 0.4_dp - 2.0e-10_dp .and. &
      constrained_orders%fit%b <= constrained_orders%fit%d + 2.0e-10_dp, &
      'inequality-constrained fractional orders')
   call check(abs(constrained_orders%fit%b - 0.4_dp) < 2.0e-7_dp, &
      'active fractional-order inequality')

   restrictions = fcvar_restrictions_t()
   allocate(restrictions%order_matrix(1, 2), restrictions%order_value(1))
   restrictions%order_matrix = reshape([1.0_dp, 0.0_dp], [1, 2])
   restrictions%order_value = 0.1_dp
   constrained_orders = fcvar_estimate_constrained_orders(x, 1, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), restrictions, &
      options, .true., 7, 200)
   call check(constrained_orders%info /= 0, &
      'infeasible fractional-order restriction')

   packed_parameters = fcvar_pack_parameters(fit, .false.)
   unpacked_fit = fcvar_unpack_parameters(packed_parameters, fit, .false.)
   call check(size(packed_parameters) == 9 .and. unpacked_fit%info == 0, &
      'identified parameter packing')
   call check(maxval(abs(unpacked_fit%alpha - fit%alpha)) < 1.0e-14_dp .and. &
      maxval(abs(unpacked_fit%beta - fit%beta)) < 1.0e-14_dp .and. &
      maxval(abs(unpacked_fit%gamma - fit%gamma)) < 1.0e-14_dp, &
      'identified parameter round trip')
   call check(abs(fcvar_full_log_likelihood(x, fit, options) - &
      fit%log_likelihood) < 2.0e-11_dp, 'full FCVAR likelihood')
   packed_parameters = fcvar_pack_parameters(restricted_fit, .false.)
   unpacked_fit = fcvar_unpack_parameters(packed_parameters, restricted_fit, &
      .false.)
   call check(size(packed_parameters) == 10 .and. &
      maxval(abs(unpacked_fit%rho - restricted_fit%rho)) < 1.0e-14_dp, &
      'restricted-constant parameter packing')
   packed_parameters = fcvar_pack_parameters(unrestricted_fit, .false.)
   unpacked_fit = fcvar_unpack_parameters(packed_parameters, unrestricted_fit, &
      .false.)
   call check(size(packed_parameters) == 11 .and. &
      maxval(abs(unpacked_fit%unrestricted_constant - &
      unrestricted_fit%unrestricted_constant)) < 1.0e-14_dp, &
      'unrestricted-constant parameter packing')

   standard_errors = fcvar_standard_errors(x, equal_estimation%fit, options, &
      .true., 1.0e-4_dp)
   call check(standard_errors%info == 0 .and. &
      size(standard_errors%standard_error) == 8, &
      'numerical FCVAR standard-error status')
   call check(maxval(abs(standard_errors%hessian - &
      transpose(standard_errors%hessian))) < 1.0e-10_dp, &
      'symmetric numerical Hessian')
   call check(all(ieee_is_finite(standard_errors%standard_error)) .and. &
      all(standard_errors%standard_error > 0.0_dp), &
      'finite positive standard errors')
   call check(maxval(abs(standard_errors%standard_error - &
      [0.0338427568473701_dp, 0.772678475310039_dp, &
      0.276944245747527_dp, 0.406613562093237_dp, &
      0.901030124139693_dp, 0.505522393725017_dp, &
      0.614694062326895_dp, 0.3228782978111_dp])) < 3.0e-4_dp, &
      'numerical standard-error reference')
   call check(abs(standard_errors%d - standard_errors%b) < 1.0e-14_dp .and. &
      abs(standard_errors%beta(1, 1)) < 1.0e-14_dp .and. &
      standard_errors%beta(2, 1) > 0.0_dp, &
      'mapped identified coefficient errors')

   roots = fcvar_characteristic_roots(fit)
   call check(roots%info == 0 .and. size(roots%roots) == 4, &
      'FCVAR characteristic-root status')
   call check(abs(roots%roots(1) - cmplx(1.0_dp, 0.0_dp, dp)) < 2.0e-11_dp, &
      'FCVAR unit root')
   call check(abs(real(roots%roots(2), dp) - 0.815524811504652_dp) < &
      2.0e-11_dp .and. abs(abs(aimag(roots%roots(2))) - &
      0.392932416982437_dp) < 2.0e-11_dp, 'FCVAR complex roots')
   call check(abs(roots%roots(4) - cmplx(-0.488774559664528_dp, &
      0.0_dp, dp)) < 2.0e-11_dp, 'FCVAR final characteristic root')
   call check(.not. roots%outside_unit_circle .and. &
      abs(roots%minimum_root_modulus - 0.488774559664528_dp) < 2.0e-11_dp, &
      'FCVAR root modulus diagnostic')

   restrictions = fcvar_restrictions_t()
   allocate(restrictions%order_matrix(1, 2), restrictions%order_value(1))
   restrictions%order_matrix = reshape([1.0_dp, 0.0_dp], [1, 2])
   restrictions%order_value = 0.4_dp
   restricted_estimation = fcvar_restricted_estimate(x, 1, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), restrictions, &
      options, .true., .true., 7, 300)
   call check(restricted_estimation%info == 0 .and. &
      restricted_estimation%restriction_rank == 1 .and. &
      restricted_estimation%test%degrees_of_freedom == 1, &
      'fractional-order restriction status')
   call check(abs(restricted_estimation%restricted_fit%d - 0.4_dp) < &
      1.0e-13_dp .and. abs(restricted_estimation%restricted_fit%b - &
      0.4_dp) < 1.0e-13_dp, 'fractional-order restriction')
   call check(abs(restricted_estimation%restricted_fit%log_likelihood - &
      37.1943624530828_dp) < 3.0e-7_dp .and. &
      abs(restricted_estimation%test%statistic - 2.93615443987596_dp) < &
      6.0e-7_dp .and. abs(restricted_estimation%test%p_value - &
      0.0866167214107244_dp) < 3.0e-8_dp, 'restricted-model LR reference')
   restrictions = fcvar_restrictions_t()
   allocate(restrictions%alpha_matrix(1, 2), restrictions%alpha_value(1))
   allocate(restrictions%beta_matrix(1, 2), restrictions%beta_value(1))
   restrictions%alpha_matrix = reshape([0.0_dp, 1.0_dp], [1, 2])
   restrictions%alpha_value = 0.0_dp
   restrictions%beta_matrix = reshape([0.0_dp, 1.0_dp], [1, 2])
   restrictions%beta_value = -1.0_dp
   switching_estimation = fcvar_estimate_restricted_fixed(x, 1, 1, 0.6_dp, &
      0.4_dp, restrictions, options, 300, 1.0e-9_dp)
   call check(switching_estimation%info == 0 .and. &
      switching_estimation%converged, 'restricted switching status')
   call check(abs(switching_estimation%fit%alpha(2, 1)) < 1.0e-11_dp .and. &
      abs(switching_estimation%fit%beta(2, 1) + 1.0_dp) < 1.0e-11_dp, &
      'restricted switching constraints')
   call check(all(ieee_is_finite(switching_estimation%fit%covariance)) .and. &
      ieee_is_finite(switching_estimation%fit%log_likelihood), &
      'restricted switching likelihood')
   call check(abs(switching_estimation%fit%log_likelihood - &
      17.6782175821406_dp) < 1.0e-8_dp, 'restricted switching reference')
   call check(maxval(abs(switching_estimation%fit%covariance - &
      matmul(transpose(switching_estimation%fit%residuals), &
      switching_estimation%fit%residuals)/ &
      real(size(switching_estimation%fit%residuals, 1), dp))) < 2.0e-13_dp, &
      'restricted switching residual covariance')
   restricted_estimation = fcvar_restricted_estimate(x, 1, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), restrictions, &
      options, .true., .true., 7, 300)
   call check(restricted_estimation%info == 0 .and. &
      restricted_estimation%restriction_rank == 2, &
      'alpha-beta restriction status')
   call check(abs(restricted_estimation%restricted_fit%alpha(2, 1)) < &
      1.0e-13_dp .and. abs(restricted_estimation%restricted_fit%beta(2, 1) + &
      1.0_dp) < 1.0e-13_dp, 'exact alpha-beta restrictions')
   call check(abs(restricted_estimation%restricted_fit%log_likelihood - &
      27.3998623378953_dp) < 5.0e-7_dp .and. &
      abs(restricted_estimation%test%statistic - 22.5251546702511_dp) < &
      1.0e-6_dp, 'alpha-beta restricted likelihood')
   restrictions = fcvar_restrictions_t()
   allocate(restrictions%level_matrix(1, 2), restrictions%level_value(1))
   restrictions%level_matrix = reshape([1.0_dp, 0.0_dp], [1, 2])
   restrictions%level_value = level_fit%level(1)
   restricted_estimation = fcvar_restricted_estimate(x, 1, 2, &
      reshape([0.3_dp, 0.7_dp, 0.2_dp, 0.6_dp], [2, 2]), restrictions, &
      level_options, .false., .true., 3, 80)
   call check(restricted_estimation%info == 0 .and. &
      restricted_estimation%restriction_rank == 1 .and. &
      abs(restricted_estimation%restricted_fit%level(1) - &
      restrictions%level_value(1)) < 2.0e-10_dp, &
      'exact FCVAR level restriction')

   forecast = fcvar_forecast(x, fit, 3, options)
   expected_path = reshape([ &
      -0.303067025318518_dp, -0.389577281021513_dp, -0.415437209301971_dp, &
      -0.186932010040259_dp, -0.205843412029856_dp, -0.282586404378645_dp], &
      [3, 2])
   call check(forecast%info == 0, 'recursive forecast status')
   call check(maxval(abs(forecast%series - expected_path)) < 3.0e-10_dp, &
      'recursive forecast reference')

   supplied_innovations = reshape([ &
      0.1_dp, 0.05_dp, -0.04_dp, -0.2_dp, 0.03_dp, 0.2_dp], [3, 2])
   simulation = fcvar_simulate_from_innovations(x, fit, supplied_innovations, &
      options)
   expected_path = reshape([ &
      -0.203067025318518_dp, -0.352709542501696_dp, -0.314842858137571_dp, &
      -0.386932010040259_dp, -0.135637323067703_dp, -0.0616899533294316_dp], &
      [3, 2])
   call check(simulation%info == 0, 'innovation-driven simulation status')
   call check(maxval(abs(simulation%series - expected_path)) < 3.0e-10_dp, &
      'innovation-driven simulation reference')

   bootstrap = fcvar_bootstrap_from_signs(x, fit, fit%residuals, &
      [1.0_dp, -1.0_dp, 1.0_dp], options)
   call check(bootstrap%info == 0, 'wild bootstrap status')
   call check(all(abs(bootstrap%innovations(1, :) - &
      (fit%residuals(1, :) - sum(fit%residuals, dim=1)/ &
      real(size(fit%residuals, 1), dp))) < 2.0e-13_dp), &
      'centered residual bootstrap innovation')

   call set_random_seed(194)
   random_path = fcvar_simulate(x, fit, 3, options)
   call check(random_path%info == 0 .and. &
      all(ieee_is_finite(random_path%series)), 'random FCVAR simulation')
   call set_random_seed(195)
   random_path = fcvar_bootstrap(x, fit, fit%residuals, 3, options)
   call check(random_path%info == 0 .and. &
      all(ieee_is_finite(random_path%series)), 'random FCVAR bootstrap')

   zero_rank = fcvar_estimate_fixed(x, 0, 0, 0.6_dp, 0.4_dp, options)
   forecast = fcvar_forecast(x, zero_rank, 2, options)
   call check(zero_rank%info == 0 .and. forecast%info == 0, &
      'zero-rank no-lag forecast')

   lr_test = fcvar_likelihood_ratio(12.1_dp, 10.0_dp, 8, 6)
   call check(lr_test%info == 0 .and. lr_test%degrees_of_freedom == 2, &
      'nested likelihood-ratio test status')
   call check(abs(lr_test%statistic - 4.2_dp) < 1.0e-14_dp, &
      'nested likelihood-ratio statistic')
   call check(abs(lr_test%p_value - 0.122456428252982_dp) < 2.0e-14_dp, &
      'nested likelihood-ratio p-value')

   rank_tests = fcvar_rank_tests(x, 1, reshape([0.2_dp, 0.8_dp, &
      0.2_dp, 0.8_dp], [2, 2]), options, .true., .true., 7, 100)
   call check(rank_tests%info == 0, 'rank testing status')
   call check(maxval(abs(rank_tests%d - [0.344872368547162_dp, &
      0.339604099274659_dp, 0.2_dp])) < 3.0e-7_dp, 'rankwise fractional orders')
   call check(maxval(abs(rank_tests%log_likelihood - &
      [25.1495238971921_dp, 38.6624396730209_dp, &
      40.2852018115273_dp])) < 3.0e-8_dp, 'rankwise log likelihoods')
   call check(all(rank_tests%free_parameters == [5, 8, 9]), &
      'rankwise free parameters')
   call check(maxval(abs(rank_tests%statistic - &
      [30.2713558286704_dp, 3.2455242770128_dp, 0.0_dp])) < 6.0e-8_dp, &
      'rank-versus-full statistics')
   call check(rank_tests%aic_rank == 2 .and. rank_tests%bic_rank == 2, &
      'information-criterion rank selection')
   call check(all(rank_tests%p_value_available .eqv. &
      [.true., .true., .false.]), &
      'fractional rank p-value availability')
   call check(maxval(abs(rank_tests%p_value(:1) - &
      [4.30968787798998e-6_dp, 0.0716187741439483_dp])) < 1.0e-12_dp, &
      'fractional rank p-values')

   white_noise = multivariate_white_noise_test(x, 3)
   call check(white_noise%info == 0, 'multivariate white-noise test status')
   call check(maxval(abs(white_noise%q_statistic - &
      [29.0718241009175_dp, 30.8016983666987_dp])) < 3.0e-11_dp, &
      'univariate Ljung-Box statistics')
   call check(maxval(abs(white_noise%q_p_value - &
      [2.1629395569045e-6_dp, 9.35821659386882e-7_dp])) < 1.0e-15_dp, &
      'univariate Ljung-Box probabilities')
   call check(maxval(abs(white_noise%lm_statistic - &
      [26.321758283943_dp, 19.896554924196_dp])) < 3.0e-11_dp, &
      'robust LM statistics')
   call check(maxval(abs(white_noise%lm_p_value - &
      [8.16673028125603e-6_dp, 1.7833064463979e-4_dp])) < 3.0e-16_dp, &
      'robust LM probabilities')
   call check(abs(white_noise%multivariate_q_statistic - &
      38.5724411695678_dp) < 4.0e-11_dp, 'multivariate Ljung-Box statistic')
   call check(abs(white_noise%multivariate_q_p_value - &
      1.23732736804882e-4_dp) < 3.0e-16_dp, &
      'multivariate Ljung-Box probability')

   lag_selection = fcvar_lag_select(x, 2, 1, reshape([0.2_dp, 0.8_dp, &
      0.2_dp, 0.8_dp], [2, 2]), 3, options, .true., .true., 7, 100)
   call check(lag_selection%info == 0, 'lag-order selection status')
   call check(maxval(abs(lag_selection%d - [0.8_dp, &
      0.339604091656892_dp, 0.2_dp])) < 3.0e-7_dp, &
      'lagwise fractional orders')
   call check(maxval(abs(lag_selection%log_likelihood - &
      [24.8344194770455_dp, 38.6624396730208_dp, &
      41.3513198414063_dp])) < 4.0e-8_dp, 'lagwise log likelihoods')
   call check(maxval(abs(lag_selection%statistic - &
      [0.0_dp, 27.6560403919507_dp, 5.37776033677102_dp])) < 8.0e-8_dp, &
      'sequential lag LR statistics')
   call check(maxval(abs(lag_selection%p_value - &
      [1.0_dp, 1.46436824341666e-5_dp, 0.250685219169742_dp])) < &
      3.0e-13_dp, 'sequential lag LR probabilities')
   call check(all(lag_selection%free_parameters == [4, 8, 12]), &
      'lagwise free parameters')
   call check(maxval(abs(lag_selection%q_p_value(:, 1) - &
      [0.0772030825209113_dp, 0.390590631980864_dp, &
      0.0431607492574439_dp])) < 5.0e-8_dp, &
      'lagwise univariate Q diagnostics')
   call check(maxval(abs(lag_selection%multivariate_q_p_value - &
      [0.0598255166322431_dp, 0.399433192310514_dp, &
      0.11868784940474_dp])) < 5.0e-8_dp, &
      'lagwise multivariate Q diagnostics')
   call check(lag_selection%aic_lag == 1 .and. &
      lag_selection%bic_lag == 1, 'lag information-criterion selection')

   rank_signs = reshape([ &
      1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, &
      -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, &
      1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, &
      1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, &
      -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, &
      -1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp], [10, 3])
   restrictions = fcvar_restrictions_t()
   allocate(restrictions%order_matrix(1, 2), restrictions%order_value(1))
   restrictions%order_matrix = reshape([1.0_dp, 0.0_dp], [1, 2])
   restrictions%order_value = 0.4_dp
   bootstrap_hypothesis = fcvar_bootstrap_hypothesis_from_signs(x, 1, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), restrictions, &
      rank_signs, options, .true., .true., 7, 300)
   call check(bootstrap_hypothesis%info == 0 .and. &
      bootstrap_hypothesis%restriction_rank == 1, &
      'supplied-sign restriction bootstrap status')
   call check(abs(bootstrap_hypothesis%observed_statistic - &
      2.93615443987596_dp) < 6.0e-7_dp, &
      'observed restriction-bootstrap statistic')
   call check(maxval(abs(bootstrap_hypothesis%statistic - &
      [0.429185952305517_dp, 0.624171751663994_dp, &
      14.8396130217157_dp])) < 2.0e-4_dp, &
      'sorted restriction-bootstrap statistics')
   call check(bootstrap_hypothesis%exceedances == 1 .and. &
      abs(bootstrap_hypothesis%p_value - 1.0_dp/3.0_dp) < 1.0e-14_dp, &
      'strict-exceedance restriction-bootstrap probability')

   call set_random_seed(197)
   random_bootstrap_hypothesis = fcvar_bootstrap_hypothesis(x, 1, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), restrictions, &
      2, options, .true., .true., 7, 300)
   call check(random_bootstrap_hypothesis%info == 0 .and. &
      all(ieee_is_finite(random_bootstrap_hypothesis%statistic)), &
      'shared-RNG restriction bootstrap')

   bootstrap_rank = fcvar_bootstrap_rank_from_signs(x, 1, 0, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), rank_signs, &
      options, .true., .true., 7, 100)
   call check(bootstrap_rank%info == 0 .and. &
      bootstrap_rank%null_fit%rank == 0 .and. &
      bootstrap_rank%alternative_fit%rank == 1, &
      'supplied-sign bootstrap rank status')
   call check(abs(bootstrap_rank%observed_statistic - &
      27.0258315516575_dp) < 2.0e-7_dp, 'observed rank LR statistic')
   call check(maxval(abs(bootstrap_rank%statistic - &
      [7.2943907181885_dp, 8.50370840889693_dp, &
      26.2566393004643_dp])) < 5.0e-6_dp, &
      'sorted bootstrap rank statistics')
   call check(bootstrap_rank%exceedances == 0 .and. &
      abs(bootstrap_rank%p_value) < 1.0e-15_dp, &
      'strict-exceedance bootstrap rank probability')

   call set_random_seed(196)
   random_bootstrap_rank = fcvar_bootstrap_rank(x, 1, 0, 1, &
      reshape([0.2_dp, 0.8_dp, 0.2_dp, 0.8_dp], [2, 2]), 2, options, &
      .true., .true., 7, 100)
   call check(random_bootstrap_rank%info == 0 .and. &
      all(ieee_is_finite(random_bootstrap_rank%statistic)), &
      'shared-RNG bootstrap rank test')

   print '(a)', 'FCVAR tests passed'

contains

   subroutine check(condition, message)
      ! Stop the test program when a condition fails.
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_fcvar
