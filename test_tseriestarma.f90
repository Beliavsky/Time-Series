! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Regression tests for the R tseriesTARMA package translation.
program test_tseriestarma
   use kind_mod, only: dp
   use random_mod, only: set_random_seed
   use tseriestarma_mod
   implicit none

   type(tseriestarma_model_t) :: model
   type(tseriestarma_model_t) :: plus_model, minus_model
   type(tseriestarma_evaluation_t) :: evaluated
   type(tseriestarma_evaluation_t) :: plus_evaluation, minus_evaluation
   type(tseriestarma_simulation_t) :: simulated
   type(tseriestarma_fit_t) :: fitted, contaminated_fit, robust_fit
   type(tseriestarma_fit2_t) :: fit2
   type(tseriestarma_tar_test_t) :: tar_test
   type(tseriestarma_tarma_test_t) :: tarma_test_fixed, tarma_test_joint
   type(tseriestarma_tar_bootstrap_t) :: bootstrap_iid, bootstrap_hansen
   type(tseriestarma_tar_bootstrap_t) :: bootstrap_rademacher, bootstrap_normal
   type(tseriestarma_tar_bootstrap_t) :: bootstrap_repeat, bootstrap_shared
   type(tseriestarma_tar_bootstrap_t) :: bootstrap_shared_repeat
   type(tseriestarma_unit_root_test_t) :: unit_root_test
   type(tseriestarma_unit_root_bootstrap_t) :: unit_bootstrap_iid
   type(tseriestarma_unit_root_bootstrap_t) :: unit_bootstrap_rademacher
   type(tseriestarma_unit_root_bootstrap_t) :: unit_bootstrap_normal
   type(tseriestarma_garch_test_t) :: garch_test, arch_test
   type(tseriestarma_derivatives_t) :: derivatives
   type(tseriestarma_forecast_t) :: forecast, repeated_forecast
   real(dp) :: standard(300), series(220), residuals(220), prediction
   real(dp) :: contaminated(220), parameters(6), numerical_gradient(6)
   real(dp) :: fit2_series(260), fit2_residuals(260), fit2_regressors(260, 1)
   real(dp) :: forecast_draws(4, 5), expected
   real(dp), allocatable :: initial(:), analytic_gradient(:), robust_weights(:)
   real(dp), allocatable :: iid_random(:, :), rademacher_random(:, :)
   real(dp), allocatable :: normal_random(:, :), hansen_random(:, :)
   real(dp), allocatable :: unit_iid_random(:, :), unit_rademacher_random(:, :)
   real(dp), allocatable :: unit_normal_random(:, :)
   real(dp), allocatable :: accessor_values(:), accessor_matrix(:, :)
   real(dp) :: unit_series(120), unit_innovations(120), critical_values(4)
   real(dp) :: garch_series(180), garch_innovations(180)
   real(dp) :: garch_variance(180), garch_standard(180)
   real(dp) :: garch_critical_values(3)
   real(dp) :: difference_step
   integer :: time, row, replication, scratch_unit, display_size

   model = tseriestarma_model([0.2_dp, 0.45_dp], [-0.1_dp, 0.65_dp], &
      [0.3_dp], [-0.2_dp], [1], [1], [1], [1], delay=1, threshold=0.0_dp, &
      innovation_sd_lower=0.5_dp, innovation_sd_upper=0.8_dp)
   call check(model%info == 0, "model construction")
   do time = 1, size(standard)
      standard(time) = 0.7_dp*sin(0.73_dp*real(time, dp)) + &
         0.3_dp*cos(0.19_dp*real(time, dp))
   end do
   simulated = tseriestarma_simulate_from_standard(model, standard, 80)
   call check(simulated%info == 0 .and. &
      size(simulated%observations) == 220 .and. &
      any(simulated%regime == 1) .and. any(simulated%regime == 2), &
      "supplied-innovation simulation")

   series = 0.0_dp
   residuals = 0.0_dp
   do time = 2, size(series)
      residuals(time) = 0.08_dp*sin(0.91_dp*real(time, dp))
      if (series(time - 1) <= model%threshold) then
         prediction = model%phi_lower(1) + &
            model%phi_lower(2)*series(time - 1) + &
            model%theta_lower(1)*residuals(time - 1)
      else
         prediction = model%phi_upper(1) + &
            model%phi_upper(2)*series(time - 1) + &
            model%theta_upper(1)*residuals(time - 1)
      end if
      series(time) = prediction + residuals(time)
   end do
   evaluated = tseriestarma_evaluate(model, series)
   call check(evaluated%info == 0 .and. evaluated%start_index == 2 .and. &
      maxval(abs(evaluated%residuals(2:) - residuals(2:))) < 1.0e-13_dp, &
      "conditional innovation recursion")
   parameters = [0.2_dp, 0.45_dp, -0.1_dp, 0.65_dp, 0.3_dp, -0.2_dp]
   derivatives = tseriestarma_residual_derivatives(model, series)
   analytic_gradient = tseriestarma_least_squares_gradient(model, series)
   difference_step = 1.0e-6_dp
   do time = 1, size(parameters)
      parameters(time) = parameters(time) + difference_step
      plus_model = tseriestarma_model(parameters(1:2), parameters(3:4), &
         parameters(5:5), parameters(6:6), [1], [1], [1], [1], &
         delay=1, threshold=0.0_dp)
      plus_evaluation = tseriestarma_evaluate(plus_model, series)
      parameters(time) = parameters(time) - 2.0_dp*difference_step
      minus_model = tseriestarma_model(parameters(1:2), parameters(3:4), &
         parameters(5:5), parameters(6:6), [1], [1], [1], [1], &
         delay=1, threshold=0.0_dp)
      minus_evaluation = tseriestarma_evaluate(minus_model, series)
      numerical_gradient(time) = (plus_evaluation%rss - &
         minus_evaluation%rss)/(2.0_dp*difference_step)
      parameters(time) = parameters(time) + difference_step
   end do
   call check(derivatives%info == 0 .and. &
      maxval(abs(analytic_gradient - numerical_gradient)) < 1.0e-6_dp, &
      "analytic residual and objective derivatives")

   tarma_test_fixed = tseriestarma_tarma_test(series, 1, 1, delay=1, &
      ma_fixed=.true., threshold_range=[-0.1_dp, 0.0_dp, 0.1_dp], &
      initial=[0.4_dp, 0.2_dp, 0.0_dp], max_iterations=150, &
      tolerance=1.0e-6_dp)
   tarma_test_joint = tseriestarma_tarma_test(series, 1, 1, delay=1, &
      ma_fixed=.false., threshold_range=[-0.1_dp, 0.0_dp, 0.1_dp], &
      initial=[0.4_dp, 0.2_dp, 0.0_dp], max_iterations=150, &
      tolerance=1.0e-6_dp)
   call check(tarma_test_fixed%info == 0 .and. &
      tarma_test_joint%info == 0 .and. &
      tarma_test_fixed%ma_fixed .and. .not. tarma_test_joint%ma_fixed .and. &
      tarma_test_fixed%degrees_of_freedom == 2 .and. &
      tarma_test_joint%degrees_of_freedom == 3 .and. &
      all(shape(tarma_test_fixed%test_values) == [3, 2]) .and. &
      all(shape(tarma_test_joint%test_values) == [3, 2]) .and. &
      all(tarma_test_fixed%statistic >= 0.0_dp) .and. &
      all(tarma_test_joint%statistic >= 0.0_dp) .and. &
      tarma_test_fixed%noise_fit%likelihood%info == 0 .and. &
      tarma_test_joint%noise_fit%likelihood%info == 0, &
      "ARMA-versus-TARMA supLM tests")

   unit_series = 0.0_dp
   do time = 1, size(unit_series)
      unit_innovations(time) = 0.08_dp*sin(0.71_dp*real(time, dp)) + &
         0.03_dp*cos(0.23_dp*real(time, dp))
   end do
   do time = 2, size(unit_series)
      unit_series(time) = unit_series(time - 1) + unit_innovations(time) + &
         0.4_dp*unit_innovations(time - 1)
   end do
   unit_root_test = tseriestarma_unit_root_test(unit_series, &
      lower_probability=0.25_dp, upper_probability=0.75_dp, &
      initial=[0.4_dp, 0.0_dp], max_iterations=150, tolerance=1.0e-6_dp)
   critical_values = tseriestarma_unit_root_critical_values(100, 0.01_dp, &
      -0.9_dp)
   call check(unit_root_test%info == 0 .and. &
      unit_root_test%effective_observations == 119 .and. &
      size(unit_root_test%test_values) == 60 .and. &
      unit_root_test%statistic >= 0.0_dp .and. &
      all(unit_root_test%critical_values > 0.0_dp) .and. &
      abs(critical_values(1) - 16.676881673042026_dp) < 1.0e-12_dp, &
      "IMA-versus-TARMA unit-root supLM test")

   allocate(unit_iid_random(121, 1))
   allocate(unit_rademacher_random(122, 1), unit_normal_random(122, 1))
   do row = 1, 121
      unit_iid_random(row, 1) = real(mod(41*row + 13, 991), dp)/991.0_dp
   end do
   unit_rademacher_random(1, 1) = 0.35_dp
   unit_normal_random(1, 1) = 0.35_dp
   do row = 2, 122
      unit_rademacher_random(row, 1) = merge(-1.0_dp, 1.0_dp, &
         mod(row, 2) == 0)
      unit_normal_random(row, 1) = sin(0.31_dp*real(row, dp))
   end do
   unit_bootstrap_iid = tseriestarma_unit_root_bootstrap_from_random( &
      unit_series, tseriestarma_bootstrap_iid, unit_iid_random, &
      threshold_range=[unit_series(30), unit_series(60), unit_series(90)], &
      max_iterations=150, tolerance=1.0e-6_dp)
   unit_bootstrap_rademacher = &
      tseriestarma_unit_root_bootstrap_from_random(unit_series, &
      tseriestarma_bootstrap_rademacher, unit_rademacher_random, &
      threshold_range=[unit_series(30), unit_series(60), unit_series(90)], &
      max_iterations=150, tolerance=1.0e-6_dp)
   unit_bootstrap_normal = tseriestarma_unit_root_bootstrap_from_random( &
      unit_series, tseriestarma_bootstrap_normal, unit_normal_random, &
      threshold_range=[unit_series(30), unit_series(60), unit_series(90)], &
      max_iterations=150, tolerance=1.0e-6_dp)
   call check(unit_bootstrap_iid%info == 0 .and. &
      unit_bootstrap_rademacher%info == 0 .and. &
      unit_bootstrap_normal%info == 0 .and. &
      size(unit_bootstrap_iid%bootstrap_statistics) == 1 .and. &
      unit_bootstrap_iid%p_value >= 0.0_dp .and. &
      unit_bootstrap_iid%p_value <= 1.0_dp, &
      "IMA-versus-TARMA bootstrap variants")

   garch_series = 0.0_dp
   garch_innovations = 0.0_dp
   do time = 1, size(garch_series)
      garch_standard(time) = sin(1.71_dp*real(time, dp)) + &
         0.5_dp*cos(0.37_dp*real(time, dp))
   end do
   garch_variance(1) = 0.08_dp/(1.0_dp - 0.1_dp - 0.8_dp)
   garch_innovations(1) = sqrt(garch_variance(1))*garch_standard(1)
   garch_series(1) = garch_innovations(1)
   do time = 2, size(garch_series)
      garch_variance(time) = 0.08_dp + &
         0.1_dp*garch_innovations(time - 1)**2 + &
         0.8_dp*garch_variance(time - 1)
      garch_innovations(time) = sqrt(garch_variance(time))* &
         garch_standard(time)
      garch_series(time) = 0.3_dp*garch_series(time - 1) + &
         0.25_dp*garch_innovations(time - 1) + garch_innovations(time)
   end do
   garch_test = tseriestarma_garch_test(garch_series, 1, 1, 1, 1, &
      delay=1, initial=[0.0_dp, 0.3_dp, 0.25_dp, 0.08_dp, 0.1_dp, &
      0.8_dp], max_iterations=80, tolerance=1.0e-6_dp)
   arch_test = tseriestarma_garch_test(garch_series, 1, 1, 1, 0, &
      delay=1, initial=[0.0_dp, 0.3_dp, 0.25_dp, 0.4_dp, 0.2_dp], &
      max_iterations=80, tolerance=1.0e-6_dp)
   garch_critical_values = tseriestarma_garch_critical_values(3, 0.25_dp)
   call check(garch_test%info == 0 .and. arch_test%info == 0 .and. &
      garch_test%null_fit%causal .and. garch_test%null_fit%invertible .and. &
      garch_test%null_fit%variance_stationary .and. &
      size(garch_test%null_fit%garch) == 1 .and. &
      size(arch_test%null_fit%garch) == 0 .and. &
      garch_test%statistic >= 0.0_dp .and. arch_test%statistic >= 0.0_dp .and. &
      abs(garch_critical_values(1) - 11.32_dp) < 1.0e-12_dp .and. &
      all(garch_test%critical_values > 0.0_dp), &
      "ARMA-GARCH versus TARMA-GARCH supLM tests")

   model = tseriestarma_model([0.15_dp, 0.35_dp], [-0.2_dp, 0.6_dp], &
      [real(dp) ::], [real(dp) ::], [1], [1], [integer ::], [integer ::], &
      delay=1, threshold=0.0_dp)
   do time = 1, size(standard)
      standard(time) = sin(1.17_dp*real(time, dp)) + &
         0.4_dp*cos(0.31_dp*real(time, dp))
   end do
   simulated = tseriestarma_simulate_from_standard(model, 0.12_dp*standard, 80)
   initial = [0.15_dp, 0.35_dp, -0.2_dp, 0.6_dp]
   fitted = tseriestarma_fit(simulated%observations, [1], [1], &
      [integer ::], [integer ::], threshold=0.0_dp, initial=initial, &
      max_iterations=200, tolerance=1.0e-7_dp)
   call check(fitted%info == 0 .and. fitted%evaluation%info == 0 .and. &
      fitted%sigma2 > 0.0_dp .and. &
      allocated(fitted%covariance) .and. &
      maxval(abs(fitted%parameters - initial)) < 0.12_dp, &
      "fixed-threshold conditional fit")

   tar_test = tseriestarma_tar_test(simulated%observations, 1, delay=1, &
      lower_probability=0.25_dp, upper_probability=0.75_dp)
   call check(tar_test%info == 0 .and. &
      tar_test%effective_observations == 219 .and. &
      tar_test%degrees_of_freedom == 2 .and. &
      all(shape(tar_test%test_values) == [110, 2]) .and. &
      size(tar_test%threshold_values) == 110 .and. &
      size(tar_test%coefficients) == 2 .and. &
      size(tar_test%residuals) == 219 .and. &
      all(tar_test%threshold_values(2:) >= &
         tar_test%threshold_values(:109)) .and. &
      all(tar_test%statistic >= 0.0_dp) .and. &
      all(tar_test%lower_regime_proportion > 0.0_dp) .and. &
      all(tar_test%lower_regime_proportion < 1.0_dp) .and. &
      tar_test%sigma2 > 0.0_dp, "AR-versus-TAR supLM tests")

   allocate(iid_random(224, 3), rademacher_random(224, 3))
   allocate(normal_random(224, 3), hansen_random(219, 3))
   do replication = 1, 3
      do row = 1, 224
         iid_random(row, replication) = real(mod(37*row + 17*replication, &
            997), dp)/997.0_dp
         rademacher_random(row, replication) = merge(-1.0_dp, 1.0_dp, &
            mod(row + replication, 2) == 0)
         normal_random(row, replication) = sin(0.37_dp*real(row, dp) + &
            0.61_dp*real(replication, dp))
      end do
      do row = 1, 219
         hansen_random(row, replication) = cos(0.29_dp*real(row, dp) + &
            0.43_dp*real(replication, dp))
      end do
   end do
   bootstrap_iid = tseriestarma_tar_bootstrap_from_random( &
      simulated%observations, 1, tseriestarma_bootstrap_iid, iid_random, &
      delay=1, burnin=5)
   bootstrap_rademacher = tseriestarma_tar_bootstrap_from_random( &
      simulated%observations, 1, tseriestarma_bootstrap_rademacher, &
      rademacher_random, delay=1, burnin=5)
   bootstrap_normal = tseriestarma_tar_bootstrap_from_random( &
      simulated%observations, 1, tseriestarma_bootstrap_normal, &
      normal_random, delay=1, burnin=5)
   bootstrap_hansen = tseriestarma_tar_bootstrap_from_random( &
      simulated%observations, 1, tseriestarma_bootstrap_hansen, &
      hansen_random, delay=1)
   bootstrap_repeat = tseriestarma_tar_bootstrap_from_random( &
      simulated%observations, 1, tseriestarma_bootstrap_normal, &
      normal_random, delay=1, burnin=5)
   call set_random_seed(2468)
   bootstrap_shared = tseriestarma_tar_bootstrap(simulated%observations, 1, &
      replications=3, bootstrap_type=tseriestarma_bootstrap_iid, delay=1, &
      burnin=5)
   call set_random_seed(2468)
   bootstrap_shared_repeat = tseriestarma_tar_bootstrap( &
      simulated%observations, 1, replications=3, &
      bootstrap_type=tseriestarma_bootstrap_iid, delay=1, burnin=5)
   call check(bootstrap_iid%info == 0 .and. &
      bootstrap_rademacher%info == 0 .and. bootstrap_normal%info == 0 .and. &
      bootstrap_hansen%info == 0 .and. bootstrap_shared%info == 0 .and. &
      all(shape(bootstrap_iid%bootstrap_statistics) == [3, 2]) .and. &
      all(shape(bootstrap_hansen%bootstrap_statistics) == [3, 2]) .and. &
      all(bootstrap_iid%p_values >= 0.0_dp) .and. &
      all(bootstrap_iid%p_values <= 1.0_dp) .and. &
      maxval(abs(bootstrap_normal%bootstrap_statistics - &
         bootstrap_repeat%bootstrap_statistics)) < 1.0e-14_dp .and. &
      maxval(abs(bootstrap_shared%bootstrap_statistics - &
         bootstrap_shared_repeat%bootstrap_statistics)) < 1.0e-14_dp, &
      "AR-versus-TAR bootstrap variants")

   contaminated = simulated%observations
   contaminated(100) = contaminated(100) + 5.0_dp
   contaminated_fit = tseriestarma_fit(contaminated, [1], [1], &
      [integer ::], [integer ::], threshold=0.0_dp, initial=initial, &
      max_iterations=200, tolerance=1.0e-7_dp)
   robust_fit = tseriestarma_robust_fit(contaminated, contaminated_fit, &
      0.5_dp, innovation_family=tseriestarma_innovation_normal, &
      trim_probabilities=[0.05_dp, 0.95_dp], max_irls_iterations=30, &
      irls_tolerance=1.0e-5_dp, max_iterations=150, tolerance=1.0e-7_dp)
   call check(robust_fit%info == 0 .and. robust_fit%robust .and. &
      robust_fit%weights(100) < 0.1_dp .and. &
      allocated(robust_fit%covariance) .and. &
      allocated(robust_fit%standard_errors), "normal robust IRLS fit")
   robust_weights = tseriestarma_robust_weights([-5.0_dp, 0.0_dp, 5.0_dp], &
      0.4_dp, tseriestarma_innovation_student, 5.0_dp)
   call check(robust_weights(2) > robust_weights(1) .and. &
      robust_weights(2) > robust_weights(3), "Student robust weights")

   fit2_series = 0.0_dp
   fit2_residuals = 0.0_dp
   do time = 1, size(fit2_series)
      fit2_regressors(time, 1) = cos(0.17_dp*real(time, dp))
   end do
   do time = 3, size(fit2_series)
      fit2_residuals(time) = 0.08_dp*sin(0.83_dp*real(time, dp)) + &
         0.03_dp*cos(0.37_dp*real(time, dp))
      if (fit2_series(time - 1) <= 0.0_dp) then
         prediction = 0.1_dp + 0.25_dp*fit2_series(time - 1) + &
            0.35_dp*fit2_series(time - 2)
      else
         prediction = -0.15_dp + 0.25_dp*fit2_series(time - 1) - &
            0.2_dp*fit2_series(time - 2)
      end if
      fit2_series(time) = prediction + 0.3_dp*fit2_residuals(time - 1) + &
         0.4_dp*fit2_regressors(time, 1) + fit2_residuals(time)
   end do
   fit2 = tseriestarma_fit2(fit2_series, [1], [2], [2], 1, &
      threshold=0.0_dp, delay=1, regressors=fit2_regressors, &
      initial=[0.3_dp, 0.25_dp, 0.1_dp, 0.35_dp, -0.15_dp, -0.2_dp, &
      0.4_dp], max_iterations=250, tolerance=1.0e-7_dp)
   call check(fit2%info == 0 .and. fit2%start_index == 3 .and. &
      size(fit2%common_ar) == 1 .and. size(fit2%lower_ar) == 1 .and. &
      size(fit2%upper_ar) == 1 .and. &
      size(fit2%regression_coefficients) == 1 .and. &
      abs(fit2%common_ar(1) - 0.25_dp) < 0.2_dp .and. &
      abs(fit2%regression_coefficients(1) - 0.4_dp) < 0.2_dp .and. &
      fit2%rss > 0.0_dp .and. fit2%forecast_model%info == 0, &
      "exact common-MA threshold regression fit")

   do time = 1, size(forecast_draws, 2)
      forecast_draws(:, time) = [0.1_dp, -0.2_dp, 0.3_dp, -0.1_dp]* &
         real(time - 3, dp)
   end do
   forecast = tseriestarma_forecast_from_standard(fitted%model, &
      simulated%observations, forecast_draws, [0.1_dp, 0.5_dp, 0.9_dp])
   repeated_forecast = tseriestarma_forecast_from_standard(fitted%model, &
      simulated%observations, forecast_draws, [0.1_dp, 0.5_dp, 0.9_dp])
   if (simulated%observations(size(simulated%observations)) <= &
      fitted%model%threshold) then
      expected = fitted%model%phi_lower(1) + fitted%model%phi_lower(2)* &
         simulated%observations(size(simulated%observations))
   else
      expected = fitted%model%phi_upper(1) + fitted%model%phi_upper(2)* &
         simulated%observations(size(simulated%observations))
   end if
   call check(forecast%info == 0 .and. &
      all(shape(forecast%paths) == [4, 5]) .and. &
      all(shape(forecast%quantiles) == [4, 3]) .and. &
      abs(forecast%point(1) - expected) < 1.0e-12_dp .and. &
      all(forecast%quantiles(:, 1) <= forecast%quantiles(:, 2)) .and. &
      all(forecast%quantiles(:, 2) <= forecast%quantiles(:, 3)) .and. &
      maxval(abs(forecast%paths - repeated_forecast%paths)) < 1.0e-14_dp, &
      "draw-driven forecast distribution")

   accessor_values = tseriestarma_coefficients(fitted)
   accessor_matrix = tseriestarma_covariance(fitted)
   call check(size(accessor_values) == size(fitted%parameters) .and. &
      all(shape(accessor_matrix) == shape(fitted%covariance)), &
      "fit coefficient and covariance accessors")
   accessor_values = tseriestarma_residuals(fit2)
   call check(size(accessor_values) == size(fit2%residuals), &
      "exact-fit residual accessor")
   accessor_values = tseriestarma_fitted_values(fitted)
   call check(size(accessor_values) == size(fitted%evaluation%fitted), &
      "conditional fitted-value accessor")
   accessor_values = tseriestarma_coefficients(garch_test%null_fit)
   call check(size(accessor_values) == 1 + size(garch_test%null_fit%ar) + &
      size(garch_test%null_fit%ma) + 1 + size(garch_test%null_fit%arch) + &
      size(garch_test%null_fit%garch), "ARMA-GARCH coefficient accessor")

   open(newunit=scratch_unit, status='scratch', action='write')
   call display(model, scratch_unit)
   call display(evaluated, scratch_unit)
   call display(simulated, scratch_unit)
   call display(fitted, scratch_unit)
   call display(forecast, scratch_unit)
   call display(fit2, scratch_unit)
   call display(tar_test, scratch_unit)
   call display(tarma_test_fixed, scratch_unit)
   call display(bootstrap_iid, scratch_unit)
   call display(unit_root_test, scratch_unit)
   call display(unit_bootstrap_iid, scratch_unit)
   call display(garch_test%null_fit, scratch_unit)
   call display(garch_test, scratch_unit)
   inquire(unit=scratch_unit, size=display_size)
   close(scratch_unit)
   call check(display_size > 0, "display methods")

   print '(a)', "tseriesTARMA tests passed"

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Test condition.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', "FAILED: "//message
         error stop 1
      end if
   end subroutine check

end program test_tseriestarma
