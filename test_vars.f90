! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Numerical tests for algorithms translated from R vars.
program test_vars
   use kind_mod, only: dp
   use vars_mod
   use urca_mod, only: johansen_result_t, johansen_test
   use random_mod, only: set_random_seed
   implicit none

   real(dp) :: series(180, 2), innovation(2)
   logical, allocatable :: restrictions(:, :)
   logical :: causes(2)
   type(vars_fit_t) :: fit, restricted
   type(vars_selection_t) :: selection
   type(vars_serial_tests_t) :: serial
   type(vars_normality_tests_t) :: normality
   type(vars_arch_tests_t) :: arch
   type(vars_stability_t) :: stability
   type(vars_test_t) :: instantaneous, granger, bootstrap_granger
   type(vars_bq_t) :: bq
   type(vars_svar_t) :: svar, scored_svar
   type(vars_svec_t) :: svec
   type(vars_vec2var_t) :: level_var
   type(vars_irf_bootstrap_t) :: irf_bootstrap
   type(vars_structural_bootstrap_t) :: structural_bootstrap
   type(vars_svec_bootstrap_t) :: svec_bootstrap
   type(vars_structural_irf_t) :: svar_irf, svec_irf
   type(vars_fevd_t) :: svar_fevd, svec_fevd
   type(johansen_result_t) :: johansen
   complex(dp), allocatable :: roots(:)
   real(dp), allocatable :: phi(:, :, :), psi(:, :, :)
   real(dp) :: a_template(2, 2), b_template(2, 2)
   logical :: estimate_a(2, 2), estimate_b(2, 2)
   logical :: short_run_zero(2, 2), long_run_zero(2, 2)
   character(len=80) :: display_heading
   integer :: t, display_unit

   series = 0.0_dp
   call set_random_seed(731)
   do t = 2, size(series, 1)
      innovation = [0.17_dp*sin(0.71_dp*real(t, dp)), &
         0.13_dp*cos(0.43_dp*real(t, dp))]
      series(t, 1) = 0.4_dp + 0.55_dp*series(t - 1, 1) - &
         0.15_dp*series(t - 1, 2) + innovation(1)
      series(t, 2) = -0.2_dp + 0.2_dp*series(t - 1, 1) + &
         0.45_dp*series(t - 1, 2) + innovation(2)
   end do

   fit = vars_fit(series, 1, vars_deterministic_constant, season=4)
   call assert_true(fit%info == 0, 'VAR fit status')
   call assert_true(all(shape(fit%ar) == [2, 2, 1]), 'VAR coefficient shape')
   call assert_true(size(fit%design, 2) == 6, 'deterministic and seasonal design')
   open(newunit=display_unit, status='scratch', action='readwrite')
   call display(fit, display_unit)
   rewind(display_unit)
   read(display_unit, '(a)') display_heading
   close(display_unit)
   call assert_true(trim(display_heading) == 'Reduced-form VAR fit', &
      'VAR display heading')

   selection = vars_select(series, 4, vars_deterministic_constant)
   call assert_true(selection%info == 0, 'VAR selection status')
   call assert_true(all(selection%selected >= 1 .and. selection%selected <= 4), &
      'VAR selected orders')

   allocate(restrictions(2, size(fit%design, 2)))
   restrictions = .true.
   restrictions(1, 2) = .false.
   restricted = vars_restrict(fit, restrictions)
   call assert_true(restricted%info == 0, 'manual restriction status')
   call assert_close(restricted%coefficients(1, 2), 0.0_dp, 1.0e-14_dp, &
      'manual restriction coefficient')

   phi = vars_phi(fit, 5)
   psi = vars_psi(fit, 5)
   roots = vars_roots(fit)
   call assert_close(phi(1, 1, 1), 1.0_dp, 1.0e-14_dp, 'Phi zero horizon')
   call assert_true(all(shape(psi) == [2, 2, 6]), 'Psi shape')
   call assert_true(size(roots) == 2, 'companion roots count')

   serial = vars_serial_test(fit, 8, 2)
   normality = vars_normality_test(fit%residuals)
   arch = vars_arch_test(fit%residuals, 2, 4)
   call assert_probability(serial%portmanteau%p_value, 'serial p-value')
   call assert_probability(normality%jarque_bera%p_value, 'normality p-value')
   call assert_true(size(normality%univariate) == 2, 'univariate normality results')
   call assert_probability(arch%multivariate%p_value, 'ARCH p-value')

   causes = [.true., .false.]
   instantaneous = vars_instantaneous_causality(fit, causes)
   call assert_probability(instantaneous%p_value, 'instantaneous causality p-value')

   stability = vars_ols_cusum(fit)
   call assert_true(stability%info == 0, 'CUSUM status')
   call assert_true(all(shape(stability%process) == [179, 2]), 'CUSUM shape')

   bq = vars_bq(fit)
   call assert_true(bq%info == 0, 'Blanchard-Quah status')
   call assert_close(bq%long_run(1, 2), 0.0_dp, 1.0e-12_dp, &
      'lower triangular long-run impact')

   a_template = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
   b_template = psi(:, :, 1)
   estimate_a = .false.
   estimate_b = reshape([.true., .true., .false., .true.], [2, 2])
   svar = vars_svar(fit, a_template, b_template, estimate_a, estimate_b, &
      max_iterations=300)
   call assert_true(svar%info == 0, 'structural VAR status')
   call assert_true(maxval(abs(svar%covariance - fit%covariance)) < 1.0e-4_dp, &
      'structural covariance reproduction')
   call assert_true(allocated(svar%b_standard_error), 'structural standard errors')
   scored_svar = vars_svar_scoring(fit, a_template, b_template, estimate_a, &
      estimate_b, max_iterations=100)
   call assert_true(scored_svar%info == 0, 'scored structural VAR status')
   call assert_true(maxval(abs(scored_svar%covariance - fit%covariance)) < 1.0e-4_dp, &
      'scored structural covariance reproduction')

   johansen = johansen_test(series, 'trace', 'none', 2, 'transitory')
   call assert_true(johansen%info == 0, 'Johansen input for SVEC')
   short_run_zero = .false.
   short_run_zero(1, 2) = .true.
   long_run_zero = .false.
   svec = vars_svec(johansen, 1, short_run_zero, long_run_zero, max_iterations=300)
   call assert_true(svec%info == 0, 'structural VEC status')
   call assert_close(svec%short_run(1, 2), 0.0_dp, 1.0e-10_dp, &
      'SVEC short-run restriction')
   call assert_true(allocated(svec%short_run_standard_error), 'SVEC standard errors')
   level_var = vars_vec2var(johansen, 1)
   call assert_true(level_var%info == 0, 'vec2var conversion status')
   call assert_true(all(shape(level_var%ar) == [2, 2, 2]), 'vec2var coefficient shape')
   svar_irf = vars_svar_irf(fit, svar, 4)
   svec_irf = vars_svec_irf(johansen, 1, svec, 4)
   svar_fevd = vars_svar_fevd(fit, svar, 4)
   svec_fevd = vars_svec_fevd(johansen, 1, svec, 4)
   call assert_true(svar_irf%info == 0 .and. svec_irf%info == 0, &
      'structural IRF status')
   call assert_close(svar_irf%response(1, 1, 1), svar%impact(1, 1), 1.0e-14_dp, &
      'SVAR impact response')
   call assert_close(svec_irf%response(1, 1, 1), svec%short_run(1, 1), 1.0e-14_dp, &
      'SVEC impact response')
   call assert_true(svar_fevd%info == 0 .and. svec_fevd%info == 0, 'FEVD status')
   call assert_close(sum(svar_fevd%share(1, :, 4)), 1.0_dp, 1.0e-12_dp, &
      'SVAR FEVD normalization')
   call assert_close(sum(svec_fevd%share(1, :, 4)), 1.0_dp, 1.0e-12_dp, &
      'SVEC FEVD normalization')

   irf_bootstrap = vars_irf_bootstrap(fit, 4, 20, level=0.9_dp)
   call assert_true(irf_bootstrap%info == 0, 'IRF bootstrap status')
   call assert_true(irf_bootstrap%successful_runs == 20, 'IRF bootstrap runs')
   call assert_true(all(irf_bootstrap%lower <= irf_bootstrap%upper), &
      'IRF bootstrap interval ordering')
   structural_bootstrap = vars_svar_bootstrap(fit, a_template, b_template, &
      estimate_a, estimate_b, 10, level=0.9_dp, max_iterations=300)
   call assert_true(structural_bootstrap%info == 0, 'structural bootstrap status')
   call assert_true(all(structural_bootstrap%lower <= structural_bootstrap%upper), &
      'structural bootstrap interval ordering')
   svec_bootstrap = vars_svec_bootstrap(series, johansen, 1, short_run_zero, &
      long_run_zero, 4, max_iterations=300)
   call assert_true(svec_bootstrap%info == 0, 'SVEC bootstrap status')
   call assert_true(svec_bootstrap%successful_runs >= 2, 'SVEC bootstrap runs')
   call assert_true(all(svec_bootstrap%short_run_standard_error >= 0.0_dp), &
      'SVEC bootstrap standard errors')

   granger = vars_granger_test(fit, causes)
   bootstrap_granger = vars_granger_bootstrap(fit, causes, 20)
   call assert_probability(granger%p_value, 'Granger p-value')
   call assert_probability(bootstrap_granger%p_value, 'bootstrap Granger p-value')

   print '(a)', 'All vars tests passed.'

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

   subroutine assert_probability(value, label)
      !! Stop when a value is outside the unit interval.
      real(dp), intent(in) :: value !! Candidate probability.
      character(len=*), intent(in) :: label !! Failure label.

      if (value < 0.0_dp .or. value > 1.0_dp) error stop 'FAIL: '//label
   end subroutine assert_probability
end program test_vars
