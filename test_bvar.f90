! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Numerical tests for algorithms translated from R BVAR.
program test_bvar
   use kind_mod, only: dp
   use bvar_mod
   implicit none

   real(dp) :: y(2, 5), x(3, 5), prior_mean(3, 2)
   real(dp), allocatable :: variance(:)
   type(bvar_dummy_observations_t) :: soc, sur
   type(bvar_evidence_t) :: evidence
   type(bvar_hyper_evaluation_t) :: evaluation
   type(bvar_metropolis_result_t) :: metropolis
   real(dp) :: hyper(6), lower(6), upper(6), shape_par(6), scale_par(6)
   real(dp) :: proposal(6, 6), normals(6, 6), uniforms(6)
   integer :: time

   y(1, :) = [1.0_dp, 1.2_dp, 1.1_dp, 1.4_dp, 1.3_dp]
   y(2, :) = [0.5_dp, 0.4_dp, 0.6_dp, 0.7_dp, 0.8_dp]
   x(1, :) = 1.0_dp
   x(2:3, :) = reshape([(0.1_dp*real(time, dp), time=1, 10)], [2, 5])
   prior_mean = 0.0_dp
   prior_mean(2, 1) = 1.0_dp
   prior_mean(3, 2) = 1.0_dp

   variance = bvar_minnesota_variance([0.5_dp, 2.0_dp], 2, 0.2_dp, &
      2.0_dp, 1.0e6_dp)
   call assert_true(size(variance) == 5, 'Minnesota variance shape')
   call assert_close(variance(2), 0.08_dp, 1.0e-14_dp, &
      'Minnesota first-lag variance')
   call assert_close(variance(5), 0.005_dp, 1.0e-14_dp, &
      'Minnesota second-lag variance')

   soc = bvar_soc_dummy(y, 2, 0.5_dp)
   call assert_true(soc%info == 0 .and. all(shape(soc%y) == [2, 2]) .and. &
      all(shape(soc%x) == [5, 2]), 'SOC dummy shapes')
   call assert_close(soc%y(1, 1), 2.2_dp, 1.0e-14_dp, 'SOC initial mean')
   call assert_true(all(abs(soc%x(1, :)) < 1.0e-14_dp), 'SOC zero intercept')

   sur = bvar_sur_dummy(y, 2, 0.5_dp)
   call assert_true(sur%info == 0 .and. all(shape(sur%y) == [2, 1]) .and. &
      all(shape(sur%x) == [5, 1]), 'SUR dummy shapes')
   call assert_close(sur%x(1, 1), 2.0_dp, 1.0e-14_dp, 'SUR intercept dummy')

   evidence = bvar_conjugate_evidence(y, x, prior_mean, &
      [10.0_dp, 0.5_dp, 0.5_dp], [1.0_dp, 1.0_dp])
   call assert_true(evidence%info == 0, 'conjugate evidence status')
   call assert_true(all(shape(evidence%beta_hat) == [3, 2]) .and. &
      all(shape(evidence%posterior_scale) == [2, 2]), &
      'conjugate evidence shapes')
   call assert_true(evidence%log_marginal_likelihood > -huge(1.0_dp) .and. &
      evidence%log_marginal_likelihood < huge(1.0_dp), &
      'finite log marginal likelihood')
   call assert_true(maxval(abs(evidence%posterior_scale - &
      transpose(evidence%posterior_scale))) < 1.0e-12_dp, &
      'symmetric posterior scale')

   call assert_close(bvar_gamma_log_density(2.0_dp, 3.0_dp, 4.0_dp), &
      2.0_dp*log(2.0_dp) - 0.5_dp - log_gamma(3.0_dp) - 3.0_dp*log(4.0_dp), &
      1.0e-14_dp, 'Gamma log density')
   call assert_close(bvar_inverse_gamma_log_density(2.0_dp, 3.0_dp, 4.0_dp), &
      3.0_dp*log(4.0_dp) - 4.0_dp*log(2.0_dp) - 2.0_dp - log_gamma(3.0_dp), &
      1.0e-14_dp, 'inverse-Gamma log density')

   hyper = [0.2_dp, 2.0_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp]
   lower = [0.01_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.01_dp, 0.01_dp]
   upper = [2.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp]
   shape_par = 2.0_dp
   scale_par = 1.0_dp
   evaluation = bvar_hyper_log_posterior(y, x, y, prior_mean, 1, 10.0_dp, &
      hyper, lower, upper, shape_par, scale_par, .true., .true.)
   call assert_true(evaluation%info == 0, 'hierarchical posterior status')
   call assert_close(evaluation%log_posterior, -12.82635219795579_dp, &
      2.0e-12_dp, 'hierarchical posterior matches BVAR reference')

   proposal = 0.0_dp
   normals = 0.0_dp
   uniforms = 0.5_dp
   metropolis = bvar_hierarchical_metropolis_from_random(y, x, y, prior_mean, &
      1, 10.0_dp, hyper, lower, upper, shape_par, scale_par, .true., .true., &
      2, 2, proposal, normals, uniforms)
   call assert_true(metropolis%info == 0, 'Metropolis status')
   call assert_true(all(shape(metropolis%hyperparameters) == [6, 2]), &
      'Metropolis retained-draw shape')
   call assert_true(maxval(abs(metropolis%hyperparameters - &
      spread(hyper, 2, 2))) < 1.0e-14_dp, 'zero-increment Metropolis draws')
   call assert_close(metropolis%acceptance_rate, 1.0_dp, 1.0e-14_dp, &
      'Metropolis acceptance rate')

   print '(a)', 'All BVAR tests passed.'

contains

   subroutine assert_true(condition, label)
      !! Stop when a logical test condition is false.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: label !! Test label.

      if (.not. condition) error stop 'FAIL: '//trim(label)
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, label)
      !! Stop when two scalar values differ beyond tolerance.
      real(dp), intent(in) :: actual !! Computed value.
      real(dp), intent(in) :: expected !! Reference value.
      real(dp), intent(in) :: tolerance !! Absolute tolerance.
      character(len=*), intent(in) :: label !! Test label.

      if (abs(actual - expected) > tolerance) then
         print '(a,2(1x,es24.16))', trim(label)//':', actual, expected
         error stop 'FAIL: '//trim(label)
      end if
   end subroutine assert_close

end program test_bvar
