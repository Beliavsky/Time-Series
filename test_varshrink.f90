! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Regression tests for the VARshrink translation.
program test_varshrink
   use kind_mod, only: dp
   use varshrink_mod
   implicit none

   type(varshrink_ridge_path_t) :: ridge
   type(varshrink_semibayes_fit_t) :: semibayes
   type(varshrink_var_fit_t) :: var_fit
   type(varshrink_covariance_t) :: shrinkage
   type(varshrink_nonparametric_fit_t) :: nonparametric
   real(dp) :: x(4, 2), y(4, 2), z(4, 4), series(8, 2)
   real(dp) :: expected(2, 2), likelihood, intensity

   x = reshape([1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, &
      0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], shape(x))
   y = matmul(x, reshape([2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp], [2, 2]))
   ridge = varshrink_multivariate_ridge(y, x, [0.0_dp, 0.5_dp])
   expected = reshape([2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp], [2, 2])
   call assert_true(ridge%info == 0, 'ridge status')
   call assert_true(maxval(abs(ridge%coefficient(:, :, 1) - expected)) < 1.0e-11_dp, &
      'unpenalized multivariate ridge')
   call assert_true(ridge%gcv(1) < ridge%gcv(2), 'ridge GCV selection')
   call assert_true(abs(ridge%gcv(2) - 2.902376033057851_dp) < 1.0e-12_dp, &
      'ridge GCV upstream value')
   call assert_true(ridge%selected == 1, 'ridge selected index')
   call assert_true(abs(ridge%effective_df(1) - 2.0_dp) < 1.0e-11_dp, &
      'ridge effective degrees of freedom')

   semibayes = varshrink_semibayes(y, x, 0.25_dp, conjugate=.true.)
   call assert_true(semibayes%info == 0, 'conjugate semibayes status')
   call assert_true(semibayes%converged, 'conjugate semibayes convergence')
   call assert_true(all(semibayes%covariance == transpose(semibayes%covariance)), &
      'conjugate semibayes symmetric covariance')
   call assert_true(sum(abs(semibayes%coefficient)) > 0.0_dp, &
      'conjugate semibayes coefficients')
   expected = reshape([1.4210526315789473_dp, -0.3157894736842106_dp, &
      0.8684210526315790_dp, 1.973684210526316_dp], [2, 2])
   call assert_true(maxval(abs(semibayes%coefficient - expected)) < 1.0e-12_dp, &
      'conjugate semibayes upstream coefficients')
   expected = reshape([0.9064327485380118_dp, -0.0263157894736841_dp, &
      -0.0263157894736841_dp, 1.2616959064327487_dp], [2, 2])
   call assert_true(maxval(abs(semibayes%covariance - expected)) < 1.0e-12_dp, &
      'conjugate semibayes upstream covariance')

   semibayes = varshrink_semibayes(y, x, 1.0_dp, conjugate=.false.)
   call assert_true(semibayes%converged, 'complete shrinkage convergence')
   call assert_true(all(semibayes%coefficient == 0.0_dp), 'complete shrinkage coefficients')

   semibayes = varshrink_semibayes(y, x, 0.2_dp, degrees_of_freedom=5.0_dp, &
      conjugate=.true.)
   call assert_true(semibayes%info == 0, 'Student-t semibayes status')
   call assert_true(semibayes%weights_estimated, 'Student-t weights estimated')
   call assert_true(all(semibayes%weights > 0.0_dp), 'Student-t positive weights')

   expected = reshape([2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp], [2, 2])
   likelihood = varshrink_log_likelihood(y - matmul(x, expected), &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))
   call assert_true(abs(likelihood + 4.0_dp*log(2.0_dp*acos(-1.0_dp))) < 1.0e-12_dp, &
      'Gaussian log likelihood')
   call assert_true(abs(varshrink_effective_df(x, 0.0_dp) - 2.0_dp) < 1.0e-11_dp, &
      'effective degrees of freedom')

   series(:, 1) = [1.0_dp, 1.2_dp, 1.5_dp, 1.7_dp, 2.0_dp, 2.2_dp, 2.5_dp, 2.7_dp]
   series(:, 2) = [0.5_dp, 0.7_dp, 0.8_dp, 1.1_dp, 1.2_dp, 1.5_dp, 1.7_dp, 1.9_dp]
   var_fit = varshrink_var_ridge(series, 1, [0.0_dp, 0.01_dp], include_intercept=.true.)
   call assert_true(var_fit%info == 0, 'VAR ridge status')
   call assert_true(all(shape(var_fit%ar) == [2, 2]), 'VAR ridge coefficient shape')
   call assert_true(all(shape(var_fit%fitted) == [7, 2]), 'VAR ridge fitted shape')
   call assert_true(maxval(abs(var_fit%covariance - transpose(var_fit%covariance))) < &
      1.0e-14_dp, 'VAR ridge covariance symmetry')

   intensity = varshrink_variance_intensity(series)
   call assert_true(intensity >= 0.0_dp .and. intensity <= 1.0_dp, &
      'variance shrinkage intensity bounds')

   z(:, 1:2) = x
   z(:, 3:4) = y
   shrinkage = varshrink_covariance_shrink(z, lambda=1.0_dp, &
      lambda_variance=1.0_dp)
   call assert_true(shrinkage%info == 0, 'covariance shrinkage status')
   call assert_true(maxval(abs(shrinkage%covariance - &
      transpose(shrinkage%covariance))) < 1.0e-14_dp, &
      'covariance shrinkage symmetry')
   call assert_true(maxval(abs(shrinkage%covariance - &
      reshape([shrinkage%variance(1), 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, shrinkage%variance(2), 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, shrinkage%variance(3), 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, shrinkage%variance(4)], [4, 4]))) < &
      1.0e-14_dp, 'complete covariance shrinkage')

   nonparametric = varshrink_nonparametric(y, x, lambda=0.0_dp, &
      lambda_variance=0.0_dp)
   call assert_true(nonparametric%info == 0, 'nonparametric status')
   call assert_true(maxval(abs(nonparametric%coefficient - &
      reshape([2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp], [2, 2]))) < 1.0e-10_dp, &
      'nonparametric unshrunk coefficients')

   print *, 'VARshrink tests passed'

contains

   subroutine assert_true(condition, label)
      !! Stop the test program when a condition is false.
      logical, intent(in) :: condition !! Flag controlling condition.
      character(len=*), intent(in) :: label !! Label.

      if (.not. condition) then
         print *, 'FAILED: ', trim(label)
         error stop 1
      end if
   end subroutine assert_true

end program test_varshrink
