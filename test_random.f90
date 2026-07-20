! SPDX-License-Identifier: MIT
! SPDX-FileComment: Regression tests for the shared random-number module.
program test_random
   use kind_mod, only: dp
   use random_mod
   implicit none
   real(dp) :: first(3, 4), second(3, 4), mean(2), covariance(2, 2), standard(2), draw(2), uniform
   integer :: info

   call set_random_seed(4217)
   call random_standard_normal_matrix(first)
   call set_random_seed(4217)
   call random_standard_normal_matrix(second)
   call check(maxval(abs(first - second)) < epsilon(1.0_dp), 'seed reproducibility')
   call check(all(first < huge(1.0_dp)) .and. all(first > -huge(1.0_dp)), 'finite normals')
   uniform = random_uniform()
   call check(uniform >= 0.0_dp .and. uniform < 1.0_dp, 'uniform range')

   mean = [2.0_dp, -1.0_dp]
   covariance = 0.0_dp
   covariance(1, 1) = 1.0_dp
   covariance(2, 2) = 1.0_dp
   standard = [0.5_dp, -0.25_dp]
   call multivariate_normal_from_standard(mean, covariance, standard, draw, info)
   call check(info == 0, 'multivariate status')
   call check(maxval(abs(draw - (mean + standard))) < 1.0e-12_dp, 'multivariate transform')
   covariance(2, 2) = -1.0_dp
   call multivariate_normal_from_standard(mean, covariance, standard, draw, info)
   call check(info == 2, 'invalid covariance')
   print '(a)', 'All random_mod tests passed.'

contains

   subroutine check(ok, name)
      !! Stop the test program when a named assertion fails.
      logical, intent(in) :: ok !! Flag controlling ok.
      character(len=*), intent(in) :: name !! Name.
      if (.not. ok) then
         print '(a)', 'FAILED: '//name
         error stop 1
      end if
   end subroutine check
end program
