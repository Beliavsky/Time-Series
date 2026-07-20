! SPDX-License-Identifier: MIT
! SPDX-FileComment: Numerical tests for reusable spline infrastructure.
program test_spline
   use kind_mod, only: dp
   use spline_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   integer, parameter :: observations = 101
   real(dp) :: x(observations), y(observations), roughness(12)
   real(dp), allocatable :: basis_matrix(:, :), penalty(:, :), prediction(:)
   type(spline_basis_t) :: basis
   type(spline_fit_t) :: fixed_fit, selected_fit
   integer :: index

   do index = 1, observations
      x(index) = real(index - 1, dp)/real(observations - 1, dp)
      y(index) = sin(2.0_dp*acos(-1.0_dp)*x(index)) + &
         0.03_dp*cos(17.0_dp*acos(-1.0_dp)*x(index))
   end do

   basis = spline_basis_create(x, 12)
   call check(basis%info == 0 .and. size(basis%knots) == 16, &
      'cubic B-spline basis construction')
   basis_matrix = spline_basis_matrix(basis, x)
   call check(all(shape(basis_matrix) == [observations, 12]) .and. &
      minval(basis_matrix) >= -1.0e-14_dp .and. &
      maxval(abs(sum(basis_matrix, dim=2) - 1.0_dp)) < 1.0e-12_dp, &
      'B-spline partition of unity')
   call check(abs(basis_matrix(1, 1) - 1.0_dp) < 1.0e-14_dp .and. &
      abs(basis_matrix(observations, 12) - 1.0_dp) < 1.0e-14_dp, &
      'B-spline endpoint behavior')

   penalty = spline_difference_penalty(12, 2)
   roughness = [(sin(0.4_dp*real(index, dp)), index=1, 12)]
   call check(all(shape(penalty) == [12, 12]) .and. &
      maxval(abs(penalty - transpose(penalty))) < 1.0e-14_dp .and. &
      dot_product(roughness, matmul(penalty, roughness)) >= -1.0e-12_dp, &
      'second-difference spline penalty')

   fixed_fit = spline_fit(x, y, 12, 0.01_dp)
   call check(fixed_fit%info == 0 .and. &
      fixed_fit%regression%effective_df > 2.0_dp .and. &
      fixed_fit%regression%effective_df < 12.1_dp .and. &
      sqrt(fixed_fit%regression%rss/real(observations, dp)) < 0.05_dp, &
      'fixed-lambda penalized spline fit')
   selected_fit = spline_gcv_fit(x, y, 12)
   call check(selected_fit%info == 0 .and. &
      selected_fit%regression%lambda >= 0.0_dp .and. &
      ieee_is_finite(selected_fit%regression%gcv), 'GCV spline selection')
   prediction = spline_predict(selected_fit, [0.0_dp, 0.25_dp, 0.5_dp, &
      0.75_dp, 1.0_dp])
   call check(size(prediction) == 5 .and. all(ieee_is_finite(prediction)) .and. &
      maxval(abs(prediction - [0.03_dp, 1.0_dp, -0.03_dp, -1.0_dp, &
      0.03_dp])) < 0.09_dp, 'penalized spline prediction')

   print '(a)', 'Spline tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_spline
