! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original random-number infrastructure for this Fortran library.
! Shared pseudo-random generation and Gaussian transformation utilities.
module random_mod
   use kind_mod, only: dp
   use linalg_mod, only: symmetric_eigen, symmetrize
   implicit none
   private

   public :: set_random_seed
   public :: random_uniform, random_standard_normal, random_standard_normal_matrix
   public :: random_gamma
   public :: random_standard_student, random_standard_johnson_su
   public :: multivariate_normal_from_standard, random_multivariate_normal
   public :: random_multivariate_normal_matrix

contains

   subroutine set_random_seed(seed)
      !! Expand one integer into the compiler's intrinsic random seed vector.
      integer, intent(in) :: seed !! Random-number seed.
      integer, allocatable :: values(:)
      integer :: i, n

      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         values(i) = modulo(abs(seed) + 104729*i + 37*i*i, huge(1) - 1) + 1
      end do
      call random_seed(put=values)
   end subroutine set_random_seed

   real(dp) function random_uniform() result(value)
      !! Return one intrinsic uniform variate on the half-open unit interval.
      call random_number(value)
   end function random_uniform

   real(dp) function random_standard_normal() result(value)
      !! Return one standard-normal variate using a Box-Muller transform.
      real(dp) :: u1, u2

      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      value = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function random_standard_normal

   recursive real(dp) function random_gamma(shape, scale) result(value)
      !! Return a gamma variate using the Marsaglia-Tsang method.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in), optional :: scale !! Scale.
      real(dp) :: factor, d, c, normal, uniform, candidate

      factor = 1.0_dp
      if (present(scale)) factor = scale
      if (shape <= 0.0_dp .or. factor <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (shape < 1.0_dp) then
         uniform = max(random_uniform(), tiny(1.0_dp))
         value = random_gamma(shape + 1.0_dp, factor)*uniform**(1.0_dp/shape)
         return
      end if
      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         normal = random_standard_normal()
         candidate = (1.0_dp + c*normal)**3
         if (candidate <= 0.0_dp) cycle
         uniform = random_uniform()
         if (uniform < 1.0_dp - 0.0331_dp*normal**4) exit
         if (log(max(uniform, tiny(1.0_dp))) < 0.5_dp*normal**2 + d*(1.0_dp - candidate + log(candidate))) exit
      end do
      value = factor*d*candidate
   end function random_gamma

   real(dp) function random_standard_student(degrees) result(value)
      !! Return a Student-t variate standardized to unit variance.
      real(dp), intent(in) :: degrees !! Degrees.
      real(dp) :: chi_square

      if (degrees <= 2.0_dp) then
         value = 0.0_dp
         return
      end if
      chi_square = 2.0_dp*random_gamma(0.5_dp*degrees)
      value = sqrt((degrees - 2.0_dp)/chi_square)*random_standard_normal()
   end function random_standard_student

   real(dp) function random_standard_johnson_su(skew, shape) result(value)
      !! Return the standardized Johnson SU variate used by tsissm.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      real(dp) :: c, omega, reciprocal_shape, w

      if (shape <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      reciprocal_shape = 1.0_dp/shape
      w = exp(reciprocal_shape**2)
      omega = -skew*reciprocal_shape
      c = sqrt(1.0_dp/(0.5_dp*(w - 1.0_dp)*(w*cosh(2.0_dp*omega) + 1.0_dp)))
      value = c*sinh(reciprocal_shape*(random_standard_normal() + skew)) + &
         c*sqrt(w)*sinh(omega)
   end function random_standard_johnson_su

   subroutine random_standard_normal_matrix(draws)
      !! Fill a matrix with independent Box-Muller standard-normal variates.
      real(dp), intent(out) :: draws(:, :) !! Draws.
      real(dp) :: u1, u2, radius, angle
      integer :: linear_index, total, row, column

      total = size(draws)
      linear_index = 1
      do while (linear_index <= total)
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         radius = sqrt(-2.0_dp*log(u1))
         angle = 2.0_dp*acos(-1.0_dp)*u2
         row = 1 + mod(linear_index - 1, size(draws, 1))
         column = 1 + (linear_index - 1)/size(draws, 1)
         draws(row, column) = radius*cos(angle)
         linear_index = linear_index + 1
         if (linear_index <= total) then
            row = 1 + mod(linear_index - 1, size(draws, 1))
            column = 1 + (linear_index - 1)/size(draws, 1)
            draws(row, column) = radius*sin(angle)
            linear_index = linear_index + 1
         end if
      end do
   end subroutine random_standard_normal_matrix

   pure subroutine multivariate_normal_from_standard(mean, covariance, standard, draw, info)
      !! Transform independent standard normals into one multivariate Gaussian draw.
      real(dp), intent(in) :: mean(:) !! Mean value or vector.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix.
      real(dp), intent(in) :: standard(:) !! Standard.
      real(dp), intent(out) :: draw(:) !! Draw.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      real(dp) :: tolerance

      if (size(covariance, 1) /= size(mean) .or. size(covariance, 2) /= size(mean) &
         .or. size(standard) /= size(mean) .or. size(draw) /= size(mean)) then
         info = 1
         return
      end if
      call symmetric_eigen(symmetrize(covariance), eigenvalues, eigenvectors, info)
      if (info /= 0) return
      tolerance = 1.0e-8_dp*max(1.0_dp, abs(eigenvalues(1)))
      if (minval(eigenvalues) < -tolerance) then
         info = 2
         return
      end if
      draw = mean + matmul(eigenvectors, sqrt(max(eigenvalues, 0.0_dp))*standard)
   end subroutine multivariate_normal_from_standard

   subroutine random_multivariate_normal(mean, covariance, draw, info)
      !! Draw one multivariate Gaussian variate using the shared intrinsic stream.
      real(dp), intent(in) :: mean(:) !! Mean value or vector.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix.
      real(dp), intent(out) :: draw(:) !! Draw.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: standard(:, :)

      allocate(standard(size(mean), 1))
      call random_standard_normal_matrix(standard)
      call multivariate_normal_from_standard(mean, covariance, standard(:, 1), draw, info)
   end subroutine random_multivariate_normal

   subroutine random_multivariate_normal_matrix(mean, covariance, draws, info)
      !! Fill rows with multivariate Gaussian draws using one factorization.
      real(dp), intent(in) :: mean(:) !! Mean value or vector.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix.
      real(dp), intent(out) :: draws(:, :) !! Draws.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: factor(:, :), standard(:, :)
      real(dp) :: tolerance

      if (size(covariance, 1) /= size(mean) .or. &
         size(covariance, 2) /= size(mean) .or. &
         size(draws, 2) /= size(mean)) then
         info = 1
         return
      end if
      call symmetric_eigen(symmetrize(covariance), eigenvalues, eigenvectors, &
         info)
      if (info /= 0) return
      tolerance = 1.0e-8_dp*max(1.0_dp, abs(eigenvalues(1)))
      if (minval(eigenvalues) < -tolerance) then
         info = 2
         return
      end if
      allocate(factor(size(mean), size(mean)))
      factor = eigenvectors*spread(sqrt(max(eigenvalues, 0.0_dp)), 1, &
         size(mean))
      allocate(standard(size(draws, 1), size(draws, 2)))
      call random_standard_normal_matrix(standard)
      draws = matmul(standard, transpose(factor)) + &
         spread(mean, 1, size(draws, 1))
      info = 0
   end subroutine random_multivariate_normal_matrix
end module random_mod
