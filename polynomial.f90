! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original polynomial infrastructure for this Fortran library.
module polynomial_mod
   use kind_mod, only: dp
   implicit none
   private

   public :: polynomial_product, polynomial_product_truncated
   public :: polynomial_value, polynomial_derivative

   type, public :: polynomial_roots_t
      !! Complex polynomial roots and numerical convergence status.
      complex(dp), allocatable :: roots(:)
      integer :: info = 0
   end type polynomial_roots_t

   public :: polynomial_roots

contains

   pure function polynomial_roots(coefficients) result(out)
      !! Return roots of an ascending-power real polynomial.
      real(dp), intent(in) :: coefficients(:) !! Ascending-power polynomial coefficients.
      type(polynomial_roots_t) :: out
      complex(dp), allocatable :: previous(:)
      complex(dp) :: numerator, denominator
      real(dp) :: tolerance, radius, angle, change
      integer :: degree, iteration, root, other, power

      if (size(coefficients) < 1) then
         out%info = 1
         return
      end if
      tolerance = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(coefficients)))
      degree = size(coefficients) - 1
      do while (degree > 0 .and. abs(coefficients(degree + 1)) <= tolerance)
         degree = degree - 1
      end do
      if (degree < 1) then
         allocate(out%roots(0))
         out%info = 0
         return
      end if

      allocate(out%roots(degree), previous(degree))
      if (degree == 1) then
         out%roots(1) = cmplx(-coefficients(1)/coefficients(2), 0.0_dp, dp)
         out%info = 0
         return
      end if
      radius = 1.0_dp + maxval(abs(coefficients(:degree)/ &
         coefficients(degree + 1)))
      do root = 1, degree
         angle = 2.0_dp*acos(-1.0_dp)*real(root - 1, dp)/ &
            real(degree, dp) + 0.1_dp
         out%roots(root) = radius*cmplx(cos(angle), sin(angle), dp)
      end do
      do iteration = 1, 2000
         previous = out%roots
         change = 0.0_dp
         do root = 1, degree
            numerator = cmplx(coefficients(degree + 1), 0.0_dp, dp)
            do power = degree, 1, -1
               numerator = numerator*previous(root) + coefficients(power)
            end do
            denominator = cmplx(1.0_dp, 0.0_dp, dp)
            do other = 1, degree
               if (other /= root) denominator = denominator* &
                  (previous(root) - previous(other))
            end do
            if (abs(denominator) <= tiny(1.0_dp)) then
               out%info = 2
               return
            end if
            out%roots(root) = previous(root) - numerator/denominator
            change = max(change, abs(out%roots(root) - previous(root)))
         end do
         if (change <= 100.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(out%roots)))) then
            out%info = 0
            return
         end if
      end do
      out%info = 3
   end function polynomial_roots

   pure function polynomial_product(first, second) result(product)
      !! Multiply ascending-power polynomial coefficient vectors.
      real(dp), intent(in) :: first(:) !! First operand.
      real(dp), intent(in) :: second(:) !! Second operand.
      real(dp) :: product(size(first) + size(second) - 1)
      integer :: i, j

      product = 0.0_dp
      do i = 1, size(first)
         do j = 1, size(second)
            product(i + j - 1) = product(i + j - 1) + first(i)*second(j)
         end do
      end do
   end function polynomial_product

   pure function polynomial_product_truncated(first, second, maximum_degree) result(product)
      !! Multiply ascending-power polynomials through a requested degree.
      real(dp), intent(in) :: first(0:) !! First operand.
      real(dp), intent(in) :: second(0:) !! Second operand.
      integer, intent(in) :: maximum_degree !! Maximum degree.
      real(dp) :: product(0:maximum_degree)
      integer :: i, j

      product = 0.0_dp
      do i = 0, min(maximum_degree, ubound(first, 1))
         do j = 0, min(maximum_degree - i, ubound(second, 1))
            product(i + j) = product(i + j) + first(i)*second(j)
         end do
      end do
   end function polynomial_product_truncated

   pure real(dp) function polynomial_value(coefficients, argument) result(value)
      !! Evaluate an ascending-power polynomial by Horner recursion.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: argument !! Argument.
      integer :: i

      value = 0.0_dp
      do i = size(coefficients), 1, -1
         value = value*argument + coefficients(i)
      end do
   end function polynomial_value

   pure function polynomial_derivative(coefficients) result(derivative)
      !! Return coefficients of the derivative of an ascending-power polynomial.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp) :: derivative(max(0, size(coefficients) - 1))
      integer :: i

      do i = 1, size(derivative)
         derivative(i) = real(i, dp)*coefficients(i + 1)
      end do
   end function polynomial_derivative

end module polynomial_mod
