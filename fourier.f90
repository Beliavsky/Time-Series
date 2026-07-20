! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original Fourier infrastructure for this Fortran library.
! Shared direct Fourier transforms for real time-series algorithms.
module fourier_mod
   use kind_mod, only: dp
   implicit none
   private
   public :: real_dft, inverse_real_dft, fft_transform

contains

   pure function fft_transform(values, inverse) result(transformed)
      !! Compute a radix-2 complex FFT, optionally with inverse normalization.
      complex(dp), intent(in) :: values(:) !! Input values.
      logical, intent(in), optional :: inverse !! Flag controlling inverse.
      complex(dp), allocatable :: transformed(:)
      complex(dp) :: temporary, factor, step
      real(dp) :: angle
      integer :: n, i, j, bit, width, half, start, offset
      logical :: backward

      n = size(values)
      backward = .false.
      if (present(inverse)) backward = inverse
      if (n < 1 .or. iand(n, n - 1) /= 0) then
         allocate(transformed(0))
         return
      end if
      transformed = values
      j = 0
      do i = 0, n - 1
         if (i < j) then
            temporary = transformed(i + 1)
            transformed(i + 1) = transformed(j + 1)
            transformed(j + 1) = temporary
         end if
         bit = n/2
         do while (bit > 0 .and. j >= bit)
            j = j - bit
            bit = bit/2
         end do
         j = j + bit
      end do
      width = 2
      do while (width <= n)
         half = width/2
         angle = -2.0_dp*acos(-1.0_dp)/real(width, dp)
         if (backward) angle = -angle
         step = cmplx(cos(angle), sin(angle), dp)
         do start = 1, n, width
            factor = cmplx(1.0_dp, 0.0_dp, dp)
            do offset = 0, half - 1
               temporary = factor*transformed(start + offset + half)
               transformed(start + offset + half) = &
                  transformed(start + offset) - temporary
               transformed(start + offset) = transformed(start + offset) + temporary
               factor = factor*step
            end do
         end do
         width = 2*width
      end do
      if (backward) transformed = transformed/real(n, dp)
   end function fft_transform

   pure function real_dft(series) result(coefficients)
      !! Return normalized DFT coefficients with array index equal to bin plus one.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      complex(dp), allocatable :: coefficients(:)
      real(dp) :: angle
      integer :: n, bin, t

      n = size(series)
      allocate(coefficients(n))
      coefficients = cmplx(0.0_dp, 0.0_dp, dp)
      if (n == 0) return
      do bin = 0, n - 1
         do t = 0, n - 1
            angle = -2.0_dp*acos(-1.0_dp)*real(bin*t, dp)/real(n, dp)
            coefficients(bin + 1) = coefficients(bin + 1) + &
               series(t + 1)*cmplx(cos(angle), sin(angle), dp)
         end do
      end do
      coefficients = coefficients/real(n, dp)
   end function real_dft

   pure function inverse_real_dft(coefficients) result(series)
      !! Reconstruct a real series from normalized conjugate-symmetric coefficients.
      complex(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), allocatable :: series(:)
      complex(dp) :: value
      real(dp) :: angle
      integer :: n, bin, t

      n = size(coefficients)
      allocate(series(n))
      do t = 0, n - 1
         value = cmplx(0.0_dp, 0.0_dp, dp)
         do bin = 0, n - 1
            angle = 2.0_dp*acos(-1.0_dp)*real(bin*t, dp)/real(n, dp)
            value = value + coefficients(bin + 1)*cmplx(cos(angle), sin(angle), dp)
         end do
         series(t + 1) = real(value, dp)
      end do
   end function inverse_real_dft
end module fourier_mod
