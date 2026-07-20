! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original utility infrastructure for this Fortran library.
! Shared scalar and option-string utilities.
module utils_mod
   use kind_mod, only: dp
   use stats_mod, only: normal_quantile
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   private

   type, public :: real_vector_t
      !! Allocatable real vector for ragged numeric collections.
      real(dp), allocatable :: values(:)
   end type real_vector_t

   type, public :: integer_vector_t
      !! Allocatable integer vector for ragged index collections.
      integer, allocatable :: values(:)
   end type integer_vector_t

   public :: quiet_nan, lowercase, inverse_standard_normal

contains

   pure function quiet_nan() result(value)
      !! Return a quiet IEEE NaN in library working precision.
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure elemental function lowercase(text) result(lower)
      !! Convert ASCII letters in an option string to lowercase.
      character(len=*), intent(in) :: text !! Text.
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(lower(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase

   pure elemental real(dp) function inverse_standard_normal(probability) result(value)
      !! Return the standard-normal quantile.
      real(dp), intent(in) :: probability !! Probability value.
      value = normal_quantile(probability)
   end function inverse_standard_normal
end module utils_mod
