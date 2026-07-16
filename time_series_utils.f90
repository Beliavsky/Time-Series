! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original utility infrastructure for this Fortran library.
! Shared scalar and option-string utilities.
module time_series_utils_mod
   use kind_mod, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: quiet_nan, lowercase, inverse_standard_normal

contains

   pure function quiet_nan() result(value)
      ! Return a quiet IEEE NaN in library working precision.
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure elemental function lowercase(text) result(lower)
      ! Convert ASCII letters in an option string to lowercase.
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(lower(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase

   pure elemental real(dp) function inverse_standard_normal(probability) result(value)
      ! Approximate the standard-normal quantile using Acklam's rational form.
      real(dp), intent(in) :: probability
      real(dp), parameter :: a(6) = [-39.69683028665376_dp, 220.9460984245205_dp, &
         -275.9285104469687_dp, 138.3577518672690_dp, -30.66479806614716_dp, 2.506628277459239_dp]
      real(dp), parameter :: b(5) = [-54.47609879822406_dp, 161.5858368580409_dp, &
         -155.6989798598866_dp, 66.80131188771972_dp, -13.28068155288572_dp]
      real(dp), parameter :: c(6) = [-0.007784894002430293_dp, -0.3223964580411365_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp]
      real(dp), parameter :: d(4) = [0.007784695709041462_dp, 0.3224671290700398_dp, &
         2.445134137142996_dp, 3.754408661907416_dp]
      real(dp) :: q, r

      if (probability <= 0.0_dp) then
         value = -huge(1.0_dp)
      else if (probability >= 1.0_dp) then
         value = huge(1.0_dp)
      else if (probability < 0.02425_dp) then
         q = sqrt(-2.0_dp*log(probability))
         value = (((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6))/ &
            ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0_dp)
      else if (probability > 0.97575_dp) then
         q = sqrt(-2.0_dp*log(1.0_dp - probability))
         value = -(((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6))/ &
            ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0_dp)
      else
         q = probability - 0.5_dp
         r = q*q
         value = (((((a(1)*r + a(2))*r + a(3))*r + a(4))*r + a(5))*r + a(6))*q/ &
            (((((b(1)*r + b(2))*r + b(3))*r + b(4))*r + b(5))*r + 1.0_dp)
      end if
   end function inverse_standard_normal
end module time_series_utils_mod
