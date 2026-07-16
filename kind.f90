! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original numeric-kind infrastructure for this Fortran library.
! Numeric kind definitions shared by the complete time-series library.
module kind_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   integer, parameter, public :: dp = real64
end module kind_mod
