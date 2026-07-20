! SPDX-License-Identifier: MIT
! SPDX-FileComment: Date utilities adapted from the DataFrame Fortran project.
! Gregorian date utilities adapted from DataFrame/date_mod.f90.
! Copyright (c) 2025 Beliavsky. Distributed under the MIT license.
module calendar_mod
   implicit none
   private

   type, public :: date_t
      ! Proleptic Gregorian calendar date.
      integer :: year = 0
      integer :: month = 0
      integer :: day = 0
   contains
      procedure :: to_string => date_to_string
   end type date_t

   public :: date_valid, date_from_iso, date_from_basic
   public :: date_is_leap_year, date_days_in_month, date_day_of_year
   public :: date_day_of_week, date_day_number, date_from_day_number
   public :: date_easter
   public :: operator(+), operator(-), operator(==), operator(/=)
   public :: operator(<), operator(<=), operator(>), operator(>=)

   interface operator(+)
      module procedure add_days_right
      module procedure add_days_left
   end interface operator(+)

   interface operator(-)
      module procedure subtract_days
      module procedure difference_days
   end interface operator(-)

   interface operator(==)
      module procedure equal_dates
   end interface operator(==)

   interface operator(/=)
      module procedure unequal_dates
   end interface operator(/=)

   interface operator(<)
      module procedure earlier_date
   end interface operator(<)

   interface operator(<=)
      module procedure earlier_or_equal_date
   end interface operator(<=)

   interface operator(>)
      module procedure later_date
   end interface operator(>)

   interface operator(>=)
      module procedure later_or_equal_date
   end interface operator(>=)

contains

   pure function date_to_string(this) result(string)
      !! Format a date as ISO yyyy-mm-dd text.
      class(date_t), intent(in) :: this !! This.
      character(len=10) :: string

      string = zero_pad_4(this%year)//'-'//zero_pad_2(this%month)//'-'// &
         zero_pad_2(this%day)
   end function date_to_string

   pure elemental logical function date_valid(value) result(valid)
      !! Report whether a value is a valid Gregorian date.
      type(date_t), intent(in) :: value !! Input value.

      valid = .false.
      if (value%month < 1 .or. value%month > 12) return
      if (value%day < 1) return
      if (value%day > date_days_in_month(value%year, value%month)) return
      valid = .true.
   end function date_valid

   pure elemental logical function date_is_leap_year(year) result(leap)
      !! Report whether a Gregorian year contains February 29.
      integer, intent(in) :: year !! Year.

      leap = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. &
         mod(year, 400) == 0
   end function date_is_leap_year

   pure elemental integer function date_days_in_month(year, month) result(days)
      !! Return the number of days in a Gregorian month.
      integer, intent(in) :: year !! Year.
      integer, intent(in) :: month !! Month.

      select case (month)
      case (1, 3, 5, 7, 8, 10, 12)
         days = 31
      case (4, 6, 9, 11)
         days = 30
      case (2)
         days = merge(29, 28, date_is_leap_year(year))
      case default
         days = 0
      end select
   end function date_days_in_month

   pure function date_from_iso(string) result(value)
      !! Parse yyyy-mm-dd text, returning a zero date on syntax failure.
      character(len=*), intent(in) :: string !! String.
      type(date_t) :: value
      character(len=len(string)) :: text
      integer :: year, month, day
      logical :: year_ok, month_ok, day_ok

      value = date_t(0, 0, 0)
      text = adjustl(string)
      if (len_trim(text) /= 10) return
      if (text(5:5) /= '-' .or. text(8:8) /= '-') return
      call parse_unsigned_integer(text(1:4), year, year_ok)
      call parse_unsigned_integer(text(6:7), month, month_ok)
      call parse_unsigned_integer(text(9:10), day, day_ok)
      if (.not. (year_ok .and. month_ok .and. day_ok)) return
      value = date_t(year, month, day)
   end function date_from_iso

   pure function date_from_basic(string) result(value)
      !! Parse yyyymmdd text, returning a zero date on syntax failure.
      character(len=*), intent(in) :: string !! String.
      type(date_t) :: value
      character(len=len(string)) :: text
      integer :: year, month, day
      logical :: year_ok, month_ok, day_ok

      value = date_t(0, 0, 0)
      text = adjustl(string)
      if (len_trim(text) /= 8) return
      call parse_unsigned_integer(text(1:4), year, year_ok)
      call parse_unsigned_integer(text(5:6), month, month_ok)
      call parse_unsigned_integer(text(7:8), day, day_ok)
      if (.not. (year_ok .and. month_ok .and. day_ok)) return
      value = date_t(year, month, day)
   end function date_from_basic

   pure elemental type(date_t) function add_days_right(value, days) result(sum_date)
      !! Add an integer number of days to a valid date.
      type(date_t), intent(in) :: value !! Input value.
      integer, intent(in) :: days !! Days.

      if (.not. date_valid(value)) then
         sum_date = date_t(0, 0, 0)
      else
         sum_date = date_from_day_number(date_day_number(value) + days)
      end if
   end function add_days_right

   pure elemental type(date_t) function add_days_left(days, value) result(sum_date)
      !! Add a valid date to an integer number of days.
      integer, intent(in) :: days !! Days.
      type(date_t), intent(in) :: value !! Input value.

      sum_date = add_days_right(value, days)
   end function add_days_left

   pure elemental type(date_t) function subtract_days(value, days) result(difference)
      !! Subtract an integer number of days from a valid date.
      type(date_t), intent(in) :: value !! Input value.
      integer, intent(in) :: days !! Days.

      difference = add_days_right(value, -days)
   end function subtract_days

   pure elemental integer function difference_days(left, right) result(difference)
      !! Return signed elapsed days between two valid dates.
      type(date_t), intent(in) :: left !! Left.
      type(date_t), intent(in) :: right !! Right.

      if (.not. date_valid(left) .or. .not. date_valid(right)) then
         difference = 0
      else
         difference = date_day_number(left) - date_day_number(right)
      end if
   end function difference_days

   pure elemental logical function equal_dates(left, right) result(equal)
      !! Test two dates for component equality.
      type(date_t), intent(in) :: left !! Left.
      type(date_t), intent(in) :: right !! Right.

      equal = left%year == right%year .and. left%month == right%month .and. &
         left%day == right%day
   end function equal_dates

   pure elemental logical function unequal_dates(left, right) result(unequal)
      !! Test two dates for inequality.
      type(date_t), intent(in) :: left !! Left.
      type(date_t), intent(in) :: right !! Right.

      unequal = .not. equal_dates(left, right)
   end function unequal_dates

   pure elemental logical function earlier_date(left, right) result(earlier)
      !! Test whether the left date precedes the right date.
      type(date_t), intent(in) :: left !! Left.
      type(date_t), intent(in) :: right !! Right.

      earlier = left%year < right%year .or. &
         (left%year == right%year .and. &
         (left%month < right%month .or. &
         (left%month == right%month .and. left%day < right%day)))
   end function earlier_date

   pure elemental logical function earlier_or_equal_date(left, right) result(earlier)
      !! Test whether the left date precedes or equals the right date.
      type(date_t), intent(in) :: left !! Left.
      type(date_t), intent(in) :: right !! Right.

      earlier = earlier_date(left, right) .or. equal_dates(left, right)
   end function earlier_or_equal_date

   pure elemental logical function later_date(left, right) result(later)
      !! Test whether the left date follows the right date.
      type(date_t), intent(in) :: left !! Left.
      type(date_t), intent(in) :: right !! Right.

      later = .not. earlier_or_equal_date(left, right)
   end function later_date

   pure elemental logical function later_or_equal_date(left, right) result(later)
      !! Test whether the left date follows or equals the right date.
      type(date_t), intent(in) :: left !! Left.
      type(date_t), intent(in) :: right !! Right.

      later = .not. earlier_date(left, right)
   end function later_or_equal_date

   pure elemental integer function date_day_number(value) result(number)
      !! Convert a Gregorian date to days relative to 1970-01-01.
      type(date_t), intent(in) :: value !! Input value.
      integer :: year, month, day, era, year_of_era, day_of_year, day_of_era
      integer :: shifted_month

      year = value%year
      month = value%month
      day = value%day
      if (month <= 2) year = year - 1
      era = floor_divide(year, 400)
      year_of_era = year - era*400
      if (month > 2) then
         shifted_month = month - 3
      else
         shifted_month = month + 9
      end if
      day_of_year = (153*shifted_month + 2)/5 + day - 1
      day_of_era = year_of_era*365 + year_of_era/4 - year_of_era/100 + day_of_year
      number = era*146097 + day_of_era - 719468
   end function date_day_number

   pure elemental type(date_t) function date_from_day_number(number) result(value)
      !! Convert days relative to 1970-01-01 to a Gregorian date.
      integer, intent(in) :: number !! Number.
      integer :: shifted, era, day_of_era, year_of_era, year, day_of_year
      integer :: shifted_month, month, day

      shifted = number + 719468
      era = floor_divide(shifted, 146097)
      day_of_era = shifted - era*146097
      year_of_era = (day_of_era - day_of_era/1460 + day_of_era/36524 - &
         day_of_era/146096)/365
      year = year_of_era + era*400
      day_of_year = day_of_era - &
         (365*year_of_era + year_of_era/4 - year_of_era/100)
      shifted_month = (5*day_of_year + 2)/153
      day = day_of_year - (153*shifted_month + 2)/5 + 1
      if (shifted_month < 10) then
         month = shifted_month + 3
      else
         month = shifted_month - 9
      end if
      if (month <= 2) year = year + 1
      value = date_t(year, month, day)
   end function date_from_day_number

   pure elemental integer function date_day_of_week(value) result(weekday)
      !! Return ISO weekday number, Monday one through Sunday seven.
      type(date_t), intent(in) :: value !! Input value.

      if (.not. date_valid(value)) then
         weekday = 0
      else
         weekday = modulo(date_day_number(value) + 3, 7) + 1
      end if
   end function date_day_of_week

   pure elemental integer function date_day_of_year(value) result(ordinal)
      !! Return a date's one-origin ordinal day within its year.
      type(date_t), intent(in) :: value !! Input value.

      if (.not. date_valid(value)) then
         ordinal = 0
      else
         ordinal = date_day_number(value) - &
            date_day_number(date_t(value%year, 1, 1)) + 1
      end if
   end function date_day_of_year

   pure elemental type(date_t) function date_easter(year) result(value)
      !! Return Gregorian Easter Sunday using the Meeus algorithm.
      integer, intent(in) :: year !! Year.
      integer :: a, b, c, d, e, f, g, h, i, k, l, m, month, day

      if (year < 1583) then
         value = date_t(0, 0, 0)
         return
      end if
      a = mod(year, 19)
      b = year/100
      c = mod(year, 100)
      d = b/4
      e = mod(b, 4)
      f = (b + 8)/25
      g = (b - f + 1)/3
      h = mod(19*a + b - d - g + 15, 30)
      i = c/4
      k = mod(c, 4)
      l = mod(32 + 2*e + 2*i - h - k, 7)
      m = (a + 11*h + 22*l)/451
      month = (h + l - 7*m + 114)/31
      day = mod(h + l - 7*m + 114, 31) + 1
      value = date_t(year, month, day)
   end function date_easter

   pure elemental integer function floor_divide(numerator, denominator) result(quotient)
      !! Return floor division for a positive denominator.
      integer, intent(in) :: numerator !! Numerator polynomial coefficients.
      integer, intent(in) :: denominator !! Denominator polynomial coefficients.

      quotient = numerator/denominator
      if (mod(numerator, denominator) < 0) quotient = quotient - 1
   end function floor_divide

   pure function zero_pad_2(number) result(string)
      !! Format a two-digit nonnegative integer with leading zeroes.
      integer, intent(in) :: number !! Number.
      character(len=2) :: string

      if (number < 0 .or. number > 99) then
         string = '**'
         return
      end if
      string(1:1) = achar(iachar('0') + number/10)
      string(2:2) = achar(iachar('0') + mod(number, 10))
   end function zero_pad_2

   pure function zero_pad_4(number) result(string)
      !! Format a four-digit nonnegative integer with leading zeroes.
      integer, intent(in) :: number !! Number.
      character(len=4) :: string
      integer :: remaining

      if (number < 0 .or. number > 9999) then
         string = '****'
         return
      end if
      remaining = number
      string(4:4) = achar(iachar('0') + mod(remaining, 10))
      remaining = remaining/10
      string(3:3) = achar(iachar('0') + mod(remaining, 10))
      remaining = remaining/10
      string(2:2) = achar(iachar('0') + mod(remaining, 10))
      remaining = remaining/10
      string(1:1) = achar(iachar('0') + mod(remaining, 10))
   end function zero_pad_4

   pure subroutine parse_unsigned_integer(string, number, valid)
      !! Parse a nonnegative decimal integer without formatted I/O.
      character(len=*), intent(in) :: string !! String.
      integer, intent(out) :: number !! Number.
      logical, intent(out) :: valid !! Flag controlling valid.
      character(len=len(string)) :: text
      integer :: i, length, digit

      number = 0
      valid = .false.
      text = adjustl(string)
      length = len_trim(text)
      if (length <= 0) return
      do i = 1, length
         if (text(i:i) < '0' .or. text(i:i) > '9') return
         digit = iachar(text(i:i)) - iachar('0')
         number = 10*number + digit
      end do
      valid = .true.
   end subroutine parse_unsigned_integer

end module calendar_mod
