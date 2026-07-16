! SPDX-License-Identifier: MIT
! SPDX-FileComment: Regression tests for the shared calendar module.
program test_calendar
   use time_series_calendar_mod
   implicit none

   type(date_t) :: value, other

   value = date_from_iso('2024-02-29')
   call check(date_valid(value), 'ISO date parsing')
   call check(value%to_string() == '2024-02-29', 'ISO date formatting')
   call check(date_from_basic('20250205') == date_t(2025, 2, 5), &
      'basic date parsing')
   call check(.not. date_valid(date_t(1900, 2, 29)) .and. &
      date_valid(date_t(2000, 2, 29)), 'Gregorian leap-year validation')
   call check(date_is_leap_year(2000) .and. .not. date_is_leap_year(1900), &
      'Gregorian leap-year rule')
   call check(date_days_in_month(2024, 2) == 29 .and. &
      date_days_in_month(2023, 2) == 28, 'days per month')

   value = date_t(1970, 1, 1)
   call check(date_day_number(value) == 0, 'Unix epoch day number')
   call check(date_day_of_week(value) == 4, 'ISO weekday calculation')
   other = value - 1
   call check(other == date_t(1969, 12, 31), 'date subtraction')
   other = date_t(2024, 2, 28) + 2
   call check(other == date_t(2024, 3, 1), 'leap-day arithmetic')
   call check(date_t(2025, 1, 10) - date_t(2025, 1, 1) == 9, &
      'signed date difference')
   call check(date_from_day_number(date_day_number(date_t(1800, 7, 14))) == &
      date_t(1800, 7, 14), 'day-number round trip')
   call check(date_day_of_year(date_t(2024, 12, 31)) == 366, &
      'ordinal day calculation')

   call check(date_easter(2024) == date_t(2024, 3, 31), 'Easter 2024')
   call check(date_easter(2025) == date_t(2025, 4, 20), 'Easter 2025')
   call check(date_easter(2026) == date_t(2026, 4, 5), 'Easter 2026')
   call check(date_t(2025, 1, 1) < date_t(2025, 1, 2) .and. &
      date_t(2025, 1, 2) >= date_t(2025, 1, 2), 'date comparisons')

contains

   subroutine check(condition, label)
      ! Stop the test program when an assertion fails.
      logical, intent(in) :: condition
      character(*), intent(in) :: label

      if (.not. condition) then
         write (*, '(a)') 'FAILED: '//label
         error stop 1
      end if
   end subroutine check

end program test_calendar
