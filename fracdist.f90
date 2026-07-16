! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R fracdist package.
! Numerical distribution functions translated from the GPL-3 fracdist package.
module fracdist_mod
   use kind_mod, only: dp
   use fracdist_tables_mod, only: fracdist_table_values
   use time_series_linalg_mod, only: invert_matrix
   use itsmr_mod, only: regularized_gamma_q
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
      ieee_value
   implicit none
   private

   type, public :: fracdist_table_t
      ! Fractional-order grid, probabilities, and simulated distribution quantiles.
      real(dp), allocatable :: b(:), probability(:), chi_square_quantile(:)
      real(dp), allocatable :: quantile(:, :)
      integer :: info = 0
   end type fracdist_table_t

   type, public :: fracdist_probability_t
      ! Upper-tail probability and status for one fracdist statistic.
      real(dp) :: p_value = 1.0_dp
      integer :: info = 0
   end type fracdist_probability_t

   type, public :: fracdist_critical_values_t
      ! Critical values for requested upper-tail significance levels.
      real(dp), allocatable :: critical_value(:)
      integer :: info = 0
   end type fracdist_critical_values_t

   public :: fracdist_get_table, fracdist_blocal
   public :: fracdist_fpval, fracdist_fpcrit
   public :: fracdist_p_value, fracdist_critical_values
   public :: fracdist_chi_square_quantile

contains

   pure function fracdist_get_table(iq, iscon) result(out)
      ! Return the embedded simulated table for one rank and constant choice.
      integer, intent(in) :: iq, iscon
      type(fracdist_table_t) :: out

      allocate(out%b(31), out%probability(221))
      allocate(out%chi_square_quantile(221), out%quantile(221, 31))
      call fracdist_table_values(iq, iscon, out%b, out%probability, &
         out%chi_square_quantile, out%quantile, out%info)
   end function fracdist_get_table

   pure function fracdist_blocal(b, estimated_quantile, b_values) result(value)
      ! Interpolate one simulated quantile locally in the fractional order.
      real(dp), intent(in) :: b, estimated_quantile(:), b_values(:)
      real(dp) :: value
      real(dp), allocatable :: design(:, :), response(:), coefficient(:)
      real(dp) :: weight
      integer :: selected, index, status

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (b < 0.51_dp .or. b > 2.0_dp .or. &
         size(estimated_quantile) /= size(b_values) .or. &
         .not. all(ieee_is_finite(estimated_quantile)) .or. &
         .not. all(ieee_is_finite(b_values))) return
      selected = count(1.0_dp - 5.0_dp*abs(b_values - b) > 0.0_dp)
      if (selected < 3) return
      allocate(design(selected, 3), response(selected))
      selected = 0
      do index = 1, size(b_values)
         weight = max(0.0_dp, 1.0_dp - 5.0_dp*abs(b_values(index) - b))
         if (weight <= 0.0_dp) cycle
         selected = selected + 1
         design(selected, :) = weight*[1.0_dp, b_values(index), &
            b_values(index)**2]
         response(selected) = weight*estimated_quantile(index)
      end do
      call quadratic_regression(design, response, coefficient, status)
      if (status == 0) value = coefficient(1) + coefficient(2)*b + &
         coefficient(3)*b*b
   end function fracdist_blocal

   pure function fracdist_fpval(iq, statistic, probabilities, quantiles, &
      chi_square_quantiles, points) result(p_value)
      ! Interpolate an upper-tail probability from one local response surface.
      integer, intent(in) :: iq
      real(dp), intent(in) :: statistic, probabilities(:), quantiles(:)
      real(dp), intent(in) :: chi_square_quantiles(:)
      integer, intent(in), optional :: points
      real(dp) :: p_value
      real(dp), allocatable :: design(:, :), response(:), coefficient(:)
      real(dp) :: transformed
      integer :: selected_points, half, nearest, first, last, status, count_points

      p_value = ieee_value(0.0_dp, ieee_quiet_nan)
      selected_points = 9
      if (present(points)) selected_points = points
      if (iq < 1 .or. iq > 12 .or. selected_points < 5 .or. &
         mod(selected_points, 2) == 0 .or. size(probabilities) /= 221 .or. &
         size(quantiles) /= 221 .or. size(chi_square_quantiles) /= 221 .or. &
         .not. ieee_is_finite(statistic)) return
      if (statistic < 0.5_dp*quantiles(1)) then
         p_value = 1.0_dp
         return
      end if
      if (statistic > 2.0_dp*quantiles(221)) then
         p_value = 0.0_dp
         return
      end if
      nearest = minloc(abs(statistic - quantiles), dim=1)
      half = selected_points/2
      if (nearest > half .and. nearest < 221 - half) then
         first = nearest - half
         last = nearest + half
      else if (nearest < half) then
         first = 1
         last = max(nearest + half, 5)
      else
         last = 221
         first = 222 - min(max(222 - nearest + half, 5), 221)
      end if
      count_points = last - first + 1
      allocate(design(count_points, 3), response(count_points))
      design(:, 1) = 1.0_dp
      design(:, 2) = quantiles(first:last)
      design(:, 3) = quantiles(first:last)**2
      response = chi_square_quantiles(first:last)
      call quadratic_regression(design, response, coefficient, status)
      if (status /= 0) return
      transformed = max(coefficient(1) + coefficient(2)*statistic + &
         coefficient(3)*statistic**2, 1.0e-6_dp)
      p_value = regularized_gamma_q(0.5_dp*real(iq*iq, dp), &
         0.5_dp*transformed)
      p_value = round_four(p_value)
   end function fracdist_fpval

   pure function fracdist_fpcrit(iq, significance, probabilities, quantiles, &
      chi_square_quantiles, points) result(critical_value)
      ! Interpolate a critical value for one upper-tail significance level.
      integer, intent(in) :: iq
      real(dp), intent(in) :: significance, probabilities(:), quantiles(:)
      real(dp), intent(in) :: chi_square_quantiles(:)
      integer, intent(in), optional :: points
      real(dp) :: critical_value
      real(dp), allocatable :: design(:, :), response(:), coefficient(:)
      real(dp) :: target_probability, target_quantile
      integer :: selected_points, half, nearest, first, last, status, count_points

      critical_value = ieee_value(0.0_dp, ieee_quiet_nan)
      selected_points = 9
      if (present(points)) selected_points = points
      if (iq < 1 .or. iq > 12 .or. significance < 0.0_dp .or. &
         significance > 1.0_dp .or. selected_points < 5 .or. &
         mod(selected_points, 2) == 0 .or. size(probabilities) /= 221 .or. &
         size(quantiles) /= 221 .or. size(chi_square_quantiles) /= 221) return
      if (significance < 1.0e-4_dp .or. significance > 0.9999_dp) then
         critical_value = quantiles(221)
         return
      end if
      target_probability = 1.0_dp - significance
      nearest = minloc(abs(target_probability - probabilities), dim=1)
      half = selected_points/2
      if (nearest > half .and. nearest < 221 - half) then
         first = nearest - half
         last = nearest + half
      else if (nearest < half) then
         first = 1
         last = max(nearest + half, 5)
      else
         last = 221
         first = 222 - max(222 - nearest + half, 5)
      end if
      count_points = last - first + 1
      allocate(design(count_points, 3), response(count_points))
      design(:, 1) = 1.0_dp
      design(:, 2) = chi_square_quantiles(first:last)
      design(:, 3) = chi_square_quantiles(first:last)**2
      response = quantiles(first:last)
      call quadratic_regression(design, response, coefficient, status)
      if (status /= 0) return
      target_quantile = fracdist_chi_square_quantile(target_probability, iq*iq)
      critical_value = coefficient(1) + coefficient(2)*target_quantile + &
         coefficient(3)*target_quantile**2
      critical_value = round_four(critical_value)
   end function fracdist_fpcrit

   pure function fracdist_p_value(iq, iscon, b, statistic) result(out)
      ! Calculate a fracdist or low-order chi-square upper-tail probability.
      integer, intent(in) :: iq, iscon
      real(dp), intent(in) :: b, statistic
      type(fracdist_probability_t) :: out
      type(fracdist_table_t) :: table
      real(dp), allocatable :: local_quantile(:)
      integer :: probability_index

      if (iq < 1 .or. iq > 12 .or. (iscon /= 0 .and. iscon /= 1) .or. &
         b < 0.0_dp .or. b > 2.0_dp .or. statistic < 0.0_dp .or. &
         .not. ieee_is_finite(b) .or. .not. ieee_is_finite(statistic)) then
         out%info = 1
         return
      end if
      if (b < 0.51_dp) then
         out%p_value = regularized_gamma_q(0.5_dp*real(iq*iq, dp), &
            0.5_dp*statistic)
         return
      end if
      table = fracdist_get_table(iq, iscon)
      if (table%info /= 0) then
         out%info = 10 + table%info
         return
      end if
      allocate(local_quantile(221))
      do probability_index = 1, 221
         local_quantile(probability_index) = fracdist_blocal(b, &
            table%quantile(probability_index, :), table%b)
      end do
      if (.not. all(ieee_is_finite(local_quantile))) then
         out%info = 20
         return
      end if
      out%p_value = fracdist_fpval(iq, statistic, table%probability, &
         local_quantile, table%chi_square_quantile)
      if (.not. ieee_is_finite(out%p_value)) out%info = 30
   end function fracdist_p_value

   pure function fracdist_critical_values(iq, iscon, b, significance) result(out)
      ! Calculate fracdist or low-order chi-square critical values.
      integer, intent(in) :: iq, iscon
      real(dp), intent(in) :: b, significance(:)
      type(fracdist_critical_values_t) :: out
      type(fracdist_table_t) :: table
      real(dp), allocatable :: local_quantile(:)
      integer :: probability_index, level

      allocate(out%critical_value(size(significance)))
      out%critical_value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (iq < 1 .or. iq > 12 .or. (iscon /= 0 .and. iscon /= 1) .or. &
         b < 0.0_dp .or. b > 2.0_dp .or. &
         .not. all(ieee_is_finite(significance)) .or. &
         any(significance < 0.0_dp) .or. any(significance > 1.0_dp)) then
         out%info = 1
         return
      end if
      if (b < 0.51_dp) then
         do level = 1, size(significance)
            out%critical_value(level) = fracdist_chi_square_quantile( &
               1.0_dp - significance(level), iq*iq)
         end do
         return
      end if
      table = fracdist_get_table(iq, iscon)
      if (table%info /= 0) then
         out%info = 10 + table%info
         return
      end if
      allocate(local_quantile(221))
      do probability_index = 1, 221
         local_quantile(probability_index) = fracdist_blocal(b, &
            table%quantile(probability_index, :), table%b)
      end do
      if (.not. all(ieee_is_finite(local_quantile))) then
         out%info = 20
         return
      end if
      do level = 1, size(significance)
         out%critical_value(level) = fracdist_fpcrit(iq, significance(level), &
            table%probability, local_quantile, table%chi_square_quantile)
      end do
      if (.not. all(ieee_is_finite(out%critical_value))) out%info = 30
   end function fracdist_critical_values

   pure function fracdist_chi_square_quantile(probability, degrees) result(value)
      ! Invert a chi-square distribution using the shared incomplete gamma.
      real(dp), intent(in) :: probability
      integer, intent(in) :: degrees
      real(dp) :: value
      real(dp) :: lower, upper, target_survival
      integer :: iteration

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (degrees < 1 .or. probability < 0.0_dp .or. probability > 1.0_dp) return
      if (probability == 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (probability == 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      target_survival = 1.0_dp - probability
      lower = 0.0_dp
      upper = max(1.0_dp, real(degrees, dp))
      do while (regularized_gamma_q(0.5_dp*real(degrees, dp), &
         0.5_dp*upper) > target_survival)
         upper = 2.0_dp*upper
         if (upper > 0.1_dp*huge(1.0_dp)) exit
      end do
      do iteration = 1, 120
         value = 0.5_dp*(lower + upper)
         if (regularized_gamma_q(0.5_dp*real(degrees, dp), 0.5_dp*value) > &
            target_survival) then
            lower = value
         else
            upper = value
         end if
      end do
      value = 0.5_dp*(lower + upper)
   end function fracdist_chi_square_quantile

   pure subroutine quadratic_regression(design, response, coefficient, info)
      ! Fit a three-term response surface by ordinary least squares.
      real(dp), intent(in) :: design(:, :), response(:)
      real(dp), allocatable, intent(out) :: coefficient(:)
      integer, intent(out) :: info
      real(dp), allocatable :: cross(:, :), inverse(:, :)

      info = 0
      if (size(design, 2) /= 3 .or. size(design, 1) /= size(response) .or. &
         size(response) < 3) then
         allocate(coefficient(0))
         info = 1
         return
      end if
      cross = matmul(transpose(design), design)
      call invert_matrix(cross, inverse, info)
      if (info /= 0) return
      coefficient = matmul(inverse, matmul(transpose(design), response))
   end subroutine quadratic_regression

   pure elemental function round_four(value) result(rounded)
      ! Round a finite value to four decimal places as in fracdist.
      real(dp), intent(in) :: value
      real(dp) :: rounded

      rounded = anint(1.0e4_dp*value)/1.0e4_dp
   end function round_four

end module fracdist_mod
