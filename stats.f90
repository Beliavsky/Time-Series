! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original statistical infrastructure for this Fortran library.
! General statistical and regression helpers.
module stats_mod
   use kind_mod, only: dp
   use linalg_mod, only: invert_matrix
   implicit none
   private

   public :: ols_fit, regression_rss
   public :: normal_quantile
   public :: sort, sorted, median, quantile
   public :: variance, standard_deviation, covariance
   public :: correlation_matrix, correlation_from_covariance

contains

   pure elemental real(dp) function normal_quantile(probability) result(value)
      !! Approximate the standard-normal quantile using Acklam's rational form.
      real(dp), intent(in) :: probability !! Probability value.
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
   end function normal_quantile

   pure subroutine sort(values)
      !! Sort a real vector in ascending order.
      real(dp), intent(inout) :: values(:) !! Input values, updated in place.
      real(dp) :: held
      integer :: i, j

      do i = 2, size(values)
         held = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= held) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = held
      end do
   end subroutine sort

   pure function sorted(values) result(ordered)
      !! Return a real vector sorted in ascending order.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp) :: ordered(size(values))

      ordered = values
      call sort(ordered)
   end function sorted

   pure real(dp) function median(values) result(value)
      !! Return the median of a real vector, or zero for an empty vector.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp), allocatable :: ordered(:)
      integer :: n

      n = size(values)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      ordered = sorted(values)
      if (mod(n, 2) == 0) then
         value = 0.5_dp*(ordered(n/2) + ordered(n/2 + 1))
      else
         value = ordered((n + 1)/2)
      end if
   end function median

   pure real(dp) function quantile(ordered, probability) result(value)
      !! Interpolate an ascending sample using R's type-7 quantile.
      real(dp), intent(in) :: ordered(:) !! Ordered.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp) :: fraction, position, selected_probability
      integer :: lower

      if (size(ordered) == 0) then
         value = 0.0_dp
         return
      end if
      selected_probability = max(0.0_dp, min(1.0_dp, probability))
      position = 1.0_dp + real(size(ordered) - 1, dp)*selected_probability
      lower = min(size(ordered), int(floor(position)))
      fraction = position - real(lower, dp)
      if (lower == size(ordered)) then
         value = ordered(lower)
      else
         value = (1.0_dp - fraction)*ordered(lower) + fraction*ordered(lower + 1)
      end if
   end function quantile

   pure real(dp) function variance(values) result(value)
      !! Return the unbiased variance of a real sample.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp) :: center

      if (size(values) < 2) then
         value = 0.0_dp
         return
      end if
      center = sum(values)/real(size(values), dp)
      value = sum((values - center)**2)/real(size(values) - 1, dp)
   end function variance

   pure real(dp) function standard_deviation(values) result(value)
      !! Return the unbiased standard deviation of a real sample.
      real(dp), intent(in) :: values(:) !! Input values.

      value = sqrt(max(0.0_dp, variance(values)))
   end function standard_deviation

   pure real(dp) function covariance(first, second) result(value)
      !! Return the unbiased covariance of two equal-length samples.
      real(dp), intent(in) :: first(:) !! First operand.
      real(dp), intent(in) :: second(:) !! Second operand.
      real(dp) :: first_mean, second_mean

      if (size(first) /= size(second) .or. size(first) < 2) then
         value = 0.0_dp
         return
      end if
      first_mean = sum(first)/real(size(first), dp)
      second_mean = sum(second)/real(size(second), dp)
      value = sum((first - first_mean)*(second - second_mean))/real(size(first) - 1, dp)
   end function covariance

   pure function correlation_matrix(values) result(correlation)
      !! Return the sample correlation matrix of a multivariate data block.
      real(dp), intent(in) :: values(:, :) !! Input values.
      real(dp) :: correlation(size(values, 2), size(values, 2))
      real(dp) :: centered(size(values, 1), size(values, 2))
      real(dp) :: scale(size(values, 2))
      integer :: row, column

      correlation = 0.0_dp
      if (size(values, 1) < 2) return
      centered = values - spread(sum(values, 1)/real(size(values, 1), dp), 1, size(values, 1))
      correlation = matmul(transpose(centered), centered)/real(size(values, 1) - 1, dp)
      do row = 1, size(values, 2)
         scale(row) = sqrt(max(correlation(row, row), tiny(1.0_dp)))
      end do
      do column = 1, size(values, 2)
         do row = 1, size(values, 2)
            correlation(row, column) = correlation(row, column)/(scale(row)*scale(column))
         end do
      end do
      do row = 1, size(values, 2)
         correlation(row, row) = 1.0_dp
      end do
   end function correlation_matrix

   pure function correlation_from_covariance(covariance_matrix) result(correlation)
      !! Convert a covariance matrix to a correlation matrix.
      real(dp), intent(in) :: covariance_matrix(:, :) !! Covariance matrix.
      real(dp), allocatable :: correlation(:, :)
      real(dp) :: scale
      integer :: i, j, n

      n = size(covariance_matrix, 1)
      if (size(covariance_matrix, 2) /= n) then
         allocate(correlation(0, 0))
         return
      end if
      allocate(correlation(n, n))
      correlation = 0.0_dp
      do j = 1, n
         do i = 1, n
            scale = covariance_matrix(i, i)*covariance_matrix(j, j)
            if (scale > 0.0_dp) correlation(i, j) = covariance_matrix(i, j)/sqrt(scale)
         end do
      end do
   end function correlation_from_covariance

   pure subroutine ols_fit(x, y, beta, standard_errors, residuals, rss, info)
      !! Fit ordinary least squares using the shared dense solver.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), allocatable, intent(out) :: beta(:) !! Regression or model coefficients.
      real(dp), allocatable, intent(out) :: standard_errors(:) !! Standard errors.
      real(dp), allocatable, intent(out) :: residuals(:) !! Model residuals.
      real(dp), intent(out) :: rss !! Rss.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: xtx(:, :), inverse(:, :)
      real(dp) :: sigma2
      integer :: i, p
      p = size(x, 2)
      if (p == 0) then
         allocate(beta(0), standard_errors(0), residuals(size(y)))
         residuals = y
         rss = sum(y*y)
         info = 0
         return
      end if
      xtx = matmul(transpose(x), x)
      call invert_matrix(xtx, inverse, info)
      if (info /= 0) return
      allocate(beta(p), standard_errors(p), residuals(size(y)))
      beta = matmul(inverse, matmul(transpose(x), y))
      residuals = y - matmul(x, beta)
      rss = sum(residuals*residuals)
      sigma2 = rss/real(size(y) - p, dp)
      do i = 1, p
         standard_errors(i) = sqrt(max(0.0_dp, sigma2*inverse(i, i)))
      end do
   end subroutine ols_fit

   pure subroutine regression_rss(x, y, rss, info)
      !! Return residual sum of squares for a possibly empty regression.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(out) :: rss !! Rss.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: beta(:), se(:), residuals(:)
      call ols_fit(x, y, beta, se, residuals, rss, info)
   end subroutine regression_rss

end module stats_mod
