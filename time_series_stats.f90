! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original statistical infrastructure for this Fortran library.
! Shared regression and statistical helpers.
module time_series_stats_mod
   use kind_mod, only: dp
   use time_series_linalg_mod, only: invert_matrix
   implicit none
   private

   type, public :: yule_walker_result_t
      ! Coefficients, innovation variance, information criterion, and status.
      real(dp), allocatable :: coefficients(:)
      real(dp) :: variance = 0.0_dp
      real(dp) :: criterion = huge(1.0_dp)
      integer :: info = 0
   end type

   type, public :: burg_result_t
      ! Burg AR coefficients, reflections, variance, covariance, and status.
      real(dp), allocatable :: coefficients(:), reflection(:), covariance(:, :)
      real(dp) :: variance = 0.0_dp
      integer :: info = 0
   end type

   type, public :: harmonic_regression_result_t
      ! Harmonic regression coefficients, inference, decomposition, and extrapolation.
      real(dp), allocatable :: coefficients(:), standard_errors(:), fitted(:), residuals(:), forecast(:)
      real(dp), allocatable :: periods(:)
      real(dp) :: rss = 0.0_dp
      integer :: polynomial_order = 0
      integer :: info = 0
   end type

   public :: ols_fit, regression_rss, yule_walker_fit, burg_fit, harmonic_regression
   public :: normal_quantile

contains

   pure real(dp) function normal_quantile(probability) result(value)
      ! Invert the standard normal distribution by monotone bisection.
      real(dp), intent(in) :: probability
      real(dp) :: lower, upper, middle
      integer :: iteration

      lower = -9.0_dp
      upper = 9.0_dp
      do iteration = 1, 100
         middle = 0.5_dp*(lower + upper)
         if (0.5_dp*erfc(-middle/sqrt(2.0_dp)) < probability) then
            lower = middle
         else
            upper = middle
         end if
      end do
      value = 0.5_dp*(lower + upper)
   end function normal_quantile

   pure subroutine ols_fit(x, y, beta, standard_errors, residuals, rss, info)
      ! Fit ordinary least squares using the shared dense solver.
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), allocatable, intent(out) :: beta(:), standard_errors(:), residuals(:)
      real(dp), intent(out) :: rss
      integer, intent(out) :: info
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
      ! Return residual sum of squares for a possibly empty regression.
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), intent(out) :: rss
      integer, intent(out) :: info
      real(dp), allocatable :: beta(:), se(:), residuals(:)
      call ols_fit(x, y, beta, se, residuals, rss, info)
   end subroutine regression_rss

   pure function yule_walker_fit(series, order) result(out)
      ! Fit a zero-mean AR model by the Yule-Walker equations.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: order
      type(yule_walker_result_t) :: out
      real(dp), allocatable :: covariance(:), toeplitz(:, :), inverse(:, :), right(:)
      real(dp) :: mean_value
      integer :: i, j, lag, n, status

      n = size(series)
      if (order < 0 .or. order >= n) then
         out%info = 1
         return
      end if
      allocate(out%coefficients(order), covariance(0:order))
      mean_value = sum(series)/real(n, dp)
      do lag = 0, order
         covariance(lag) = dot_product(series(:n - lag) - mean_value, &
            series(lag + 1:) - mean_value)/real(n, dp)
      end do
      if (covariance(0) <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      if (order == 0) then
         out%variance = covariance(0)
      else
         allocate(toeplitz(order, order), right(order))
         do i = 1, order
            right(i) = covariance(i)
            do j = 1, order
               toeplitz(i, j) = covariance(abs(i - j))
            end do
         end do
         call invert_matrix(toeplitz, inverse, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%coefficients = matmul(inverse, right)
         out%variance = covariance(0) - dot_product(out%coefficients, right)
      end if
      if (out%variance <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      out%criterion = log(out%variance) + 2.0_dp*real(order, dp)/real(n, dp)
   end function yule_walker_fit

   pure function burg_fit(series, order) result(out)
      ! Fit a stable AR model by Burg forward-backward error recursion.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: order
      type(burg_result_t) :: out
      real(dp), allocatable :: centered(:), forward(:), backward(:), next_forward(:), next_backward(:)
      real(dp), allocatable :: previous(:), sample_covariance(:), toeplitz(:, :), inverse(:, :)
      real(dp) :: mean_value, denominator, reflection
      integer :: n, stage, j, lag, status

      n = size(series)
      if (order < 1 .or. order >= n) then
         out%info = 1
         return
      end if
      mean_value = sum(series)/real(n, dp)
      centered = series - mean_value
      out%variance = sum(centered**2)/real(n, dp)
      if (out%variance <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      allocate(out%coefficients(order), out%reflection(order), previous(order))
      out%coefficients = 0.0_dp
      forward = centered(2:)
      backward = centered(:n - 1)
      do stage = 1, order
         denominator = sum(forward**2 + backward**2)
         if (denominator <= tiny(1.0_dp)) then
            out%info = 3
            return
         end if
         reflection = 2.0_dp*dot_product(forward, backward)/denominator
         reflection = max(-1.0_dp, min(1.0_dp, reflection))
         out%reflection(stage) = reflection
         previous = out%coefficients
         out%coefficients(stage) = reflection
         do j = 1, stage - 1
            out%coefficients(j) = previous(j) - reflection*previous(stage - j)
         end do
         out%variance = out%variance*(1.0_dp - reflection**2)
         if (stage < order) then
            next_forward = forward - reflection*backward
            next_backward = backward - reflection*forward
            forward = next_forward(2:)
            backward = next_backward(:size(next_backward) - 1)
         end if
      end do
      if (out%variance <= tiny(1.0_dp)) then
         out%info = 4
         return
      end if
      allocate(sample_covariance(0:order - 1), toeplitz(order, order))
      do lag = 0, order - 1
         sample_covariance(lag) = dot_product(centered(:n - lag), centered(lag + 1:))/real(n, dp)
      end do
      do stage = 1, order
         do j = 1, order
            toeplitz(stage, j) = sample_covariance(abs(stage - j))
         end do
      end do
      call invert_matrix(toeplitz, inverse, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%covariance = out%variance*inverse/real(n, dp)
   end function burg_fit

   pure function harmonic_regression(series, periods, horizon, polynomial_order) result(out)
      ! Fit polynomial trend and Fourier harmonics and extrapolate their signal.
      real(dp), intent(in) :: series(:), periods(:)
      integer, intent(in), optional :: horizon, polynomial_order
      type(harmonic_regression_result_t) :: out
      real(dp), allocatable :: design(:, :), extended_design(:, :), residuals(:)
      real(dp), allocatable :: beta(:), standard_errors(:)
      real(dp) :: angle
      integer :: n, future, degree, columns, t, power, harmonic, column, status

      n = size(series)
      future = 0
      if (present(horizon)) future = horizon
      degree = 0
      if (present(polynomial_order)) degree = polynomial_order
      columns = degree + 1 + 2*size(periods)
      if (n < 1 .or. future < 0 .or. degree < 0 .or. columns >= n .or. &
         any(periods <= 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(extended_design(n + future, columns))
      do t = 1, n + future
         do power = 0, degree
            extended_design(t, power + 1) = real(t, dp)**power
         end do
         column = degree + 1
         do harmonic = 1, size(periods)
            angle = 2.0_dp*acos(-1.0_dp)*real(t - 1, dp)/periods(harmonic)
            extended_design(t, column + 1) = cos(angle)
            extended_design(t, column + 2) = sin(angle)
            column = column + 2
         end do
      end do
      design = extended_design(:n, :)
      call ols_fit(design, series, beta, standard_errors, residuals, out%rss, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%coefficients = beta
      out%standard_errors = standard_errors
      out%fitted = matmul(design, beta)
      out%residuals = residuals
      allocate(out%forecast(future))
      if (future > 0) out%forecast = matmul(extended_design(n + 1:, :), beta)
      out%periods = periods
      out%polynomial_order = degree
   end function harmonic_regression
end module time_series_stats_mod
