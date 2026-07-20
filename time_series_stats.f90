! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original time-series statistical infrastructure for this Fortran library.
! Shared autocorrelation, autoregression, and harmonic-regression helpers.
module time_series_stats_mod
   use kind_mod, only: dp
   use linalg_mod, only: invert_matrix
   use stats_mod, only: ols_fit
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

   public :: yule_walker_fit, burg_fit, harmonic_regression
   public :: acf_values, pacf_values, ccf_values

contains

   pure function acf_values(y, lag_max, covariance_values) result(a)
      !! Compute biased autocorrelations or autocovariances through lag_max.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: lag_max !! Lag max.
      logical, intent(in), optional :: covariance_values !! Flag controlling covariance values.
      real(dp), allocatable :: a(:)
      real(dp) :: mu, c0
      integer :: k, n
      logical :: return_covariance

      n = size(y)
      allocate(a(lag_max + 1))
      if (n == 0 .or. lag_max < 0 .or. lag_max >= n) then
         a = 0.0_dp
         return
      end if
      mu = sum(y)/real(n, dp)
      do k = 0, lag_max
         a(k + 1) = sum((y(1:n - k) - mu)*(y(1 + k:n) - mu))/real(n, dp)
      end do
      return_covariance = .false.
      if (present(covariance_values)) return_covariance = covariance_values
      if (.not. return_covariance) then
         c0 = a(1)
         if (c0 > tiny(1.0_dp)) then
            a = a/c0
         else
            a = 0.0_dp
            a(1) = 1.0_dp
         end if
      end if
   end function acf_values

   pure function pacf_values(y, lag_max) result(pacf)
      !! Compute partial autocorrelations using Durbin-Levinson recursion.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: lag_max !! Lag max.
      real(dp), allocatable :: pacf(:)
      real(dp), allocatable :: rho(:), previous(:), current(:)
      real(dp) :: reflection, denominator
      integer :: k, j

      allocate(pacf(max(0, lag_max)))
      if (lag_max <= 0) return
      rho = acf_values(y, lag_max)
      allocate(previous(lag_max), current(lag_max))
      previous = 0.0_dp
      previous(1) = rho(2)
      pacf(1) = previous(1)
      do k = 2, lag_max
         denominator = 1.0_dp - sum([(previous(j)*rho(j + 1), j=1, k - 1)])
         if (abs(denominator) <= tiny(1.0_dp)) then
            pacf(k:) = 0.0_dp
            return
         end if
         reflection = (rho(k + 1) - &
            sum([(previous(j)*rho(k - j + 1), j=1, k - 1)]))/denominator
         current = 0.0_dp
         do j = 1, k - 1
            current(j) = previous(j) - reflection*previous(k - j)
         end do
         current(k) = reflection
         pacf(k) = reflection
         previous = current
      end do
   end function pacf_values

   pure function ccf_values(x, y, lag_max) result(correlation)
      !! Compute normalized cross-correlations at positive and negative lags.
      real(dp), intent(in) :: x(:) !! Input data or predictor values.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: lag_max !! Lag max.
      real(dp), allocatable :: correlation(:)
      real(dp) :: x_mean, y_mean, denominator
      integer :: k, n

      n = min(size(x), size(y))
      allocate(correlation(-lag_max:lag_max))
      correlation = 0.0_dp
      if (n == 0 .or. lag_max < 0 .or. lag_max >= n) return
      x_mean = sum(x(:n))/real(n, dp)
      y_mean = sum(y(:n))/real(n, dp)
      denominator = sqrt(sum((x(:n) - x_mean)**2)*sum((y(:n) - y_mean)**2))
      if (denominator <= tiny(1.0_dp)) return
      do k = 0, lag_max
         correlation(k) = sum((x(1:n - k) - x_mean)*(y(1 + k:n) - y_mean))/denominator
         correlation(-k) = sum((x(1 + k:n) - x_mean)*(y(1:n - k) - y_mean))/denominator
      end do
   end function ccf_values

   pure function yule_walker_fit(series, order) result(out)
      !! Fit a zero-mean AR model by the Yule-Walker equations.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
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
      !! Fit a stable AR model by Burg forward-backward error recursion.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
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
      !! Fit polynomial trend and Fourier harmonics and extrapolate their signal.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: periods(:) !! Periods.
      integer, intent(in), optional :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: polynomial_order !! Polynomial order.
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
