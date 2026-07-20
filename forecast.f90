! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R forecast package.
! Numerical time-series algorithms modelled after the R forecast package.
module forecast_mod
   use kind_mod, only: dp
   use rolling_forecast_mod, only: rolling_forecast_result_t
   use time_series_stats_mod, only: acf_values, pacf_values, ccf_values
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)

   type, public :: forecast_result_t
      real(dp), allocatable :: mean(:), fitted(:), residuals(:), level(:), trend(:), seasonal(:)
      real(dp) :: sigma = 0.0_dp
   end type

   type, public :: accuracy_result_t
      real(dp) :: me, rmse, mae, mpe, mape, mase, acf1
   end type

   type, public :: dm_result_t
      real(dp) :: statistic
      integer :: horizon, power
   end type

   public :: meanf, naive, rwf, snaive, ses, holt, holt_winters, croston, thetaf
   public :: box_cox, inv_box_cox, fourier, moving_average
   public :: acf_values, pacf_values, ccf_values, forecast_accuracy, dm_test
   public :: rolling_forecast_accuracy

contains

   pure function meanf(y, h) result(out)
      !! Forecast every horizon from the sample mean.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: h !! H.
      type(forecast_result_t) :: out
      real(dp) :: mu
      integer :: n
      n = size(y)
      call allocate_result(out, n, h)
      mu = sum(y)/real(n, dp)
      out%mean = mu
      out%fitted = mu
      out%residuals = y - mu
      out%sigma = residual_sd(out%residuals, 1)
   end function meanf

   pure function naive(y, h) result(out)
      !! Forecast every horizon from the final observation.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: h !! H.
      type(forecast_result_t) :: out
      integer :: n
      n = size(y)
      call allocate_result(out, n, h)
      out%mean = y(n)
      out%fitted(1) = ieee_value(0.0_dp, ieee_quiet_nan)
      if (n > 1) out%fitted(2:) = y(:n - 1)
      out%residuals = y - out%fitted
      out%sigma = residual_sd(out%residuals, 0)
   end function naive

   pure function rwf(y, h, drift) result(out)
      !! Fit a random-walk forecast with optional linear drift.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: h !! H.
      logical, intent(in), optional :: drift !! Flag controlling drift.
      type(forecast_result_t) :: out
      real(dp) :: slope
      integer :: i, n
      logical :: use_drift
      n = size(y)
      use_drift = .false.
      if (present(drift)) use_drift = drift
      out = naive(y, h)
      if (use_drift .and. n > 1) then
         slope = (y(n) - y(1))/real(n - 1, dp)
         do i = 1, h
            out%mean(i) = y(n) + real(i, dp)*slope
         end do
         do i = 2, n
            out%fitted(i) = y(i - 1) + slope
         end do
         out%residuals = y - out%fitted
         out%sigma = residual_sd(out%residuals, 1)
      end if
   end function rwf

   pure function snaive(y, period, h) result(out)
      !! Repeat the most recent seasonal cycle over the forecast horizon.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: h !! H.
      type(forecast_result_t) :: out
      integer :: i, n
      n = size(y)
      call allocate_result(out, n, h)
      out%fitted = ieee_value(0.0_dp, ieee_quiet_nan)
      if (period < 1 .or. period > n) return
      do i = period + 1, n
         out%fitted(i) = y(i - period)
      end do
      do i = 1, h
         out%mean(i) = y(n - period + 1 + mod(i - 1, period))
      end do
      out%residuals = y - out%fitted
      out%sigma = residual_sd(out%residuals, 0)
   end function snaive

   pure function ses(y, h, alpha, initial) result(out)
      !! Apply simple exponential smoothing with a fixed smoothing weight.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: h !! H.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: initial !! Initial value.
      type(forecast_result_t) :: out
      real(dp) :: level
      integer :: i, n
      n = size(y)
      call allocate_result(out, n, h)
      level = y(1)
      if (present(initial)) level = initial
      out%fitted(1) = level
      allocate (out%level(n))
      out%level(1) = level
      do i = 2, n
         level = alpha*y(i - 1) + (1.0_dp - alpha)*level
         out%fitted(i) = level
         out%level(i) = level
      end do
      level = alpha*y(n) + (1.0_dp - alpha)*level
      out%mean = level
      out%residuals = y - out%fitted
      out%sigma = residual_sd(out%residuals, 1)
   end function ses

   pure function holt(y, h, alpha, beta, phi) result(out)
      !! Apply Holt's level-trend method with optional trend damping.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: h !! H.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: beta !! Regression or model coefficients.
      real(dp), intent(in), optional :: phi !! Autoregressive or model coefficient.
      type(forecast_result_t) :: out
      real(dp) :: level, trend, old_level, damp, phisum
      integer :: i, n
      n = size(y)
      call allocate_result(out, n, h)
      allocate (out%level(n), out%trend(n))
      damp = 1.0_dp
      if (present(phi)) damp = phi
      level = y(1)
      trend = merge(y(2) - y(1), 0.0_dp, n > 1)
      out%level(1) = level
      out%trend(1) = trend
      out%fitted(1) = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = 2, n
         out%fitted(i) = level + damp*trend
         old_level = level
         level = alpha*y(i) + (1.0_dp - alpha)*(level + damp*trend)
         trend = beta*(level - old_level) + (1.0_dp - beta)*damp*trend
         out%level(i) = level
         out%trend(i) = trend
      end do
      phisum = 0.0_dp
      do i = 1, h
         phisum = phisum + damp**i
         out%mean(i) = level + phisum*trend
      end do
      out%residuals = y - out%fitted
      out%sigma = residual_sd(out%residuals, 0)
   end function holt

   pure function holt_winters(y, period, h, alpha, beta, gamma, multiplicative, phi) result(out)
      !! Apply seasonal Holt-Winters smoothing with additive or multiplicative seasonality.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: h !! H.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: beta !! Regression or model coefficients.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      logical, intent(in), optional :: multiplicative !! Flag controlling multiplicative.
      real(dp), intent(in), optional :: phi !! Autoregressive or model coefficient.
      type(forecast_result_t) :: out
      real(dp) :: level, trend, old_level, damp, phisum, season
      real(dp), allocatable :: s(:)
      integer :: i, j, n
      logical :: mult
      n = size(y)
      call allocate_result(out, n, h)
      allocate (out%level(n), out%trend(n), out%seasonal(n), s(period))
      mult = .false.
      if (present(multiplicative)) mult = multiplicative
      damp = 1.0_dp
      if (present(phi)) damp = phi
      level = sum(y(1:period))/real(period, dp)
      trend = 0.0_dp
      if (n >= 2*period) trend = (sum(y(period + 1:2*period)) - sum(y(1:period)))/real(period*period, dp)
      if (mult) then
         s = y(1:period)/level
      else
         s = y(1:period) - level
      end if
      out%fitted = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = 1, n
         j = 1 + mod(i - 1, period)
         season = s(j)
         if (i > period) then
            if (mult) then
               out%fitted(i) = (level + damp*trend)*season
            else
               out%fitted(i) = level + damp*trend + season
            end if
         end if
         old_level = level
         if (mult) then
            level = alpha*y(i)/season + (1 - alpha)*(level + damp*trend)
            s(j) = gamma*y(i)/level + (1 - gamma)*season
         else
            level = alpha*(y(i) - season) + (1 - alpha)*(level + damp*trend)
            s(j) = gamma*(y(i) - level) + (1 - gamma)*season
         end if
         trend = beta*(level - old_level) + (1 - beta)*damp*trend
         out%level(i) = level
         out%trend(i) = trend
         out%seasonal(i) = s(j)
      end do
      phisum = 0
      do i = 1, h
         phisum = phisum + damp**i
         season = s(1 + mod(n + i - 1, period))
         if (mult) then
            out%mean(i) = (level + phisum*trend)*season
         else
            out%mean(i) = level + phisum*trend + season
         end if
      end do
      out%residuals = y - out%fitted
      out%sigma = residual_sd(out%residuals, 0)
   end function holt_winters

   pure function croston(y, h, alpha, method) result(out)
      !! Forecast intermittent demand using Croston, SBA, or SBJ updates.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: h !! H.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      character(len=*), intent(in), optional :: method !! Algorithm or estimation method.
      type(forecast_result_t) :: out
      real(dp) :: a, demand, interval, ratio, coeff
      integer :: i, n, last, gap
      character(len=16) :: kind
      n = size(y)
      call allocate_result(out, n, h)
      a = .1_dp
      if (present(alpha)) a = alpha
      kind = 'croston'
      if (present(method)) kind = adjustl(method)
      coeff = 1
      if (trim(kind) == 'sba') coeff = 1 - a/2
      if (trim(kind) == 'sbj') coeff = 1 - a/(2 - a)
      out%fitted = ieee_value(0.0_dp, ieee_quiet_nan)
      last = 0
      demand = 0
      interval = 0
      ratio = 0
      do i = 1, n
         if (abs(y(i)) > tiny(1.0_dp)) then
            if (last == 0) then
               demand = y(i)
               interval = real(i, dp)
            else
               gap = i - last
               demand = demand + a*(y(i) - demand)
               interval = interval + a*(real(gap, dp) - interval)
            end if
            ratio = coeff*demand/interval
            last = i
         end if
         if (i < n) out%fitted(i + 1) = ratio
      end do
      out%mean = ratio
      out%residuals = y - out%fitted
      out%sigma = residual_sd(out%residuals, 0)
   end function croston

   pure function thetaf(y, h, alpha) result(out)
      !! Forecast with the classical two-line Theta method.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: h !! H.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      type(forecast_result_t) :: out
      type(forecast_result_t) :: smooth
      real(dp) :: a, slope, intercept
      real(dp), allocatable :: theta2(:)
      integer :: i, n
      n = size(y)
      a = .2_dp
      if (present(alpha)) a = alpha
      slope = 12.0_dp*sum([((real(i, dp) - (n + 1.0_dp)/2)*y(i), i=1, n)])/real(n*(n*n - 1), dp)
      intercept = sum(y)/real(n, dp) - slope*(n + 1.0_dp)/2
      theta2 = 2.0_dp*y - [(intercept + slope*real(i, dp), i=1, n)]
      smooth = ses(theta2, h, a)
      call allocate_result(out, n, h)
      do i = 1, h
         out%mean(i) = 0.5_dp*(smooth%mean(i) + intercept + slope*real(n + i, dp))
      end do
      out%fitted = 0.5_dp*(smooth%fitted + [(intercept + slope*real(i, dp), i=1, n)])
      out%residuals = y - out%fitted
      out%sigma = residual_sd(out%residuals, 2)
   end function thetaf

   pure elemental function box_cox(x, lambda) result(y)
      !! Apply the Box-Cox power transform to one value.
      real(dp), intent(in) :: x !! Input data or predictor values.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp) :: y
      if (abs(lambda) < epsilon(1.0_dp)) then
         y = log(x)
      else
         y = (x**lambda - 1)/lambda
      end if
   end function box_cox

   pure elemental function inv_box_cox(x, lambda) result(y)
      !! Invert the Box-Cox power transform for one value.
      real(dp), intent(in) :: x !! Input data or predictor values.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp) :: y
      if (abs(lambda) < epsilon(1.0_dp)) then
         y = exp(x)
      else
         y = (lambda*x + 1)**(1/lambda)
      end if
   end function inv_box_cox

   pure function fourier(n, period, order, start) result(x)
      !! Generate sine and cosine regressors for a seasonal period.
      integer, intent(in) :: n !! Number of observations or elements.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: order !! Model or polynomial order.
      integer, intent(in), optional :: start !! Start.
      real(dp), allocatable::x(:, :)
      integer::i, k, s
      s = 1
      if (present(start)) s = start
      allocate (x(n, 2*order))
      do k = 1, order
         do i = 1, n
            x(i, 2*k - 1) = sin(2*pi*k*real(s + i - 1, dp)/period)
            x(i, 2*k) = cos(2*pi*k*real(s + i - 1, dp)/period)
         end do
      end do
   end function fourier

   pure function moving_average(y, order, centre) result(x)
      !! Compute an odd or centered-even simple moving average.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      logical, intent(in), optional :: centre !! Flag controlling centre.
      real(dp), allocatable::x(:)
      real(dp), allocatable::w(:)
      integer::i, j, m, left
      logical::cent
      cent = .true.
      if (present(centre)) cent = centre
      m = order
      if (mod(order, 2) == 0 .and. cent) then
         allocate (w(order + 1))
         w = 1.0_dp/order
         w(1) = .5_dp/order
         w(order + 1) = w(1)
      else
         allocate (w(order))
         w = 1.0_dp/order
      end if
      allocate (x(size(y)))
      x = ieee_value(0.0_dp, ieee_quiet_nan)
      left = size(w)/2
      do i = 1 + left, size(y) - (size(w) - left - 1)
         x(i) = sum([(w(j)*y(i - left + j - 1), j=1, size(w))])
      end do
   end function moving_average

   pure function forecast_accuracy(actual, predicted, training, period) result(a)
      !! Calculate common point-forecast error and scale-free accuracy measures.
      real(dp), intent(in) :: actual(:) !! Observed values used for evaluation.
      real(dp), intent(in) :: predicted(:) !! Predicted values.
      real(dp), intent(in), optional :: training(:) !! Training observations.
      integer, intent(in), optional :: period !! Seasonal period.
      type(accuracy_result_t)::a
      real(dp), allocatable::e(:), pe(:)
      real(dp)::scale
      integer::n, m
      n = min(size(actual), size(predicted))
      e = actual(:n) - predicted(:n)
      pe = 100*e/actual(:n)
      a%me = sum(e)/n
      a%rmse = sqrt(sum(e*e)/n)
      a%mae = sum(abs(e))/n
      a%mpe = sum(pe)/n
      a%mape = sum(abs(pe))/n
      a%mase = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(training)) then
         m = 1
         if (present(period)) m = period
         scale = sum(abs(training(1 + m:) - training(:size(training) - m)))/real(size(training) - m, dp)
         a%mase = a%mae/scale
      end if
      if (n > 1) then
         a%acf1 = sum((e(:n - 1) - sum(e)/n)*(e(2:) - sum(e)/n))/sum((e - sum(e)/n)**2)
      else
         a%acf1 = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
   end function forecast_accuracy

   pure function rolling_forecast_accuracy(rolling, training, period) &
      result(accuracy)
      !! Summarize every valid rolling-forecast horizon with common accuracy measures.
      type(rolling_forecast_result_t), intent(in) :: rolling !! Rolling forecast table.
      real(dp), intent(in), optional :: training(:) !! Training observations for MASE scaling.
      integer, intent(in), optional :: period !! Seasonal period for MASE scaling.
      type(accuracy_result_t), allocatable :: accuracy(:)
      real(dp), allocatable :: actual(:), predicted(:)
      integer :: horizon_index

      if (rolling%info /= 0 .or. .not. allocated(rolling%horizon)) then
         allocate(accuracy(0))
         return
      end if
      allocate(accuracy(size(rolling%horizon)))
      do horizon_index = 1, size(rolling%horizon)
         actual = pack(rolling%actual(:, horizon_index), &
            rolling%valid(:, horizon_index))
         predicted = pack(rolling%forecast(:, horizon_index), &
            rolling%valid(:, horizon_index))
         if (size(actual) == 0) cycle
         if (present(training)) then
            if (present(period)) then
               accuracy(horizon_index) = forecast_accuracy(actual, predicted, &
                  training, period)
            else
               accuracy(horizon_index) = forecast_accuracy(actual, predicted, &
                  training)
            end if
         else
            accuracy(horizon_index) = forecast_accuracy(actual, predicted)
         end if
      end do
   end function rolling_forecast_accuracy

   pure function dm_test(e1, e2, h, power, bartlett) result(out)
      !! Calculate the horizon-corrected Diebold-Mariano test statistic.
      real(dp), intent(in) :: e1(:) !! E1.
      real(dp), intent(in) :: e2(:) !! E2.
      integer, intent(in), optional :: h !! H.
      integer, intent(in), optional :: power !! Power.
      logical, intent(in), optional :: bartlett !! Flag controlling bartlett.
      type(dm_result_t)::out
      real(dp), allocatable::d(:), cv(:)
      real(dp)::dv, kfac, w
      integer::i, n, hh, pp
      logical::bw
      hh = 1
      if (present(h)) hh = h
      pp = 2
      if (present(power)) pp = power
      bw = .false.
      if (present(bartlett)) bw = bartlett
      n = min(size(e1), size(e2))
      d = abs(e1(:n))**pp - abs(e2(:n))**pp
      cv = acf_values(d, hh - 1, .true.)
      dv = cv(1)
      do i = 1, hh - 1
         w = 1
         if (bw) w = 1 - real(i, dp)/hh
         dv = dv + 2*w*cv(i + 1)
      end do
      dv = dv/n
      kfac = sqrt(real(n + 1 - 2*hh, dp)/n + real(hh*(hh - 1), dp)/(n*n))
      out%statistic = (sum(d)/n)/sqrt(dv)*kfac
      out%horizon = hh
      out%power = pp
   end function dm_test

   pure subroutine allocate_result(out, n, h)
      !! Allocate and initialize the common forecast result arrays.
      type(forecast_result_t), intent(out) :: out !! Procedure result.
      integer, intent(in) :: n !! Number of observations or elements.
      integer, intent(in) :: h !! H.
      allocate (out%mean(h), out%fitted(n), out%residuals(n))
      out%mean = ieee_value(0.0_dp, ieee_quiet_nan)
      out%fitted = ieee_value(0.0_dp, ieee_quiet_nan)
      out%residuals = ieee_value(0.0_dp, ieee_quiet_nan)
   end subroutine allocate_result

   pure function residual_sd(e, df_used) result(s)
      !! Estimate residual standard deviation while ignoring non-finite values.
      real(dp), intent(in) :: e(:) !! E.
      integer, intent(in) :: df_used !! Degrees of freedom used.
      real(dp)::s, mu
      integer::i, n
      n = 0
      mu = 0
      do i = 1, size(e)
      if (ieee_is_finite(e(i))) then
         n = n + 1
         mu = mu + e(i)
      end if
      end do
      if (n <= df_used + 1) then
         s = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      mu = mu/n
      s = 0
      do i = 1, size(e)
         if (ieee_is_finite(e(i))) s = s + (e(i) - mu)**2
      end do
      s = sqrt(s/real(n - df_used, dp))
   end function residual_sd
end module forecast_mod
