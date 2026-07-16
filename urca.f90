! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R urca package.
! Unit-root and stationarity tests translated from the R package urca.
module urca_mod
   use kind_mod, only: dp
   use time_series_utils_mod, only: lowercase
   use time_series_linalg_mod, only: invert_matrix, cholesky_lower, symmetric_eigen
   use time_series_stats_mod, only: ols_fit, regression_rss
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   private

   type, public :: adf_result_t
      ! Augmented Dickey-Fuller statistics, critical values, and test regression.
      real(dp), allocatable :: statistic(:), critical_values(:, :)
      real(dp), allocatable :: coefficients(:), standard_errors(:), residuals(:)
      integer :: lags = 0
      integer :: info = 0
      character(len=8) :: model = 'none'
   end type

   type, public :: kpss_result_t
      ! KPSS stationarity statistic, critical values, and detrended residuals.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: critical_values(4) = 0.0_dp
      real(dp), allocatable :: residuals(:)
      integer :: lags = 0
      integer :: info = 0
      character(len=3) :: model = 'mu'
   end type

   type, public :: pp_result_t
      ! Phillips-Perron statistic, auxiliary statistics, and test regression.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: critical_values(3) = 0.0_dp
      real(dp), allocatable :: auxiliary(:), coefficients(:), standard_errors(:), residuals(:)
      integer :: lags = 0
      integer :: info = 0
      character(len=8) :: statistic_type = 'Z-tau'
      character(len=8) :: model = 'constant'
   end type

   type, public :: ers_result_t
      ! ERS DF-GLS or point-optimal statistic and GLS-detrended series.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: critical_values(3) = 0.0_dp
      real(dp), allocatable :: detrended(:), coefficients(:), standard_errors(:), residuals(:)
      integer :: lags = 0
      integer :: info = 0
      character(len=8) :: statistic_type = 'DF-GLS'
      character(len=8) :: model = 'constant'
   end type

   type, public :: za_result_t
      ! Zivot-Andrews minimum statistic, break location, and test regression.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: critical_values(3) = 0.0_dp
      real(dp), allocatable :: break_statistics(:), coefficients(:), standard_errors(:), residuals(:)
      integer :: break_point = 0
      integer :: lags = 0
      integer :: info = 0
      character(len=9) :: model = 'intercept'
   end type

   type, public :: johansen_result_t
      ! Johansen eigenvalues, rank statistics, and VECM coefficient matrices.
      real(dp), allocatable :: eigenvalues(:), statistic(:), critical_values(:, :)
      real(dp), allocatable :: cointegration(:, :), loading(:, :), pi(:, :), delta(:, :), gamma(:, :)
      real(dp), allocatable :: r0(:, :), rk(:, :)
      integer :: lag = 0
      integer :: variables = 0
      integer :: info = 0
      character(len=6) :: test_type = 'trace'
      character(len=10) :: deterministic = 'none'
      character(len=10) :: specification = 'longrun'
   end type

   public :: adf_test, kpss_test, pp_test, ers_test, za_test, johansen_test

contains

   pure function adf_test(y, model, lags, select_lags) result(out)
      ! Perform the urca augmented Dickey-Fuller test.
      real(dp), intent(in) :: y(:)
      character(len=*), intent(in), optional :: model, select_lags
      integer, intent(in), optional :: lags
      type(adf_result_t) :: out
      real(dp), allocatable :: dy(:), dependent(:), x(:, :), beta(:), se(:), residuals(:)
      real(dp), allocatable :: xr(:, :), criteria(:)
      real(dp) :: rss, rss_restricted, criterion, penalty
      integer :: n, max_lag, chosen_lag, q, nobs, info, zlag_column
      character(len=8) :: kind, selection

      kind = 'none'
      if (present(model)) kind = adjustl(model)
      selection = 'fixed'
      if (present(select_lags)) selection = lowercase(adjustl(select_lags))
      max_lag = 1
      if (present(lags)) max_lag = lags
      out%model = kind
      out%lags = max_lag
      if (size(y) < max_lag + 3 .or. max_lag < 0) then
         out%info = 1
         return
      end if
      if (trim(kind) /= 'none' .and. trim(kind) /= 'drift' .and. trim(kind) /= 'trend') then
         out%info = 2
         return
      end if
      if (trim(selection) /= 'fixed' .and. trim(selection) /= 'aic' .and. trim(selection) /= 'bic') then
         out%info = 3
         return
      end if
      n = size(y) - 1
      allocate (dy(n))
      dy = y(2:) - y(:size(y) - 1)
      nobs = n - max_lag
      allocate (dependent(nobs))
      dependent = dy(max_lag + 1:)
      chosen_lag = max_lag
      if (max_lag > 0 .and. trim(selection) /= 'fixed') then
         allocate (criteria(max_lag))
         criteria = huge(1.0_dp)
         penalty = 2.0_dp
         if (trim(selection) == 'bic') penalty = log(real(nobs, dp))
         do q = 1, max_lag
            call build_adf_matrix(y, dy, kind, max_lag, q, x, zlag_column)
            call ols_fit(x, dependent, beta, se, residuals, rss, info)
            if (info == 0) then
               criterion = real(nobs, dp)*log(rss/real(nobs, dp)) + &
                           penalty*real(size(beta) + 1, dp)
               criteria(q) = criterion
            end if
         end do
         chosen_lag = minloc(criteria, dim=1)
      end if
      call build_adf_matrix(y, dy, kind, max_lag, chosen_lag, x, zlag_column)
      call ols_fit(x, dependent, beta, se, residuals, rss, info)
      if (info /= 0) then
         out%info = 10 + info
         return
      end if
      out%lags = chosen_lag
      out%coefficients = beta
      out%standard_errors = se
      out%residuals = residuals
      select case (trim(kind))
      case ('none')
         allocate (out%statistic(1))
         out%statistic(1) = beta(zlag_column)/se(zlag_column)
      case ('drift')
         allocate (out%statistic(2))
         out%statistic(1) = beta(zlag_column)/se(zlag_column)
         call restricted_diff_matrix(dy, max_lag, chosen_lag, .false., xr)
         call regression_rss(xr, dependent, rss_restricted, info)
         out%statistic(2) = f_statistic(rss_restricted, rss, size(x, 2), size(xr, 2), nobs)
      case ('trend')
         allocate (out%statistic(3))
         out%statistic(1) = beta(zlag_column)/se(zlag_column)
         call restricted_diff_matrix(dy, max_lag, chosen_lag, .false., xr)
         call regression_rss(xr, dependent, rss_restricted, info)
         out%statistic(2) = f_statistic(rss_restricted, rss, size(x, 2), size(xr, 2), nobs)
         call restricted_diff_matrix(dy, max_lag, chosen_lag, .true., xr)
         call regression_rss(xr, dependent, rss_restricted, info)
         out%statistic(3) = f_statistic(rss_restricted, rss, size(x, 2), size(xr, 2), nobs)
      end select
      call adf_critical_values(kind, n, out%critical_values)
   end function adf_test

   pure function kpss_test(y, model, lag_rule, use_lag) result(out)
      ! Perform the urca KPSS level- or trend-stationarity test.
      real(dp), intent(in) :: y(:)
      character(len=*), intent(in), optional :: model, lag_rule
      integer, intent(in), optional :: use_lag
      type(kpss_result_t) :: out
      real(dp), allocatable :: x(:, :), beta(:), se(:), partial(:), residuals(:)
      real(dp) :: rss, numerator, denominator, covariance, weight
      integer :: i, n, lmax, info
      character(len=8) :: kind, rule

      kind = 'mu'
      if (present(model)) kind = adjustl(model)
      rule = 'short'
      if (present(lag_rule)) rule = lowercase(adjustl(lag_rule))
      n = size(y)
      if (n < 3) then
         out%info = 1
         return
      end if
      if (trim(kind) /= 'mu' .and. trim(kind) /= 'tau') then
         out%info = 2
         return
      end if
      if (present(use_lag)) then
         lmax = use_lag
      else if (trim(rule) == 'short') then
         lmax = int(4.0_dp*(real(n, dp)/100.0_dp)**0.25_dp)
      else if (trim(rule) == 'long') then
         lmax = int(12.0_dp*(real(n, dp)/100.0_dp)**0.25_dp)
      else if (trim(rule) == 'nil') then
         lmax = 0
      else
         out%info = 3
         return
      end if
      if (lmax < 0 .or. lmax >= n) then
         out%info = 4
         return
      end if
      if (trim(kind) == 'mu') then
         allocate (x(n, 1))
         x = 1.0_dp
         out%critical_values = [0.347_dp, 0.463_dp, 0.574_dp, 0.739_dp]
      else
         allocate (x(n, 2))
         x(:, 1) = 1.0_dp
         x(:, 2) = [(real(i, dp), i=1, n)]
         out%critical_values = [0.119_dp, 0.146_dp, 0.176_dp, 0.216_dp]
      end if
      call ols_fit(x, y, beta, se, residuals, rss, info)
      if (info /= 0) then
         out%info = 10 + info
         return
      end if
      allocate (partial(n))
      partial(1) = residuals(1)
      do i = 2, n
         partial(i) = partial(i - 1) + residuals(i)
      end do
      numerator = sum(partial*partial)/real(n*n, dp)
      denominator = sum(residuals*residuals)/real(n, dp)
      do i = 1, lmax
         covariance = dot_product(residuals(i + 1:), residuals(:n - i))
         weight = 1.0_dp - real(i, dp)/real(lmax + 1, dp)
         denominator = denominator + 2.0_dp*weight*covariance/real(n, dp)
      end do
      out%statistic = numerator/denominator
      out%residuals = residuals
      out%lags = lmax
      out%model = kind(1:min(3, len_trim(kind)))
   end function kpss_test

   pure function pp_test(x, statistic_type, model, lag_rule, use_lag) result(out)
      ! Perform the urca Phillips-Perron Z-tau or Z-alpha test.
      real(dp), intent(in) :: x(:)
      character(len=*), intent(in), optional :: statistic_type, model, lag_rule
      integer, intent(in), optional :: use_lag
      type(pp_result_t) :: out
      real(dp), allocatable :: y(:), ylag(:), design(:, :), beta(:), se(:), residuals(:)
      real(dp) :: rss, s, sigma, lambda, lambda_prime, covariance, weight
      real(dp) :: my, myy, myybar, mty, moment, tstat, my_stat, beta_stat
      integer :: i, n, lmax, info
      character(len=8) :: kind, regression, rule

      kind = 'Z-tau'
      if (present(statistic_type)) kind = adjustl(statistic_type)
      regression = 'constant'
      if (present(model)) regression = adjustl(model)
      rule = 'short'
      if (present(lag_rule)) rule = lowercase(adjustl(lag_rule))
      out%statistic_type = kind
      out%model = regression
      n = size(x) - 1
      if (n < 3) then
         out%info = 1
         return
      end if
      if (trim(kind) /= 'Z-tau' .and. trim(kind) /= 'Z-alpha') then
         out%info = 2
         return
      end if
      if (trim(regression) /= 'constant' .and. trim(regression) /= 'trend') then
         out%info = 3
         return
      end if
      if (present(use_lag)) then
         lmax = use_lag
      else if (trim(rule) == 'short') then
         lmax = int(4.0_dp*(real(n, dp)/100.0_dp)**0.25_dp)
      else if (trim(rule) == 'long') then
         lmax = int(12.0_dp*(real(n, dp)/100.0_dp)**0.25_dp)
      else
         out%info = 4
         return
      end if
      if (lmax < 0 .or. lmax >= n) then
         out%info = 5
         return
      end if
      allocate (y(n), ylag(n))
      y = x(2:)
      ylag = x(:size(x) - 1)
      if (trim(regression) == 'constant') then
         allocate (design(n, 2))
         design(:, 1) = 1.0_dp
         design(:, 2) = ylag
      else
         allocate (design(n, 3))
         design(:, 1) = 1.0_dp
         design(:, 2) = ylag
         design(:, 3) = [(real(i, dp) - real(n, dp)/2.0_dp, i=1, n)]
      end if
      call ols_fit(design, y, beta, se, residuals, rss, info)
      if (info /= 0) then
         out%info = 10 + info
         return
      end if
      s = rss/real(n, dp)
      sigma = s
      do i = 1, lmax
         covariance = dot_product(residuals(i + 1:), residuals(:n - i))
         weight = 1.0_dp - real(i, dp)/real(lmax + 1, dp)
         sigma = sigma + 2.0_dp*weight*covariance/real(n, dp)
      end do
      lambda = 0.5_dp*(sigma - s)
      lambda_prime = lambda/sigma
      my = sum(y)/real(n, dp)**1.5_dp
      myy = sum(y*y)/real(n*n, dp)
      myybar = sum((y - sum(y)/real(n, dp))**2)/real(n*n, dp)
      if (trim(regression) == 'constant') then
         my_stat = sqrt(s/sigma)*(beta(1)/se(1)) + lambda_prime*sqrt(sigma)*my/(sqrt(myy)*sqrt(myybar))
         allocate (out%auxiliary(1))
         out%auxiliary(1) = round_four(my_stat)
         if (trim(kind) == 'Z-tau') then
            tstat = (beta(2) - 1.0_dp)/se(2)
            out%statistic = sqrt(s/sigma)*tstat - lambda_prime*sqrt(sigma)/sqrt(myybar)
         else
            out%statistic = real(n, dp)*(beta(2) - 1.0_dp) - lambda/myybar
         end if
         out%critical_values = [ &
                               -3.4335_dp - 5.999_dp/n - 29.25_dp/(n*n), &
                               -2.8621_dp - 2.738_dp/n - 8.36_dp/(n*n), &
                               -2.5671_dp - 1.438_dp/n - 4.48_dp/(n*n)]
      else
         mty = sum([(real(i, dp)*y(i), i=1, n)])/real(n, dp)**2.5_dp
         moment = (1.0_dp - 1.0_dp/real(n*n, dp))*myy - 12.0_dp*mty*mty &
                  + 12.0_dp*(1.0_dp + 1.0_dp/n)*mty*my &
                  - (4.0_dp + 6.0_dp/n + 2.0_dp/(n*n))*my*my
         my_stat = sqrt(s/sigma)*(beta(1)/se(1)) - lambda_prime*sqrt(sigma)*my &
                   /(sqrt(moment)*sqrt(moment + my*my))
         beta_stat = sqrt(s/sigma)*(beta(3)/se(3)) - lambda_prime*sqrt(sigma)*(0.5_dp*my - mty) &
                     /(sqrt(moment/12.0_dp)*sqrt(myybar))
         allocate (out%auxiliary(2))
         out%auxiliary = [round_four(my_stat), round_four(beta_stat)]
         if (trim(kind) == 'Z-tau') then
            tstat = (beta(2) - 1.0_dp)/se(2)
            out%statistic = sqrt(s/sigma)*tstat - lambda_prime*sqrt(sigma)/sqrt(moment)
         else
            out%statistic = real(n, dp)*(beta(2) - 1.0_dp) - lambda/moment
         end if
         out%critical_values = [ &
                               -3.9638_dp - 8.353_dp/n - 47.44_dp/(n*n), &
                               -3.4126_dp - 4.039_dp/n - 17.83_dp/(n*n), &
                               -3.1279_dp - 2.418_dp/n - 7.58_dp/(n*n)]
      end if
      if (trim(kind) == 'Z-alpha') out%critical_values = ieee_value(0.0_dp, ieee_quiet_nan)
      out%coefficients = beta
      out%standard_errors = se
      out%residuals = residuals
      out%lags = lmax
   end function pp_test

   pure function ers_test(y, statistic_type, model, lag_max) result(out)
      ! Perform the Elliott-Rothenberg-Stock DF-GLS or point-optimal test.
      real(dp), intent(in) :: y(:)
      character(len=*), intent(in), optional :: statistic_type, model
      integer, intent(in), optional :: lag_max
      type(ers_result_t) :: out
      real(dp), allocatable :: ya(:), zdet(:, :), beta(:), se(:), res(:), design(:, :), dependent(:)
      real(dp), allocatable :: null_res(:), criteria(:)
      real(dp) :: ahat, rss, sig_null, sig_res, sum_lag, scale
      integer :: i, n, max_lag, q, selected, info, row
      character(len=8) :: kind, regression

      kind = 'DF-GLS'
      if (present(statistic_type)) kind = adjustl(statistic_type)
      regression = 'constant'
      if (present(model)) regression = adjustl(model)
      max_lag = 4
      if (present(lag_max)) max_lag = lag_max
      out%statistic_type = kind
      out%model = regression
      if (max_lag < 0 .or. size(y) < max_lag + 4) then
         out%info = 1
         return
      end if
      if (trim(kind) /= 'DF-GLS' .and. trim(kind) /= 'P-test') then
         out%info = 2
         return
      end if
      if (trim(regression) /= 'constant' .and. trim(regression) /= 'trend') then
         out%info = 3
         return
      end if
      n = size(y)
      if (n < 50) then
         row = 1
      else if (n < 100) then
         row = 2
      else if (n <= 200) then
         row = 3
      else
         row = 4
      end if
      if (trim(regression) == 'constant') then
         ahat = 1.0_dp - 7.0_dp/real(n, dp)
         allocate (ya(n), zdet(n, 1))
         ya(1) = y(1)
         ya(2:) = y(2:) - ahat*y(:n - 1)
         zdet(1, 1) = 1.0_dp
         zdet(2:, 1) = 1.0_dp - ahat
      else
         ahat = 1.0_dp - 13.5_dp/real(n, dp)
         allocate (ya(n), zdet(n, 2))
         ya(1) = y(1)
         ya(2:) = y(2:) - ahat*y(:n - 1)
         zdet(1, :) = 1.0_dp
         zdet(2:, 1) = 1.0_dp - ahat
         do i = 2, n
            zdet(i, 2) = real(i, dp) - ahat*real(i - 1, dp)
         end do
      end if
      call ols_fit(zdet, ya, beta, se, res, rss, info)
      if (info /= 0) then
         out%info = 10 + info
         return
      end if
      allocate (out%detrended(n))
      if (trim(regression) == 'constant') then
         out%detrended = y - beta(1)
      else
         out%detrended = y - beta(1) - beta(2)*[(real(i, dp), i=1, n)]
      end if
      if (trim(kind) == 'DF-GLS') then
         call build_ers_regression(out%detrended, max_lag, design, dependent)
         call ols_fit(design, dependent, beta, se, res, rss, info)
         if (info /= 0) then
            out%info = 20 + info
            return
         end if
         out%statistic = beta(1)/se(1)
         out%lags = max_lag
         if (trim(regression) == 'constant') then
            out%critical_values = [ &
                                  -2.5658_dp - 1.960_dp/n - 10.04_dp/(n*n), &
                                  -1.9393_dp - 0.398_dp/n, &
                                  -1.6156_dp - 0.181_dp/n]
         else
            call ers_trend_critical(row, out%critical_values)
         end if
         out%critical_values = anint(100.0_dp*out%critical_values)/100.0_dp
         out%coefficients = beta
         out%standard_errors = se
         out%residuals = res
      else
         allocate (null_res(n))
         null_res(1) = 0.0_dp
         null_res(2:) = y(2:) - y(:n - 1)
         if (trim(regression) == 'trend') null_res = null_res - sum(null_res)/real(n, dp)
         sig_null = sum(null_res*null_res)
         sig_res = sum(res*res)
         selected = 0
         if (max_lag > 0) then
            allocate (criteria(max_lag))
            do q = 1, max_lag
               call build_ers_level_regression(y, q, design, dependent)
               call ols_fit(design, dependent, beta, se, res, rss, info)
               criteria(q) = huge(1.0_dp)
               if (info == 0) criteria(q) = real(size(dependent), dp)*log(rss/size(dependent)) &
                                            + log(real(size(dependent), dp))*real(size(beta) + 1, dp)
            end do
            selected = minloc(criteria, dim=1)
         end if
         call build_ers_level_regression(y, selected, design, dependent)
         call ols_fit(design, dependent, beta, se, res, rss, info)
         if (info /= 0) then
            out%info = 30 + info
            return
         end if
         sum_lag = 0.0_dp
         if (selected > 0) sum_lag = sum(beta(3:))
         scale = (rss/real(size(dependent) - size(beta), dp))/(1.0_dp - sum_lag)**2
         out%statistic = (sig_res - ahat*sig_null)/scale
         out%lags = selected
         call ers_p_critical(row, regression, out%critical_values)
         out%critical_values = anint(100.0_dp*out%critical_values)/100.0_dp
      end if
   end function ers_test

   pure function za_test(y, model, lag) result(out)
      ! Perform the Zivot-Andrews unit-root test with one endogenous break.
      real(dp), intent(in) :: y(:)
      character(len=*), intent(in), optional :: model
      integer, intent(in), optional :: lag
      type(za_result_t) :: out
      real(dp), allocatable :: design(:, :), dependent(:), beta(:), se(:), res(:)
      real(dp) :: rss, value
      integer :: n, q, break_at, info
      character(len=9) :: kind

      kind = 'intercept'
      if (present(model)) kind = adjustl(model)
      q = 0
      if (present(lag)) q = lag
      out%model = kind
      out%lags = q
      n = size(y)
      if (q < 0 .or. n < q + 5) then
         out%info = 1
         return
      end if
      if (trim(kind) /= 'intercept' .and. trim(kind) /= 'trend' .and. trim(kind) /= 'both') then
         out%info = 2
         return
      end if
      allocate (out%break_statistics(n - 1))
      out%break_statistics = huge(1.0_dp)
      do break_at = 1, n - 1
         call build_za_regression(y, q, kind, break_at, design, dependent)
         call ols_fit(design, dependent, beta, se, res, rss, info)
         if (info /= 0) cycle
         value = (beta(2) - 1.0_dp)/se(2)
         out%break_statistics(break_at) = value
      end do
      out%break_point = minloc(out%break_statistics, dim=1)
      if (out%break_statistics(out%break_point) >= huge(1.0_dp)) then
         out%info = 3
         return
      end if
      out%statistic = out%break_statistics(out%break_point)
      call build_za_regression(y, q, kind, out%break_point, design, dependent)
      call ols_fit(design, dependent, beta, se, res, rss, info)
      out%coefficients = beta
      out%standard_errors = se
      out%residuals = res
      select case (trim(kind))
      case ('intercept')
         out%critical_values = [-5.34_dp, -4.80_dp, -4.58_dp]
      case ('trend')
         out%critical_values = [-4.93_dp, -4.42_dp, -4.11_dp]
      case default
         out%critical_values = [-5.57_dp, -5.08_dp, -4.82_dp]
      end select
   end function za_test

   pure function johansen_test(x, test_type, deterministic, lag, specification, season, exogenous) result(out)
      ! Estimate the Johansen cointegration system used by urca ca.jo.
      real(dp), intent(in) :: x(:, :)
      character(len=*), intent(in), optional :: test_type, deterministic, specification
      integer, intent(in), optional :: lag
      integer, intent(in), optional :: season
      real(dp), intent(in), optional :: exogenous(:, :)
      type(johansen_result_t) :: out
      real(dp), allocatable :: dx(:, :), z0(:, :), z1(:, :), z1base(:, :), zk(:, :)
      real(dp), allocatable :: seasonal(:, :)
      real(dp), allocatable :: m00(:, :), m11(:, :), mkk(:, :), m01(:, :), m0k(:, :), mk1(:, :)
      real(dp), allocatable :: m11inv(:, :), s00(:, :), s0k(:, :), skk(:, :), s00inv(:, :), skkinv(:, :)
      real(dp), allocatable :: chol(:, :), cholinv(:, :), eigmat(:, :), eigvec(:, :), raw_v(:, :)
      real(dp), allocatable :: middle(:, :), middle_inv(:, :)
      real(dp) :: sample
      integer :: n, p, k, rows, pe, i, j, block, info
      integer :: seasonal_columns, exogenous_columns, column
      character(len=10) :: test, det, spec

      test = 'trace'
      if (present(test_type)) test = lowercase(adjustl(test_type))
      det = 'none'
      if (present(deterministic)) det = lowercase(adjustl(deterministic))
      spec = 'longrun'
      if (present(specification)) spec = lowercase(adjustl(specification))
      k = 2
      if (present(lag)) k = lag
      if (trim(test) == 'trace') then
         out%test_type = 'trace'
      else
         out%test_type = 'eigen'
      end if
      out%deterministic = det
      out%specification = spec
      out%lag = k
      n = size(x, 1)
      p = size(x, 2)
      out%variables = p
      if (k < 2 .or. n <= k + 2 .or. p < 2) then
         out%info = 1
         return
      end if
      if (trim(test) /= 'trace' .and. trim(test) /= 'eigen') then
         out%info = 2
         return
      end if
      if (trim(det) /= 'none' .and. trim(det) /= 'const' .and. trim(det) /= 'trend') then
         out%info = 3
         return
      end if
      if (trim(spec) /= 'longrun' .and. trim(spec) /= 'transitory') then
         out%info = 4
         return
      end if
      seasonal_columns = 0
      if (present(season)) then
         if (season < 2) then
            out%info = 5
            return
         end if
         seasonal_columns = season - 1
      end if
      exogenous_columns = 0
      if (present(exogenous)) then
         if (size(exogenous, 1) /= n) then
            out%info = 6
            return
         end if
         exogenous_columns = size(exogenous, 2)
      end if
      rows = n - k
      allocate (dx(n - 1, p), z0(rows, p))
      dx = x(2:, :) - x(:n - 1, :)
      z0 = dx(k:, :)
      if (trim(det) == 'const') then
         allocate (z1base(rows, p*(k - 1)))
      else
         allocate (z1base(rows, 1 + p*(k - 1)))
         z1base(:, 1) = 1.0_dp
      end if
      do block = 1, k - 1
         if (trim(det) == 'const') then
            z1base(:, (block - 1)*p + 1:block*p) = dx(k - block:n - 1 - block, :)
         else
            z1base(:, 2 + (block - 1)*p:1 + block*p) = dx(k - block:n - 1 - block, :)
         end if
      end do
      if (seasonal_columns > 0) then
         allocate (seasonal(rows, seasonal_columns))
         seasonal = -1.0_dp/real(season, dp)
         do i = 1, rows
            j = 1 + mod(k + i - 1, season)
            if (j <= seasonal_columns) seasonal(i, j) = seasonal(i, j) + 1.0_dp
         end do
      end if
      allocate (z1(rows, size(z1base, 2) + seasonal_columns + exogenous_columns))
      if (trim(det) == 'const') then
         column = 1
         if (seasonal_columns > 0) then
            z1(:, column:column + seasonal_columns - 1) = seasonal
            column = column + seasonal_columns
         end if
         if (exogenous_columns > 0) then
            z1(:, column:column + exogenous_columns - 1) = exogenous(k + 1:, :)
            column = column + exogenous_columns
         end if
         z1(:, column:) = z1base
      else
         z1(:, 1) = z1base(:, 1)
         column = 2
         if (seasonal_columns > 0) then
            z1(:, column:column + seasonal_columns - 1) = seasonal
            column = column + seasonal_columns
         end if
         if (exogenous_columns > 0) then
            z1(:, column:column + exogenous_columns - 1) = exogenous(k + 1:, :)
            column = column + exogenous_columns
         end if
         z1(:, column:) = z1base(:, 2:)
      end if
      pe = p
      if (trim(det) == 'const' .or. trim(det) == 'trend') pe = p + 1
      allocate (zk(rows, pe))
      if (trim(spec) == 'longrun') then
         zk(:, :p) = x(1:n - k, :)
         if (trim(det) == 'trend') zk(:, pe) = [(real(i, dp), i=1, n - k)]
      else
         zk(:, :p) = x(k:n - 1, :)
         if (trim(det) == 'trend') zk(:, pe) = [(real(k + i - 1, dp), i=1, rows)]
      end if
      if (trim(det) == 'const') zk(:, pe) = 1.0_dp
      sample = real(rows, dp)
      m00 = matmul(transpose(z0), z0)/sample
      m11 = matmul(transpose(z1), z1)/sample
      mkk = matmul(transpose(zk), zk)/sample
      m01 = matmul(transpose(z0), z1)/sample
      m0k = matmul(transpose(z0), zk)/sample
      mk1 = matmul(transpose(zk), z1)/sample
      call invert_matrix(m11, m11inv, info)
      if (info /= 0) then
         out%info = 10 + info
         return
      end if
      out%r0 = z0 - matmul(z1, matmul(m11inv, transpose(m01)))
      out%rk = zk - matmul(z1, matmul(m11inv, transpose(mk1)))
      s00 = m00 - matmul(m01, matmul(m11inv, transpose(m01)))
      s0k = m0k - matmul(m01, matmul(m11inv, transpose(mk1)))
      skk = mkk - matmul(mk1, matmul(m11inv, transpose(mk1)))
      call invert_matrix(s00, s00inv, info)
      if (info /= 0) then
         out%info = 20 + info
         return
      end if
      call cholesky_lower(skk, chol, info)
      if (info /= 0) then
         out%info = 30 + info
         return
      end if
      call invert_matrix(chol, cholinv, info)
      eigmat = matmul(cholinv, matmul(transpose(s0k), matmul(s00inv, matmul(s0k, transpose(cholinv)))))
      eigmat = 0.5_dp*(eigmat + transpose(eigmat))
      call symmetric_eigen(eigmat, out%eigenvalues, eigvec, info)
      if (info /= 0) then
         out%info = 40 + info
         return
      end if
      raw_v = matmul(transpose(cholinv), eigvec)
      allocate (out%cointegration(pe, pe))
      out%cointegration = raw_v
      do j = 1, pe
         if (abs(out%cointegration(1, j)) <= epsilon(1.0_dp)) then
            out%info = 50 + j
            return
         end if
         out%cointegration(:, j) = out%cointegration(:, j)/out%cointegration(1, j)
      end do
      middle = matmul(transpose(out%cointegration), matmul(skk, out%cointegration))
      call invert_matrix(middle, middle_inv, info)
      out%loading = matmul(s0k, matmul(out%cointegration, middle_inv))
      call invert_matrix(skk, skkinv, info)
      out%pi = matmul(s0k, skkinv)
      out%delta = s00 - matmul(s0k, matmul(out%cointegration, &
                                           matmul(middle_inv, matmul(transpose(out%cointegration), transpose(s0k)))))
      out%gamma = matmul(m01, m11inv) - matmul(out%pi, matmul(mk1, m11inv))
      allocate (out%statistic(p), out%critical_values(p, 3))
      do i = 1, p
         if (trim(test) == 'trace') then
            out%statistic(p - i + 1) = -sample*sum(log(1.0_dp - out%eigenvalues(i:pe)))
         else
            out%statistic(p - i + 1) = -sample*log(1.0_dp - out%eigenvalues(i))
         end if
      end do
      call johansen_critical_values(p, det, test, out%critical_values)
   end function johansen_test

   pure subroutine build_adf_matrix(y, dy, kind, max_lag, q, x, zlag_column)
      ! Build an ADF regression on the common maximum-lag sample.
      real(dp), intent(in) :: y(:), dy(:)
      character(len=*), intent(in) :: kind
      integer, intent(in) :: max_lag, q
      real(dp), allocatable, intent(out) :: x(:, :)
      integer, intent(out) :: zlag_column
      integer :: i, nobs, column

      nobs = size(dy) - max_lag
      select case (trim(kind))
      case ('none')
         allocate (x(nobs, 1 + q))
         zlag_column = 1
         x(:, 1) = y(max_lag + 1:size(y) - 1)
         column = 2
      case ('drift')
         allocate (x(nobs, 2 + q))
         x(:, 1) = 1.0_dp
         zlag_column = 2
         x(:, 2) = y(max_lag + 1:size(y) - 1)
         column = 3
      case default
         allocate (x(nobs, 3 + q))
         x(:, 1) = 1.0_dp
         zlag_column = 2
         x(:, 2) = y(max_lag + 1:size(y) - 1)
         x(:, 3) = [(real(max_lag + i, dp), i=1, nobs)]
         column = 4
      end select
      do i = 1, q
         x(:, column + i - 1) = dy(max_lag + 1 - i:size(dy) - i)
      end do
   end subroutine build_adf_matrix

   pure subroutine restricted_diff_matrix(dy, max_lag, q, intercept, x)
      ! Build the restricted ADF regression used by the auxiliary phi tests.
      real(dp), intent(in) :: dy(:)
      integer, intent(in) :: max_lag, q
      logical, intent(in) :: intercept
      real(dp), allocatable, intent(out) :: x(:, :)
      integer :: i, offset, nobs

      nobs = size(dy) - max_lag
      offset = merge(1, 0, intercept)
      allocate (x(nobs, q + offset))
      if (intercept) x(:, 1) = 1.0_dp
      do i = 1, q
         x(:, offset + i) = dy(max_lag + 1 - i:size(dy) - i)
      end do
   end subroutine restricted_diff_matrix

   pure subroutine build_ers_regression(y, q, design, dependent)
      ! Build the no-intercept DF-GLS regression with q lagged differences.
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: q
      real(dp), allocatable, intent(out) :: design(:, :), dependent(:)
      real(dp), allocatable :: dy(:)
      integer :: i, nobs

      allocate (dy(size(y) - 1))
      dy = y(2:) - y(:size(y) - 1)
      nobs = size(dy) - q
      allocate (dependent(nobs), design(nobs, 1 + q))
      dependent = dy(q + 1:)
      design(:, 1) = y(q + 1:size(y) - 1)
      do i = 1, q
         design(:, i + 1) = dy(q + 1 - i:size(dy) - i)
      end do
   end subroutine build_ers_regression

   pure subroutine build_ers_level_regression(y, q, design, dependent)
      ! Build the intercept, lagged-level, and lagged-difference ERS scale regression.
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: q
      real(dp), allocatable, intent(out) :: design(:, :), dependent(:)
      real(dp), allocatable :: dy(:)
      integer :: i, nobs

      allocate (dy(size(y) - 1))
      dy = y(2:) - y(:size(y) - 1)
      nobs = size(dy) - q
      allocate (dependent(nobs), design(nobs, 2 + q))
      dependent = dy(q + 1:)
      design(:, 1) = 1.0_dp
      design(:, 2) = y(q + 1:size(y) - 1)
      do i = 1, q
         design(:, i + 2) = dy(q + 1 - i:size(dy) - i)
      end do
   end subroutine build_ers_level_regression

   pure subroutine build_za_regression(y, q, kind, break_at, design, dependent)
      ! Build one complete-case Zivot-Andrews candidate-break regression.
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: q, break_at
      character(len=*), intent(in) :: kind
      real(dp), allocatable, intent(out) :: design(:, :), dependent(:)
      real(dp), allocatable :: dy(:)
      integer :: first, nobs, columns, i, t, column

      allocate (dy(size(y) - 1))
      dy = y(2:) - y(:size(y) - 1)
      first = q + 2
      nobs = size(y) - first + 1
      columns = 4 + q
      if (trim(kind) == 'both') columns = columns + 1
      allocate (dependent(nobs), design(nobs, columns))
      dependent = y(first:)
      design(:, 1) = 1.0_dp
      design(:, 2) = y(first - 1:size(y) - 1)
      design(:, 3) = [(real(t, dp), t=first, size(y))]
      do i = 1, q
         design(:, 3 + i) = dy(first - 1 - i:size(dy) - i)
      end do
      column = 4 + q
      select case (trim(kind))
      case ('intercept')
         do t = first, size(y)
            design(t - first + 1, column) = merge(1.0_dp, 0.0_dp, t > break_at)
         end do
      case ('trend')
         do t = first, size(y)
            design(t - first + 1, column) = merge(real(t - break_at, dp), 0.0_dp, t > break_at)
         end do
      case default
         do t = first, size(y)
            design(t - first + 1, column) = merge(1.0_dp, 0.0_dp, t > break_at)
            design(t - first + 1, column + 1) = merge(real(t - break_at, dp), 0.0_dp, t > break_at)
         end do
      end select
   end subroutine build_za_regression

   pure subroutine ers_trend_critical(row, values)
      ! Return ERS DF-GLS trend critical values for one sample-size row.
      integer, intent(in) :: row
      real(dp), intent(out) :: values(3)
      real(dp), parameter :: table(4, 3) = reshape([ &
                                                        -3.77_dp, -3.58_dp, -3.46_dp, -3.48_dp, &
                                                        -3.19_dp, -3.03_dp, -2.93_dp, -2.89_dp, &
                                                    -2.89_dp, -2.74_dp, -2.64_dp, -2.57_dp], [4, 3])
      values = table(row, :)
   end subroutine ers_trend_critical

   pure subroutine ers_p_critical(row, model, values)
      ! Return ERS point-optimal critical values for model and sample size.
      integer, intent(in) :: row
      character(len=*), intent(in) :: model
      real(dp), intent(out) :: values(3)
      real(dp), parameter :: constant_table(4, 3) = reshape([ &
                                                               1.87_dp, 1.95_dp, 1.91_dp, 1.99_dp, &
                                                               2.97_dp, 3.11_dp, 3.17_dp, 3.26_dp, &
                                                        3.91_dp, 4.17_dp, 4.33_dp, 4.48_dp], [4, 3])
      real(dp), parameter :: trend_table(4, 3) = reshape([ &
                                                              4.22_dp, 4.26_dp, 4.05_dp, 3.96_dp, &
                                                              5.72_dp, 5.64_dp, 5.66_dp, 5.62_dp, &
                                                        6.77_dp, 6.79_dp, 6.86_dp, 6.89_dp], [4, 3])
      if (trim(model) == 'constant') then
         values = constant_table(row, :)
      else
         values = trend_table(row, :)
      end if
   end subroutine ers_p_critical



   pure subroutine johansen_critical_values(p, deterministic, test_type, values)
      ! Return the Osterwald-Lenum critical values embedded in urca ca.jo.
      integer, intent(in) :: p
      character(len=*), intent(in) :: deterministic, test_type
      real(dp), intent(out) :: values(p, 3)
      real(dp), parameter :: none_eigen(11, 3) = reshape([ &
                                           6.50, 12.91, 18.90, 24.78, 30.84, 36.25, 42.06, 48.43, 54.01, 59.00, 65.07, &
                                           8.18, 14.90, 21.07, 27.14, 33.32, 39.43, 44.91, 51.07, 57.00, 62.42, 68.27, &
                                  11.65, 19.19, 25.75, 32.14, 38.78, 44.59, 51.30, 57.07, 63.37, 68.61, 74.36], [11, 3])
      real(dp), parameter :: none_trace(11, 3) = reshape([ &
                                      6.50, 15.66, 28.71, 45.23, 66.49, 85.18, 118.99, 151.38, 186.54, 226.34, 269.53, &
                                      8.18, 17.95, 31.52, 48.28, 70.60, 90.39, 124.25, 157.11, 192.84, 232.49, 277.39, &
                            11.65, 23.52, 37.22, 55.43, 78.87, 104.20, 136.06, 168.92, 204.79, 246.27, 292.65], [11, 3])
      real(dp), parameter :: const_eigen(11, 3) = reshape([ &
                                           7.52, 13.75, 19.77, 25.56, 31.66, 37.45, 43.25, 48.91, 54.35, 60.25, 66.02, &
                                           9.24, 15.67, 22.00, 28.14, 34.40, 40.30, 46.45, 52.00, 57.42, 63.57, 69.74, &
                                  12.97, 20.20, 26.81, 33.24, 39.79, 46.82, 51.91, 57.95, 63.71, 69.94, 76.63], [11, 3])
      real(dp), parameter :: const_trace(11, 3) = reshape([ &
                                      7.52, 17.85, 32.00, 49.65, 71.86, 97.18, 126.58, 159.48, 196.37, 236.54, 282.45, &
                                     9.24, 19.96, 34.91, 53.12, 76.07, 102.14, 131.70, 165.58, 202.92, 244.15, 291.40, &
                            12.97, 24.60, 41.07, 60.16, 84.45, 111.01, 143.09, 177.20, 215.74, 257.68, 307.64], [11, 3])
      real(dp), parameter :: trend_eigen(11, 3) = reshape([ &
                                          10.49, 16.85, 23.11, 29.12, 34.75, 40.91, 46.32, 52.16, 57.87, 63.18, 69.26, &
                                          12.25, 18.96, 25.54, 31.46, 37.52, 43.97, 49.42, 55.50, 61.29, 66.23, 72.72, &
                                  16.26, 23.65, 30.34, 36.65, 42.36, 49.51, 54.71, 62.46, 67.88, 73.73, 79.23], [11, 3])
      real(dp), parameter :: trend_trace(11, 3) = reshape([ &
                                    10.49, 22.76, 39.06, 59.14, 83.20, 110.42, 141.01, 176.67, 215.17, 256.72, 303.13, &
                                    12.25, 25.32, 42.44, 62.99, 87.31, 114.90, 146.76, 182.82, 222.21, 263.42, 310.81, &
                            16.26, 30.45, 48.45, 70.05, 96.58, 124.75, 158.49, 196.08, 234.41, 279.07, 327.45], [11, 3])

      if (p > 11) then
         values = ieee_value(0.0, ieee_quiet_nan)
      else if (trim(deterministic) == 'none' .and. trim(test_type) == 'eigen') then
         values = none_eigen(:p, :)
      else if (trim(deterministic) == 'none') then
         values = none_trace(:p, :)
      else if (trim(deterministic) == 'const' .and. trim(test_type) == 'eigen') then
         values = const_eigen(:p, :)
      else if (trim(deterministic) == 'const') then
         values = const_trace(:p, :)
      else if (trim(test_type) == 'eigen') then
         values = trend_eigen(:p, :)
      else
         values = trend_trace(:p, :)
      end if
   end subroutine johansen_critical_values



   pure elemental real(dp) function f_statistic(rss_restricted, rss_full, p_full, p_restricted, nobs) result(value)
      ! Compute the nested-regression F statistic used by ur.df.
      real(dp), intent(in) :: rss_restricted, rss_full
      integer, intent(in) :: p_full, p_restricted, nobs
      value = ((rss_restricted - rss_full)/real(p_full - p_restricted, dp)) &
              /(rss_full/real(nobs - p_full, dp))
   end function f_statistic

   pure subroutine adf_critical_values(kind, n, values)
      ! Select urca ADF critical values by model and effective sample size.
      character(len=*), intent(in) :: kind
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: values(:, :)
      real(dp), parameter :: tau1(6, 3) = reshape([ &
                             -2.66_dp, -2.62_dp, -2.60_dp, -2.58_dp, -2.58_dp, -2.58_dp, &
                             -1.95_dp, -1.95_dp, -1.95_dp, -1.95_dp, -1.95_dp, -1.95_dp, &
                      -1.60_dp, -1.61_dp, -1.61_dp, -1.62_dp, -1.62_dp, -1.62_dp], [6, 3])
      real(dp), parameter :: drift(6, 6) = reshape([ &
                             -3.75_dp, -3.58_dp, -3.51_dp, -3.46_dp, -3.44_dp, -3.43_dp, &
                             -3.00_dp, -2.93_dp, -2.89_dp, -2.88_dp, -2.87_dp, -2.86_dp, &
                             -2.63_dp, -2.60_dp, -2.58_dp, -2.57_dp, -2.57_dp, -2.57_dp, &
                                   7.88_dp, 7.06_dp, 6.70_dp, 6.52_dp, 6.47_dp, 6.43_dp, &
                                   5.18_dp, 4.86_dp, 4.71_dp, 4.63_dp, 4.61_dp, 4.59_dp, &
                            4.12_dp, 3.94_dp, 3.86_dp, 3.81_dp, 3.79_dp, 3.78_dp], [6, 6])
      real(dp), parameter :: trend(6, 9) = reshape([ &
                             -4.38_dp, -4.15_dp, -4.04_dp, -3.99_dp, -3.98_dp, -3.96_dp, &
                             -3.60_dp, -3.50_dp, -3.45_dp, -3.43_dp, -3.42_dp, -3.41_dp, &
                             -3.24_dp, -3.18_dp, -3.15_dp, -3.13_dp, -3.13_dp, -3.12_dp, &
                                   8.21_dp, 7.02_dp, 6.50_dp, 6.22_dp, 6.15_dp, 6.09_dp, &
                                   5.68_dp, 5.13_dp, 4.88_dp, 4.75_dp, 4.71_dp, 4.68_dp, &
                                   4.67_dp, 4.31_dp, 4.16_dp, 4.07_dp, 4.05_dp, 4.03_dp, &
                                  10.61_dp, 9.31_dp, 8.73_dp, 8.43_dp, 8.34_dp, 8.27_dp, &
                                   7.24_dp, 6.73_dp, 6.49_dp, 6.49_dp, 6.30_dp, 6.25_dp, &
                            5.91_dp, 5.61_dp, 5.47_dp, 5.47_dp, 5.36_dp, 5.34_dp], [6, 9])
      integer :: row

      if (n < 25) then
         row = 1
      else if (n < 50) then
         row = 2
      else if (n < 100) then
         row = 3
      else if (n < 250) then
         row = 4
      else if (n < 500) then
         row = 5
      else
         row = 6
      end if
      select case (trim(kind))
      case ('none')
         allocate (values(1, 3))
         values(1, :) = tau1(row, :)
      case ('drift')
         allocate (values(2, 3))
         values(1, :) = drift(row, 1:3)
         values(2, :) = drift(row, 4:6)
      case default
         allocate (values(3, 3))
         values(1, :) = trend(row, 1:3)
         values(2, :) = trend(row, 4:6)
         values(3, :) = trend(row, 7:9)
      end select
   end subroutine adf_critical_values



   pure elemental real(dp) function round_four(value) result(rounded)
      ! Round a value to four decimal places as in ur.pp auxiliary output.
      real(dp), intent(in) :: value
      rounded = anint(10000.0_dp*value)/10000.0_dp
   end function round_four
end module urca_mod
