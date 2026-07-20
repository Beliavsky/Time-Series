! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Distinct algorithms translated from the R VARshrink package.
module varshrink_mod
   use kind_mod, only: dp
   use linalg_mod, only: symmetric_eigen, symmetric_pseudoinverse, identity_matrix
   use linalg_mod, only: inverse_logdet
   use stats_mod, only: median
   implicit none
   private

   type, public :: varshrink_ridge_path_t
      ! Multivariate ridge estimates and generalized cross-validation scores.
      real(dp), allocatable :: coefficient(:, :, :)
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: gcv(:)
      real(dp), allocatable :: effective_df(:)
      real(dp), allocatable :: predictor_mean(:)
      real(dp), allocatable :: predictor_scale(:)
      real(dp), allocatable :: response_mean(:)
      integer :: selected = 0
      integer :: info = 0
      logical :: scaled = .false.
   end type varshrink_ridge_path_t

   type, public :: varshrink_semibayes_fit_t
      ! Semiparametric Bayesian coefficient and noise-covariance estimates.
      real(dp), allocatable :: coefficient(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: degrees_of_freedom = huge(1.0_dp)
      integer :: iterations = 0
      integer :: info = 0
      logical :: conjugate = .false.
      logical :: converged = .false.
      logical :: weights_estimated = .false.
   end type varshrink_semibayes_fit_t

   type, public :: varshrink_var_fit_t
      ! VAR-form ridge estimates with lag matrices and in-sample diagnostics.
      real(dp), allocatable :: ar(:, :)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: gcv = huge(1.0_dp)
      real(dp) :: effective_df = 0.0_dp
      integer :: order = 0
      integer :: info = 0
      logical :: include_intercept = .true.
   end type varshrink_var_fit_t

   type, public :: varshrink_covariance_t
      ! James-Stein covariance estimate and its two shrinkage intensities.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: variance(:)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: lambda_variance = 0.0_dp
      integer :: info = 0
      logical :: lambda_estimated = .false.
      logical :: lambda_variance_estimated = .false.
   end type varshrink_covariance_t

   type, public :: varshrink_nonparametric_fit_t
      ! Centered nonparametric shrinkage regression with recovered intercept.
      real(dp), allocatable :: coefficient(:, :)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: lambda_variance = 0.0_dp
      integer :: info = 0
   end type varshrink_nonparametric_fit_t

   public :: varshrink_multivariate_ridge
   public :: varshrink_semibayes
   public :: varshrink_var_ridge
   public :: varshrink_effective_df
   public :: varshrink_log_likelihood
   public :: varshrink_variance_intensity
   public :: varshrink_covariance_shrink
   public :: varshrink_nonparametric

contains

   pure function varshrink_multivariate_ridge(response, predictors, lambda, do_scale) result(out)
      !! Fit multivariate ridge regressions and select the minimum-GCV estimate.
      real(dp), intent(in) :: response(:, :) !! N-by-K matrix of dependent variables.
      real(dp), intent(in) :: predictors(:, :) !! N-by-M matrix of regressors.
      real(dp), intent(in) :: lambda(:) !! Candidate nonnegative ridge parameters.
      logical, intent(in), optional :: do_scale !! Center and scale predictors when true.
      type(varshrink_ridge_path_t) :: out
      real(dp), allocatable :: x(:, :), y(:, :), cross(:, :), inverse(:, :)
      real(dp), allocatable :: values(:), vectors(:, :), residuals(:, :)
      real(dp) :: denominator, score, scale_value
      integer :: n, m, d, candidate, column, info
      logical :: scale_data

      n = size(predictors, 1)
      m = size(predictors, 2)
      d = size(response, 2)
      allocate(out%coefficient(m, d, size(lambda)))
      allocate(out%lambda(size(lambda)), out%gcv(size(lambda)))
      allocate(out%effective_df(size(lambda)))
      allocate(out%predictor_mean(m), out%predictor_scale(m))
      allocate(out%response_mean(d))
      out%lambda = lambda
      out%gcv = huge(1.0_dp)
      out%effective_df = 0.0_dp
      out%coefficient = 0.0_dp
      out%predictor_mean = 0.0_dp
      out%predictor_scale = 1.0_dp
      out%response_mean = 0.0_dp
      out%selected = 0
      out%info = 0
      scale_data = .false.
      if (present(do_scale)) scale_data = do_scale
      out%scaled = scale_data

      if (n <= 0 .or. m <= 0 .or. d <= 0 .or. size(response, 1) /= n .or. &
         size(lambda) <= 0 .or. any(lambda < 0.0_dp)) then
         out%info = 1
         return
      end if

      allocate(x(n, m), y(n, d))
      x = predictors
      y = response
      if (scale_data) then
         do column = 1, m
            out%predictor_mean(column) = sum(x(:, column))/real(n, dp)
            x(:, column) = x(:, column) - out%predictor_mean(column)
            if (n > 1) then
               scale_value = sqrt(sum(x(:, column)**2)/real(n - 1, dp))
            else
               scale_value = 1.0_dp
            end if
            if (scale_value <= sqrt(tiny(1.0_dp))) scale_value = 1.0_dp
            out%predictor_scale(column) = scale_value
            x(:, column) = x(:, column)/scale_value
         end do
         do column = 1, d
            out%response_mean(column) = sum(y(:, column))/real(n, dp)
            y(:, column) = y(:, column) - out%response_mean(column)
         end do
      end if

      allocate(cross(m, m), values(m), vectors(m, m))
      cross = matmul(transpose(x), x)
      call symmetric_eigen(cross, values, vectors, info)
      if (info /= 0) then
         out%info = 2
         return
      end if

      allocate(inverse(m, m), residuals(n, d))
      do candidate = 1, size(lambda)
         call symmetric_pseudoinverse(cross + real(n, dp)*lambda(candidate)*identity_matrix(m), &
            inverse, info)
         if (info /= 0) then
            out%info = 3
            return
         end if
         out%coefficient(:, :, candidate) = matmul(inverse, matmul(transpose(x), y))
         residuals = y - matmul(x, out%coefficient(:, :, candidate))
         out%effective_df(candidate) = sum(max(values, 0.0_dp)/ &
            (max(values, 0.0_dp) + real(n, dp)*lambda(candidate)))
         denominator = real(n, dp) - out%effective_df(candidate)
         if (denominator > sqrt(tiny(1.0_dp))) then
            score = real(n, dp)*sum(residuals**2)/(denominator**2)
            out%gcv(candidate) = score
         end if
      end do
      out%selected = minloc(out%gcv, dim=1)
   end function varshrink_multivariate_ridge

   pure function varshrink_semibayes(response, predictors, lambda, degrees_of_freedom, &
      conjugate, tolerance, m0, max_iterations) result(out)
      !! Fit VARshrink's fixed-lambda semiparametric Bayesian regression.
      real(dp), intent(in) :: response(:, :) !! N-by-K matrix of dependent variables.
      real(dp), intent(in) :: predictors(:, :) !! N-by-M matrix of regressors.
      real(dp), intent(in) :: lambda !! Shrinkage intensity in the interval [0, 1].
      real(dp), intent(in), optional :: degrees_of_freedom !! Student-t degrees of freedom.
      logical, intent(in), optional :: conjugate !! Use the conjugate covariance prior.
      real(dp), intent(in), optional :: tolerance !! Relative convergence tolerance.
      real(dp), intent(in), optional :: m0 !! Inverse-Wishart prior degrees of freedom.
      integer, intent(in), optional :: max_iterations !! Maximum update iterations.
      type(varshrink_semibayes_fit_t) :: out
      real(dp), allocatable :: coefficient(:, :), covariance(:, :), previous_covariance(:, :)
      real(dp), allocatable :: weights(:), previous_weights(:), inverse_covariance(:, :)
      real(dp), allocatable :: residuals(:, :), quadratic(:)
      real(dp) :: nu, tol, prior_df, coefficient_change, weight_change
      integer :: n, m, d, iteration, limit, info
      logical :: use_conjugate, normal_noise, ncj_phase

      n = size(predictors, 1)
      m = size(predictors, 2)
      d = size(response, 2)
      allocate(out%coefficient(m, d), out%covariance(d, d), out%weights(n))
      allocate(out%fitted(n, d), out%residuals(n, d))
      out%coefficient = 0.0_dp
      out%covariance = 0.0_dp
      out%weights = 1.0_dp
      out%fitted = 0.0_dp
      out%residuals = response
      out%lambda = lambda
      out%iterations = 0
      out%info = 0
      out%converged = .false.
      out%weights_estimated = .false.

      nu = huge(1.0_dp)
      if (present(degrees_of_freedom)) nu = degrees_of_freedom
      use_conjugate = .false.
      if (present(conjugate)) use_conjugate = conjugate
      tol = 1.0e-4_dp
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      prior_df = real(d, dp)
      if (present(m0)) prior_df = m0
      limit = 200
      if (present(max_iterations)) limit = max(1, max_iterations)
      out%degrees_of_freedom = nu
      out%conjugate = use_conjugate

      if (n <= 1 .or. m <= 0 .or. d <= 0 .or. size(response, 1) /= n .or. &
         lambda < 0.0_dp .or. lambda > 1.0_dp .or. prior_df <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (lambda >= 1.0_dp) then
         out%converged = .true.
         return
      end if

      normal_noise = nu >= huge(1.0_dp)/2.0_dp
      if (.not. normal_noise .and. nu <= 0.0_dp) then
         out%info = 2
         return
      end if

      allocate(coefficient(m, d), covariance(d, d), previous_covariance(d, d))
      allocate(weights(n), previous_weights(n), inverse_covariance(d, d))
      allocate(residuals(n, d), quadratic(n))
      weights = 1.0_dp
      call conjugate_update(response, predictors, weights, lambda, prior_df, &
         coefficient, covariance, info)
      if (info /= 0) then
         out%info = 3
         return
      end if

      if (normal_noise .and. use_conjugate) then
         out%coefficient = coefficient
         out%covariance = covariance
         out%weights = weights
         out%fitted = matmul(predictors, coefficient)
         out%residuals = response - out%fitted
         out%converged = .true.
         return
      end if

      ncj_phase = normal_noise
      do iteration = 1, limit
         previous_covariance = covariance
         previous_weights = weights
         call symmetric_pseudoinverse(covariance, inverse_covariance, info)
         if (info /= 0) then
            out%info = 4
            return
         end if

         if (.not. normal_noise .and. .not. ncj_phase) then
            residuals = response - matmul(predictors, coefficient)
            quadratic = row_quadratic(residuals, inverse_covariance)
            weights = (nu + real(d, dp))/(nu + quadratic)
            out%weights_estimated = .true.
         end if

         if (use_conjugate .or. .not. ncj_phase) then
            call conjugate_update(response, predictors, weights, lambda, prior_df, &
               coefficient, covariance, info)
         else
            call nonconjugate_update(response, predictors, weights, lambda, &
               previous_covariance, coefficient, covariance, prior_df, info)
         end if
         if (info /= 0) then
            out%info = 5
            return
         end if

         coefficient_change = sum((covariance - previous_covariance)**2)
         weight_change = sum((weights - previous_weights)**2)
         out%iterations = iteration
         if (normal_noise .or. ncj_phase) then
            if (coefficient_change <= tol*max(sum(previous_covariance**2), &
               tiny(1.0_dp))) then
               out%converged = .true.
               exit
            end if
         else if (weight_change <= tol**2*max(sum(previous_weights**2), tiny(1.0_dp))) then
            if (use_conjugate) then
               out%converged = .true.
               exit
            else
               ncj_phase = .true.
            end if
         end if
      end do

      out%coefficient = coefficient
      out%covariance = covariance
      out%weights = weights
      out%fitted = matmul(predictors, coefficient)
      out%residuals = response - out%fitted
   end function varshrink_semibayes

   pure function varshrink_var_ridge(series, order, lambda, include_intercept) result(out)
      !! Fit a VAR with VARshrink's multivariate ridge and GCV selection.
      real(dp), intent(in) :: series(:, :) !! T-by-K matrix of endogenous variables.
      integer, intent(in) :: order !! VAR lag order.
      real(dp), intent(in) :: lambda(:) !! Candidate nonnegative ridge parameters.
      logical, intent(in), optional :: include_intercept !! Include a constant when true.
      type(varshrink_var_fit_t) :: out
      type(varshrink_ridge_path_t) :: path
      real(dp), allocatable :: x(:, :), y(:, :), coefficient(:, :), residuals(:, :)
      logical :: use_intercept
      integer :: total, n, d, m, lag, row, selected

      total = size(series, 1)
      d = size(series, 2)
      use_intercept = .true.
      if (present(include_intercept)) use_intercept = include_intercept
      out%order = order
      out%include_intercept = use_intercept
      out%info = 0
      if (order < 1 .or. total <= order .or. d < 2 .or. size(lambda) < 1) then
         out%info = 1
         return
      end if

      n = total - order
      m = d*order
      if (use_intercept) m = m + 1
      allocate(x(n, m), y(n, d))
      y = series(order + 1:total, :)
      x = 0.0_dp
      do row = 1, n
         do lag = 1, order
            x(row, (lag - 1)*d + 1:lag*d) = series(order + row - lag, :)
         end do
         if (use_intercept) x(row, m) = 1.0_dp
      end do

      path = varshrink_multivariate_ridge(y, x, lambda)
      if (path%info /= 0) then
         out%info = path%info + 1
         return
      end if
      selected = path%selected
      coefficient = path%coefficient(:, :, selected)
      allocate(out%ar(d, d*order), out%intercept(d))
      allocate(out%fitted(n, d), out%residuals(n, d), out%covariance(d, d))
      out%ar = transpose(coefficient(1:d*order, :))
      out%intercept = 0.0_dp
      if (use_intercept) out%intercept = coefficient(m, :)
      out%fitted = matmul(x, coefficient)
      residuals = y - out%fitted
      out%residuals = residuals
      out%covariance = matmul(transpose(residuals), residuals)/real(n, dp)
      out%lambda = lambda(selected)
      out%gcv = path%gcv(selected)
      out%effective_df = path%effective_df(selected)
   end function varshrink_var_ridge

   pure function varshrink_effective_df(predictors, lambda0, noise_variance, weights) result(value)
      !! Compute the shrinkage-adjusted regression degrees of freedom.
      real(dp), intent(in) :: predictors(:, :) !! N-by-M matrix of regressors.
      real(dp), intent(in) :: lambda0 !! Rescaled nonnegative shrinkage parameter.
      real(dp), intent(in), optional :: noise_variance !! Equation innovation variance.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights.
      real(dp) :: value
      real(dp), allocatable :: cross(:, :), values(:), vectors(:, :)
      real(dp) :: variance
      integer :: n, m, info

      n = size(predictors, 1)
      m = size(predictors, 2)
      value = 0.0_dp
      if (n <= 0 .or. m <= 0 .or. lambda0 < 0.0_dp) return
      allocate(cross(m, m), values(m), vectors(m, m))
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) return
         cross = matmul(transpose(predictors), predictors*spread(weights, 2, m))
      else
         cross = matmul(transpose(predictors), predictors)
      end if
      call symmetric_eigen(cross, values, vectors, info)
      if (info /= 0) return
      variance = 1.0_dp
      if (present(noise_variance)) variance = max(noise_variance, 0.0_dp)
      where (values < 0.0_dp) values = 0.0_dp
      value = sum(values/(values + lambda0*variance), mask=values > sqrt(tiny(1.0_dp)))
   end function varshrink_effective_df

   pure function varshrink_log_likelihood(residuals, covariance, degrees_of_freedom) result(value)
      !! Evaluate the multivariate Gaussian or Student-t residual log likelihood.
      real(dp), intent(in) :: residuals(:, :) !! N-by-K matrix of model residuals.
      real(dp), intent(in) :: covariance(:, :) !! K-by-K innovation scale matrix.
      real(dp), intent(in), optional :: degrees_of_freedom !! Student-t degrees of freedom.
      real(dp) :: value
      real(dp), allocatable :: inverse(:, :), quadratic(:)
      real(dp) :: logdet, nu, pi
      integer :: n, d, info

      n = size(residuals, 1)
      d = size(residuals, 2)
      value = -huge(1.0_dp)
      if (n <= 0 .or. d <= 0 .or. size(covariance, 1) /= d .or. &
         size(covariance, 2) /= d) return
      allocate(inverse(d, d))
      call inverse_logdet(covariance, inverse, logdet, info, &
         100.0_dp*epsilon(1.0_dp))
      if (info /= 0) return
      quadratic = row_quadratic(residuals, inverse)
      pi = acos(-1.0_dp)
      nu = huge(1.0_dp)
      if (present(degrees_of_freedom)) nu = degrees_of_freedom
      if (nu >= huge(1.0_dp)/2.0_dp) then
         value = -0.5_dp*real(n*d, dp)*log(2.0_dp*pi) - &
            0.5_dp*real(n, dp)*logdet - 0.5_dp*sum(quadratic)
      else if (nu > 0.0_dp) then
         value = real(n, dp)*(log_gamma(0.5_dp*(nu + real(d, dp))) - &
            log_gamma(0.5_dp*nu)) - 0.5_dp*real(n*d, dp)*log(nu*pi) - &
            0.5_dp*real(n, dp)*logdet - 0.5_dp*(nu + real(d, dp))* &
            sum(log(1.0_dp + quadratic/nu))
      end if
   end function varshrink_log_likelihood

   pure function varshrink_variance_intensity(series) result(lambda)
      !! Estimate VARshrink's Stein-type intensity for marginal variances.
      real(dp), intent(in) :: series(:, :) !! T-by-K matrix of time-series observations.
      real(dp) :: lambda
      real(dp), allocatable :: centered(:, :), squares(:, :), variances(:)
      real(dp) :: var_r, cov_r, denominator
      integer :: n, d, column, lag

      n = size(series, 1)
      d = size(series, 2)
      lambda = 0.0_dp
      if (n <= 2 .or. d <= 0) return
      allocate(centered(n, d), squares(n, d), variances(d))
      centered = series - spread(sum(series, dim=1)/real(n, dp), 1, n)
      squares = centered**2
      variances = sum(squares, dim=1)/real(n - 1, dp)
      do column = 1, d
         squares(:, column) = squares(:, column) - &
            sum(squares(:, column))/real(n, dp)
      end do
      var_r = sum(squares**2)/real((n - 1)**2, dp)
      cov_r = 0.0_dp
      do lag = 1, n - 1
         cov_r = cov_r + 2.0_dp*real(n - lag, dp)* &
            sum(squares(1:n - lag, :)*squares(lag + 1:n, :))/ &
            real((n - 1)**2*n, dp)
      end do
      denominator = sum((variances - median(variances))**2)
      if (denominator > tiny(1.0_dp)) lambda = min(1.0_dp, max(0.0_dp, &
         (var_r + cov_r)/denominator))
   end function varshrink_variance_intensity

   pure function varshrink_covariance_shrink(data, lambda, lambda_variance) result(out)
      !! Estimate covariance by separate correlation and variance shrinkage.
      real(dp), intent(in) :: data(:, :) !! N-by-P matrix of observations.
      real(dp), intent(in), optional :: lambda !! Correlation shrinkage intensity.
      real(dp), intent(in), optional :: lambda_variance !! Variance shrinkage intensity.
      type(varshrink_covariance_t) :: out
      real(dp), allocatable :: centered(:, :), standardized(:, :)
      real(dp), allocatable :: empirical_variance(:), correlation(:, :)
      real(dp) :: target, numerator, denominator, scale_product
      integer :: n, d, row, column, other

      n = size(data, 1)
      d = size(data, 2)
      allocate(out%covariance(d, d), out%variance(d))
      out%covariance = 0.0_dp
      out%variance = 0.0_dp
      out%info = 0
      if (n < 3 .or. d <= 0) then
         out%info = 1
         return
      end if

      allocate(centered(n, d), standardized(n, d), empirical_variance(d))
      allocate(correlation(d, d))
      centered = data - spread(sum(data, dim=1)/real(n, dp), 1, n)
      empirical_variance = sum(centered**2, dim=1)/real(n - 1, dp)
      if (any(empirical_variance <= tiny(1.0_dp))) then
         out%info = 2
         return
      end if
      target = median(empirical_variance)

      if (present(lambda_variance)) then
         out%lambda_variance = min(1.0_dp, max(0.0_dp, lambda_variance))
      else
         numerator = 0.0_dp
         denominator = 0.0_dp
         do column = 1, d
            numerator = numerator + sum((centered(:, column)**2 - &
               sum(centered(:, column)**2)/real(n, dp))**2)/real(n, dp)
            denominator = denominator + (sum(centered(:, column)**2)/ &
               real(n, dp) - target*real(n - 1, dp)/real(n, dp))**2
         end do
         if (denominator <= tiny(1.0_dp)) then
            out%lambda_variance = 1.0_dp
         else
            out%lambda_variance = min(1.0_dp, max(0.0_dp, &
               numerator/denominator/real(n - 1, dp)))
         end if
         out%lambda_variance_estimated = .true.
      end if
      out%variance = out%lambda_variance*target + &
         (1.0_dp - out%lambda_variance)*empirical_variance

      do column = 1, d
         standardized(:, column) = centered(:, column)/sqrt(empirical_variance(column))
      end do
      correlation = matmul(transpose(standardized), standardized)/real(n - 1, dp)
      if (present(lambda)) then
         out%lambda = min(1.0_dp, max(0.0_dp, lambda))
      else if (d == 1) then
         out%lambda = 1.0_dp
         out%lambda_estimated = .true.
      else
         denominator = 0.0_dp
         numerator = 0.0_dp
         do column = 1, d
            do other = 1, d
               if (column == other) cycle
               denominator = denominator + correlation(column, other)**2
               do row = 1, n
                  numerator = numerator + standardized(row, column)**2* &
                     standardized(row, other)**2/real(n, dp)
               end do
            end do
         end do
         numerator = numerator - denominator
         if (denominator <= tiny(1.0_dp)) then
            out%lambda = 1.0_dp
         else
            out%lambda = min(1.0_dp, max(0.0_dp, numerator/denominator/ &
               real(n - 1, dp)))
         end if
         out%lambda_estimated = .true.
      end if

      correlation = (1.0_dp - out%lambda)*correlation
      do column = 1, d
         correlation(column, column) = 1.0_dp
      end do
      do column = 1, d
         do other = 1, d
            scale_product = sqrt(out%variance(column)*out%variance(other))
            out%covariance(column, other) = correlation(column, other)*scale_product
         end do
      end do
   end function varshrink_covariance_shrink

   pure function varshrink_nonparametric(response, predictors, lambda, &
      lambda_variance) result(out)
      !! Fit VARshrink's covariance-based nonparametric shrinkage regression.
      real(dp), intent(in) :: response(:, :) !! N-by-K matrix of dependent variables.
      real(dp), intent(in) :: predictors(:, :) !! N-by-M matrix of regressors.
      real(dp), intent(in), optional :: lambda !! Correlation shrinkage intensity.
      real(dp), intent(in), optional :: lambda_variance !! Variance shrinkage intensity.
      type(varshrink_nonparametric_fit_t) :: out
      type(varshrink_covariance_t) :: shrinkage
      real(dp), allocatable :: centered_x(:, :), centered_y(:, :), joined(:, :)
      real(dp), allocatable :: inverse(:, :)
      real(dp), allocatable :: x_mean(:), y_mean(:)
      integer :: n, m, d, info

      n = size(predictors, 1)
      m = size(predictors, 2)
      d = size(response, 2)
      allocate(out%coefficient(m, d), out%intercept(d))
      allocate(out%covariance(m + d, m + d))
      allocate(out%fitted(n, d), out%residuals(n, d))
      out%coefficient = 0.0_dp
      out%intercept = 0.0_dp
      out%covariance = 0.0_dp
      out%fitted = 0.0_dp
      out%residuals = response
      out%info = 0
      if (n < 3 .or. m <= 0 .or. d <= 0 .or. size(response, 1) /= n) then
         out%info = 1
         return
      end if

      x_mean = sum(predictors, dim=1)/real(n, dp)
      y_mean = sum(response, dim=1)/real(n, dp)
      centered_x = predictors - spread(x_mean, 1, n)
      centered_y = response - spread(y_mean, 1, n)
      allocate(joined(n, m + d))
      joined(:, 1:m) = centered_x
      joined(:, m + 1:m + d) = centered_y
      if (present(lambda) .and. present(lambda_variance)) then
         shrinkage = varshrink_covariance_shrink(joined, lambda, lambda_variance)
      else if (present(lambda)) then
         shrinkage = varshrink_covariance_shrink(joined, lambda=lambda)
      else if (present(lambda_variance)) then
         shrinkage = varshrink_covariance_shrink(joined, &
            lambda_variance=lambda_variance)
      else
         shrinkage = varshrink_covariance_shrink(joined)
      end if
      if (shrinkage%info /= 0) then
         out%info = shrinkage%info + 1
         return
      end if
      allocate(inverse(m, m))
      call symmetric_pseudoinverse(shrinkage%covariance(1:m, 1:m), inverse, info)
      if (info /= 0) then
         out%info = 4
         return
      end if
      out%coefficient = matmul(inverse, shrinkage%covariance(1:m, m + 1:m + d))
      out%intercept = y_mean - matmul(x_mean, out%coefficient)
      out%fitted = matmul(predictors, out%coefficient) + &
         spread(out%intercept, 1, n)
      out%residuals = response - out%fitted
      out%covariance = shrinkage%covariance
      out%lambda = shrinkage%lambda
      out%lambda_variance = shrinkage%lambda_variance
   end function varshrink_nonparametric

   pure subroutine conjugate_update(response, predictors, weights, lambda, prior_df, &
      coefficient, covariance, info)
      !! Update coefficients and covariance under the conjugate shrinkage form.
      real(dp), intent(in) :: response(:, :) !! N-by-K matrix of dependent variables.
      real(dp), intent(in) :: predictors(:, :) !! N-by-M matrix of regressors.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights.
      real(dp), intent(in) :: lambda !! Shrinkage intensity in the interval [0, 1].
      real(dp), intent(in) :: prior_df !! Inverse-Wishart prior degrees of freedom.
      real(dp), intent(out) :: coefficient(:, :) !! Updated M-by-K coefficient matrix.
      real(dp), intent(out) :: covariance(:, :) !! Updated K-by-K noise covariance.
      integer, intent(out) :: info !! Zero on success; nonzero on numerical failure.
      real(dp), allocatable :: weighted_x(:, :), weighted_y(:, :)
      real(dp), allocatable :: system(:, :), inverse(:, :), residuals(:, :)
      integer :: n, m, d, diagonal

      n = size(predictors, 1)
      m = size(predictors, 2)
      d = size(response, 2)
      weighted_x = predictors*spread(weights, 2, m)
      weighted_y = response*spread(weights, 2, d)
      system = (1.0_dp - lambda)*matmul(transpose(predictors), weighted_x)/ &
         real(n - 1, dp) + lambda*identity_matrix(m)
      allocate(inverse(m, m))
      call symmetric_pseudoinverse(system, inverse, info)
      if (info /= 0) return
      coefficient = matmul(inverse, (1.0_dp - lambda)* &
         matmul(transpose(predictors), weighted_y)/real(n - 1, dp))
      residuals = response - matmul(predictors, coefficient)
      covariance = matmul(transpose(response), residuals*spread(weights, 2, d))
      do diagonal = 1, d
         covariance(diagonal, diagonal) = covariance(diagonal, diagonal) + &
            prior_df + real(d + 1, dp)
      end do
      covariance = 0.5_dp*(covariance + transpose(covariance))/ &
         (prior_df + real(n + d + 1, dp))
   end subroutine conjugate_update

   pure subroutine nonconjugate_update(response, predictors, weights, lambda, &
      previous_covariance, coefficient, covariance, prior_df, info)
      !! Update coefficients using VARshrink's non-conjugate eigenvalue formula.
      real(dp), intent(in) :: response(:, :) !! N-by-K matrix of dependent variables.
      real(dp), intent(in) :: predictors(:, :) !! N-by-M matrix of regressors.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights.
      real(dp), intent(in) :: lambda !! Shrinkage intensity in the interval [0, 1].
      real(dp), intent(in) :: previous_covariance(:, :) !! Previous K-by-K covariance.
      real(dp), intent(out) :: coefficient(:, :) !! Updated M-by-K coefficient matrix.
      real(dp), intent(out) :: covariance(:, :) !! Updated K-by-K noise covariance.
      real(dp), intent(in) :: prior_df !! Inverse-Wishart prior degrees of freedom.
      integer, intent(out) :: info !! Zero on success; nonzero on numerical failure.
      real(dp), allocatable :: x_values(:), x_vectors(:, :)
      real(dp), allocatable :: v_values(:), v_vectors(:, :), inverse_covariance(:, :)
      real(dp), allocatable :: cross(:, :), transformed(:, :), residuals(:, :)
      real(dp) :: theta
      integer :: n, m, d, row, column, diagonal

      n = size(predictors, 1)
      m = size(predictors, 2)
      d = size(response, 2)
      cross = matmul(transpose(predictors), predictors*spread(weights, 2, m))
      allocate(x_values(m), x_vectors(m, m), v_values(d), v_vectors(d, d))
      allocate(inverse_covariance(d, d))
      call symmetric_eigen(cross, x_values, x_vectors, info)
      if (info /= 0) return
      call symmetric_eigen(previous_covariance, v_values, v_vectors, info)
      if (info /= 0) return
      call symmetric_pseudoinverse(previous_covariance, inverse_covariance, info)
      if (info /= 0) return
      transformed = matmul(transpose(x_vectors), matmul(transpose(predictors), &
         matmul(response*spread(weights, 2, d), inverse_covariance)))
      transformed = matmul(transformed, v_vectors)
      theta = lambda*real(n - 1, dp)/max(1.0_dp - lambda, tiny(1.0_dp))
      do column = 1, d
         do row = 1, m
            transformed(row, column) = transformed(row, column)/ &
               (max(x_values(row), 0.0_dp)/max(v_values(column), 1.0e-18_dp) + theta)
         end do
      end do
      coefficient = matmul(x_vectors, matmul(transformed, transpose(v_vectors)))
      residuals = response - matmul(predictors, coefficient)
      covariance = matmul(transpose(response), residuals*spread(weights, 2, d))
      do diagonal = 1, d
         covariance(diagonal, diagonal) = covariance(diagonal, diagonal) + &
            prior_df + real(d + 1, dp)
      end do
      covariance = 0.5_dp*(covariance + transpose(covariance))/ &
         (prior_df + real(n + d + 1, dp))
   end subroutine nonconjugate_update

   pure function row_quadratic(matrix, metric) result(values)
      !! Return row-wise quadratic forms against a symmetric metric.
      real(dp), intent(in) :: matrix(:, :) !! Matrix whose rows define the vectors.
      real(dp), intent(in) :: metric(:, :) !! Symmetric quadratic-form matrix.
      real(dp), allocatable :: values(:)

      values = sum(matrix*matmul(matrix, metric), dim=2)
   end function row_quadratic

end module varshrink_mod
