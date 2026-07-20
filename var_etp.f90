! SPDX-License-Identifier: GPL-2.0-only
! SPDX-FileComment: Distinct algorithms translated from the R VAR.etp package.
module var_etp_mod
   !! Bias correction and prediction methods for stationary vector autoregressions.
   use kind_mod, only: dp
   use linalg_mod, only: invert_matrix, general_eigenvalues, solve_complex_system, &
      kronecker_product
   use stats_mod, only: normal_quantile, sort, quantile, ols_fit
   use random_mod, only: random_uniform
   use mts_mod, only: regularized_beta_mts
   use special_functions_mod, only: regularized_gamma_q
   use vars_mod, only: vars_fit_t, vars_fit, vars_phi, vars_roots, &
      vars_deterministic_none, vars_deterministic_constant, vars_deterministic_both
   implicit none
   private

   type, public :: var_etp_bias_result_t
      !! Bias estimate and stationarity-adjusted VAR fit.
      type(vars_fit_t) :: model
      real(dp), allocatable :: bias(:, :)
      real(dp) :: adjustment = 1.0_dp
      integer :: successful_runs = 0
      integer :: info = 0
   end type var_etp_bias_result_t

   type, public :: var_etp_forecast_result_t
      !! Point forecasts, prediction bounds, and horizon-specific MSE matrices.
      real(dp), allocatable :: forecast(:, :)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
      real(dp), allocatable :: mse(:, :, :)
      real(dp) :: level = 0.95_dp
      integer :: successful_runs = 0
      integer :: info = 0
   end type var_etp_forecast_result_t

   type, public :: var_etp_test_result_t
      !! Asymptotic and bootstrap results for a system restriction test.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: bootstrap_p_value = 1.0_dp
      real(dp) :: df1 = 0.0_dp
      real(dp) :: df2 = 0.0_dp
      integer :: successful_runs = 0
      integer :: info = 0
   end type var_etp_test_result_t

   type, public :: var_etp_predictive_result_t
      !! Ordinary and improved augmented predictive-regression estimates.
      real(dp), allocatable :: ordinary_coefficients(:)
      real(dp), allocatable :: improved_coefficients(:)
      real(dp), allocatable :: ordinary_covariance(:, :)
      real(dp), allocatable :: improved_slope_covariance(:, :)
      real(dp), allocatable :: predictor_ar(:, :)
      real(dp), allocatable :: corrected_predictor_ar(:, :)
      type(var_etp_test_result_t) :: ordinary_test
      type(var_etp_test_result_t) :: improved_test
      integer :: order = 0
      integer :: info = 0
   end type var_etp_predictive_result_t

   type, public :: var_etp_order_result_t
      !! AIC and BIC values and selected predictive-regression orders.
      real(dp), allocatable :: aic(:)
      real(dp), allocatable :: bic(:)
      integer :: aic_order = 0
      integer :: bic_order = 0
      integer :: info = 0
   end type var_etp_order_result_t

   public :: var_etp_pope
   public :: var_etp_bootstrap_bias
   public :: var_etp_forecast
   public :: var_etp_bootstrap_prediction
   public :: var_etp_restrict
   public :: var_etp_wald_test
   public :: var_etp_lr_test
   public :: var_etp_predictive_regression
   public :: var_etp_predictive_order
   public :: var_etp_predictive_forecast

contains

   pure function var_etp_pope(series, order, deterministic_type) result(out)
      !! Apply the Nicholls-Pope asymptotic bias formula and Kilian adjustment.
      real(dp), intent(in) :: series(:, :) !! Time-ordered multivariate observations.
      integer, intent(in) :: order !! Vector autoregressive lag order.
      integer, intent(in), optional :: deterministic_type !! Deterministic-term selector.
      type(var_etp_bias_result_t) :: out
      type(vars_fit_t) :: model
      real(dp), allocatable :: companion(:, :), innovation(:, :), lag_moment(:, :)
      real(dp), allocatable :: lag_moment_inverse(:, :), sum_matrix(:, :), inverse(:, :)
      real(dp), allocatable :: raw_bias(:, :), corrected(:, :)
      complex(dp), allocatable :: roots(:), complex_matrix(:, :), complex_inverse(:, :)
      integer :: selected_type, n, k, p, kp, root_index, status

      selected_type = vars_deterministic_constant
      if (present(deterministic_type)) selected_type = deterministic_type
      model = vars_fit(series, order, selected_type)
      if (model%info /= 0 .or. order < 1) then
         out%info = 1
         return
      end if
      n = size(series, 1)
      k = size(series, 2)
      p = order
      kp = k*p
      companion = companion_matrix(model%ar)
      allocate(innovation(kp, kp))
      innovation = 0.0_dp
      innovation(:k, :k) = model%covariance
      lag_moment = matmul(transpose(model%design(:, :kp)), model%design(:, :kp))/ &
         real(n - p, dp)
      call invert_matrix(lag_moment, lag_moment_inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      call general_eigenvalues(companion, roots, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      allocate(sum_matrix(kp, kp))
      sum_matrix = 0.0_dp
      call invert_matrix(identity(kp) - transpose(companion), inverse, status)
      if (status /= 0) then
         out%info = 4
         return
      end if
      sum_matrix = inverse
      call invert_matrix(identity(kp) - matmul(transpose(companion), &
         transpose(companion)), inverse, status)
      if (status /= 0) then
         out%info = 5
         return
      end if
      sum_matrix = sum_matrix + matmul(transpose(companion), inverse)
      allocate(complex_matrix(kp, kp), complex_inverse(kp, kp))
      do root_index = 1, size(roots)
         complex_matrix = cmplx(identity(kp), 0.0_dp, dp) - &
            roots(root_index)*cmplx(transpose(companion), 0.0_dp, dp)
         call invert_complex_matrix(complex_matrix, complex_inverse, status)
         if (status /= 0) then
            out%info = 6
            return
         end if
         sum_matrix = sum_matrix + real(roots(root_index)*complex_inverse, dp)
      end do
      raw_bias = -matmul(innovation, matmul(sum_matrix, lag_moment_inverse))/real(n, dp)
      allocate(out%bias(size(model%coefficients, 1), size(model%coefficients, 2)))
      out%bias = 0.0_dp
      out%bias(:, :kp) = raw_bias(:k, :)
      call stationarity_adjusted_coefficients(model, out%bias, corrected, out%adjustment)
      out%model = model_with_coefficients(model, corrected)
      if (out%model%info /= 0) out%info = 7
   end function var_etp_pope

   function var_etp_bootstrap_bias(series, order, runs, deterministic_type) result(out)
      !! Estimate VAR coefficient bias by residual bootstrap with stationarity correction.
      real(dp), intent(in) :: series(:, :) !! Time-ordered multivariate observations.
      integer, intent(in) :: order !! Vector autoregressive lag order.
      integer, intent(in) :: runs !! Number of first-stage bootstrap replications.
      integer, intent(in), optional :: deterministic_type !! Deterministic-term selector.
      type(var_etp_bias_result_t) :: out
      type(vars_fit_t) :: model, bootstrap_model
      real(dp), allocatable :: simulated(:, :), mean_coefficients(:, :), corrected(:, :)
      integer :: selected_type, run, successful

      selected_type = vars_deterministic_constant
      if (present(deterministic_type)) selected_type = deterministic_type
      model = vars_fit(series, order, selected_type)
      if (model%info /= 0 .or. runs < 1) then
         out%info = 1
         return
      end if
      allocate(mean_coefficients(size(model%coefficients, 1), size(model%coefficients, 2)))
      mean_coefficients = 0.0_dp
      successful = 0
      do run = 1, runs
         simulated = simulate_var(series, model)
         bootstrap_model = vars_fit(simulated, order, selected_type)
         if (bootstrap_model%info /= 0) cycle
         successful = successful + 1
         mean_coefficients = mean_coefficients + bootstrap_model%coefficients
      end do
      if (successful < 1) then
         out%info = 2
         return
      end if
      mean_coefficients = mean_coefficients/real(successful, dp)
      out%bias = mean_coefficients - model%coefficients
      call stationarity_adjusted_coefficients(model, out%bias, corrected, out%adjustment)
      out%model = model_with_coefficients(model, corrected)
      out%successful_runs = successful
   end function var_etp_bootstrap_bias

   pure function var_etp_forecast(series, order, horizon, level, &
      deterministic_type) result(out)
      !! Form normal prediction intervals using VAR.etp parameter-adjusted forecast MSE.
      real(dp), intent(in) :: series(:, :) !! Time-ordered multivariate observations.
      integer, intent(in) :: order !! Vector autoregressive lag order.
      integer, intent(in) :: horizon !! Number of forecast steps.
      real(dp), intent(in), optional :: level !! Central prediction probability.
      integer, intent(in), optional :: deterministic_type !! Deterministic-term selector.
      type(var_etp_forecast_result_t) :: out
      type(vars_fit_t) :: model
      real(dp) :: selected_level, cutoff
      integer :: selected_type, step

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      selected_type = vars_deterministic_constant
      if (present(deterministic_type)) selected_type = deterministic_type
      if (horizon < 1 .or. selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      model = vars_fit(series, order, selected_type)
      if (model%info /= 0) then
         out%info = 2
         return
      end if
      out%forecast = forecast_from_model(series, model, horizon)
      out%mse = forecast_mse(model, horizon)
      allocate(out%lower(horizon, size(series, 2)), out%upper(horizon, size(series, 2)))
      cutoff = normal_quantile(0.5_dp + 0.5_dp*selected_level)
      do step = 1, horizon
         out%lower(step, :) = out%forecast(step, :) - cutoff*sqrt(max(0.0_dp, &
            diagonal(out%mse(:, :, step))))
         out%upper(step, :) = out%forecast(step, :) + cutoff*sqrt(max(0.0_dp, &
            diagonal(out%mse(:, :, step))))
      end do
      out%level = selected_level
   end function var_etp_forecast

   function var_etp_bootstrap_prediction(series, order, horizon, runs, level, &
      deterministic_type, bias_corrected, bias_runs) result(out)
      !! Compute forward/backward bootstrap or bootstrap-after-bootstrap intervals.
      real(dp), intent(in) :: series(:, :) !! Time-ordered multivariate observations.
      integer, intent(in) :: order !! Vector autoregressive lag order.
      integer, intent(in) :: horizon !! Number of forecast steps.
      integer, intent(in) :: runs !! Number of prediction-bootstrap replications.
      real(dp), intent(in), optional :: level !! Central prediction probability.
      integer, intent(in), optional :: deterministic_type !! Deterministic-term selector.
      logical, intent(in), optional :: bias_corrected !! Apply first-stage bootstrap bias correction.
      integer, intent(in), optional :: bias_runs !! Number of first-stage bias replications.
      type(var_etp_forecast_result_t) :: out
      type(vars_fit_t) :: forward_model, backward_model, bootstrap_model
      type(var_etp_bias_result_t) :: forward_bias, backward_bias
      real(dp), allocatable :: reversed(:, :), synthetic_reversed(:, :), synthetic(:, :)
      real(dp), allocatable :: draws(:, :, :), path(:, :), corrected(:, :)
      real(dp) :: selected_level, tail
      integer :: selected_type, selected_bias_runs, run, successful, step, variable
      logical :: apply_bias

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      selected_type = vars_deterministic_constant
      if (present(deterministic_type)) selected_type = deterministic_type
      apply_bias = .false.
      if (present(bias_corrected)) apply_bias = bias_corrected
      selected_bias_runs = 200
      if (present(bias_runs)) selected_bias_runs = bias_runs
      if (horizon < 1 .or. runs < 1 .or. selected_level <= 0.0_dp .or. &
         selected_level >= 1.0_dp .or. (apply_bias .and. selected_bias_runs < 1)) then
         out%info = 1
         return
      end if
      forward_model = vars_fit(series, order, selected_type)
      reversed = reverse_rows(series)
      backward_model = vars_fit(reversed, order, selected_type)
      if (forward_model%info /= 0 .or. backward_model%info /= 0) then
         out%info = 2
         return
      end if
      if (apply_bias) then
         forward_bias = var_etp_bootstrap_bias(series, order, selected_bias_runs, selected_type)
         backward_bias = var_etp_bootstrap_bias(reversed, order, selected_bias_runs, selected_type)
         if (forward_bias%info /= 0 .or. backward_bias%info /= 0) then
            out%info = 3
            return
         end if
         forward_model = forward_bias%model
         backward_model = backward_bias%model
      end if
      out%forecast = forecast_from_model(series, forward_model, horizon)
      allocate(draws(horizon, size(series, 2), runs))
      successful = 0
      do run = 1, runs
         synthetic_reversed = simulate_var(reversed, backward_model)
         synthetic = reverse_rows(synthetic_reversed)
         bootstrap_model = vars_fit(synthetic, order, selected_type)
         if (bootstrap_model%info /= 0) cycle
         if (apply_bias) then
            call stationarity_adjusted_coefficients(bootstrap_model, forward_bias%bias, &
               corrected, tail)
            bootstrap_model = model_with_coefficients(bootstrap_model, corrected)
         end if
         path = forecast_with_resampled_errors(series, bootstrap_model, &
            forward_model%residuals, horizon)
         successful = successful + 1
         draws(:, :, successful) = path
      end do
      if (successful < 1) then
         out%info = 4
         return
      end if
      allocate(out%lower(horizon, size(series, 2)), out%upper(horizon, size(series, 2)))
      tail = 0.5_dp*(1.0_dp - selected_level)
      do variable = 1, size(series, 2)
         do step = 1, horizon
            path = reshape(draws(step, variable, :successful), [successful, 1])
            call sort(path(:, 1))
            out%lower(step, variable) = quantile(path(:, 1), tail)
            out%upper(step, variable) = quantile(path(:, 1), 1.0_dp - tail)
         end do
      end do
      out%level = selected_level
      out%successful_runs = successful
   end function var_etp_bootstrap_prediction

   pure function var_etp_restrict(series, order, restrictions, values, &
      deterministic_type) result(model)
      !! Estimate a VAR under lag-equation-predictor linear equality restrictions.
      real(dp), intent(in) :: series(:, :) !! Time-ordered multivariate observations.
      integer, intent(in) :: order !! Vector autoregressive lag order.
      integer, intent(in) :: restrictions(:, :) !! Rows containing lag, equation, and predictor.
      real(dp), intent(in), optional :: values(:) !! Restricted coefficient values.
      integer, intent(in), optional :: deterministic_type !! Deterministic-term selector.
      type(vars_fit_t) :: model
      type(vars_fit_t) :: unrestricted
      real(dp), allocatable :: constraint(:, :), target(:), covariance(:, :)
      real(dp), allocatable :: middle(:, :), inverse(:, :), coefficients(:), difference(:)
      integer :: selected_type, status

      selected_type = vars_deterministic_constant
      if (present(deterministic_type)) selected_type = deterministic_type
      unrestricted = vars_fit(series, order, selected_type)
      if (unrestricted%info /= 0 .or. size(restrictions, 2) < 3) then
         model = unrestricted
         model%info = 1
         return
      end if
      call restriction_system(unrestricted, restrictions, values, constraint, target, status)
      if (status /= 0) then
         model = unrestricted
         model%info = 2
         return
      end if
      covariance = coefficient_covariance(unrestricted)
      middle = matmul(constraint, matmul(covariance, transpose(constraint)))
      call invert_matrix(middle, inverse, status)
      if (status /= 0) then
         model = unrestricted
         model%info = 3
         return
      end if
      coefficients = reshape(unrestricted%coefficients, [size(unrestricted%coefficients)])
      difference = matmul(constraint, coefficients) - target
      coefficients = coefficients - matmul(covariance, &
         matmul(transpose(constraint), matmul(inverse, difference)))
      model = model_with_coefficients(unrestricted, reshape(coefficients, &
         shape(unrestricted%coefficients)))
   end function var_etp_restrict

   function var_etp_wald_test(series, order, restrictions, values, &
      deterministic_type, bootstrap_method, runs) result(out)
      !! Test coefficient restrictions with asymptotic or null-imposed bootstrap inference.
      real(dp), intent(in) :: series(:, :) !! Time-ordered multivariate observations.
      integer, intent(in) :: order !! Vector autoregressive lag order.
      integer, intent(in) :: restrictions(:, :) !! Rows containing lag, equation, and predictor.
      real(dp), intent(in), optional :: values(:) !! Restricted coefficient values.
      integer, intent(in), optional :: deterministic_type !! Deterministic-term selector.
      integer, intent(in), optional :: bootstrap_method !! Zero, iid resampling, or Mammen wild bootstrap.
      integer, intent(in), optional :: runs !! Number of bootstrap replications.
      type(var_etp_test_result_t) :: out
      type(vars_fit_t) :: unrestricted, restricted, bootstrap_fit
      real(dp), allocatable :: constraint(:, :), target(:), covariance(:, :)
      real(dp), allocatable :: middle(:, :), inverse(:, :), coefficients(:), difference(:)
      real(dp), allocatable :: simulated(:, :)
      integer :: selected_type, method, selected_runs, status, run, successful, exceedances

      selected_type = vars_deterministic_constant
      if (present(deterministic_type)) selected_type = deterministic_type
      method = 0
      if (present(bootstrap_method)) method = bootstrap_method
      selected_runs = 500
      if (present(runs)) selected_runs = runs
      unrestricted = vars_fit(series, order, selected_type)
      if (unrestricted%info /= 0 .or. method < 0 .or. method > 2 .or. &
         (method > 0 .and. selected_runs < 1)) then
         out%info = 1
         return
      end if
      call restriction_system(unrestricted, restrictions, values, constraint, target, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      covariance = coefficient_covariance(unrestricted)
      middle = matmul(constraint, matmul(covariance, transpose(constraint)))
      call invert_matrix(middle, inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      coefficients = reshape(unrestricted%coefficients, [size(unrestricted%coefficients)])
      difference = matmul(constraint, coefficients) - target
      out%df1 = real(size(constraint, 1), dp)
      out%df2 = real(size(unrestricted%response, 1) - size(unrestricted%design, 2), dp)
      out%statistic = dot_product(difference, matmul(inverse, difference))/out%df1
      out%p_value = f_upper_tail(out%statistic, out%df1, out%df2)
      if (method == 0) return
      restricted = var_etp_restrict(series, order, restrictions, values, selected_type)
      if (restricted%info /= 0) then
         out%info = 4
         return
      end if
      successful = 0
      exceedances = 0
      do run = 1, selected_runs
         simulated = simulate_null(series, restricted, method)
         bootstrap_fit = vars_fit(simulated, order, selected_type)
         if (bootstrap_fit%info /= 0) cycle
         call wald_statistic(bootstrap_fit, constraint, target, difference(1), status)
         if (status /= 0) cycle
         successful = successful + 1
         if (difference(1) > out%statistic) exceedances = exceedances + 1
      end do
      if (successful < 1) then
         out%info = 5
         return
      end if
      out%bootstrap_p_value = real(exceedances, dp)/real(successful, dp)
      out%successful_runs = successful
   end function var_etp_wald_test

   function var_etp_lr_test(series, order, null_restrictions, &
      alternative_restrictions, null_values, alternative_values, &
      deterministic_type, bootstrap_method, runs) result(out)
      !! Test nested VAR restrictions by likelihood ratio with optional null bootstrap.
      real(dp), intent(in) :: series(:, :) !! Time-ordered multivariate observations.
      integer, intent(in) :: order !! Vector autoregressive lag order.
      integer, intent(in) :: null_restrictions(:, :) !! Restrictions imposed under the null.
      integer, intent(in) :: alternative_restrictions(:, :) !! Restrictions under the alternative.
      real(dp), intent(in), optional :: null_values(:) !! Null restricted coefficient values.
      real(dp), intent(in), optional :: alternative_values(:) !! Alternative restricted values.
      integer, intent(in), optional :: deterministic_type !! Deterministic-term selector.
      integer, intent(in), optional :: bootstrap_method !! Zero, iid resampling, or Mammen wild bootstrap.
      integer, intent(in), optional :: runs !! Number of bootstrap replications.
      type(var_etp_test_result_t) :: out
      type(vars_fit_t) :: null_model, alternative_model, bootstrap_null, bootstrap_alternative
      real(dp), allocatable :: simulated(:, :)
      real(dp) :: bootstrap_statistic
      integer :: selected_type, method, selected_runs, restrictions_count
      integer :: run, successful, exceedances

      selected_type = vars_deterministic_constant
      if (present(deterministic_type)) selected_type = deterministic_type
      method = 0
      if (present(bootstrap_method)) method = bootstrap_method
      selected_runs = 500
      if (present(runs)) selected_runs = runs
      null_model = var_etp_restrict(series, order, null_restrictions, &
         null_values, selected_type)
      if (size(alternative_restrictions, 1) == 0) then
         alternative_model = vars_fit(series, order, selected_type)
      else
         alternative_model = var_etp_restrict(series, order, alternative_restrictions, &
            alternative_values, selected_type)
      end if
      restrictions_count = size(null_restrictions, 1) - size(alternative_restrictions, 1)
      if (null_model%info /= 0 .or. alternative_model%info /= 0 .or. &
         restrictions_count < 1 .or. method < 0 .or. method > 2 .or. &
         (method > 0 .and. selected_runs < 1)) then
         out%info = 1
         return
      end if
      out%statistic = real(size(null_model%response, 1), dp)*abs( &
         log_determinant(null_model%covariance) - &
         log_determinant(alternative_model%covariance))
      out%df1 = real(restrictions_count, dp)
      out%df2 = real(size(null_model%response, 1) - size(alternative_model%design, 2), dp)
      out%p_value = regularized_gamma_q(0.5_dp*out%df1, 0.5_dp*out%statistic)
      if (method == 0) return
      successful = 0
      exceedances = 0
      do run = 1, selected_runs
         simulated = simulate_null(series, null_model, method)
         bootstrap_null = var_etp_restrict(simulated, order, null_restrictions, &
            null_values, selected_type)
         if (size(alternative_restrictions, 1) == 0) then
            bootstrap_alternative = vars_fit(simulated, order, selected_type)
         else
            bootstrap_alternative = var_etp_restrict(simulated, order, &
               alternative_restrictions, alternative_values, selected_type)
         end if
         if (bootstrap_null%info /= 0 .or. bootstrap_alternative%info /= 0) cycle
         bootstrap_statistic = real(size(bootstrap_null%response, 1), dp)*abs( &
            log_determinant(bootstrap_null%covariance) - &
            log_determinant(bootstrap_alternative%covariance))
         successful = successful + 1
         if (bootstrap_statistic > out%statistic) exceedances = exceedances + 1
      end do
      if (successful < 1) then
         out%info = 2
         return
      end if
      out%bootstrap_p_value = real(exceedances, dp)/real(successful, dp)
      out%successful_runs = successful
   end function var_etp_lr_test

   pure function var_etp_predictive_regression(predictors, response, order, &
      restriction, target) result(out)
      !! Estimate Kim's improved augmented regression for persistent predictors.
      real(dp), intent(in) :: predictors(:, :) !! Predictor observations in columns.
      real(dp), intent(in) :: response(:) !! Response observations.
      integer, intent(in) :: order !! Common predictor autoregressive order.
      real(dp), intent(in), optional :: restriction(:, :) !! Linear restrictions on lag slopes.
      real(dp), intent(in), optional :: target(:) !! Restriction target values.
      type(var_etp_predictive_result_t) :: out
      real(dp), allocatable :: design(:, :), augmented(:, :), dependent(:)
      real(dp), allocatable :: beta(:), se(:), residuals(:), covariance(:, :)
      real(dp), allocatable :: ar_design(:, :), ar_beta(:), ar_se(:), ar_residual(:)
      real(dp), allocatable :: ar_covariance(:, :, :), corrected_covariance(:, :, :)
      real(dp), allocatable :: correction(:, :), inverse(:, :), innovation_projection(:)
      real(dp), allocatable :: corrected_ar(:)
      real(dp), allocatable :: combined_ar_covariance(:, :), slope_covariance(:, :)
      real(dp), allocatable :: phi(:), phi_outer(:, :), multiplier(:, :)
      real(dp) :: rss, variance
      integer :: n, k, p, rows, variable, status, first, last

      n = size(predictors, 1)
      k = size(predictors, 2)
      p = order
      rows = n - p
      out%order = p
      if (size(response) /= n .or. p < 1 .or. rows <= 1 + k*(p + 1)) then
         out%info = 1
         return
      end if
      design = predictive_design(predictors, p)
      dependent = response(p + 1:)
      call regression_fit(design, dependent, beta, residuals, covariance, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%ordinary_coefficients = beta
      out%ordinary_covariance = covariance
      allocate(out%predictor_ar(p + 1, k), out%corrected_predictor_ar(p + 1, k))
      allocate(ar_covariance(p, p, k), corrected_covariance(p, p, k))
      allocate(augmented(rows, 1 + k*p + k))
      augmented(:, :1 + k*p) = design
      do variable = 1, k
         ar_design = univariate_ar_design(predictors(:, variable), p)
         call regression_fit(ar_design, predictors(p + 1:, variable), ar_beta, &
            ar_residual, covariance, status)
         if (status /= 0) then
            out%info = 3
            return
         end if
         out%predictor_ar(:p, variable) = ar_beta(:p)
         out%predictor_ar(p + 1, variable) = ar_beta(p + 1)
         ar_covariance(:, :, variable) = covariance(:p, :p)
         call shaman_stine_correct(predictors(:, variable), ar_beta, corrected_ar, &
            correction, ar_residual, status)
         if (status /= 0) then
            out%info = 4
            return
         end if
         out%corrected_predictor_ar(:, variable) = corrected_ar
         call invert_matrix(correction, inverse, status)
         corrected_covariance(:, :, variable) = matmul(inverse, &
            matmul(ar_covariance(:, :, variable), transpose(inverse)))
         augmented(:, 1 + k*p + variable) = ar_residual
      end do
      call regression_fit(augmented, dependent, beta, residuals, covariance, status)
      if (status /= 0) then
         out%info = 5
         return
      end if
      out%improved_coefficients = beta
      allocate(combined_ar_covariance(k*p, k*p))
      combined_ar_covariance = 0.0_dp
      do variable = 1, k
         first = (variable - 1)*p + 1
         last = variable*p
         combined_ar_covariance(first:last, first:last) = &
            corrected_covariance(:, :, variable)
      end do
      phi = beta(2 + k*p:1 + k*p + k)
      phi_outer = spread(phi, 2, k)*spread(phi, 1, k)
      multiplier = kronecker_product(phi_outer, ones(p, p))
      slope_covariance = multiplier*combined_ar_covariance + covariance(2:1 + k*p, 2:1 + k*p)
      out%improved_slope_covariance = slope_covariance
      if (present(restriction)) then
         if (size(restriction, 2) /= k*p .or. &
            (present(target) .and. size(target) /= size(restriction, 1))) then
            out%info = 6
            return
         end if
         call predictive_wald(out%ordinary_coefficients(2:1 + k*p), &
            out%ordinary_covariance(2:1 + k*p, 2:1 + k*p), restriction, target, &
            real(rows - size(design, 2), dp), out%ordinary_test)
         call predictive_wald(out%improved_coefficients(2:1 + k*p), &
            out%improved_slope_covariance, restriction, target, &
            real(rows - size(augmented, 2), dp), out%improved_test)
      end if
   end function var_etp_predictive_regression

   pure function var_etp_predictive_order(predictors, response, maximum_order) result(out)
      !! Select the predictive-regression order by joint residual AIC and BIC.
      real(dp), intent(in) :: predictors(:, :) !! Predictor observations in columns.
      real(dp), intent(in) :: response(:) !! Response observations.
      integer, intent(in) :: maximum_order !! Largest candidate order.
      type(var_etp_order_result_t) :: out
      real(dp), allocatable :: design(:, :), dependent(:), beta(:), residual(:), covariance(:, :)
      real(dp), allocatable :: residual_matrix(:, :), ar_design(:, :), ar_beta(:)
      real(dp), allocatable :: ar_residual(:), residual_covariance(:, :)
      real(dp) :: logdet
      integer :: order, variable, rows, k, status

      k = size(predictors, 2)
      if (size(response) /= size(predictors, 1) .or. maximum_order < 1) then
         out%info = 1
         return
      end if
      allocate(out%aic(maximum_order), out%bic(maximum_order))
      do order = 1, maximum_order
         rows = size(response) - order
         design = predictive_design(predictors, order)
         dependent = response(order + 1:)
         call regression_fit(design, dependent, beta, residual, covariance, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         allocate(residual_matrix(rows, k + 1))
         residual_matrix(:, 1) = residual
         do variable = 1, k
            ar_design = univariate_ar_design(predictors(:, variable), order)
            call regression_fit(ar_design, predictors(order + 1:, variable), ar_beta, &
               ar_residual, covariance, status)
            residual_matrix(:, variable + 1) = ar_residual
         end do
         residual_covariance = matmul(transpose(residual_matrix), residual_matrix)/real(rows, dp)
         logdet = log_determinant(residual_covariance)
         out%aic(order) = logdet + 2.0_dp*real((k + 1)*order, dp)/real(rows, dp)
         out%bic(order) = logdet + log(real(rows, dp))*real((k + 1)*order, dp)/real(rows, dp)
         deallocate(residual_matrix)
      end do
      out%aic_order = minloc(out%aic, dim=1)
      out%bic_order = minloc(out%bic, dim=1)
   end function var_etp_predictive_order

   pure function var_etp_predictive_forecast(predictors, response, model, horizon) &
      result(forecast)
      !! Forecast the response and predictors from an improved augmented regression.
      real(dp), intent(in) :: predictors(:, :) !! Predictor observations in columns.
      real(dp), intent(in) :: response(:) !! Response observations.
      type(var_etp_predictive_result_t), intent(in) :: model !! Improved regression fit.
      integer, intent(in) :: horizon !! Number of dynamic forecast steps.
      real(dp), allocatable :: forecast(:, :)
      real(dp), allocatable :: x_history(:, :), y_history(:), row(:)
      integer :: n, k, p, step, variable, lag, column

      n = size(predictors, 1)
      k = size(predictors, 2)
      p = model%order
      allocate(forecast(horizon, k + 1), x_history(n + horizon, k), &
         y_history(n + horizon), row(1 + k*p))
      x_history(:n, :) = predictors
      y_history(:n) = response
      do step = 1, horizon
         row = 0.0_dp
         row(1) = 1.0_dp
         column = 2
         do variable = 1, k
            do lag = 1, p
               row(column) = x_history(n + step - lag, variable)
               column = column + 1
            end do
         end do
         y_history(n + step) = dot_product(row, model%improved_coefficients(:1 + k*p))
         do variable = 1, k
            x_history(n + step, variable) = model%corrected_predictor_ar(p + 1, variable)
            do lag = 1, p
               x_history(n + step, variable) = x_history(n + step, variable) + &
                  model%corrected_predictor_ar(lag, variable)* &
                  x_history(n + step - lag, variable)
            end do
         end do
         forecast(step, 1) = y_history(n + step)
         forecast(step, 2:) = x_history(n + step, :)
      end do
   end function var_etp_predictive_forecast

   pure function predictive_design(predictors, order) result(design)
      !! Construct intercept and predictor-major lag columns for predictive regression.
      real(dp), intent(in) :: predictors(:, :) !! Predictor observations in columns.
      integer, intent(in) :: order !! Common lag order.
      real(dp), allocatable :: design(:, :)
      integer :: n, k, variable, lag, column

      n = size(predictors, 1)
      k = size(predictors, 2)
      allocate(design(n - order, 1 + k*order))
      design(:, 1) = 1.0_dp
      column = 2
      do variable = 1, k
         do lag = 1, order
            design(:, column) = predictors(order + 1 - lag:n - lag, variable)
            column = column + 1
         end do
      end do
   end function predictive_design

   pure function univariate_ar_design(series, order) result(design)
      !! Construct lag columns followed by an intercept for an AR regression.
      real(dp), intent(in) :: series(:) !! Univariate predictor series.
      integer, intent(in) :: order !! Autoregressive order.
      real(dp), allocatable :: design(:, :)
      integer :: lag, n

      n = size(series)
      allocate(design(n - order, order + 1))
      do lag = 1, order
         design(:, lag) = series(order + 1 - lag:n - lag)
      end do
      design(:, order + 1) = 1.0_dp
   end function univariate_ar_design

   pure subroutine regression_fit(design, response, coefficients, residuals, covariance, info)
      !! Fit a linear regression and return its coefficient covariance matrix.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      real(dp), intent(in) :: response(:) !! Regression response vector.
      real(dp), allocatable, intent(out) :: coefficients(:) !! OLS coefficients.
      real(dp), allocatable, intent(out) :: residuals(:) !! OLS residuals.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Coefficient covariance matrix.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: standard_errors(:), inverse(:, :)
      real(dp) :: rss, variance

      call ols_fit(design, response, coefficients, standard_errors, residuals, rss, info)
      if (info /= 0) return
      call invert_matrix(matmul(transpose(design), design), inverse, info)
      if (info /= 0) return
      variance = rss/real(max(1, size(design, 1) - size(design, 2)), dp)
      covariance = variance*inverse
   end subroutine regression_fit

   pure subroutine shaman_stine_correct(series, coefficients, corrected, &
      correction, residuals, info)
      !! Apply the Shaman-Stine finite-sample AR coefficient correction.
      real(dp), intent(in) :: series(:) !! Univariate predictor series.
      real(dp), intent(in) :: coefficients(:) !! AR coefficients followed by intercept.
      real(dp), allocatable, intent(out) :: corrected(:) !! Bias-corrected AR coefficients and intercept.
      real(dp), allocatable, intent(out) :: correction(:, :) !! Linear correction matrix.
      real(dp), allocatable, intent(out) :: residuals(:) !! Corrected AR residuals.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: bias_matrix(:, :), right_hand(:), inverse(:, :), design(:, :)
      integer :: p, n, i, j, m, position

      p = size(coefficients) - 1
      n = size(series)
      allocate(bias_matrix(p + 1, p + 1))
      bias_matrix = 0.0_dp
      if (p == 1) then
         bias_matrix = reshape([0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp], [2, 2])
      else
         do i = 1, p + 1
            bias_matrix(i, i) = real(i - 1, dp)
         end do
         if (mod(p, 2) == 0) then
            m = p/2
            do i = 0, m - 1
               do j = 0, m - i - 1
                  position = i + 3 + 2*j
                  bias_matrix(position, i + 1) = bias_matrix(position, i + 1) - 1.0_dp
                  bias_matrix(position, p + 1 - i) = &
                     bias_matrix(position, p + 1 - i) + 1.0_dp
               end do
            end do
         else
            m = (p - 1)/2
            do i = 0, m
               do j = 0, m - i
                  position = i + 2 + 2*j
                  if (i > 0) bias_matrix(position, i) = bias_matrix(position, i) - 1.0_dp
                  bias_matrix(position, p + 1 - i) = &
                     bias_matrix(position, p + 1 - i) + 1.0_dp
               end do
            end do
         end if
         do i = 1, p + 1
            do j = 1, p + 1
               if (j < i .and. i <= p - j + 2) bias_matrix(i, j) = bias_matrix(i, j) - 1.0_dp
               if (p - j + 2 < i .and. i <= j) bias_matrix(i, j) = bias_matrix(i, j) + 1.0_dp
            end do
         end do
         bias_matrix(:, 1) = -bias_matrix(:, 1)
      end if
      correction = identity(p) - bias_matrix(2:, 2:)/real(n, dp)
      right_hand = coefficients(:p) + bias_matrix(2:, 1)/real(n, dp)
      call invert_matrix(correction, inverse, info)
      if (info /= 0) return
      allocate(corrected(p + 1))
      corrected(:p) = matmul(inverse, right_hand)
      corrected(p + 1) = sum(series)/real(n, dp)*(1.0_dp - sum(corrected(:p)))
      design = univariate_ar_design(series, p)
      residuals = series(p + 1:) - matmul(design, corrected)
   end subroutine shaman_stine_correct

   pure subroutine predictive_wald(slopes, covariance, restriction, target, df2, out)
      !! Evaluate a predictive-regression slope restriction F test.
      real(dp), intent(in) :: slopes(:) !! Predictive-regression lag slopes.
      real(dp), intent(in) :: covariance(:, :) !! Slope covariance matrix.
      real(dp), intent(in) :: restriction(:, :) !! Linear restriction matrix.
      real(dp), intent(in), optional :: target(:) !! Restriction target values.
      real(dp), intent(in) :: df2 !! Denominator degrees of freedom.
      type(var_etp_test_result_t), intent(out) :: out !! Completed test result.
      real(dp), allocatable :: difference(:), middle(:, :), inverse(:, :)
      integer :: status

      difference = matmul(restriction, slopes)
      if (present(target)) difference = difference - target
      middle = matmul(restriction, matmul(covariance, transpose(restriction)))
      call invert_matrix(middle, inverse, status)
      if (status /= 0) then
         out%info = 1
         return
      end if
      out%df1 = real(size(restriction, 1), dp)
      out%df2 = df2
      out%statistic = dot_product(difference, matmul(inverse, difference))/out%df1
      out%p_value = f_upper_tail(out%statistic, out%df1, out%df2)
   end subroutine predictive_wald

   pure subroutine restriction_system(model, restrictions, values, constraint, target, info)
      !! Convert package lag-equation-predictor rows to a coefficient constraint system.
      type(vars_fit_t), intent(in) :: model !! Unrestricted fitted VAR.
      integer, intent(in) :: restrictions(:, :) !! Lag, equation, and predictor rows.
      real(dp), intent(in), optional :: values(:) !! Restricted coefficient values.
      real(dp), allocatable, intent(out) :: constraint(:, :) !! Constraint coefficient matrix.
      real(dp), allocatable, intent(out) :: target(:) !! Constraint target vector.
      integer, intent(out) :: info !! Zero on success.
      integer :: row, lag, equation, predictor, position, k

      k = size(model%response, 2)
      allocate(constraint(size(restrictions, 1), size(model%coefficients)))
      allocate(target(size(restrictions, 1)))
      constraint = 0.0_dp
      target = 0.0_dp
      if (size(restrictions, 2) < 3 .or. &
         (present(values) .and. size(values) /= size(restrictions, 1))) then
         info = 1
         return
      end if
      do row = 1, size(restrictions, 1)
         lag = restrictions(row, 1)
         equation = restrictions(row, 2)
         predictor = restrictions(row, 3)
         if (lag < 1 .or. lag > model%order .or. equation < 1 .or. equation > k .or. &
            predictor < 1 .or. predictor > k) then
            info = 2
            return
         end if
         position = equation + ((lag - 1)*k + predictor - 1)*k
         constraint(row, position) = 1.0_dp
      end do
      if (present(values)) target = values
      info = 0
   end subroutine restriction_system

   pure function coefficient_covariance(model) result(covariance)
      !! Form the Kronecker covariance of vectorized unrestricted VAR coefficients.
      type(vars_fit_t), intent(in) :: model !! Unrestricted fitted VAR.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: inverse(:, :)
      integer :: status

      call invert_matrix(matmul(transpose(model%design), model%design), inverse, status)
      if (status /= 0) then
         allocate(covariance(size(model%coefficients), size(model%coefficients)))
         covariance = huge(1.0_dp)
         return
      end if
      covariance = kronecker_product(inverse, model%covariance)
   end function coefficient_covariance

   pure subroutine wald_statistic(model, constraint, target, statistic, info)
      !! Evaluate a system Wald F statistic for a fitted bootstrap model.
      type(vars_fit_t), intent(in) :: model !! Fitted unrestricted VAR.
      real(dp), intent(in) :: constraint(:, :) !! Coefficient constraint matrix.
      real(dp), intent(in) :: target(:) !! Constraint target vector.
      real(dp), intent(out) :: statistic !! Wald F statistic.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: covariance(:, :), middle(:, :), inverse(:, :)
      real(dp), allocatable :: coefficients(:), difference(:)

      covariance = coefficient_covariance(model)
      middle = matmul(constraint, matmul(covariance, transpose(constraint)))
      call invert_matrix(middle, inverse, info)
      if (info /= 0) return
      coefficients = reshape(model%coefficients, [size(model%coefficients)])
      difference = matmul(constraint, coefficients) - target
      statistic = dot_product(difference, matmul(inverse, difference))/ &
         real(size(constraint, 1), dp)
   end subroutine wald_statistic

   function simulate_null(series, model, method) result(simulated)
      !! Simulate under a restricted null using iid or Mammen wild innovations.
      real(dp), intent(in) :: series(:, :) !! Source observations supplying initial values.
      type(vars_fit_t), intent(in) :: model !! Restricted data-generating VAR.
      integer, intent(in) :: method !! One for iid or two for Mammen wild bootstrap.
      real(dp), allocatable :: simulated(:, :)
      real(dp), allocatable :: centered(:, :), prediction(:), innovation(:)
      real(dp) :: probability, weight
      integer :: n, k, p, rows, t, lag, index

      n = size(series, 1)
      k = size(series, 2)
      p = model%order
      rows = size(model%residuals, 1)
      centered = model%residuals - spread(sum(model%residuals, dim=1)/real(rows, dp), 1, rows)
      allocate(simulated(n, k), prediction(k), innovation(k))
      simulated(:p, :) = series(:p, :)
      probability = (sqrt(5.0_dp) + 1.0_dp)/(2.0_dp*sqrt(5.0_dp))
      do t = p + 1, n
         prediction = deterministic_prediction(model, t)
         do lag = 1, p
            prediction = prediction + matmul(model%ar(:, :, lag), simulated(t - lag, :))
         end do
         if (method == 1) then
            index = min(rows, 1 + int(random_uniform()*real(rows, dp)))
            innovation = centered(index, :)
         else
            weight = -(sqrt(5.0_dp) - 1.0_dp)/2.0_dp
            if (random_uniform() > probability) weight = (sqrt(5.0_dp) + 1.0_dp)/2.0_dp
            innovation = weight*centered(t - p, :)
         end if
         simulated(t, :) = prediction + innovation
      end do
   end function simulate_null

   pure function forecast_mse(model, horizon) result(mse)
      !! Compute Lutkepohl forecast MSE including VAR coefficient uncertainty.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      integer, intent(in) :: horizon !! Number of forecast steps.
      real(dp), allocatable :: mse(:, :, :)
      real(dp), allocatable :: phi(:, :, :), transition(:, :), powers(:, :, :)
      real(dp), allocatable :: moment(:, :), moment_inverse(:, :), innovation(:, :)
      real(dp), allocatable :: estimation(:, :), left(:, :), right(:, :)
      real(dp) :: multiplier
      integer :: k, p, q, rows, h, i, j, status

      k = size(model%ar, 1)
      p = model%order
      q = size(model%design, 2)
      rows = size(model%design, 1)
      phi = vars_phi(model, horizon - 1)
      moment = matmul(transpose(model%design), model%design)/real(rows, dp)
      call invert_matrix(moment, moment_inverse, status)
      allocate(mse(k, k, horizon))
      if (status /= 0) then
         mse = huge(1.0_dp)
         return
      end if
      transition = design_transition(model)
      allocate(powers(q, q, horizon))
      powers(:, :, 1) = identity(q)
      do h = 2, horizon
         powers(:, :, h) = matmul(powers(:, :, h - 1), transition)
      end do
      do h = 1, horizon
         allocate(innovation(k, k), estimation(k, k))
         innovation = 0.0_dp
         estimation = 0.0_dp
         do i = 0, h - 1
            innovation = innovation + matmul(phi(:, :, i + 1), &
               matmul(model%covariance, transpose(phi(:, :, i + 1))))
            do j = 0, h - 1
               left = powers(:, :, h - i)
               right = powers(:, :, h - j)
               multiplier = matrix_trace(matmul(transpose(left), &
                  matmul(moment_inverse, matmul(right, moment))))
               estimation = estimation + multiplier*matmul(phi(:, :, i + 1), &
                  matmul(model%covariance, transpose(phi(:, :, j + 1))))
            end do
         end do
         mse(:, :, h) = innovation + estimation/real(rows, dp)
         deallocate(innovation, estimation)
      end do
   end function forecast_mse

   function simulate_var(series, model) result(simulated)
      !! Simulate a fitted VAR over its sample by centered residual resampling.
      real(dp), intent(in) :: series(:, :) !! Source observations supplying initial values.
      type(vars_fit_t), intent(in) :: model !! Data-generating VAR model.
      real(dp), allocatable :: simulated(:, :)
      real(dp), allocatable :: centered(:, :), prediction(:)
      integer :: n, k, p, rows, t, lag, index

      n = size(series, 1)
      k = size(series, 2)
      p = model%order
      rows = size(model%residuals, 1)
      centered = model%residuals - spread(sum(model%residuals, dim=1)/real(rows, dp), 1, rows)
      allocate(simulated(n, k), prediction(k))
      simulated(:p, :) = series(:p, :)
      do t = p + 1, n
         prediction = deterministic_prediction(model, t)
         do lag = 1, p
            prediction = prediction + matmul(model%ar(:, :, lag), simulated(t - lag, :))
         end do
         index = min(rows, 1 + int(random_uniform()*real(rows, dp)))
         simulated(t, :) = prediction + centered(index, :)
      end do
   end function simulate_var

   function forecast_with_resampled_errors(series, model, residuals, horizon) result(forecast)
      !! Generate dynamic forecasts with independently resampled future innovations.
      real(dp), intent(in) :: series(:, :) !! Observed history.
      type(vars_fit_t), intent(in) :: model !! VAR coefficient model.
      real(dp), intent(in) :: residuals(:, :) !! Innovation pool.
      integer, intent(in) :: horizon !! Number of forecast steps.
      real(dp), allocatable :: forecast(:, :)
      real(dp), allocatable :: history(:, :), centered(:, :), prediction(:)
      integer :: n, k, p, rows, step, lag, index

      n = size(series, 1)
      k = size(series, 2)
      p = model%order
      rows = size(residuals, 1)
      centered = residuals - spread(sum(residuals, dim=1)/real(rows, dp), 1, rows)
      allocate(history(n + horizon, k), forecast(horizon, k), prediction(k))
      history(:n, :) = series
      do step = 1, horizon
         prediction = deterministic_prediction(model, n + step)
         do lag = 1, p
            prediction = prediction + matmul(model%ar(:, :, lag), history(n + step - lag, :))
         end do
         index = min(rows, 1 + int(random_uniform()*real(rows, dp)))
         history(n + step, :) = prediction + centered(index, :)
         forecast(step, :) = history(n + step, :)
      end do
   end function forecast_with_resampled_errors

   pure function forecast_from_model(series, model, horizon) result(forecast)
      !! Generate deterministic dynamic forecasts from a fitted VAR.
      real(dp), intent(in) :: series(:, :) !! Observed history.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      integer, intent(in) :: horizon !! Number of forecast steps.
      real(dp), allocatable :: forecast(:, :)
      real(dp), allocatable :: history(:, :), prediction(:)
      integer :: n, k, p, step, lag

      n = size(series, 1)
      k = size(series, 2)
      p = model%order
      allocate(history(n + horizon, k), forecast(horizon, k), prediction(k))
      history(:n, :) = series
      do step = 1, horizon
         prediction = deterministic_prediction(model, n + step)
         do lag = 1, p
            prediction = prediction + matmul(model%ar(:, :, lag), history(n + step - lag, :))
         end do
         history(n + step, :) = prediction
         forecast(step, :) = prediction
      end do
   end function forecast_from_model

   pure function deterministic_prediction(model, time_index) result(value)
      !! Evaluate VAR.etp constant and trend terms at a time index.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      integer, intent(in) :: time_index !! One-based observation or forecast index.
      real(dp), allocatable :: value(:)

      allocate(value(size(model%ar, 1)))
      value = 0.0_dp
      if (size(model%deterministic, 2) >= 1) value = model%deterministic(:, 1)
      if (size(model%deterministic, 2) >= 2) value = value + &
         model%deterministic(:, 2)*real(time_index, dp)
   end function deterministic_prediction

   pure subroutine stationarity_adjusted_coefficients(model, bias, coefficients, adjustment)
      !! Shrink a bias correction until the corrected VAR is stationary.
      type(vars_fit_t), intent(in) :: model !! Uncorrected fitted VAR.
      real(dp), intent(in) :: bias(:, :) !! Estimated coefficient bias.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Adjusted coefficients.
      real(dp), intent(out) :: adjustment !! Retained fraction of the bias correction.
      type(vars_fit_t) :: candidate
      complex(dp), allocatable :: roots(:)

      adjustment = 1.0_dp
      do
         coefficients = model%coefficients - adjustment*bias
         candidate = model_with_coefficients(model, coefficients)
         roots = vars_roots(candidate)
         if (maxval(abs(roots)) < 1.0_dp) exit
         adjustment = adjustment - 0.01_dp
         if (adjustment <= 0.0_dp) then
            adjustment = 0.0_dp
            coefficients = model%coefficients
            exit
         end if
      end do
   end subroutine stationarity_adjusted_coefficients

   pure function model_with_coefficients(model, coefficients) result(out)
      !! Replace VAR coefficients and recompute fitted moments.
      type(vars_fit_t), intent(in) :: model !! Source fitted VAR.
      real(dp), intent(in) :: coefficients(:, :) !! Replacement coefficient matrix.
      type(vars_fit_t) :: out
      integer :: k, rows

      out = model
      if (any(shape(coefficients) /= shape(model%coefficients))) then
         out%info = 1
         return
      end if
      k = size(model%response, 2)
      rows = size(model%response, 1)
      out%coefficients = coefficients
      out%ar = reshape(coefficients(:, :k*model%order), [k, k, model%order])
      if (size(coefficients, 2) > k*model%order) &
         out%deterministic = coefficients(:, k*model%order + 1:)
      out%fitted = matmul(model%design, transpose(coefficients))
      out%residuals = model%response - out%fitted
      out%covariance = matmul(transpose(out%residuals), out%residuals)/ &
         real(max(1, rows - size(model%design, 2)), dp)
   end function model_with_coefficients

   pure function design_transition(model) result(transition)
      !! Construct the augmented companion transition for the VAR design state.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      real(dp), allocatable :: transition(:, :)
      integer :: k, p, q, deterministic_columns

      k = size(model%ar, 1)
      p = model%order
      q = size(model%design, 2)
      deterministic_columns = q - k*p
      allocate(transition(q, q))
      transition = 0.0_dp
      transition(:k, :) = model%coefficients
      if (p > 1) transition(k + 1:k*p, :k*(p - 1)) = identity(k*(p - 1))
      if (deterministic_columns >= 1) transition(k*p + 1, k*p + 1) = 1.0_dp
      if (deterministic_columns >= 2) then
         transition(k*p + 2, k*p + 1) = 1.0_dp
         transition(k*p + 2, k*p + 2) = 1.0_dp
      end if
   end function design_transition

   pure function companion_matrix(ar) result(companion)
      !! Construct the standard first-order companion matrix of a VAR.
      real(dp), intent(in) :: ar(:, :, :) !! Autoregressive coefficient matrices.
      real(dp), allocatable :: companion(:, :)
      integer :: k, p, lag

      k = size(ar, 1)
      p = size(ar, 3)
      allocate(companion(k*p, k*p))
      companion = 0.0_dp
      do lag = 1, p
         companion(:k, (lag - 1)*k + 1:lag*k) = ar(:, :, lag)
      end do
      if (p > 1) companion(k + 1:, :k*(p - 1)) = identity(k*(p - 1))
   end function companion_matrix

   pure subroutine invert_complex_matrix(matrix, inverse, info)
      !! Invert a complex matrix by repeated shared linear solves.
      complex(dp), intent(in) :: matrix(:, :) !! Square complex matrix.
      complex(dp), allocatable, intent(out) :: inverse(:, :) !! Matrix inverse.
      integer, intent(out) :: info !! Zero on success.
      complex(dp), allocatable :: right_hand(:), solution(:)
      integer :: column, n

      n = size(matrix, 1)
      allocate(inverse(n, n), right_hand(n), solution(n))
      inverse = cmplx(0.0_dp, 0.0_dp, dp)
      do column = 1, n
         right_hand = cmplx(0.0_dp, 0.0_dp, dp)
         right_hand(column) = cmplx(1.0_dp, 0.0_dp, dp)
         call solve_complex_system(matrix, right_hand, solution, info)
         if (info /= 0) return
         inverse(:, column) = solution
      end do
      info = 0
   end subroutine invert_complex_matrix

   pure function reverse_rows(matrix) result(reversed)
      !! Return a matrix with its observation order reversed.
      real(dp), intent(in) :: matrix(:, :) !! Input observation matrix.
      real(dp), allocatable :: reversed(:, :)
      integer :: row, n

      n = size(matrix, 1)
      allocate(reversed(n, size(matrix, 2)))
      do row = 1, n
         reversed(row, :) = matrix(n - row + 1, :)
      end do
   end function reverse_rows

   pure function diagonal(matrix) result(values)
      !! Extract the diagonal of a rectangular matrix.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), allocatable :: values(:)
      integer :: index, n

      n = min(size(matrix, 1), size(matrix, 2))
      allocate(values(n))
      do index = 1, n
         values(index) = matrix(index, index)
      end do
   end function diagonal

   pure real(dp) function matrix_trace(matrix) result(value)
      !! Return the trace of a square matrix.
      real(dp), intent(in) :: matrix(:, :) !! Input square matrix.
      integer :: index

      value = 0.0_dp
      do index = 1, min(size(matrix, 1), size(matrix, 2))
         value = value + matrix(index, index)
      end do
   end function matrix_trace

   pure real(dp) function log_determinant(matrix) result(value)
      !! Return the log absolute determinant by pivoted elimination.
      real(dp), intent(in) :: matrix(:, :) !! Input square matrix.
      real(dp), allocatable :: work(:, :), row(:)
      real(dp) :: pivot_value, factor
      integer :: n, column, pivot, index

      n = size(matrix, 1)
      work = matrix
      value = 0.0_dp
      do column = 1, n
         pivot = column - 1 + maxloc(abs(work(column:, column)), dim=1)
         pivot_value = work(pivot, column)
         if (abs(pivot_value) <= tiny(1.0_dp)) then
            value = log(tiny(1.0_dp))
            return
         end if
         if (pivot /= column) then
            row = work(column, :)
            work(column, :) = work(pivot, :)
            work(pivot, :) = row
         end if
         value = value + log(abs(work(column, column)))
         do index = column + 1, n
            factor = work(index, column)/work(column, column)
            work(index, column:) = work(index, column:) - factor*work(column, column:)
         end do
      end do
   end function log_determinant

   pure real(dp) function f_upper_tail(statistic, df1, df2) result(probability)
      !! Return the upper-tail probability of an F statistic.
      real(dp), intent(in) :: statistic !! Nonnegative F statistic.
      real(dp), intent(in) :: df1 !! Numerator degrees of freedom.
      real(dp), intent(in) :: df2 !! Denominator degrees of freedom.
      real(dp) :: argument

      if (statistic < 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
         probability = 1.0_dp
         return
      end if
      argument = df2/(df2 + df1*statistic)
      probability = regularized_beta_mts(argument, 0.5_dp*df2, 0.5_dp*df1)
   end function f_upper_tail

   pure function identity(n) result(matrix)
      !! Construct an identity matrix.
      integer, intent(in) :: n !! Matrix dimension.
      real(dp) :: matrix(n, n)
      integer :: index

      matrix = 0.0_dp
      do index = 1, n
         matrix(index, index) = 1.0_dp
      end do
   end function identity

   pure function ones(rows, columns) result(matrix)
      !! Construct a matrix of ones.
      integer, intent(in) :: rows !! Number of rows.
      integer, intent(in) :: columns !! Number of columns.
      real(dp) :: matrix(rows, columns)

      matrix = 1.0_dp
   end function ones

end module var_etp_mod
