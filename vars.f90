! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Distinct algorithms translated from the R vars package.
module vars_mod
   !! Reduced-form VAR specification, diagnostics, restrictions, and identification.
   use kind_mod, only: dp
   use linalg_mod, only: invert_matrix, cholesky_lower, symmetric_eigen, &
      general_eigenvalues, kronecker_product
   use stats_mod, only: ols_fit, sort, quantile
   use special_functions_mod, only: regularized_gamma_q
   use mts_mod, only: regularized_beta_mts
   use optimization_mod, only: optimization_result_t, nelder_mead_minimize, &
      finite_difference_hessian
   use urca_mod, only: johansen_result_t, johansen_test
   use random_mod, only: random_uniform
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: vars_deterministic_none = 0
   integer, parameter, public :: vars_deterministic_constant = 1
   integer, parameter, public :: vars_deterministic_trend = 2
   integer, parameter, public :: vars_deterministic_both = 3

   type, public :: vars_fit_t
      !! Reduced-form VAR coefficients, aligned design, and residual covariance.
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: ar(:, :, :)
      real(dp), allocatable :: deterministic(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: response(:, :)
      logical, allocatable :: restrictions(:, :)
      integer :: order = 0
      integer :: deterministic_type = vars_deterministic_constant
      integer :: season = 0
      integer :: exogenous_count = 0
      integer :: info = 0
   end type vars_fit_t

   type, public :: vars_selection_t
      !! Lag-selection criteria and minimizing orders in AIC, HQ, SC, and FPE order.
      real(dp), allocatable :: criteria(:, :)
      integer :: selected(4) = 0
      integer :: info = 0
   end type vars_selection_t

   type, public :: vars_test_t
      !! Scalar test statistic, reference degrees of freedom, and p-value.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: df1 = 0.0_dp
      real(dp) :: df2 = 0.0_dp
      integer :: info = 0
   end type vars_test_t

   type, public :: vars_serial_tests_t
      !! Multivariate portmanteau, Breusch-Godfrey, and Edgerton-Shukur tests.
      type(vars_test_t) :: portmanteau
      type(vars_test_t) :: adjusted_portmanteau
      type(vars_test_t) :: breusch_godfrey
      type(vars_test_t) :: edgerton_shukur
   end type vars_serial_tests_t

   type, public :: vars_normality_tests_t
      !! Multivariate Jarque-Bera, skewness, and kurtosis tests.
      type(vars_test_t) :: jarque_bera
      type(vars_test_t) :: skewness
      type(vars_test_t) :: kurtosis
      type(vars_test_t), allocatable :: univariate(:)
   end type vars_normality_tests_t

   type, public :: vars_arch_tests_t
      !! Multivariate and equationwise ARCH LM tests.
      type(vars_test_t) :: multivariate
      type(vars_test_t), allocatable :: univariate(:)
      integer :: info = 0
   end type vars_arch_tests_t

   type, public :: vars_stability_t
      !! Equationwise OLS-CUSUM fluctuation processes.
      real(dp), allocatable :: process(:, :)
      real(dp) :: boundary_95 = 1.358101516_dp
      integer :: info = 0
   end type vars_stability_t

   type, public :: vars_bq_t
      !! Blanchard-Quah short-run and long-run impact matrices.
      real(dp), allocatable :: short_run(:, :)
      real(dp), allocatable :: long_run(:, :)
      real(dp), allocatable :: structural_covariance(:, :)
      integer :: info = 0
   end type vars_bq_t

   type, public :: vars_svar_t
      !! Restricted structural A/B estimates and reduced-form covariance reproduction.
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: b(:, :)
      real(dp), allocatable :: impact(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: a_standard_error(:, :)
      real(dp), allocatable :: b_standard_error(:, :)
      real(dp) :: objective = huge(1.0_dp)
      type(vars_test_t) :: overidentification
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type vars_svar_t

   type, public :: vars_svec_t
      !! Short-run and long-run structural impacts identified from a Johansen model.
      real(dp), allocatable :: short_run(:, :)
      real(dp), allocatable :: long_run(:, :)
      real(dp), allocatable :: long_run_multiplier(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: short_run_standard_error(:, :)
      real(dp), allocatable :: long_run_standard_error(:, :)
      real(dp) :: objective = huge(1.0_dp)
      type(vars_test_t) :: overidentification
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type vars_svec_t

   type, public :: vars_vec2var_t
      !! Level-VAR coefficient matrices converted from a reduced-rank VECM.
      real(dp), allocatable :: ar(:, :, :)
      real(dp), allocatable :: pi(:, :)
      integer :: order = 0
      integer :: rank = 0
      integer :: info = 0
   end type vars_vec2var_t

   type, public :: vars_irf_bootstrap_t
      !! Point impulse responses and residual-bootstrap confidence limits.
      real(dp), allocatable :: point(:, :, :)
      real(dp), allocatable :: lower(:, :, :)
      real(dp), allocatable :: upper(:, :, :)
      real(dp) :: level = 0.95_dp
      integer :: successful_runs = 0
      integer :: info = 0
   end type vars_irf_bootstrap_t

   type, public :: vars_structural_bootstrap_t
      !! Point structural impacts and residual-bootstrap confidence limits.
      real(dp), allocatable :: point(:, :)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
      real(dp) :: level = 0.95_dp
      integer :: successful_runs = 0
      integer :: info = 0
   end type vars_structural_bootstrap_t

   type, public :: vars_svec_bootstrap_t
      !! Bootstrap standard errors for SVEC short- and long-run impact matrices.
      real(dp), allocatable :: short_run(:, :)
      real(dp), allocatable :: long_run(:, :)
      real(dp), allocatable :: short_run_standard_error(:, :)
      real(dp), allocatable :: long_run_standard_error(:, :)
      integer :: successful_runs = 0
      integer :: info = 0
   end type vars_svec_bootstrap_t

   type, public :: vars_structural_irf_t
      !! Structural impulse responses from horizon zero through the requested horizon.
      real(dp), allocatable :: response(:, :, :)
      integer :: info = 0
   end type vars_structural_irf_t

   type, public :: vars_fevd_t
      !! Forecast-error variance shares by response, shock, and forecast horizon.
      real(dp), allocatable :: share(:, :, :)
      integer :: info = 0
   end type vars_fevd_t

   interface display
      module procedure display_vars_fit
   end interface display

   public :: display, display_vars_fit
   public :: vars_fit, vars_select, vars_restrict, vars_restrict_ser
   public :: vars_phi, vars_psi, vars_roots
   public :: vars_serial_test, vars_normality_test, vars_arch_test
   public :: vars_instantaneous_causality, vars_ols_cusum, vars_bq, vars_svar
   public :: vars_svar_scoring
   public :: vars_svec, vars_vec2var
   public :: vars_irf_bootstrap, vars_svar_bootstrap, vars_granger_test
   public :: vars_granger_bootstrap
   public :: vars_svec_bootstrap
   public :: vars_svar_irf, vars_svec_irf, vars_svar_fevd, vars_svec_fevd

contains

   subroutine display_vars_fit(model, unit)
      !! Display a concise summary of a fitted reduced-form VAR model.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      character(len=24) :: deterministic_name
      integer :: destination, lag, row

      destination = output_unit
      if (present(unit)) destination = unit

      write(destination, '(a)') 'Reduced-form VAR fit'
      write(destination, '(a, i0)') 'Status: ', model%info
      if (model%info /= 0) return
      if (.not. allocated(model%ar)) then
         write(destination, '(a)') 'Model coefficients are not allocated.'
         return
      end if

      select case (model%deterministic_type)
      case (vars_deterministic_none)
         deterministic_name = 'none'
      case (vars_deterministic_constant)
         deterministic_name = 'constant'
      case (vars_deterministic_trend)
         deterministic_name = 'trend'
      case (vars_deterministic_both)
         deterministic_name = 'constant and trend'
      case default
         deterministic_name = 'unknown'
      end select

      write(destination, '(a, i0)') 'Variables: ', size(model%ar, 1)
      write(destination, '(a, i0)') 'Observations used: ', size(model%response, 1)
      write(destination, '(a, i0)') 'Lag order: ', model%order
      write(destination, '(a, a)') 'Deterministic terms: ', trim(deterministic_name)
      write(destination, '(a, i0)') 'Seasonal period: ', model%season
      write(destination, '(a, i0)') 'Exogenous regressors: ', model%exogenous_count

      if (allocated(model%deterministic)) then
         if (size(model%deterministic, 2) > 0) then
            write(destination, '(a)') 'Deterministic and exogenous coefficients:'
            do row = 1, size(model%deterministic, 1)
               write(destination, '(*(es14.6, 1x))') model%deterministic(row, :)
            end do
         end if
      end if

      do lag = 1, model%order
         write(destination, '(a, i0, a)') 'AR(', lag, ') coefficients:'
         do row = 1, size(model%ar, 1)
            write(destination, '(*(es14.6, 1x))') model%ar(row, :, lag)
         end do
      end do

      if (allocated(model%covariance)) then
         write(destination, '(a)') 'Residual covariance:'
         do row = 1, size(model%covariance, 1)
            write(destination, '(*(es14.6, 1x))') model%covariance(row, :)
         end do
      end if
   end subroutine display_vars_fit

   pure function vars_fit(series, order, deterministic_type, season, exogenous) result(out)
      !! Fit a VAR with vars-style deterministic, centered seasonal, and exogenous terms.
      real(dp), intent(in) :: series(:, :) !! Time-ordered endogenous observations.
      integer, intent(in) :: order !! Autoregressive lag order.
      integer, intent(in), optional :: deterministic_type !! Deterministic-term selector.
      integer, intent(in), optional :: season !! Seasonal period; zero omits seasonal dummies.
      real(dp), intent(in), optional :: exogenous(:, :) !! Contemporaneous exogenous regressors.
      type(vars_fit_t) :: out
      real(dp), allocatable :: design(:, :), response(:, :)
      integer :: selected_type, selected_season, n, k, rows, columns, t, lag, offset

      n = size(series, 1)
      k = size(series, 2)
      selected_type = vars_deterministic_constant
      if (present(deterministic_type)) selected_type = deterministic_type
      selected_season = 0
      if (present(season)) selected_season = abs(season)
      if (order < 1 .or. order >= n .or. k < 2 .or. selected_type < 0 .or. &
         selected_type > 3 .or. selected_season == 1) then
         out%info = 1
         return
      end if
      if (present(exogenous)) then
         if (size(exogenous, 1) /= n) then
            out%info = 2
            return
         end if
         out%exogenous_count = size(exogenous, 2)
      end if
      rows = n - order
      columns = k*order + deterministic_column_count(selected_type) + &
         max(0, selected_season - 1) + out%exogenous_count
      allocate(design(rows, columns), response(rows, k))
      design = 0.0_dp
      response = series(order + 1:, :)
      do t = 1, rows
         offset = 0
         do lag = 1, order
            design(t, offset + 1:offset + k) = series(order + t - lag, :)
            offset = offset + k
         end do
         call set_deterministic_row(design(t, :), offset, selected_type, &
            order + t, selected_season, out%exogenous_count, &
            exogenous_row(exogenous, order + t, out%exogenous_count))
      end do
      out = fit_design(design, response, order, selected_type, selected_season, &
         out%exogenous_count)
   end function vars_fit

   pure function vars_select(series, lag_max, deterministic_type, season, exogenous) result(out)
      !! Select VAR order using vars common-sample AIC, HQ, SC, and FPE criteria.
      real(dp), intent(in) :: series(:, :) !! Time-ordered endogenous observations.
      integer, intent(in) :: lag_max !! Largest candidate lag order.
      integer, intent(in), optional :: deterministic_type !! Deterministic-term selector.
      integer, intent(in), optional :: season !! Seasonal period.
      real(dp), intent(in), optional :: exogenous(:, :) !! Contemporaneous exogenous regressors.
      type(vars_selection_t) :: out
      type(vars_fit_t) :: candidate, full_candidate
      real(dp), allocatable :: common_series(:, :), common_exogenous(:, :)
      real(dp) :: logdet, determinant_value, penalty, sample
      integer :: selected_type, selected_season, k, lag, det_count, info, row, drop

      selected_type = vars_deterministic_constant
      if (present(deterministic_type)) selected_type = deterministic_type
      selected_season = 0
      if (present(season)) selected_season = abs(season)
      k = size(series, 2)
      if (lag_max < 1 .or. size(series, 1) <= 2*lag_max .or. k < 2) then
         out%info = 1
         return
      end if
      allocate(out%criteria(4, lag_max))
      out%criteria = huge(1.0_dp)
      common_series = series
      sample = real(size(series, 1) - lag_max, dp)
      det_count = deterministic_column_count(selected_type) + max(0, selected_season - 1)
      if (present(exogenous)) det_count = det_count + size(exogenous, 2)
      do lag = 1, lag_max
         if (present(exogenous)) then
            common_exogenous = exogenous
            full_candidate = vars_fit(common_series, lag, selected_type, selected_season, common_exogenous)
         else
            full_candidate = vars_fit(common_series, lag, selected_type, selected_season)
         end if
         if (full_candidate%info /= 0) cycle
         drop = lag_max - lag
         candidate = fit_design(full_candidate%design(drop + 1:, :), &
            full_candidate%response(drop + 1:, :), lag, selected_type, selected_season, &
            full_candidate%exogenous_count)
         if (candidate%info /= 0) cycle
         call covariance_on_common_sample(candidate, 0, determinant_value, info)
         if (info /= 0 .or. determinant_value <= tiny(1.0_dp)) cycle
         logdet = log(determinant_value)
         penalty = real(lag*k*k + k*det_count, dp)
         out%criteria(1, lag) = logdet + 2.0_dp*penalty/sample
         out%criteria(2, lag) = logdet + 2.0_dp*log(log(sample))*penalty/sample
         out%criteria(3, lag) = logdet + log(sample)*penalty/sample
         out%criteria(4, lag) = ((sample + real(lag*k + det_count, dp))/ &
            (sample - real(lag*k + det_count, dp)))**k*exp(logdet)
      end do
      do row = 1, 4
         out%selected(row) = minloc(out%criteria(row, :), dim=1)
      end do
   end function vars_select

   pure function vars_restrict(model, restrictions) result(out)
      !! Refit each VAR equation under a logical coefficient-inclusion matrix.
      type(vars_fit_t), intent(in) :: model !! Unrestricted fitted VAR.
      logical, intent(in) :: restrictions(:, :) !! Equation-by-regressor inclusion mask.
      type(vars_fit_t) :: out
      integer :: equation, selected_count, status

      out = model
      if (model%info /= 0 .or. size(restrictions, 1) /= size(model%response, 2) .or. &
         size(restrictions, 2) /= size(model%design, 2)) then
         out%info = 1
         return
      end if
      out%restrictions = restrictions
      out%coefficients = 0.0_dp
      do equation = 1, size(model%response, 2)
         selected_count = count(restrictions(equation, :))
         if (selected_count < 1) then
            out%info = 2
            return
         end if
         call restricted_equation(model%design, model%response(:, equation), &
            restrictions(equation, :), out%coefficients(equation, :), status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
      end do
      call update_fit_moments(out)
   end function vars_restrict

   pure function vars_restrict_ser(model, threshold) result(out)
      !! Apply vars sequential elimination of the smallest absolute coefficient t-ratio.
      type(vars_fit_t), intent(in) :: model !! Unrestricted fitted VAR.
      real(dp), intent(in), optional :: threshold !! Absolute t-ratio retention threshold.
      type(vars_fit_t) :: out
      logical, allocatable :: mask(:, :)
      real(dp), allocatable :: beta(:), se(:), residuals(:), selected(:, :)
      real(dp) :: cutoff, rss
      integer :: equation, weakest, status

      cutoff = 2.0_dp
      if (present(threshold)) cutoff = abs(threshold)
      if (model%info /= 0) then
         out = model
         return
      end if
      allocate(mask(size(model%response, 2), size(model%design, 2)))
      mask = .true.
      do equation = 1, size(model%response, 2)
         do
            selected = pack_design(model%design, mask(equation, :))
            call ols_fit(selected, model%response(:, equation), beta, se, residuals, rss, status)
            if (status /= 0 .or. size(beta) < 1) then
               out = model
               out%info = 2
               return
            end if
            if (minval(abs(beta)/max(se, tiny(1.0_dp))) >= cutoff) exit
            if (count(mask(equation, :)) == 1) then
               out = model
               out%info = 3
               return
            end if
            weakest = packed_index(mask(equation, :), &
               minloc(abs(beta)/max(se, tiny(1.0_dp)), dim=1))
            mask(equation, weakest) = .false.
         end do
      end do
      out = vars_restrict(model, mask)
   end function vars_restrict_ser

   pure function vars_phi(model, horizon) result(phi)
      !! Return reduced-form moving-average coefficient matrices through a horizon.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      integer, intent(in) :: horizon !! Largest response horizon.
      real(dp), allocatable :: phi(:, :, :)
      integer :: k, h, lag

      k = size(model%ar, 1)
      allocate(phi(k, k, max(1, horizon + 1)))
      phi = 0.0_dp
      if (horizon < 0 .or. model%info /= 0) return
      phi(:, :, 1) = identity(k)
      do h = 1, horizon
         do lag = 1, min(h, model%order)
            phi(:, :, h + 1) = phi(:, :, h + 1) + matmul(phi(:, :, h - lag + 1), &
               model%ar(:, :, lag))
         end do
      end do
   end function vars_phi

   pure function vars_psi(model, horizon) result(psi)
      !! Return orthogonalized impulse coefficients using the residual Cholesky factor.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      integer, intent(in) :: horizon !! Largest response horizon.
      real(dp), allocatable :: psi(:, :, :)
      real(dp), allocatable :: phi(:, :, :), lower(:, :)
      integer :: h, status

      phi = vars_phi(model, horizon)
      allocate(psi(size(phi, 1), size(phi, 2), lbound(phi, 3):ubound(phi, 3)))
      call cholesky_lower(model%covariance, lower, status)
      if (status /= 0) then
         psi = 0.0_dp
         return
      end if
      do h = lbound(phi, 3), ubound(phi, 3)
         psi(:, :, h) = matmul(phi(:, :, h), lower)
      end do
   end function vars_psi

   pure function vars_roots(model) result(roots)
      !! Return complex roots of the VAR companion matrix.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      complex(dp), allocatable :: roots(:)
      real(dp), allocatable :: companion(:, :)
      integer :: k, p, i, status

      k = size(model%ar, 1)
      p = model%order
      allocate(companion(k*p, k*p))
      companion = 0.0_dp
      do i = 1, p
         companion(1:k, (i - 1)*k + 1:i*k) = model%ar(:, :, i)
      end do
      if (p > 1) companion(k + 1:, :k*(p - 1)) = identity(k*(p - 1))
      call general_eigenvalues(companion, roots, status)
      if (status /= 0) roots = cmplx(0.0_dp, 0.0_dp, dp)
   end function vars_roots

   pure function vars_serial_test(model, portmanteau_lags, bg_lags) result(out)
      !! Compute vars asymptotic/adjusted portmanteau, BG, and ES serial tests.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      integer, intent(in), optional :: portmanteau_lags !! Portmanteau lag count.
      integer, intent(in), optional :: bg_lags !! BG and ES lag count.
      type(vars_serial_tests_t) :: out
      real(dp), allocatable :: c0(:, :), c0_inverse(:, :), ci(:, :)
      real(dp), allocatable :: lagged(:, :), unrestricted(:, :), restricted(:, :)
      real(dp), allocatable :: sigma0(:, :), sigma1(:, :)
      real(dp) :: trace_sum, adjusted_sum, determinant0, determinant1
      real(dp) :: r_squared, m, q, effective_n, exponent
      integer :: h, lag, n, k, status, predictors

      n = size(model%residuals, 1)
      k = size(model%residuals, 2)
      h = min(16, n - 1)
      if (present(portmanteau_lags)) h = min(abs(portmanteau_lags), n - 1)
      lag = min(5, n - 1)
      if (present(bg_lags)) lag = min(abs(bg_lags), n - 1)
      if (model%info /= 0 .or. h < 1 .or. lag < 1) then
         out%portmanteau%info = 1
         out%adjusted_portmanteau%info = 1
         out%breusch_godfrey%info = 1
         out%edgerton_shukur%info = 1
         return
      end if
      c0 = matmul(transpose(model%residuals), model%residuals)/real(n, dp)
      call invert_matrix(c0, c0_inverse, status)
      if (status /= 0) then
         out%portmanteau%info = 2
         return
      end if
      trace_sum = 0.0_dp
      adjusted_sum = 0.0_dp
      do predictors = 1, h
         ci = matmul(transpose(model%residuals(predictors + 1:, :)), &
            model%residuals(:n - predictors, :))/real(n, dp)
         r_squared = matrix_trace(matmul(transpose(ci), &
            matmul(c0_inverse, matmul(ci, c0_inverse))))
         trace_sum = trace_sum + r_squared
         adjusted_sum = adjusted_sum + r_squared/real(n - predictors, dp)
      end do
      out%portmanteau = chi_square_test(real(n, dp)*trace_sum, &
         real(k*k*(h - model%order), dp))
      out%adjusted_portmanteau = chi_square_test(real(n*n, dp)*adjusted_sum, &
         real(k*k*(h - model%order), dp))
      lagged = lagged_residual_design(model%residuals, lag)
      allocate(unrestricted(n, size(model%design, 2) + size(lagged, 2)))
      unrestricted(:, :size(model%design, 2)) = model%design
      unrestricted(:, size(model%design, 2) + 1:) = lagged
      call residual_covariance(unrestricted, model%residuals, sigma0, status)
      if (status /= 0) then
         out%breusch_godfrey%info = 2
         return
      end if
      call residual_covariance(model%design, model%residuals, sigma1, status)
      determinant0 = matrix_determinant(sigma0)
      determinant1 = matrix_determinant(sigma1)
      call invert_matrix(sigma1, c0_inverse, status)
      out%breusch_godfrey = chi_square_test(real(n, dp)*(real(k, dp) - &
         matrix_trace(matmul(c0_inverse, sigma0))), real(lag*k*k, dp))
      r_squared = 1.0_dp - determinant0/max(determinant1, tiny(1.0_dp))
      m = real(k*lag, dp)
      q = 0.5_dp*real(k, dp)*m - 1.0_dp
      predictors = size(model%design, 2)
      effective_n = real(n - predictors, dp) - m - 0.5_dp*real(k, dp) + 0.5_dp*m - 0.5_dp
      exponent = sqrt((real(k*k, dp)*m*m - 4.0_dp)/(real(k*k, dp) + m*m - 5.0_dp))
      out%edgerton_shukur%statistic = (1.0_dp - (1.0_dp - r_squared)**(1.0_dp/exponent))/ &
         max((1.0_dp - r_squared)**(1.0_dp/exponent), tiny(1.0_dp))* &
         (effective_n*exponent - q)/(real(k, dp)*m)
      out%edgerton_shukur%df1 = real(lag*k*k, dp)
      out%edgerton_shukur%df2 = floor(effective_n*exponent - q)
      out%edgerton_shukur%p_value = f_upper_tail(out%edgerton_shukur%statistic, &
         out%edgerton_shukur%df1, out%edgerton_shukur%df2)
   end function vars_serial_test

   pure function vars_normality_test(residuals) result(out)
      !! Compute the vars componentwise multivariate Jarque-Bera decomposition.
      real(dp), intent(in) :: residuals(:, :) !! VAR residual matrix.
      type(vars_normality_tests_t) :: out
      real(dp), allocatable :: centered(:, :), covariance(:, :), lower(:, :), inverse(:, :), z(:, :)
      real(dp), allocatable :: skewness(:), kurtosis(:)
      real(dp) :: skew_statistic, kurt_statistic, second, third, fourth
      integer :: n, k, status, j

      n = size(residuals, 1)
      k = size(residuals, 2)
      if (n < 2 .or. k < 1) then
         out%jarque_bera%info = 1
         return
      end if
      centered = residuals - spread(sum(residuals, dim=1)/real(n, dp), 1, n)
      covariance = matmul(transpose(centered), centered)/real(n, dp)
      call cholesky_lower(covariance, lower, status)
      if (status /= 0) then
         out%jarque_bera%info = 2
         return
      end if
      call invert_matrix(lower, inverse, status)
      z = matmul(centered, transpose(inverse))
      allocate(skewness(k), kurtosis(k))
      skewness = sum(z**3, dim=1)/real(n, dp)
      kurtosis = sum(z**4, dim=1)/real(n, dp)
      skew_statistic = real(n, dp)*dot_product(skewness, skewness)/6.0_dp
      kurt_statistic = real(n, dp)*dot_product(kurtosis - 3.0_dp, &
         kurtosis - 3.0_dp)/24.0_dp
      out%skewness = chi_square_test(skew_statistic, real(k, dp))
      out%kurtosis = chi_square_test(kurt_statistic, real(k, dp))
      out%jarque_bera = chi_square_test(skew_statistic + kurt_statistic, real(2*k, dp))
      allocate(out%univariate(k))
      do j = 1, k
         second = sum(centered(:, j)**2)/real(n, dp)
         third = sum(centered(:, j)**3)/real(n, dp)
         fourth = sum(centered(:, j)**4)/real(n, dp)
         if (second <= tiny(1.0_dp)) then
            out%univariate(j)%info = 1
         else
            out%univariate(j) = chi_square_test(real(n, dp)*(third*third/second**3)/6.0_dp + &
               real(n, dp)*(fourth/second**2 - 3.0_dp)**2/24.0_dp, 2.0_dp)
         end if
      end do
   end function vars_normality_test

   pure function vars_arch_test(residuals, multivariate_lags, univariate_lags) result(out)
      !! Compute vars multivariate and equationwise ARCH LM tests.
      real(dp), intent(in) :: residuals(:, :) !! VAR residual matrix.
      integer, intent(in), optional :: multivariate_lags !! Multivariate ARCH lag count.
      integer, intent(in), optional :: univariate_lags !! Equationwise ARCH lag count.
      type(vars_arch_tests_t) :: out
      real(dp), allocatable :: standardized(:, :), vech_values(:, :), design(:, :)
      real(dp), allocatable :: omega0(:, :), omega1(:, :), inverse(:, :), squared(:, :)
      real(dp) :: r_squared, rss0, rss1
      integer :: n, k, q, multi_lag, single_lag, t, i, j, column, status

      n = size(residuals, 1)
      k = size(residuals, 2)
      multi_lag = min(5, n - 2)
      if (present(multivariate_lags)) multi_lag = min(abs(multivariate_lags), n - 2)
      single_lag = min(16, n - 2)
      if (present(univariate_lags)) single_lag = min(abs(univariate_lags), n - 2)
      if (n < 4 .or. k < 1 .or. multi_lag < 1 .or. single_lag < 1) then
         out%info = 1
         return
      end if
      standardized = residuals - spread(sum(residuals, dim=1)/real(n, dp), 1, n)
      do j = 1, k
         standardized(:, j) = standardized(:, j)/max(sqrt(sum(standardized(:, j)**2)/ &
            real(n - 1, dp)), tiny(1.0_dp))
      end do
      q = k*(k + 1)/2
      allocate(vech_values(n, q))
      do t = 1, n
         column = 0
         do j = 1, k
            do i = j, k
               column = column + 1
               vech_values(t, column) = standardized(t, i)*standardized(t, j)
            end do
         end do
      end do
      design = lagged_with_intercept(vech_values, multi_lag)
      call residual_covariance(design, vech_values(multi_lag + 1:, :), omega1, status)
      omega0 = covariance_matrix_centered(vech_values(multi_lag + 1:, :))
      call invert_matrix(omega0, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      r_squared = 1.0_dp - matrix_trace(matmul(omega1, inverse))/real(q, dp)
      out%multivariate = chi_square_test(real(n - multi_lag, dp)*real(q, dp)*r_squared, &
         real(multi_lag*q*q, dp))
      allocate(out%univariate(k), squared(n, 1))
      do j = 1, k
         squared(:, 1) = standardized(:, j)**2
         design = lagged_with_intercept(squared, single_lag)
         rss0 = sum((squared(single_lag + 1:, 1) - &
            sum(squared(single_lag + 1:, 1))/real(n - single_lag, dp))**2)
         call scalar_regression_rss(design, squared(single_lag + 1:, 1), rss1, status)
         r_squared = 1.0_dp - rss1/max(rss0, tiny(1.0_dp))
         out%univariate(j) = chi_square_test(real(n - single_lag, dp)*r_squared, &
            real(single_lag, dp))
      end do
   end function vars_arch_test

   pure function vars_instantaneous_causality(model, causes) result(out)
      !! Test zero contemporaneous covariance between selected and remaining variables.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      logical, intent(in) :: causes(:) !! Variables assigned to the cause block.
      type(vars_test_t) :: out
      real(dp), allocatable :: selected(:), covariance(:, :), inverse(:, :)
      integer, allocatable :: first(:), second(:)
      integer :: k, pairs, i, j, a, b, status

      k = size(model%covariance, 1)
      if (size(causes) /= k .or. count(causes) < 1 .or. count(.not. causes) < 1) then
         out%info = 1
         return
      end if
      pairs = count(causes)*count(.not. causes)
      allocate(selected(pairs), covariance(pairs, pairs), first(pairs), second(pairs))
      a = 0
      do i = 1, k
         if (.not. causes(i)) cycle
         do j = 1, k
            if (causes(j)) cycle
            a = a + 1
            first(a) = i
            second(a) = j
            selected(a) = model%covariance(i, j)
         end do
      end do
      do a = 1, pairs
         do b = 1, pairs
            covariance(a, b) = model%covariance(first(a), first(b))* &
               model%covariance(second(a), second(b)) + &
               model%covariance(first(a), second(b))* &
               model%covariance(second(a), first(b))
         end do
      end do
      call invert_matrix(covariance, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out = chi_square_test(real(size(model%residuals, 1), dp)* &
         dot_product(selected, matmul(inverse, selected)), real(pairs, dp))
   end function vars_instantaneous_causality

   pure function vars_ols_cusum(model) result(out)
      !! Return equationwise standardized cumulative OLS residual processes.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      type(vars_stability_t) :: out
      real(dp) :: scale
      integer :: n, k, equation, t

      n = size(model%residuals, 1)
      k = size(model%residuals, 2)
      if (model%info /= 0 .or. n < 2) then
         out%info = 1
         return
      end if
      allocate(out%process(n, k))
      do equation = 1, k
         scale = sqrt(sum(model%residuals(:, equation)**2)/real(n, dp))*sqrt(real(n, dp))
         do t = 1, n
            out%process(t, equation) = sum(model%residuals(:t, equation))/ &
               max(scale, tiny(1.0_dp))
         end do
      end do
   end function vars_ols_cusum

   pure function vars_bq(model) result(out)
      !! Identify Blanchard-Quah shocks by a lower-triangular long-run impact matrix.
      type(vars_fit_t), intent(in) :: model !! Fitted stable VAR model.
      type(vars_bq_t) :: out
      real(dp), allocatable :: long_run_multiplier(:, :), inverse(:, :), covariance(:, :)
      integer :: lag, k, status

      k = size(model%ar, 1)
      long_run_multiplier = identity(k)
      do lag = 1, model%order
         long_run_multiplier = long_run_multiplier - model%ar(:, :, lag)
      end do
      call invert_matrix(long_run_multiplier, inverse, status)
      if (status /= 0) then
         out%info = 1
         return
      end if
      covariance = matmul(inverse, matmul(model%covariance, transpose(inverse)))
      call cholesky_lower(covariance, out%long_run, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%short_run = matmul(long_run_multiplier, out%long_run)
      out%structural_covariance = matmul(out%short_run, transpose(out%short_run))
   end function vars_bq

   pure function vars_svar(model, a_template, b_template, estimate_a, estimate_b, &
      initial, max_iterations, tolerance) result(out)
      !! Estimate a restricted structural A/B model by direct Gaussian likelihood.
      type(vars_fit_t), intent(in) :: model !! Fitted reduced-form VAR.
      real(dp), intent(in) :: a_template(:, :) !! Fixed A entries and starting values.
      real(dp), intent(in) :: b_template(:, :) !! Fixed B entries and starting values.
      logical, intent(in) :: estimate_a(:, :) !! Mask selecting free A entries.
      logical, intent(in) :: estimate_b(:, :) !! Mask selecting free B entries.
      real(dp), intent(in), optional :: initial(:) !! Free parameters in A-then-B order.
      integer, intent(in), optional :: max_iterations !! Optimizer iteration limit.
      real(dp), intent(in), optional :: tolerance !! Optimizer convergence tolerance.
      type(vars_svar_t) :: out
      type(optimization_result_t) :: optimized
      real(dp), allocatable :: starting(:), inverse_a(:, :), hessian(:, :), parameter_covariance(:, :)
      real(dp) :: selected_tolerance
      integer :: k, free_count, limit, status, degrees

      k = size(model%covariance, 1)
      free_count = count(estimate_a) + count(estimate_b)
      limit = 1000
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-7_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (model%info /= 0 .or. any(shape(a_template) /= [k, k]) .or. &
         any(shape(b_template) /= [k, k]) .or. &
         any(shape(estimate_a) /= [k, k]) .or. any(shape(estimate_b) /= [k, k]) .or. &
         free_count < 1 .or. free_count > k*(k + 1)/2) then
         out%info = 1
         return
      end if
      if (present(initial)) then
         if (size(initial) /= free_count) then
            out%info = 2
            return
         end if
         starting = initial
      else
         allocate(starting(free_count))
         call pack_structural_parameters(a_template, b_template, estimate_a, estimate_b, starting)
      end if
      optimized = nelder_mead_minimize(objective, starting, limit, selected_tolerance, 0.1_dp)
      call unpack_structural_parameters(optimized%parameters, a_template, b_template, &
         estimate_a, estimate_b, out%a, out%b)
      call invert_matrix(out%a, inverse_a, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%impact = matmul(inverse_a, out%b)
      out%covariance = matmul(out%impact, transpose(out%impact))
      out%objective = optimized%objective
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%info = optimized%info
      hessian = finite_difference_hessian(objective, optimized%parameters)
      call invert_matrix(hessian, parameter_covariance, status)
      call structural_standard_errors(parameter_covariance, estimate_a, estimate_b, &
         out%a_standard_error, out%b_standard_error, status)
      if (status /= 0 .and. out%info == 0) out%info = 20 + status
      degrees = k*(k + 1)/2 - free_count
      if (degrees > 0) out%overidentification = chi_square_test(real(size(model%residuals, 1), dp)* &
         log(max(matrix_determinant(out%covariance), tiny(1.0_dp))/ &
         max(matrix_determinant(model%covariance), tiny(1.0_dp))), real(degrees, dp))

   contains

      pure function objective(parameters) result(value)
         !! Evaluate the concentrated structural Gaussian objective.
         real(dp), intent(in) :: parameters(:) !! Free structural parameters.
         real(dp) :: value
         real(dp), allocatable :: a(:, :), b(:, :), inverse_b(:, :), transformed(:, :)
         real(dp) :: determinant_a, determinant_b, quadratic
         integer :: evaluation_status

         call unpack_structural_parameters(parameters, a_template, b_template, &
            estimate_a, estimate_b, a, b)
         determinant_a = abs(matrix_determinant(a))
         determinant_b = abs(matrix_determinant(b))
         if (determinant_a <= tiny(1.0_dp) .or. determinant_b <= tiny(1.0_dp)) then
            value = huge(1.0_dp)
            return
         end if
         call invert_matrix(b, inverse_b, evaluation_status)
         if (evaluation_status /= 0) then
            value = huge(1.0_dp)
            return
         end if
         transformed = matmul(inverse_b, a)
         quadratic = matrix_trace(matmul(transformed, &
            matmul(model%covariance, transpose(transformed))))
         value = real(size(model%residuals, 1), dp)*(-log(determinant_a) + &
            log(determinant_b) + 0.5_dp*quadratic)
      end function objective
   end function vars_svar

   pure function vars_svar_scoring(model, a_template, b_template, estimate_a, &
      estimate_b, initial, max_iterations, tolerance, max_step) result(out)
      !! Estimate a restricted structural A/B model by the vars Fisher-scoring recursion.
      type(vars_fit_t), intent(in) :: model !! Fitted reduced-form VAR.
      real(dp), intent(in) :: a_template(:, :) !! Fixed A entries and starting values.
      real(dp), intent(in) :: b_template(:, :) !! Fixed B entries and starting values.
      logical, intent(in) :: estimate_a(:, :) !! Mask selecting free A entries.
      logical, intent(in) :: estimate_b(:, :) !! Mask selecting free B entries.
      real(dp), intent(in), optional :: initial(:) !! Free parameters in A-then-B order.
      integer, intent(in), optional :: max_iterations !! Scoring iteration limit.
      real(dp), intent(in), optional :: tolerance !! Parameter convergence tolerance.
      real(dp), intent(in), optional :: max_step !! Largest absolute scoring step.
      type(vars_svar_t) :: out
      real(dp), allocatable :: parameters(:), mapping(:, :)
      real(dp), allocatable :: a(:, :), b(:, :), inverse_a(:, :), inverse_b(:, :)
      real(dp), allocatable :: binva(:, :), inverse_binva(:, :), identity_k(:, :)
      real(dp), allocatable :: information_left(:, :), information_middle(:, :)
      real(dp), allocatable :: information_right(:, :), information(:, :), information_inverse(:, :)
      real(dp), allocatable :: score_core(:), score_ab(:), score(:), direction(:)
      real(dp), allocatable :: parameter_covariance(:, :)
      real(dp) :: selected_tolerance, selected_max_step, scale
      integer :: k, k_squared, free_count, limit, iteration, status, offset, i, j, degrees

      k = size(model%covariance, 1)
      k_squared = k*k
      free_count = count(estimate_a) + count(estimate_b)
      limit = 100
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-7_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_step = 1.0_dp
      if (present(max_step)) selected_max_step = max_step
      if (model%info /= 0 .or. any(shape(a_template) /= [k, k]) .or. &
         any(shape(b_template) /= [k, k]) .or. &
         any(shape(estimate_a) /= [k, k]) .or. any(shape(estimate_b) /= [k, k]) .or. &
         free_count < 1 .or. free_count > k*(k + 1)/2 .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. selected_max_step <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(initial)) then
         if (size(initial) /= free_count) then
            out%info = 2
            return
         end if
         parameters = initial
      else
         allocate(parameters(free_count))
         call pack_structural_parameters(a_template, b_template, estimate_a, estimate_b, parameters)
      end if
      allocate(mapping(2*k_squared, free_count))
      mapping = 0.0_dp
      offset = 0
      do j = 1, k
         do i = 1, k
            if (.not. estimate_a(i, j)) cycle
            offset = offset + 1
            mapping(i + (j - 1)*k, offset) = 1.0_dp
         end do
      end do
      do j = 1, k
         do i = 1, k
            if (.not. estimate_b(i, j)) cycle
            offset = offset + 1
            mapping(k_squared + i + (j - 1)*k, offset) = 1.0_dp
         end do
      end do
      identity_k = identity(k)
      information_middle = identity(k_squared) + commutation_matrix(k)
      allocate(information_left(2*k_squared, k_squared))
      allocate(information_right(k_squared, 2*k_squared))
      allocate(score_core(k_squared), score_ab(2*k_squared), score(free_count))
      allocate(direction(free_count))
      do iteration = 1, limit
         call unpack_structural_parameters(parameters, a_template, b_template, &
            estimate_a, estimate_b, a, b)
         call invert_matrix(a, inverse_a, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         call invert_matrix(b, inverse_b, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         binva = matmul(inverse_b, a)
         call invert_matrix(binva, inverse_binva, status)
         if (status /= 0) then
            out%info = 30 + status
            return
         end if
         information_left(:k_squared, :) = kronecker_product(inverse_binva, transpose(inverse_b))
         information_left(k_squared + 1:, :) = &
            -kronecker_product(identity_k, transpose(inverse_b))
         information_right(:, :k_squared) = &
            kronecker_product(transpose(inverse_binva), inverse_b)
         information_right(:, k_squared + 1:) = -kronecker_product(identity_k, inverse_b)
         information = real(size(model%residuals, 1), dp)*matmul(transpose(mapping), &
            matmul(information_left, matmul(information_middle, &
            matmul(information_right, mapping))))
         call invert_matrix(information, information_inverse, status)
         if (status /= 0) then
            out%info = 40 + status
            return
         end if
         score_core = real(size(model%residuals, 1), dp)*( &
            reshape(transpose(inverse_binva), [k_squared]) - &
            matmul(kronecker_product(model%covariance, identity_k), &
            reshape(binva, [k_squared])))
         score_ab(:k_squared) = matmul(kronecker_product(identity_k, &
            transpose(inverse_b)), score_core)
         score_ab(k_squared + 1:) = -matmul(kronecker_product(binva, &
            transpose(inverse_b)), score_core)
         score = matmul(transpose(mapping), score_ab)
         direction = matmul(information_inverse, score)
         scale = maxval(abs(direction))
         if (scale > selected_max_step) direction = direction*selected_max_step/scale
         parameters = parameters + direction
         out%iterations = iteration
         if (maxval(abs(direction)) <= selected_tolerance) then
            out%converged = .true.
            exit
         end if
      end do
      if (.not. out%converged) out%info = 4
      call unpack_structural_parameters(parameters, a_template, b_template, &
         estimate_a, estimate_b, out%a, out%b)
      call invert_matrix(out%a, inverse_a, status)
      if (status /= 0) then
         out%info = 50 + status
         return
      end if
      out%impact = matmul(inverse_a, out%b)
      out%covariance = matmul(out%impact, transpose(out%impact))
      call invert_matrix(out%b, inverse_b, status)
      if (status /= 0) then
         out%info = 55 + status
         return
      end if
      out%objective = real(size(model%residuals, 1), dp)*( &
         -log(max(abs(matrix_determinant(out%a)), tiny(1.0_dp))) + &
         log(max(abs(matrix_determinant(out%b)), tiny(1.0_dp))) + &
         0.5_dp*matrix_trace(matmul(matmul(inverse_b, out%a), &
         matmul(model%covariance, transpose(matmul(inverse_b, out%a))))))
      parameter_covariance = information_inverse
      call structural_standard_errors(parameter_covariance, estimate_a, estimate_b, &
         out%a_standard_error, out%b_standard_error, status)
      if (status /= 0 .and. out%info == 0) out%info = 60 + status
      degrees = k*(k + 1)/2 - free_count
      if (degrees > 0) out%overidentification = chi_square_test(real(size(model%residuals, 1), dp)* &
         log(max(matrix_determinant(out%covariance), tiny(1.0_dp))/ &
         max(matrix_determinant(model%covariance), tiny(1.0_dp))), real(degrees, dp))
   end function vars_svar_scoring

   pure subroutine pack_structural_parameters(a, b, estimate_a, estimate_b, parameters)
      !! Pack free structural parameters in A-then-B column-major order.
      real(dp), intent(in) :: a(:, :) !! Structural A template.
      real(dp), intent(in) :: b(:, :) !! Structural B template.
      logical, intent(in) :: estimate_a(:, :) !! Free A mask.
      logical, intent(in) :: estimate_b(:, :) !! Free B mask.
      real(dp), intent(out) :: parameters(:) !! Packed parameters.
      integer :: i, j, offset

      offset = 0
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            if (.not. estimate_a(i, j)) cycle
            offset = offset + 1
            parameters(offset) = a(i, j)
         end do
      end do
      do j = 1, size(b, 2)
         do i = 1, size(b, 1)
            if (.not. estimate_b(i, j)) cycle
            offset = offset + 1
            parameters(offset) = b(i, j)
         end do
      end do
   end subroutine pack_structural_parameters

   pure subroutine unpack_structural_parameters(parameters, a_template, b_template, &
      estimate_a, estimate_b, a, b)
      !! Expand packed free parameters into structural A and B matrices.
      real(dp), intent(in) :: parameters(:) !! Packed parameters.
      real(dp), intent(in) :: a_template(:, :) !! Structural A template.
      real(dp), intent(in) :: b_template(:, :) !! Structural B template.
      logical, intent(in) :: estimate_a(:, :) !! Free A mask.
      logical, intent(in) :: estimate_b(:, :) !! Free B mask.
      real(dp), allocatable, intent(out) :: a(:, :) !! Expanded A matrix.
      real(dp), allocatable, intent(out) :: b(:, :) !! Expanded B matrix.
      integer :: i, j, offset

      a = a_template
      b = b_template
      offset = 0
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            if (.not. estimate_a(i, j)) cycle
            offset = offset + 1
            a(i, j) = parameters(offset)
         end do
      end do
      do j = 1, size(b, 2)
         do i = 1, size(b, 1)
            if (.not. estimate_b(i, j)) cycle
            offset = offset + 1
            b(i, j) = parameters(offset)
         end do
      end do
   end subroutine unpack_structural_parameters

   pure subroutine structural_standard_errors(covariance, estimate_a, estimate_b, &
      a_standard_error, b_standard_error, info)
      !! Expand free-parameter covariance diagonal into A and B standard-error matrices.
      real(dp), intent(in) :: covariance(:, :) !! Free-parameter covariance matrix.
      logical, intent(in) :: estimate_a(:, :) !! Free A mask.
      logical, intent(in) :: estimate_b(:, :) !! Free B mask.
      real(dp), allocatable, intent(out) :: a_standard_error(:, :) !! A standard errors.
      real(dp), allocatable, intent(out) :: b_standard_error(:, :) !! B standard errors.
      integer, intent(out) :: info !! Status code.
      integer :: i, j, offset

      allocate(a_standard_error(size(estimate_a, 1), size(estimate_a, 2)))
      allocate(b_standard_error(size(estimate_b, 1), size(estimate_b, 2)))
      a_standard_error = 0.0_dp
      b_standard_error = 0.0_dp
      if (size(covariance, 1) /= count(estimate_a) + count(estimate_b)) then
         info = 1
         return
      end if
      offset = 0
      do j = 1, size(estimate_a, 2)
         do i = 1, size(estimate_a, 1)
            if (.not. estimate_a(i, j)) cycle
            offset = offset + 1
            a_standard_error(i, j) = sqrt(max(0.0_dp, covariance(offset, offset)))
         end do
      end do
      do j = 1, size(estimate_b, 2)
         do i = 1, size(estimate_b, 1)
            if (.not. estimate_b(i, j)) cycle
            offset = offset + 1
            b_standard_error(i, j) = sqrt(max(0.0_dp, covariance(offset, offset)))
         end do
      end do
      info = 0
   end subroutine structural_standard_errors

   pure function vars_svec(johansen, rank, short_run_zero, long_run_zero, &
      max_iterations, tolerance) result(out)
      !! Estimate vars SVEC short-run impacts under linear short- and long-run zeros.
      type(johansen_result_t), intent(in) :: johansen !! Fitted Johansen model.
      integer, intent(in) :: rank !! Cointegration rank.
      logical, intent(in) :: short_run_zero(:, :) !! Zero restrictions on short-run impacts.
      logical, intent(in) :: long_run_zero(:, :) !! Zero restrictions on long-run impacts.
      integer, intent(in), optional :: max_iterations !! Optimizer iteration limit.
      real(dp), intent(in), optional :: tolerance !! Optimizer convergence tolerance.
      type(vars_svec_t) :: out
      type(optimization_result_t) :: optimized
      real(dp), allocatable :: alpha(:, :), beta(:, :), alpha_orthogonal(:, :)
      real(dp), allocatable :: beta_orthogonal(:, :), gamma_sum(:, :), middle(:, :), inverse(:, :)
      real(dp), allocatable :: restrictions(:, :), basis(:, :), lower(:, :), starting(:)
      real(dp), allocatable :: hessian(:, :), coordinate_covariance(:, :), vector_covariance(:, :)
      real(dp), allocatable :: long_transformation(:, :), long_covariance(:, :)
      real(dp) :: selected_tolerance
      integer :: k, lag_count, free_count, restriction_count, limit, status, i, j, l, row

      k = johansen%variables
      limit = 1000
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-7_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (johansen%info /= 0 .or. rank < 1 .or. rank >= k .or. &
         any(shape(short_run_zero) /= [k, k]) .or. &
         any(shape(long_run_zero) /= [k, k])) then
         out%info = 1
         return
      end if
      alpha = johansen%loading(:, :rank)
      beta = johansen%cointegration(:k, :rank)
      alpha_orthogonal = null_space(transpose(alpha), k - rank)
      if (size(alpha_orthogonal, 2) /= k - rank) then
         out%info = 10
         return
      end if
      beta_orthogonal = null_space(transpose(beta), k - rank)
      if (size(beta_orthogonal, 2) /= k - rank) then
         out%info = 20
         return
      end if
      allocate(gamma_sum(k, k))
      gamma_sum = 0.0_dp
      lag_count = max(0, johansen%lag - 1)
      if (lag_count > 0 .and. size(johansen%gamma, 2) >= k*lag_count) then
         do l = 1, lag_count
            gamma_sum = gamma_sum + johansen%gamma(:, size(johansen%gamma, 2) - &
               k*lag_count + (l - 1)*k + 1:size(johansen%gamma, 2) - k*lag_count + l*k)
         end do
      end if
      middle = matmul(transpose(alpha_orthogonal), &
         matmul(identity(k) - gamma_sum, beta_orthogonal))
      call invert_matrix(middle, inverse, status)
      if (status /= 0) then
         out%info = 30 + status
         return
      end if
      out%long_run_multiplier = matmul(beta_orthogonal, &
         matmul(inverse, transpose(alpha_orthogonal)))
      restriction_count = count(short_run_zero) + count(long_run_zero)
      allocate(restrictions(restriction_count, k*k))
      restrictions = 0.0_dp
      row = 0
      do j = 1, k
         do i = 1, k
            if (.not. short_run_zero(i, j)) cycle
            row = row + 1
            restrictions(row, i + (j - 1)*k) = 1.0_dp
         end do
      end do
      do j = 1, k
         do i = 1, k
            if (.not. long_run_zero(i, j)) cycle
            row = row + 1
            do l = 1, k
               restrictions(row, l + (j - 1)*k) = out%long_run_multiplier(i, l)
            end do
         end do
      end do
      free_count = k*k - matrix_rank_from_eigen(restrictions)
      if (free_count < 1 .or. free_count > k*(k + 1)/2) then
         out%info = 2
         return
      end if
      basis = null_space(restrictions, free_count)
      if (size(basis, 2) /= free_count) then
         out%info = 40
         return
      end if
      call cholesky_lower(johansen%delta, lower, status)
      if (status /= 0) then
         out%info = 50 + status
         return
      end if
      starting = matmul(transpose(basis), reshape(lower, [k*k]))
      optimized = nelder_mead_minimize(objective, starting, limit, selected_tolerance, 0.1_dp)
      out%short_run = reshape(matmul(basis, optimized%parameters), [k, k])
      do j = 1, k
         if (out%short_run(j, j) < 0.0_dp) out%short_run(:, j) = -out%short_run(:, j)
      end do
      out%long_run = matmul(out%long_run_multiplier, out%short_run)
      out%covariance = matmul(out%short_run, transpose(out%short_run))
      out%objective = optimized%objective
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%info = optimized%info
      hessian = finite_difference_hessian(objective, optimized%parameters)
      call invert_matrix(hessian, coordinate_covariance, status)
      if (status == 0) then
         vector_covariance = matmul(basis, matmul(coordinate_covariance, transpose(basis)))
         out%short_run_standard_error = reshape(sqrt(max(0.0_dp, &
            [(vector_covariance(i, i), i=1, k*k)])), [k, k])
         long_transformation = kronecker_product(identity(k), out%long_run_multiplier)
         long_covariance = matmul(long_transformation, &
            matmul(vector_covariance, transpose(long_transformation)))
         out%long_run_standard_error = reshape(sqrt(max(0.0_dp, &
            [(long_covariance(i, i), i=1, k*k)])), [k, k])
      else if (out%info == 0) then
         out%info = 60 + status
      end if
      if (k*(k + 1)/2 > free_count) out%overidentification = chi_square_test( &
         real(size(johansen%r0, 1), dp)*log(max(matrix_determinant(out%covariance), &
         tiny(1.0_dp))/max(matrix_determinant(johansen%delta), tiny(1.0_dp))), &
         real(k*(k + 1)/2 - free_count, dp))

   contains

      pure function objective(parameters) result(value)
         !! Evaluate the structural covariance likelihood on the restriction null space.
         real(dp), intent(in) :: parameters(:) !! Free null-space coordinates.
         real(dp) :: value
         real(dp), allocatable :: b(:, :), inverse_b(:, :)
         real(dp) :: determinant_b, quadratic
         integer :: evaluation_status

         b = reshape(matmul(basis, parameters), [k, k])
         determinant_b = abs(matrix_determinant(b))
         if (determinant_b <= tiny(1.0_dp)) then
            value = huge(1.0_dp)
            return
         end if
         call invert_matrix(b, inverse_b, evaluation_status)
         if (evaluation_status /= 0) then
            value = huge(1.0_dp)
            return
         end if
         quadratic = matrix_trace(matmul(inverse_b, &
            matmul(johansen%delta, transpose(inverse_b))))
         value = real(size(johansen%r0, 1), dp)*(log(determinant_b) + 0.5_dp*quadratic)
      end function objective
   end function vars_svec

   pure function null_space(matrix, dimension) result(basis)
      !! Return an orthonormal basis for a requested matrix null-space dimension.
      real(dp), intent(in) :: matrix(:, :) !! Constraint matrix.
      integer, intent(in) :: dimension !! Expected null-space dimension.
      real(dp), allocatable :: basis(:, :)
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :), gram(:, :)
      integer :: columns, status

      columns = size(matrix, 2)
      if (dimension < 0 .or. dimension > columns) then
         allocate(basis(columns, 0))
         return
      end if
      gram = matmul(transpose(matrix), matrix)
      call symmetric_eigen(gram, eigenvalues, eigenvectors, status)
      if (status /= 0) then
         allocate(basis(columns, 0))
         return
      end if
      allocate(basis(columns, dimension))
      if (dimension > 0) basis = eigenvectors(:, columns - dimension + 1:)
   end function null_space

   pure integer function matrix_rank_from_eigen(matrix) result(rank_value)
      !! Return numerical row rank from eigenvalues of the normal matrix.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      real(dp) :: threshold
      integer :: status

      if (size(matrix, 1) == 0 .or. size(matrix, 2) == 0) then
         rank_value = 0
         return
      end if
      call symmetric_eigen(matmul(transpose(matrix), matrix), eigenvalues, eigenvectors, status)
      if (status /= 0) then
         rank_value = size(matrix, 1)
         return
      end if
      threshold = sqrt(epsilon(1.0_dp))*max(1.0_dp, maxval(abs(eigenvalues)))
      rank_value = count(eigenvalues > threshold)
   end function matrix_rank_from_eigen

   pure function vars_vec2var(johansen, rank) result(out)
      !! Convert a Johansen VECM into level-VAR coefficient matrices as vars vec2var.
      type(johansen_result_t), intent(in) :: johansen !! Fitted Johansen model.
      integer, intent(in) :: rank !! Cointegration rank.
      type(vars_vec2var_t) :: out
      real(dp), allocatable :: gamma(:, :, :)
      integer :: k, p, lag, start_column

      k = johansen%variables
      p = johansen%lag
      if (johansen%info /= 0 .or. rank < 1 .or. rank >= k .or. p < 1) then
         out%info = 1
         return
      end if
      out%order = p
      out%rank = rank
      out%pi = matmul(johansen%loading(:, :rank), &
         transpose(johansen%cointegration(:k, :rank)))
      allocate(out%ar(k, k, p))
      out%ar = 0.0_dp
      if (p == 1) then
         out%ar(:, :, 1) = identity(k) + out%pi
         return
      end if
      if (size(johansen%gamma, 2) < k*(p - 1)) then
         out%info = 2
         return
      end if
      allocate(gamma(k, k, p - 1))
      start_column = size(johansen%gamma, 2) - k*(p - 1)
      do lag = 1, p - 1
         gamma(:, :, lag) = johansen%gamma(:, start_column + (lag - 1)*k + 1: &
            start_column + lag*k)
      end do
      if (trim(johansen%specification) == 'transitory') then
         out%ar(:, :, 1) = identity(k) + out%pi + gamma(:, :, 1)
         do lag = 2, p - 1
            out%ar(:, :, lag) = gamma(:, :, lag) - gamma(:, :, lag - 1)
         end do
         out%ar(:, :, p) = -gamma(:, :, p - 1)
      else
         out%ar(:, :, 1) = identity(k) + gamma(:, :, 1)
         do lag = 2, p - 1
            out%ar(:, :, lag) = gamma(:, :, lag) - gamma(:, :, lag - 1)
         end do
         out%ar(:, :, p) = out%pi - gamma(:, :, p - 1)
      end if
   end function vars_vec2var

   pure function vars_svar_irf(model, structural, horizon, cumulative) result(out)
      !! Compute structural impulse responses for an identified stationary VAR.
      type(vars_fit_t), intent(in) :: model !! Fitted reduced-form VAR.
      type(vars_svar_t), intent(in) :: structural !! Identified structural VAR.
      integer, intent(in) :: horizon !! Largest response horizon.
      logical, intent(in), optional :: cumulative !! Cumulate responses over horizons.
      type(vars_structural_irf_t) :: out
      logical :: use_cumulative

      use_cumulative = .false.
      if (present(cumulative)) use_cumulative = cumulative
      if (model%info /= 0 .or. structural%info /= 0 .or. horizon < 0 .or. &
         size(model%ar, 1) /= size(structural%impact, 1)) then
         out%info = 1
         return
      end if
      out%response = structural_responses(model%ar, structural%impact, horizon)
      if (use_cumulative) call cumulative_responses(out%response)
   end function vars_svar_irf

   pure function vars_svec_irf(johansen, rank, structural, horizon, cumulative) result(out)
      !! Compute structural impulse responses for an identified cointegrated VAR.
      type(johansen_result_t), intent(in) :: johansen !! Fitted Johansen model.
      integer, intent(in) :: rank !! Cointegration rank.
      type(vars_svec_t), intent(in) :: structural !! Identified structural VEC model.
      integer, intent(in) :: horizon !! Largest response horizon.
      logical, intent(in), optional :: cumulative !! Cumulate responses over horizons.
      type(vars_structural_irf_t) :: out
      type(vars_vec2var_t) :: level_var
      logical :: use_cumulative

      use_cumulative = .false.
      if (present(cumulative)) use_cumulative = cumulative
      level_var = vars_vec2var(johansen, rank)
      if (level_var%info /= 0 .or. structural%info /= 0 .or. horizon < 0 .or. &
         size(level_var%ar, 1) /= size(structural%short_run, 1)) then
         out%info = 1
         return
      end if
      out%response = structural_responses(level_var%ar, structural%short_run, horizon)
      if (use_cumulative) call cumulative_responses(out%response)
   end function vars_svec_irf

   pure function vars_svar_fevd(model, structural, horizon) result(out)
      !! Compute structural forecast-error variance shares for a stationary VAR.
      type(vars_fit_t), intent(in) :: model !! Fitted reduced-form VAR.
      type(vars_svar_t), intent(in) :: structural !! Identified structural VAR.
      integer, intent(in) :: horizon !! Largest positive forecast horizon.
      type(vars_fevd_t) :: out
      real(dp), allocatable :: response(:, :, :)

      if (model%info /= 0 .or. structural%info /= 0 .or. horizon < 1 .or. &
         size(model%ar, 1) /= size(structural%impact, 1)) then
         out%info = 1
         return
      end if
      response = structural_responses(model%ar, structural%impact, horizon - 1)
      out%share = variance_decomposition(response)
   end function vars_svar_fevd

   pure function vars_svec_fevd(johansen, rank, structural, horizon) result(out)
      !! Compute structural forecast-error variance shares for a cointegrated VAR.
      type(johansen_result_t), intent(in) :: johansen !! Fitted Johansen model.
      integer, intent(in) :: rank !! Cointegration rank.
      type(vars_svec_t), intent(in) :: structural !! Identified structural VEC model.
      integer, intent(in) :: horizon !! Largest positive forecast horizon.
      type(vars_fevd_t) :: out
      type(vars_vec2var_t) :: level_var
      real(dp), allocatable :: response(:, :, :)

      level_var = vars_vec2var(johansen, rank)
      if (level_var%info /= 0 .or. structural%info /= 0 .or. horizon < 1 .or. &
         size(level_var%ar, 1) /= size(structural%short_run, 1)) then
         out%info = 1
         return
      end if
      response = structural_responses(level_var%ar, structural%short_run, horizon - 1)
      out%share = variance_decomposition(response)
   end function vars_svec_fevd

   pure function structural_responses(ar, impact, horizon) result(response)
      !! Form moving-average responses from level-VAR coefficients and an impact matrix.
      real(dp), intent(in) :: ar(:, :, :) !! Level-VAR autoregressive coefficient matrices.
      real(dp), intent(in) :: impact(:, :) !! Contemporaneous structural impact matrix.
      integer, intent(in) :: horizon !! Largest response horizon.
      real(dp), allocatable :: response(:, :, :)
      real(dp), allocatable :: phi(:, :, :)
      integer :: k, h, lag

      k = size(ar, 1)
      allocate(phi(k, k, horizon + 1), response(k, k, horizon + 1))
      phi = 0.0_dp
      phi(:, :, 1) = identity(k)
      do h = 1, horizon
         do lag = 1, min(size(ar, 3), h)
            phi(:, :, h + 1) = phi(:, :, h + 1) + &
               matmul(ar(:, :, lag), phi(:, :, h - lag + 1))
         end do
      end do
      do h = 1, horizon + 1
         response(:, :, h) = matmul(phi(:, :, h), impact)
      end do
   end function structural_responses

   pure function variance_decomposition(response) result(share)
      !! Normalize cumulative squared structural responses into variance shares.
      real(dp), intent(in) :: response(:, :, :) !! Structural responses from horizon zero.
      real(dp), allocatable :: share(:, :, :)
      real(dp), allocatable :: cumulative(:, :)
      real(dp) :: total
      integer :: h, equation

      allocate(share(size(response, 1), size(response, 2), size(response, 3)))
      allocate(cumulative(size(response, 1), size(response, 2)))
      cumulative = 0.0_dp
      do h = 1, size(response, 3)
         cumulative = cumulative + response(:, :, h)**2
         do equation = 1, size(response, 1)
            total = sum(cumulative(equation, :))
            if (total > tiny(1.0_dp)) then
               share(equation, :, h) = cumulative(equation, :)/total
            else
               share(equation, :, h) = 0.0_dp
            end if
         end do
      end do
   end function variance_decomposition

   function vars_irf_bootstrap(model, horizon, runs, level, orthogonal, cumulative) result(out)
      !! Compute residual-bootstrap confidence limits for reduced-form or orthogonal IRFs.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      integer, intent(in) :: horizon !! Largest impulse-response horizon.
      integer, intent(in) :: runs !! Number of bootstrap replications.
      real(dp), intent(in), optional :: level !! Central confidence level.
      logical, intent(in), optional :: orthogonal !! Use Cholesky-orthogonalized responses.
      logical, intent(in), optional :: cumulative !! Cumulate responses over horizons.
      type(vars_irf_bootstrap_t) :: out
      type(vars_fit_t) :: bootstrap_model
      real(dp), allocatable :: draws(:, :, :, :), response(:, :, :)
      real(dp) :: selected_level
      integer :: run, successful
      logical :: use_orthogonal, use_cumulative

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      use_orthogonal = .true.
      if (present(orthogonal)) use_orthogonal = orthogonal
      use_cumulative = .false.
      if (present(cumulative)) use_cumulative = cumulative
      if (model%info /= 0 .or. horizon < 0 .or. runs < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      if (use_orthogonal) then
         out%point = vars_psi(model, horizon)
      else
         out%point = vars_phi(model, horizon)
      end if
      if (use_cumulative) call cumulative_responses(out%point)
      allocate(draws(size(out%point, 1), size(out%point, 2), size(out%point, 3), runs))
      successful = 0
      do run = 1, runs
         bootstrap_model = residual_bootstrap_fit(model)
         if (bootstrap_model%info /= 0) cycle
         if (use_orthogonal) then
            response = vars_psi(bootstrap_model, horizon)
         else
            response = vars_phi(bootstrap_model, horizon)
         end if
         if (use_cumulative) call cumulative_responses(response)
         successful = successful + 1
         draws(:, :, :, successful) = response
      end do
      if (successful < 1) then
         out%info = 2
         return
      end if
      call interval_bounds_4d(draws(:, :, :, :successful), selected_level, &
         out%lower, out%upper)
      out%level = selected_level
      out%successful_runs = successful
   end function vars_irf_bootstrap

   function vars_svar_bootstrap(model, a_template, b_template, estimate_a, estimate_b, &
      runs, level, max_iterations, tolerance) result(out)
      !! Compute residual-bootstrap confidence limits for structural impact matrices.
      type(vars_fit_t), intent(in) :: model !! Fitted reduced-form VAR.
      real(dp), intent(in) :: a_template(:, :) !! Fixed A entries and starting values.
      real(dp), intent(in) :: b_template(:, :) !! Fixed B entries and starting values.
      logical, intent(in) :: estimate_a(:, :) !! Free A mask.
      logical, intent(in) :: estimate_b(:, :) !! Free B mask.
      integer, intent(in) :: runs !! Number of bootstrap replications.
      real(dp), intent(in), optional :: level !! Central confidence level.
      integer, intent(in), optional :: max_iterations !! Structural optimizer iteration limit.
      real(dp), intent(in), optional :: tolerance !! Structural optimizer tolerance.
      type(vars_structural_bootstrap_t) :: out
      type(vars_fit_t) :: bootstrap_model
      type(vars_svar_t) :: point_fit, bootstrap_fit
      real(dp), allocatable :: draws(:, :, :)
      real(dp) :: selected_level, selected_tolerance
      integer :: run, successful, limit

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      selected_tolerance = 1.0e-7_dp
      if (present(tolerance)) selected_tolerance = tolerance
      limit = 1000
      if (present(max_iterations)) limit = max_iterations
      if (runs < 1 .or. selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      point_fit = vars_svar(model, a_template, b_template, estimate_a, estimate_b, &
         max_iterations=limit, tolerance=selected_tolerance)
      if (point_fit%info /= 0) then
         out%info = 2
         return
      end if
      out%point = point_fit%impact
      allocate(draws(size(out%point, 1), size(out%point, 2), runs))
      successful = 0
      do run = 1, runs
         bootstrap_model = residual_bootstrap_fit(model)
         if (bootstrap_model%info /= 0) cycle
         bootstrap_fit = vars_svar(bootstrap_model, a_template, b_template, &
            estimate_a, estimate_b, max_iterations=limit, tolerance=selected_tolerance)
         if (bootstrap_fit%info /= 0) cycle
         successful = successful + 1
         draws(:, :, successful) = bootstrap_fit%impact
      end do
      if (successful < 1) then
         out%info = 3
         return
      end if
      call interval_bounds_3d(draws(:, :, :successful), selected_level, &
         out%lower, out%upper)
      out%level = selected_level
      out%successful_runs = successful
   end function vars_svar_bootstrap

   function vars_svec_bootstrap(series, johansen, rank, short_run_zero, &
      long_run_zero, runs, season, exogenous, max_iterations, tolerance) result(out)
      !! Estimate SVEC impact-matrix standard errors by residual bootstrap refitting.
      real(dp), intent(in) :: series(:, :) !! Original time-ordered level observations.
      type(johansen_result_t), intent(in) :: johansen !! Johansen fit for the original series.
      integer, intent(in) :: rank !! Cointegration rank.
      logical, intent(in) :: short_run_zero(:, :) !! Zero restrictions on short-run impacts.
      logical, intent(in) :: long_run_zero(:, :) !! Zero restrictions on long-run impacts.
      integer, intent(in) :: runs !! Number of bootstrap replications.
      integer, intent(in), optional :: season !! Seasonal period used in the Johansen fit.
      real(dp), intent(in), optional :: exogenous(:, :) !! Exogenous regressors used in the fit.
      integer, intent(in), optional :: max_iterations !! Structural optimizer iteration limit.
      real(dp), intent(in), optional :: tolerance !! Structural optimizer tolerance.
      type(vars_svec_bootstrap_t) :: out
      type(vars_svec_t) :: point_fit, bootstrap_fit
      type(vars_vec2var_t) :: level_var
      type(johansen_result_t) :: bootstrap_johansen
      real(dp), allocatable :: simulated(:, :), short_draws(:, :, :), long_draws(:, :, :)
      real(dp) :: selected_tolerance
      integer :: limit, run, successful, shock

      limit = 1000
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-7_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (johansen%info /= 0 .or. runs < 2 .or. size(series, 2) /= johansen%variables .or. &
         size(series, 1) <= johansen%lag + 2) then
         out%info = 1
         return
      end if
      point_fit = vars_svec(johansen, rank, short_run_zero, long_run_zero, &
         limit, selected_tolerance)
      level_var = vars_vec2var(johansen, rank)
      if (point_fit%info /= 0 .or. level_var%info /= 0) then
         out%info = 2
         return
      end if
      out%short_run = point_fit%short_run
      out%long_run = point_fit%long_run
      allocate(short_draws(johansen%variables, johansen%variables, runs))
      allocate(long_draws(johansen%variables, johansen%variables, runs))
      successful = 0
      do run = 1, runs
         simulated = simulate_level_var_bootstrap(series, level_var%ar)
         bootstrap_johansen = johansen_test(simulated, johansen%test_type, &
            johansen%deterministic, johansen%lag, johansen%specification, &
            season=season, exogenous=exogenous)
         if (bootstrap_johansen%info /= 0) cycle
         bootstrap_fit = vars_svec(bootstrap_johansen, rank, short_run_zero, &
            long_run_zero, limit, selected_tolerance)
         if (bootstrap_fit%info /= 0) cycle
         do shock = 1, johansen%variables
            if (dot_product(bootstrap_fit%short_run(:, shock), &
               point_fit%short_run(:, shock)) < 0.0_dp) then
               bootstrap_fit%short_run(:, shock) = -bootstrap_fit%short_run(:, shock)
               bootstrap_fit%long_run(:, shock) = -bootstrap_fit%long_run(:, shock)
            end if
         end do
         successful = successful + 1
         short_draws(:, :, successful) = bootstrap_fit%short_run
         long_draws(:, :, successful) = bootstrap_fit%long_run
      end do
      if (successful < 2) then
         out%info = 3
         return
      end if
      out%short_run_standard_error = bootstrap_standard_error( &
         short_draws(:, :, :successful))
      out%long_run_standard_error = bootstrap_standard_error( &
         long_draws(:, :, :successful))
      out%successful_runs = successful
   end function vars_svec_bootstrap

   function simulate_level_var_bootstrap(series, ar) result(simulated)
      !! Simulate a level VAR by resampling its centered empirical innovations.
      real(dp), intent(in) :: series(:, :) !! Original level observations.
      real(dp), intent(in) :: ar(:, :, :) !! Level-VAR autoregressive coefficients.
      real(dp), allocatable :: simulated(:, :)
      real(dp), allocatable :: residuals(:, :), centered(:, :), prediction(:), mean_residual(:)
      integer :: n, k, p, rows, t, lag, index

      n = size(series, 1)
      k = size(series, 2)
      p = size(ar, 3)
      rows = n - p
      allocate(residuals(rows, k), centered(rows, k), mean_residual(k))
      do t = p + 1, n
         residuals(t - p, :) = series(t, :)
         do lag = 1, p
            residuals(t - p, :) = residuals(t - p, :) - &
               matmul(ar(:, :, lag), series(t - lag, :))
         end do
      end do
      mean_residual = sum(residuals, dim=1)/real(rows, dp)
      centered = residuals - spread(mean_residual, 1, rows)
      allocate(simulated(n, k), prediction(k))
      simulated(:p, :) = series(:p, :)
      do t = p + 1, n
         prediction = mean_residual
         do lag = 1, p
            prediction = prediction + matmul(ar(:, :, lag), simulated(t - lag, :))
         end do
         index = min(rows, 1 + int(random_uniform()*real(rows, dp)))
         simulated(t, :) = prediction + centered(index, :)
      end do
   end function simulate_level_var_bootstrap

   pure function bootstrap_standard_error(draws) result(standard_error)
      !! Compute elementwise sample standard deviations across bootstrap draws.
      real(dp), intent(in) :: draws(:, :, :) !! Bootstrap matrices indexed by replication.
      real(dp), allocatable :: standard_error(:, :)
      real(dp), allocatable :: mean_draw(:, :)
      integer :: runs

      runs = size(draws, 3)
      mean_draw = sum(draws, dim=3)/real(runs, dp)
      standard_error = sqrt(sum((draws - spread(mean_draw, 3, runs))**2, dim=3)/ &
         real(runs - 1, dp))
   end function bootstrap_standard_error

   pure function vars_granger_test(model, causes) result(out)
      !! Compute the vars system Wald F test for Granger non-causality.
      type(vars_fit_t), intent(in) :: model !! Fitted VAR model.
      logical, intent(in) :: causes(:) !! Variables whose lag coefficients are tested.
      type(vars_test_t) :: out
      real(dp), allocatable :: gram_inverse(:, :), selected(:), covariance(:, :), inverse(:, :)
      integer, allocatable :: equation(:), predictor(:)
      integer :: k, n, restriction_count, response_index, cause_index, lag, item, other, status

      k = size(model%response, 2)
      n = size(model%response, 1)
      if (model%info /= 0 .or. size(causes) /= k .or. count(causes) < 1 .or. &
         count(.not. causes) < 1) then
         out%info = 1
         return
      end if
      restriction_count = count(causes)*count(.not. causes)*model%order
      allocate(selected(restriction_count), covariance(restriction_count, restriction_count))
      allocate(equation(restriction_count), predictor(restriction_count))
      item = 0
      do response_index = 1, k
         if (causes(response_index)) cycle
         do lag = 1, model%order
            do cause_index = 1, k
               if (.not. causes(cause_index)) cycle
               item = item + 1
               equation(item) = response_index
               predictor(item) = (lag - 1)*k + cause_index
               selected(item) = model%coefficients(response_index, predictor(item))
            end do
         end do
      end do
      call invert_matrix(matmul(transpose(model%design), model%design), gram_inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      do item = 1, restriction_count
         do other = 1, restriction_count
            covariance(item, other) = model%covariance(equation(item), equation(other))* &
               gram_inverse(predictor(item), predictor(other))
         end do
      end do
      call invert_matrix(covariance, inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      out%statistic = dot_product(selected, matmul(inverse, selected))/real(restriction_count, dp)
      out%df1 = real(restriction_count, dp)
      out%df2 = real(k*(n - size(model%design, 2)), dp)
      out%p_value = f_upper_tail(out%statistic, out%df1, out%df2)
   end function vars_granger_test

   function vars_granger_bootstrap(model, causes, runs) result(out)
      !! Estimate a Granger non-causality p-value by restricted residual bootstrap.
      type(vars_fit_t), intent(in) :: model !! Fitted unrestricted VAR.
      logical, intent(in) :: causes(:) !! Variables whose lag coefficients are tested.
      integer, intent(in) :: runs !! Number of bootstrap replications.
      type(vars_test_t) :: out
      type(vars_fit_t) :: restricted_model, bootstrap_model
      type(vars_test_t) :: bootstrap_test
      logical, allocatable :: mask(:, :)
      integer :: k, response_index, cause_index, lag, run, exceedances, successful

      out = vars_granger_test(model, causes)
      k = size(model%response, 2)
      if (out%info /= 0 .or. runs < 1) then
         out%info = 1
         return
      end if
      mask = model%restrictions
      do response_index = 1, k
         if (causes(response_index)) cycle
         do lag = 1, model%order
            do cause_index = 1, k
               if (causes(cause_index)) mask(response_index, (lag - 1)*k + cause_index) = .false.
            end do
         end do
      end do
      restricted_model = vars_restrict(model, mask)
      if (restricted_model%info /= 0) then
         out%info = 2
         return
      end if
      exceedances = 0
      successful = 0
      do run = 1, runs
         bootstrap_model = residual_bootstrap_fit(restricted_model)
         if (bootstrap_model%info /= 0) cycle
         bootstrap_test = vars_granger_test(bootstrap_model, causes)
         if (bootstrap_test%info /= 0) cycle
         successful = successful + 1
         if (bootstrap_test%statistic > out%statistic) exceedances = exceedances + 1
      end do
      if (successful < 1) then
         out%info = 3
         return
      end if
      out%p_value = real(exceedances, dp)/real(successful, dp)
      out%df2 = real(successful, dp)
   end function vars_granger_bootstrap

   function residual_bootstrap_fit(model) result(out)
      !! Simulate and refit a VAR using centered residual resampling.
      type(vars_fit_t), intent(in) :: model !! Source fitted VAR.
      type(vars_fit_t) :: out
      real(dp), allocatable :: simulated(:, :), centered(:, :), design(:, :), response(:, :)
      real(dp), allocatable :: prediction(:)
      integer :: n, k, p, t, lag, index, deterministic_start

      n = size(model%response, 1)
      k = size(model%response, 2)
      p = model%order
      allocate(simulated(n + p, k), centered(n, k), prediction(k))
      centered = model%residuals - spread(sum(model%residuals, dim=1)/real(n, dp), 1, n)
      do lag = 1, p
         simulated(p + 1 - lag, :) = model%design(1, (lag - 1)*k + 1:lag*k)
      end do
      deterministic_start = k*p + 1
      do t = 1, n
         prediction = 0.0_dp
         do lag = 1, p
            prediction = prediction + matmul(model%ar(:, :, lag), simulated(p + t - lag, :))
         end do
         if (deterministic_start <= size(model%design, 2)) prediction = prediction + &
            matmul(model%deterministic, model%design(t, deterministic_start:))
         index = min(n, 1 + int(random_uniform()*real(n, dp)))
         simulated(p + t, :) = prediction + centered(index, :)
      end do
      allocate(design(n, size(model%design, 2)), response(n, k))
      response = simulated(p + 1:, :)
      do t = 1, n
         do lag = 1, p
            design(t, (lag - 1)*k + 1:lag*k) = simulated(p + t - lag, :)
         end do
      end do
      if (deterministic_start <= size(model%design, 2)) &
         design(:, deterministic_start:) = model%design(:, deterministic_start:)
      out = fit_design(design, response, p, model%deterministic_type, model%season, &
         model%exogenous_count)
   end function residual_bootstrap_fit

   pure subroutine cumulative_responses(response)
      !! Cumulate impulse responses over the horizon dimension.
      real(dp), intent(inout) :: response(:, :, :) !! Responses updated in place.
      integer :: horizon

      do horizon = 2, size(response, 3)
         response(:, :, horizon) = response(:, :, horizon) + response(:, :, horizon - 1)
      end do
   end subroutine cumulative_responses

   pure subroutine interval_bounds_4d(draws, level, lower, upper)
      !! Compute elementwise central bounds from four-dimensional bootstrap draws.
      real(dp), intent(in) :: draws(:, :, :, :) !! Bootstrap draws in the last dimension.
      real(dp), intent(in) :: level !! Central interval level.
      real(dp), allocatable, intent(out) :: lower(:, :, :) !! Lower bounds.
      real(dp), allocatable, intent(out) :: upper(:, :, :) !! Upper bounds.
      real(dp), allocatable :: values(:)
      real(dp) :: tail
      integer :: i, j, h

      allocate(lower(size(draws, 1), size(draws, 2), size(draws, 3)))
      allocate(upper(size(draws, 1), size(draws, 2), size(draws, 3)))
      tail = 0.5_dp*(1.0_dp - level)
      do h = 1, size(draws, 3)
         do j = 1, size(draws, 2)
            do i = 1, size(draws, 1)
               values = draws(i, j, h, :)
               call sort(values)
               lower(i, j, h) = quantile(values, tail)
               upper(i, j, h) = quantile(values, 1.0_dp - tail)
            end do
         end do
      end do
   end subroutine interval_bounds_4d

   pure subroutine interval_bounds_3d(draws, level, lower, upper)
      !! Compute elementwise central bounds from three-dimensional bootstrap draws.
      real(dp), intent(in) :: draws(:, :, :) !! Bootstrap draws in the last dimension.
      real(dp), intent(in) :: level !! Central interval level.
      real(dp), allocatable, intent(out) :: lower(:, :) !! Lower bounds.
      real(dp), allocatable, intent(out) :: upper(:, :) !! Upper bounds.
      real(dp), allocatable :: values(:)
      real(dp) :: tail
      integer :: i, j

      allocate(lower(size(draws, 1), size(draws, 2)))
      allocate(upper(size(draws, 1), size(draws, 2)))
      tail = 0.5_dp*(1.0_dp - level)
      do j = 1, size(draws, 2)
         do i = 1, size(draws, 1)
            values = draws(i, j, :)
            call sort(values)
            lower(i, j) = quantile(values, tail)
            upper(i, j) = quantile(values, 1.0_dp - tail)
         end do
      end do
   end subroutine interval_bounds_3d

   pure function fit_design(design, response, order, deterministic_type, season, &
      exogenous_count) result(out)
      !! Fit all VAR equations to an already aligned design matrix.
      real(dp), intent(in) :: design(:, :) !! Aligned predictor matrix.
      real(dp), intent(in) :: response(:, :) !! Aligned response matrix.
      integer, intent(in) :: order !! Autoregressive order.
      integer, intent(in) :: deterministic_type !! Deterministic selector.
      integer, intent(in) :: season !! Seasonal period.
      integer, intent(in) :: exogenous_count !! Exogenous column count.
      type(vars_fit_t) :: out
      real(dp), allocatable :: beta(:), se(:), residuals(:)
      real(dp) :: rss
      integer :: equation, status, k, deterministic_columns

      k = size(response, 2)
      allocate(out%coefficients(k, size(design, 2)))
      do equation = 1, k
         call ols_fit(design, response(:, equation), beta, se, residuals, rss, status)
         if (status /= 0) then
            out%info = status
            return
         end if
         out%coefficients(equation, :) = beta
      end do
      out%design = design
      out%response = response
      out%order = order
      out%deterministic_type = deterministic_type
      out%season = season
      out%exogenous_count = exogenous_count
      out%ar = reshape(out%coefficients(:, :k*order), [k, k, order])
      deterministic_columns = size(design, 2) - k*order
      allocate(out%deterministic(k, deterministic_columns))
      if (deterministic_columns > 0) out%deterministic = out%coefficients(:, k*order + 1:)
      allocate(out%restrictions(k, size(design, 2)))
      out%restrictions = .true.
      call update_fit_moments(out)
   end function fit_design

   pure subroutine update_fit_moments(model)
      !! Recompute fitted values, residuals, and covariance after coefficient changes.
      type(vars_fit_t), intent(inout) :: model !! Model updated in place.
      integer :: n, k

      n = size(model%response, 1)
      k = size(model%response, 2)
      model%fitted = matmul(model%design, transpose(model%coefficients))
      model%residuals = model%response - model%fitted
      model%covariance = matmul(transpose(model%residuals), model%residuals)/ &
         real(max(1, n - size(model%design, 2)), dp)
      model%ar = reshape(model%coefficients(:, :k*model%order), [k, k, model%order])
      if (size(model%coefficients, 2) > k*model%order) &
         model%deterministic = model%coefficients(:, k*model%order + 1:)
   end subroutine update_fit_moments

   pure subroutine set_deterministic_row(row, offset, deterministic_type, time_index, &
      season, exogenous_count, exogenous_values)
      !! Fill deterministic, centered seasonal, and exogenous columns in one design row.
      real(dp), intent(inout) :: row(:) !! Design row updated in place.
      integer, intent(inout) :: offset !! Current column offset, updated in place.
      integer, intent(in) :: deterministic_type !! Deterministic selector.
      integer, intent(in) :: time_index !! One-based observation index.
      integer, intent(in) :: season !! Seasonal period.
      integer, intent(in) :: exogenous_count !! Exogenous column count.
      real(dp), intent(in) :: exogenous_values(:) !! Current exogenous row.
      integer :: category, j

      if (deterministic_type == vars_deterministic_constant .or. &
         deterministic_type == vars_deterministic_both) then
         offset = offset + 1
         row(offset) = 1.0_dp
      end if
      if (deterministic_type == vars_deterministic_trend .or. &
         deterministic_type == vars_deterministic_both) then
         offset = offset + 1
         row(offset) = real(time_index, dp)
      end if
      if (season > 1) then
         category = modulo(time_index - 1, season) + 1
         do j = 1, season - 1
            offset = offset + 1
            row(offset) = -1.0_dp/real(season, dp)
            if (category == j) row(offset) = row(offset) + 1.0_dp
         end do
      end if
      if (exogenous_count > 0) then
         row(offset + 1:offset + exogenous_count) = exogenous_values
         offset = offset + exogenous_count
      end if
   end subroutine set_deterministic_row

   pure function exogenous_row(exogenous, row, count_values) result(values)
      !! Return an exogenous row or an empty vector when exogenous data are absent.
      real(dp), intent(in), optional :: exogenous(:, :) !! Exogenous observations.
      integer, intent(in) :: row !! Requested row.
      integer, intent(in) :: count_values !! Expected value count.
      real(dp), allocatable :: values(:)

      allocate(values(count_values))
      if (count_values > 0 .and. present(exogenous)) values = exogenous(row, :)
   end function exogenous_row

   pure integer function deterministic_column_count(deterministic_type) result(columns)
      !! Return the number of constant and trend columns for a selector.
      integer, intent(in) :: deterministic_type !! Deterministic selector.

      select case (deterministic_type)
      case (vars_deterministic_none)
         columns = 0
      case (vars_deterministic_constant, vars_deterministic_trend)
         columns = 1
      case (vars_deterministic_both)
         columns = 2
      case default
         columns = -1
      end select
   end function deterministic_column_count

   pure subroutine covariance_on_common_sample(model, drop, determinant_value, info)
      !! Return the residual-covariance determinant after dropping early aligned rows.
      type(vars_fit_t), intent(in) :: model !! Candidate model.
      integer, intent(in) :: drop !! Number of early rows to omit.
      real(dp), intent(out) :: determinant_value !! Covariance determinant.
      integer, intent(out) :: info !! Status code.
      real(dp), allocatable :: covariance(:, :)
      integer :: n

      n = size(model%residuals, 1) - drop
      if (n < 1) then
         determinant_value = 0.0_dp
         info = 1
         return
      end if
      covariance = matmul(transpose(model%residuals(drop + 1:, :)), &
         model%residuals(drop + 1:, :))/real(n, dp)
      determinant_value = matrix_determinant(covariance)
      if (determinant_value > 0.0_dp) then
         info = 0
      else
         determinant_value = 0.0_dp
         info = 2
      end if
   end subroutine covariance_on_common_sample

   pure subroutine restricted_equation(design, response, mask, coefficients, info)
      !! Fit one equation to selected design columns and expand its coefficient vector.
      real(dp), intent(in) :: design(:, :) !! Full design matrix.
      real(dp), intent(in) :: response(:) !! Equation response.
      logical, intent(in) :: mask(:) !! Included design columns.
      real(dp), intent(out) :: coefficients(:) !! Expanded coefficient vector.
      integer, intent(out) :: info !! Status code.
      real(dp), allocatable :: selected(:, :), beta(:), se(:), residuals(:)
      real(dp) :: rss
      integer :: j, packed

      selected = pack_design(design, mask)
      call ols_fit(selected, response, beta, se, residuals, rss, info)
      if (info /= 0) return
      coefficients = 0.0_dp
      packed = 0
      do j = 1, size(mask)
         if (.not. mask(j)) cycle
         packed = packed + 1
         coefficients(j) = beta(packed)
      end do
   end subroutine restricted_equation

   pure function pack_design(design, mask) result(selected)
      !! Pack selected columns from a design matrix.
      real(dp), intent(in) :: design(:, :) !! Full design matrix.
      logical, intent(in) :: mask(:) !! Selected columns.
      real(dp), allocatable :: selected(:, :)
      integer :: j, column

      allocate(selected(size(design, 1), count(mask)))
      column = 0
      do j = 1, size(mask)
         if (.not. mask(j)) cycle
         column = column + 1
         selected(:, column) = design(:, j)
      end do
   end function pack_design

   pure integer function packed_index(mask, packed_position) result(index_value)
      !! Map a packed true-element position back to its original logical index.
      logical, intent(in) :: mask(:) !! Logical selection mask.
      integer, intent(in) :: packed_position !! One-based packed position.
      integer :: j, count_true

      index_value = 0
      count_true = 0
      do j = 1, size(mask)
         if (.not. mask(j)) cycle
         count_true = count_true + 1
         if (count_true == packed_position) then
            index_value = j
            return
         end if
      end do
   end function packed_index

   pure function identity(n) result(matrix)
      !! Return a real identity matrix.
      integer, intent(in) :: n !! Matrix order.
      real(dp) :: matrix(n, n)
      integer :: i

      matrix = 0.0_dp
      do i = 1, n
         matrix(i, i) = 1.0_dp
      end do
   end function identity

   pure function commutation_matrix(n) result(matrix)
      !! Return the matrix mapping column-major vec(A) to vec(transpose(A)).
      integer, intent(in) :: n !! Square matrix order.
      real(dp) :: matrix(n*n, n*n)
      integer :: i, j

      matrix = 0.0_dp
      do j = 1, n
         do i = 1, n
            matrix(j + (i - 1)*n, i + (j - 1)*n) = 1.0_dp
         end do
      end do
   end function commutation_matrix

   pure real(dp) function matrix_trace(matrix) result(value)
      !! Return the trace of a square matrix.
      real(dp), intent(in) :: matrix(:, :) !! Square matrix.
      integer :: i

      value = 0.0_dp
      do i = 1, min(size(matrix, 1), size(matrix, 2))
         value = value + matrix(i, i)
      end do
   end function matrix_trace

   pure real(dp) function matrix_determinant(matrix) result(value)
      !! Return a square matrix determinant by pivoted elimination.
      real(dp), intent(in) :: matrix(:, :) !! Square matrix.
      real(dp), allocatable :: work(:, :), row_values(:)
      real(dp) :: pivot
      integer :: i, pivot_row, n

      n = size(matrix, 1)
      if (size(matrix, 2) /= n) then
         value = 0.0_dp
         return
      end if
      work = matrix
      value = 1.0_dp
      do i = 1, n
         pivot_row = i - 1 + maxloc(abs(work(i:, i)), dim=1)
         if (abs(work(pivot_row, i)) <= tiny(1.0_dp)) then
            value = 0.0_dp
            return
         end if
         if (pivot_row /= i) then
            row_values = work(i, :)
            work(i, :) = work(pivot_row, :)
            work(pivot_row, :) = row_values
            value = -value
         end if
         pivot = work(i, i)
         value = value*pivot
         if (i < n) work(i + 1:, i + 1:) = work(i + 1:, i + 1:) - &
            spread(work(i + 1:, i)/pivot, 2, n - i)*spread(work(i, i + 1:), 1, n - i)
      end do
   end function matrix_determinant

   pure function chi_square_test(statistic, degrees) result(out)
      !! Construct a chi-square upper-tail test result.
      real(dp), intent(in) :: statistic !! Test statistic.
      real(dp), intent(in) :: degrees !! Reference degrees of freedom.
      type(vars_test_t) :: out

      out%statistic = max(0.0_dp, statistic)
      out%df1 = degrees
      if (degrees > 0.0_dp) then
         out%p_value = regularized_gamma_q(0.5_dp*degrees, 0.5_dp*out%statistic)
      else
         out%info = 1
      end if
   end function chi_square_test

   pure real(dp) function f_upper_tail(statistic, df1, df2) result(probability)
      !! Return the upper-tail probability of an F statistic.
      real(dp), intent(in) :: statistic !! F statistic.
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

   pure function lagged_residual_design(residuals, lags) result(design)
      !! Return zero-padded lagged residual regressors matching the input rows.
      real(dp), intent(in) :: residuals(:, :) !! Residual matrix.
      integer, intent(in) :: lags !! Lag count.
      real(dp), allocatable :: design(:, :)
      integer :: n, k, lag

      n = size(residuals, 1)
      k = size(residuals, 2)
      allocate(design(n, k*lags))
      design = 0.0_dp
      do lag = 1, lags
         design(lag + 1:, (lag - 1)*k + 1:lag*k) = residuals(:n - lag, :)
      end do
   end function lagged_residual_design

   pure subroutine residual_covariance(design, response, covariance, info)
      !! Fit a multivariate regression and return its maximum-likelihood covariance.
      real(dp), intent(in) :: design(:, :) !! Predictor matrix.
      real(dp), intent(in) :: response(:, :) !! Multivariate response.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Residual covariance.
      integer, intent(out) :: info !! Status code.
      real(dp), allocatable :: beta(:), se(:), residuals(:), all_residuals(:, :)
      real(dp) :: rss
      integer :: equation

      allocate(all_residuals(size(response, 1), size(response, 2)))
      do equation = 1, size(response, 2)
         call ols_fit(design, response(:, equation), beta, se, residuals, rss, info)
         if (info /= 0) return
         all_residuals(:, equation) = residuals
      end do
      covariance = matmul(transpose(all_residuals), all_residuals)/real(size(response, 1), dp)
   end subroutine residual_covariance

   pure function lagged_with_intercept(values, lags) result(design)
      !! Return an intercept and lagged copies aligned after the requested lag count.
      real(dp), intent(in) :: values(:, :) !! Input multivariate observations.
      integer, intent(in) :: lags !! Lag count.
      real(dp), allocatable :: design(:, :)
      integer :: n, q, lag

      n = size(values, 1)
      q = size(values, 2)
      allocate(design(n - lags, 1 + q*lags))
      design(:, 1) = 1.0_dp
      do lag = 1, lags
         design(:, 2 + (lag - 1)*q:1 + lag*q) = values(lags + 1 - lag:n - lag, :)
      end do
   end function lagged_with_intercept

   pure function covariance_matrix_centered(values) result(covariance)
      !! Return covariance about column sample means.
      real(dp), intent(in) :: values(:, :) !! Multivariate observations.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: centered(:, :)

      centered = values - spread(sum(values, dim=1)/real(size(values, 1), dp), 1, size(values, 1))
      covariance = matmul(transpose(centered), centered)/real(max(1, size(values, 1) - 1), dp)
   end function covariance_matrix_centered

   pure subroutine scalar_regression_rss(design, response, rss, info)
      !! Fit a scalar regression and return its residual sum of squares.
      real(dp), intent(in) :: design(:, :) !! Predictor matrix.
      real(dp), intent(in) :: response(:) !! Response vector.
      real(dp), intent(out) :: rss !! Residual sum of squares.
      integer, intent(out) :: info !! Status code.
      real(dp), allocatable :: beta(:), se(:), residuals(:)

      call ols_fit(design, response, beta, se, residuals, rss, info)
   end subroutine scalar_regression_rss
end module vars_mod
