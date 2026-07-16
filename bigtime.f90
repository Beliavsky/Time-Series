! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R bigtime package.
! Sparse large-VAR algorithms translated from the R bigtime package.
module bigtime_mod
   use kind_mod, only: dp
   use time_series_linalg_mod, only: symmetric_eigen, inverse_logdet
   use time_series_linalg_mod, only: general_eigenvalues
   use time_series_random_mod, only: random_uniform
   use time_series_random_mod, only: random_standard_normal_matrix
   implicit none
   private

   integer, parameter, public :: bigtime_penalty_l1 = 1
   integer, parameter, public :: bigtime_penalty_hlag = 2
   integer, parameter, public :: bigtime_ic_aic = 1
   integer, parameter, public :: bigtime_ic_bic = 2
   integer, parameter, public :: bigtime_ic_hq = 3
   integer, parameter, public :: bigtime_sparsity_dense = 0
   integer, parameter, public :: bigtime_sparsity_l1 = 1
   integer, parameter, public :: bigtime_sparsity_hlag = 2

   type, public :: bigtime_var_fit_t
      ! Sparse VAR coefficient estimates and in-sample diagnostics.
      real(dp), allocatable :: phi(:, :)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: lag_order = 0
      integer :: iterations = 0
      integer :: nonzero = 0
      integer :: info = 0
      logical :: converged = .false.
   end type bigtime_var_fit_t

   type, public :: bigtime_var_path_t
      ! Warm-started sparse VAR estimates over a regularization path.
      real(dp), allocatable :: phi(:, :, :)
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: objective(:)
      integer, allocatable :: iterations(:)
      integer, allocatable :: nonzero(:)
      logical, allocatable :: converged(:)
      real(dp), allocatable :: lambda(:)
      integer :: lag_order = 0
      integer :: penalty = 0
      integer :: info = 0
   end type bigtime_var_path_t

   type, public :: bigtime_varx_fit_t
      ! Sparse VARX coefficient estimates and in-sample diagnostics.
      real(dp), allocatable :: phi(:, :)
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp) :: lambda_phi = 0.0_dp
      real(dp) :: lambda_beta = 0.0_dp
      real(dp) :: alpha = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: ar_order = 0
      integer :: exogenous_order = 0
      integer :: iterations = 0
      integer :: nonzero = 0
      integer :: info = 0
      logical :: converged = .false.
   end type bigtime_varx_fit_t

   type, public :: bigtime_varx_path_t
      ! Sparse VARX estimates over paired endogenous and exogenous penalties.
      real(dp), allocatable :: phi(:, :, :)
      real(dp), allocatable :: beta(:, :, :)
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: objective(:)
      integer, allocatable :: iterations(:)
      integer, allocatable :: nonzero(:)
      logical, allocatable :: converged(:)
      real(dp), allocatable :: lambda_phi(:)
      real(dp), allocatable :: lambda_beta(:)
      real(dp) :: alpha = 0.0_dp
      integer :: ar_order = 0
      integer :: exogenous_order = 0
      integer :: penalty = 0
      integer :: info = 0
   end type bigtime_varx_path_t

   type, public :: bigtime_varx_grid_t
      ! Separate descending regularization grids for VARX coefficient blocks.
      real(dp), allocatable :: lambda_phi(:)
      real(dp), allocatable :: lambda_beta(:)
      integer :: info = 0
   end type bigtime_varx_grid_t

   type, public :: bigtime_varma_fit_t
      ! Two-stage sparse VARMA estimates and innovation-proxy diagnostics.
      type(bigtime_var_fit_t) :: phase1
      type(bigtime_varx_fit_t) :: phase2
      real(dp), allocatable :: phi(:, :)
      real(dp), allocatable :: theta(:, :)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: innovations(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      integer :: phase1_order = 0
      integer :: ar_order = 0
      integer :: ma_order = 0
      integer :: info = 0
      logical :: estimated_innovations = .false.
      logical :: converged = .false.
   end type bigtime_varma_fit_t

   type, public :: bigtime_varma_path_t
      ! Two-stage sparse VARMA path sharing one Phase I innovation proxy.
      type(bigtime_var_fit_t) :: phase1
      type(bigtime_varx_path_t) :: phase2
      real(dp), allocatable :: innovations(:, :)
      integer :: phase1_order = 0
      integer :: ar_order = 0
      integer :: ma_order = 0
      integer :: info = 0
      logical :: estimated_innovations = .false.
   end type bigtime_varma_path_t

   type, public :: bigtime_cv_result_t
      ! Expanding-window forecast errors and selected regularization indices.
      real(dp), allocatable :: squared_error(:, :)
      real(dp), allocatable :: mean_squared_error(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: mean_nonzero(:)
      integer :: best_index = 0
      integer :: one_se_index = 0
      integer :: info = 0
   end type bigtime_cv_result_t

   type, public :: bigtime_varma_cv_t
      ! Phase I innovation proxy and Phase II expanding-window selection.
      type(bigtime_var_fit_t) :: phase1
      type(bigtime_cv_result_t) :: phase2
      real(dp), allocatable :: innovations(:, :)
      integer :: info = 0
      logical :: estimated_innovations = .false.
   end type bigtime_varma_cv_t

   type, public :: bigtime_ic_result_t
      ! Information criteria and minimizing path indices.
      real(dp), allocatable :: value(:, :)
      integer :: selected(3) = 0
      integer :: info = 0
   end type bigtime_ic_result_t

   type, public :: bigtime_forecast_t
      ! Recursive point forecasts for one fitted sparse model.
      real(dp), allocatable :: mean(:, :)
      integer :: info = 0
   end type bigtime_forecast_t

   type, public :: bigtime_path_forecast_t
      ! Recursive point forecasts for every regularization-path slice.
      real(dp), allocatable :: mean(:, :, :)
      integer :: info = 0
   end type bigtime_path_forecast_t

   type, public :: bigtime_stability_t
      ! Companion roots and stability diagnostics for sparse AR coefficients.
      complex(dp), allocatable :: roots(:)
      real(dp) :: maximum_modulus = huge(1.0_dp)
      integer :: info = 0
      logical :: stable = .false.
   end type bigtime_stability_t

   type, public :: bigtime_lag_order_t
      ! Largest active lag for primary and optional secondary coefficient blocks.
      integer, allocatable :: primary(:, :)
      integer, allocatable :: secondary(:, :)
      integer :: info = 0
   end type bigtime_lag_order_t

   type, public :: bigtime_coefficient_t
      ! Sparse AR and companion coefficients with stability diagnostics.
      real(dp), allocatable :: phi(:, :)
      real(dp), allocatable :: companion(:, :)
      real(dp) :: maximum_modulus = huge(1.0_dp)
      integer :: scaling_iterations = 0
      integer :: info = 0
   end type bigtime_coefficient_t

   type, public :: bigtime_simulation_t
      ! Simulated VAR observations and the supplied recursion inputs.
      real(dp), allocatable :: series(:, :)
      real(dp), allocatable :: innovations(:, :)
      real(dp), allocatable :: initial_state(:)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: phi(:, :)
      integer :: burnin = 0
      integer :: info = 0
   end type bigtime_simulation_t

   public :: bigtime_soft_threshold
   public :: bigtime_hlag_prox
   public :: bigtime_sparse_var
   public :: bigtime_sparse_var_path
   public :: bigtime_var_lambda_grid
   public :: bigtime_sparse_varx
   public :: bigtime_sparse_varx_path
   public :: bigtime_varx_lambda_grid
   public :: bigtime_varma_innovation_proxy
   public :: bigtime_sparse_varma
   public :: bigtime_sparse_varma_path
   public :: bigtime_var_cv
   public :: bigtime_varx_cv
   public :: bigtime_varma_cv
   public :: bigtime_var_path_ic
   public :: bigtime_varx_path_ic
   public :: bigtime_varma_path_ic
   public :: bigtime_select_var_path
   public :: bigtime_select_varx_path
   public :: bigtime_var_forecast
   public :: bigtime_var_path_forecast
   public :: bigtime_varx_forecast
   public :: bigtime_varx_path_forecast
   public :: bigtime_varma_forecast
   public :: bigtime_var_stability
   public :: bigtime_active_lags
   public :: bigtime_var_coefficients_from_draws
   public :: bigtime_random_var_coefficients
   public :: bigtime_var_simulate_from_innovations
   public :: bigtime_var_simulate
   public :: bigtime_companion_matrix

contains

   pure elemental real(dp) function bigtime_soft_threshold(value, threshold) &
      result(out)
      ! Apply the scalar lasso proximal operator.
      real(dp), intent(in) :: value, threshold

      if (value > threshold) then
         out = value - threshold
      else if (value < -threshold) then
         out = value + threshold
      else
         out = 0.0_dp
      end if
   end function bigtime_soft_threshold

   pure function bigtime_hlag_prox(values, threshold, variables, lag_order) &
      result(out)
      ! Apply bigtime's elementwise hierarchical-lag proximal operator.
      real(dp), intent(in) :: values(:), threshold
      integer, intent(in) :: variables, lag_order
      real(dp) :: out(size(values))
      real(dp) :: group_norm, shrinkage
      integer :: predictor, first_lag, lag

      out = values
      if (threshold <= 0.0_dp) return
      if (variables < 1 .or. lag_order < 1) then
         out = 0.0_dp
         return
      end if
      if (size(values) /= variables*lag_order) then
         out = 0.0_dp
         return
      end if

      do predictor = 1, variables
         do first_lag = lag_order, 1, -1
            group_norm = 0.0_dp
            do lag = first_lag, lag_order
               group_norm = group_norm + &
                  out(predictor + (lag - 1)*variables)**2
            end do
            group_norm = sqrt(group_norm)
            if (group_norm < threshold*(1.0_dp + 1.0e-8_dp)) then
               do lag = first_lag, lag_order
                  out(predictor + (lag - 1)*variables) = 0.0_dp
               end do
            else
               shrinkage = 1.0_dp - threshold/group_norm
               do lag = first_lag, lag_order
                  out(predictor + (lag - 1)*variables) = shrinkage* &
                     out(predictor + (lag - 1)*variables)
               end do
            end if
         end do
      end do
   end function bigtime_hlag_prox

   pure function bigtime_sparse_var(series, lag_order, lambda, penalty, &
      tolerance, max_iterations, initial_phi) result(out)
      ! Estimate one sparse VAR by bigtime's row-wise accelerated proximal method.
      real(dp), intent(in) :: series(:, :), lambda
      integer, intent(in) :: lag_order, penalty
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: initial_phi(:, :)
      type(bigtime_var_fit_t) :: out
      real(dp), allocatable :: response(:, :), design(:, :)
      real(dp), allocatable :: centered_response(:, :), centered_design(:, :)
      real(dp), allocatable :: response_mean(:), design_mean(:)
      real(dp), allocatable :: gram(:, :), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: coefficient(:), row_response(:)
      real(dp) :: selected_tolerance, step
      integer :: selected_max_iterations, variables, observations, predictors
      integer :: equation, eigen_info, row_iterations
      logical :: row_converged

      out%lag_order = lag_order
      out%lambda = lambda
      observations = size(series, 1)
      variables = size(series, 2)
      if (variables < 1 .or. lag_order < 1 .or. &
         observations <= lag_order .or. lambda < 0.0_dp) then
         out%info = 1
         return
      end if
      if (penalty /= bigtime_penalty_l1 .and. &
         penalty /= bigtime_penalty_hlag) then
         out%info = 2
         return
      end if
      predictors = variables*lag_order
      if (present(initial_phi)) then
         if (any(shape(initial_phi) /= [variables, predictors])) then
            out%info = 3
            return
         end if
      end if

      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      if (selected_tolerance <= 0.0_dp .or. selected_max_iterations < 1) then
         out%info = 1
         return
      end if

      call build_var_data(series, lag_order, response, design)
      allocate(response_mean(variables), design_mean(predictors))
      response_mean = sum(response, dim=1)/real(size(response, 1), dp)
      design_mean = sum(design, dim=1)/real(size(design, 1), dp)
      allocate(centered_response(size(response, 1), variables))
      allocate(centered_design(size(design, 1), predictors))
      centered_response = response - &
         spread(response_mean, 1, size(response, 1))
      centered_design = design - spread(design_mean, 1, size(design, 1))

      allocate(gram(predictors, predictors))
      gram = matmul(transpose(centered_design), centered_design)
      call symmetric_eigen(gram, eigenvalues, eigenvectors, eigen_info)
      if (eigen_info /= 0) then
         out%info = 10 + eigen_info
         return
      end if

      allocate(out%phi(variables, predictors))
      out%phi = 0.0_dp
      if (present(initial_phi)) out%phi = initial_phi
      out%converged = .true.
      if (maxval(abs(eigenvalues)) <= tiny(1.0_dp)) then
         out%iterations = 0
      else
         step = 1.0_dp/maxval(eigenvalues)
         allocate(coefficient(predictors))
         allocate(row_response(size(centered_response, 1)))
         do equation = 1, variables
            row_response = centered_response(:, equation)
            call sparse_var_row(row_response, centered_design, &
               out%phi(equation, :), lambda, penalty, variables, lag_order, &
               step, selected_tolerance, selected_max_iterations, &
               coefficient, row_iterations, row_converged)
            out%phi(equation, :) = coefficient
            out%iterations = max(out%iterations, row_iterations)
            out%converged = out%converged .and. row_converged
         end do
      end if

      allocate(out%intercept(variables))
      out%intercept = response_mean - matmul(out%phi, design_mean)
      allocate(out%fitted(size(response, 1), variables))
      out%fitted = matmul(design, transpose(out%phi)) + &
         spread(out%intercept, 1, size(response, 1))
      allocate(out%residuals(size(response, 1), variables))
      out%residuals = response - out%fitted
      out%objective = 0.5_dp*sum(out%residuals**2) + &
         lambda*penalty_value(out%phi, penalty, variables, lag_order)
      out%nonzero = count(abs(out%phi) > &
         100.0_dp*epsilon(1.0_dp))
      if (out%converged) then
         out%info = 0
      else
         out%info = 4
      end if
   end function bigtime_sparse_var

   pure function bigtime_sparse_var_path(series, lag_order, lambdas, penalty, &
      tolerance, max_iterations, warm_start) result(out)
      ! Estimate a sequence of sparse VARs, optionally using adjacent warm starts.
      real(dp), intent(in) :: series(:, :), lambdas(:)
      integer, intent(in) :: lag_order, penalty
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      logical, intent(in), optional :: warm_start
      type(bigtime_var_path_t) :: out
      type(bigtime_var_fit_t) :: fit
      real(dp), allocatable :: starting_phi(:, :)
      real(dp) :: selected_tolerance
      integer :: selected_max_iterations, variables, predictors, index
      logical :: use_warm_start

      out%lag_order = lag_order
      variables = size(series, 2)
      predictors = variables*lag_order
      if (size(lambdas) < 1 .or. any(lambdas < 0.0_dp) .or. &
         variables < 1 .or. lag_order < 1) then
         out%info = 1
         return
      end if
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      use_warm_start = .true.
      if (present(warm_start)) use_warm_start = warm_start

      allocate(out%phi(variables, predictors, size(lambdas)))
      allocate(out%intercept(variables, size(lambdas)))
      allocate(out%objective(size(lambdas)))
      allocate(out%iterations(size(lambdas)))
      allocate(out%nonzero(size(lambdas)))
      allocate(out%converged(size(lambdas)))
      allocate(out%lambda(size(lambdas)))
      out%lambda = lambdas
      out%penalty = penalty
      allocate(starting_phi(variables, predictors))
      starting_phi = 0.0_dp
      do index = 1, size(lambdas)
         fit = bigtime_sparse_var(series, lag_order, lambdas(index), penalty, &
            selected_tolerance, selected_max_iterations, starting_phi)
         if (fit%info /= 0 .and. fit%info /= 4) then
            out%info = 100*index + fit%info
            return
         end if
         out%phi(:, :, index) = fit%phi
         out%intercept(:, index) = fit%intercept
         out%objective(index) = fit%objective
         out%iterations(index) = fit%iterations
         out%nonzero(index) = fit%nonzero
         out%converged(index) = fit%converged
         if (use_warm_start) then
            starting_phi = fit%phi
         else
            starting_phi = 0.0_dp
         end if
      end do
      if (all(out%converged)) then
         out%info = 0
      else
         out%info = 4
      end if
   end function bigtime_sparse_var_path

   pure function bigtime_var_lambda_grid(series, lag_order, penalty, &
      grid_ratio, grid_size) result(lambdas)
      ! Construct bigtime's descending geometric grid from an all-zero bound.
      real(dp), intent(in) :: series(:, :), grid_ratio
      integer, intent(in) :: lag_order, penalty, grid_size
      real(dp), allocatable :: lambdas(:)
      real(dp), allocatable :: response(:, :), design(:, :)
      real(dp), allocatable :: centered_response(:, :), centered_design(:, :)
      real(dp), allocatable :: gradient(:, :), trial(:)
      real(dp) :: upper, lower, middle, exponent
      integer :: variables, predictors, equation, iteration
      logical :: zero_solution

      allocate(lambdas(max(0, grid_size)))
      if (grid_size < 1) return
      lambdas = 0.0_dp
      variables = size(series, 2)
      predictors = variables*lag_order
      if (variables < 1 .or. lag_order < 1 .or. &
         size(series, 1) <= lag_order .or. grid_ratio < 1.0_dp .or. &
         (penalty /= bigtime_penalty_l1 .and. &
         penalty /= bigtime_penalty_hlag)) return

      call build_var_data(series, lag_order, response, design)
      allocate(centered_response(size(response, 1), variables))
      allocate(centered_design(size(design, 1), predictors))
      centered_response = response - spread( &
         sum(response, dim=1)/real(size(response, 1), dp), &
         1, size(response, 1))
      centered_design = design - spread( &
         sum(design, dim=1)/real(size(design, 1), dp), &
         1, size(design, 1))
      allocate(gradient(variables, predictors))
      gradient = matmul(transpose(centered_response), centered_design)
      upper = max(1.0_dp, maxval(abs(gradient)))
      allocate(trial(predictors))
      do
         zero_solution = .true.
         do equation = 1, variables
            trial = apply_penalty_prox(gradient(equation, :), upper, &
               penalty, variables, lag_order)
            zero_solution = zero_solution .and. &
               maxval(abs(trial)) <= 100.0_dp*epsilon(1.0_dp)
         end do
         if (zero_solution) exit
         upper = 2.0_dp*upper
      end do
      lower = 0.0_dp
      do iteration = 1, 60
         middle = 0.5_dp*(lower + upper)
         zero_solution = .true.
         do equation = 1, variables
            trial = apply_penalty_prox(gradient(equation, :), middle, &
               penalty, variables, lag_order)
            zero_solution = zero_solution .and. &
               maxval(abs(trial)) <= 100.0_dp*epsilon(1.0_dp)
         end do
         if (zero_solution) then
            upper = middle
         else
            lower = middle
         end if
      end do
      upper = max(upper*(1.0_dp + 1.0e-8_dp), tiny(1.0_dp))
      if (grid_size == 1) then
         lambdas(1) = upper
         return
      end if
      do iteration = 1, grid_size
         exponent = real(iteration - 1, dp)/real(grid_size - 1, dp)
         lambdas(iteration) = upper/grid_ratio**exponent
      end do
   end function bigtime_var_lambda_grid

   pure function bigtime_sparse_varx(endogenous, exogenous, ar_order, &
      exogenous_order, lambda_phi, lambda_beta, penalty, alpha, tolerance, &
      max_iterations, initial_phi, initial_beta) result(out)
      ! Estimate one sparse VARX by bigtime's joint accelerated proximal method.
      real(dp), intent(in) :: endogenous(:, :), exogenous(:, :)
      integer, intent(in) :: ar_order, exogenous_order, penalty
      real(dp), intent(in) :: lambda_phi, lambda_beta
      real(dp), intent(in), optional :: alpha, tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: initial_phi(:, :), initial_beta(:, :)
      type(bigtime_varx_fit_t) :: out
      real(dp), allocatable :: response(:, :), ar_design(:, :), x_design(:, :)
      real(dp), allocatable :: centered_response(:, :)
      real(dp), allocatable :: centered_ar(:, :), centered_x(:, :)
      real(dp), allocatable :: response_mean(:), ar_mean(:), x_mean(:)
      real(dp), allocatable :: combined(:, :), gram(:, :)
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: row_response(:), phi_row(:), beta_row(:)
      real(dp) :: selected_alpha, selected_tolerance, step
      integer :: selected_max_iterations, variables, exogenous_variables
      integer :: ar_predictors, x_predictors, equation, eigen_info
      integer :: row_iterations
      logical :: row_converged

      out%ar_order = ar_order
      out%exogenous_order = exogenous_order
      out%lambda_phi = lambda_phi
      out%lambda_beta = lambda_beta
      variables = size(endogenous, 2)
      exogenous_variables = size(exogenous, 2)
      if (size(endogenous, 1) /= size(exogenous, 1) .or. &
         variables < 1 .or. exogenous_variables < 1 .or. ar_order < 1 .or. &
         exogenous_order < 1 .or. size(endogenous, 1) <= &
         max(ar_order, exogenous_order) .or. lambda_phi < 0.0_dp .or. &
         lambda_beta < 0.0_dp) then
         out%info = 1
         return
      end if
      if (penalty /= bigtime_penalty_l1 .and. &
         penalty /= bigtime_penalty_hlag) then
         out%info = 2
         return
      end if
      ar_predictors = variables*ar_order
      x_predictors = exogenous_variables*exogenous_order
      if (present(initial_phi)) then
         if (any(shape(initial_phi) /= [variables, ar_predictors])) then
            out%info = 3
            return
         end if
      end if
      if (present(initial_beta)) then
         if (any(shape(initial_beta) /= [variables, x_predictors])) then
            out%info = 3
            return
         end if
      end if

      selected_alpha = 0.0_dp
      if (present(alpha)) selected_alpha = alpha
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 100
      if (present(max_iterations)) selected_max_iterations = max_iterations
      out%alpha = selected_alpha
      if (selected_alpha < 0.0_dp .or. selected_tolerance <= 0.0_dp .or. &
         selected_max_iterations < 4) then
         out%info = 1
         return
      end if

      call build_varx_data(endogenous, exogenous, ar_order, exogenous_order, &
         response, ar_design, x_design)
      allocate(response_mean(variables), ar_mean(ar_predictors))
      allocate(x_mean(x_predictors))
      response_mean = sum(response, dim=1)/real(size(response, 1), dp)
      ar_mean = sum(ar_design, dim=1)/real(size(ar_design, 1), dp)
      x_mean = sum(x_design, dim=1)/real(size(x_design, 1), dp)
      allocate(centered_response(size(response, 1), variables))
      allocate(centered_ar(size(ar_design, 1), ar_predictors))
      allocate(centered_x(size(x_design, 1), x_predictors))
      centered_response = response - &
         spread(response_mean, 1, size(response, 1))
      centered_ar = ar_design - spread(ar_mean, 1, size(ar_design, 1))
      centered_x = x_design - spread(x_mean, 1, size(x_design, 1))
      allocate(combined(size(response, 1), ar_predictors + x_predictors))
      combined(:, :ar_predictors) = centered_ar
      combined(:, ar_predictors + 1:) = centered_x
      allocate(gram(size(combined, 2), size(combined, 2)))
      gram = matmul(transpose(combined), combined)
      call symmetric_eigen(gram, eigenvalues, eigenvectors, eigen_info)
      if (eigen_info /= 0) then
         out%info = 10 + eigen_info
         return
      end if

      allocate(out%phi(variables, ar_predictors))
      allocate(out%beta(variables, x_predictors))
      out%phi = 0.0_dp
      out%beta = 0.0_dp
      if (present(initial_phi)) out%phi = initial_phi
      if (present(initial_beta)) out%beta = initial_beta
      out%converged = .true.
      if (maxval(abs(eigenvalues)) <= tiny(1.0_dp)) then
         out%iterations = 0
      else
         step = 1.0_dp/maxval(eigenvalues)
         allocate(row_response(size(response, 1)))
         allocate(phi_row(ar_predictors), beta_row(x_predictors))
         do equation = 1, variables
            row_response = centered_response(:, equation)
            call sparse_varx_row(row_response, centered_ar, centered_x, &
               out%phi(equation, :), out%beta(equation, :), lambda_phi, &
               lambda_beta, penalty, variables, ar_order, &
               exogenous_variables, exogenous_order, step, selected_alpha, &
               selected_tolerance, selected_max_iterations, phi_row, beta_row, &
               row_iterations, row_converged)
            out%phi(equation, :) = phi_row
            out%beta(equation, :) = beta_row
            out%iterations = max(out%iterations, row_iterations)
            out%converged = out%converged .and. row_converged
         end do
      end if

      allocate(out%intercept(variables))
      out%intercept = response_mean - matmul(out%phi, ar_mean) - &
         matmul(out%beta, x_mean)
      allocate(out%fitted(size(response, 1), variables))
      out%fitted = matmul(ar_design, transpose(out%phi)) + &
         matmul(x_design, transpose(out%beta)) + &
         spread(out%intercept, 1, size(response, 1))
      allocate(out%residuals(size(response, 1), variables))
      out%residuals = response - out%fitted
      out%objective = 0.5_dp*sum(out%residuals**2) + &
         lambda_phi*penalty_value(out%phi, penalty, variables, ar_order) + &
         lambda_beta*penalty_value(out%beta, penalty, &
         exogenous_variables, exogenous_order)
      if (maxval(abs(eigenvalues)) > tiny(1.0_dp)) then
         out%objective = out%objective + &
            0.5_dp*selected_alpha/step*(sum(out%phi**2) + sum(out%beta**2))
      end if
      out%nonzero = count(abs(out%phi) > 100.0_dp*epsilon(1.0_dp)) + &
         count(abs(out%beta) > 100.0_dp*epsilon(1.0_dp))
      if (out%converged) then
         out%info = 0
      else
         out%info = 4
      end if
   end function bigtime_sparse_varx

   pure function bigtime_sparse_varx_path(endogenous, exogenous, ar_order, &
      exogenous_order, lambda_phi, lambda_beta, penalty, alpha, tolerance, &
      max_iterations, warm_start) result(out)
      ! Estimate sparse VARX models over paired penalty values.
      real(dp), intent(in) :: endogenous(:, :), exogenous(:, :)
      integer, intent(in) :: ar_order, exogenous_order, penalty
      real(dp), intent(in) :: lambda_phi(:), lambda_beta(:)
      real(dp), intent(in), optional :: alpha, tolerance
      integer, intent(in), optional :: max_iterations
      logical, intent(in), optional :: warm_start
      type(bigtime_varx_path_t) :: out
      type(bigtime_varx_fit_t) :: fit
      real(dp), allocatable :: starting_phi(:, :), starting_beta(:, :)
      real(dp) :: selected_alpha, selected_tolerance
      integer :: selected_max_iterations, variables, exogenous_variables
      integer :: ar_predictors, x_predictors, index
      logical :: use_warm_start

      out%ar_order = ar_order
      out%exogenous_order = exogenous_order
      if (size(lambda_phi) < 1 .or. size(lambda_phi) /= size(lambda_beta) .or. &
         any(lambda_phi < 0.0_dp) .or. any(lambda_beta < 0.0_dp)) then
         out%info = 1
         return
      end if
      selected_alpha = 0.0_dp
      if (present(alpha)) selected_alpha = alpha
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 100
      if (present(max_iterations)) selected_max_iterations = max_iterations
      use_warm_start = .true.
      if (present(warm_start)) use_warm_start = warm_start
      variables = size(endogenous, 2)
      exogenous_variables = size(exogenous, 2)
      ar_predictors = variables*ar_order
      x_predictors = exogenous_variables*exogenous_order
      if (variables < 1 .or. exogenous_variables < 1 .or. &
         ar_predictors < 1 .or. x_predictors < 1) then
         out%info = 1
         return
      end if

      allocate(out%phi(variables, ar_predictors, size(lambda_phi)))
      allocate(out%beta(variables, x_predictors, size(lambda_phi)))
      allocate(out%intercept(variables, size(lambda_phi)))
      allocate(out%objective(size(lambda_phi)))
      allocate(out%iterations(size(lambda_phi)))
      allocate(out%nonzero(size(lambda_phi)))
      allocate(out%converged(size(lambda_phi)))
      allocate(out%lambda_phi(size(lambda_phi)))
      allocate(out%lambda_beta(size(lambda_beta)))
      out%lambda_phi = lambda_phi
      out%lambda_beta = lambda_beta
      out%penalty = penalty
      out%alpha = selected_alpha
      allocate(starting_phi(variables, ar_predictors))
      allocate(starting_beta(variables, x_predictors))
      starting_phi = 0.0_dp
      starting_beta = 0.0_dp
      do index = 1, size(lambda_phi)
         fit = bigtime_sparse_varx(endogenous, exogenous, ar_order, &
            exogenous_order, lambda_phi(index), lambda_beta(index), penalty, &
            selected_alpha, selected_tolerance, selected_max_iterations, &
            starting_phi, starting_beta)
         if (fit%info /= 0 .and. fit%info /= 4) then
            out%info = 100*index + fit%info
            return
         end if
         out%phi(:, :, index) = fit%phi
         out%beta(:, :, index) = fit%beta
         out%intercept(:, index) = fit%intercept
         out%objective(index) = fit%objective
         out%iterations(index) = fit%iterations
         out%nonzero(index) = fit%nonzero
         out%converged(index) = fit%converged
         if (use_warm_start) then
            starting_phi = fit%phi
            starting_beta = fit%beta
         else
            starting_phi = 0.0_dp
            starting_beta = 0.0_dp
         end if
      end do
      if (all(out%converged)) then
         out%info = 0
      else
         out%info = 4
      end if
   end function bigtime_sparse_varx_path

   pure function bigtime_varx_lambda_grid(endogenous, exogenous, ar_order, &
      exogenous_order, penalty, phi_ratio, phi_size, beta_ratio, beta_size) &
      result(out)
      ! Construct separate bigtime penalty grids by joint zero-model searches.
      real(dp), intent(in) :: endogenous(:, :), exogenous(:, :)
      integer, intent(in) :: ar_order, exogenous_order, penalty
      real(dp), intent(in) :: phi_ratio, beta_ratio
      integer, intent(in) :: phi_size, beta_size
      type(bigtime_varx_grid_t) :: out
      real(dp), allocatable :: response(:, :), ar_design(:, :), x_design(:, :)
      real(dp), allocatable :: centered_response(:, :)
      real(dp), allocatable :: centered_ar(:, :), centered_x(:, :)
      real(dp), allocatable :: phi_gradient(:, :), beta_gradient(:, :)
      real(dp) :: phi_upper, beta_upper, phi_lower, beta_lower
      real(dp) :: phi_middle, beta_middle, exponent
      integer :: variables, x_variables, equation, iteration
      logical :: phi_zero, beta_zero

      allocate(out%lambda_phi(max(0, phi_size)))
      allocate(out%lambda_beta(max(0, beta_size)))
      if (phi_size < 1 .or. beta_size < 1 .or. phi_ratio < 1.0_dp .or. &
         beta_ratio < 1.0_dp .or. size(endogenous, 1) /= size(exogenous, 1) .or. &
         size(endogenous, 2) < 1 .or. size(exogenous, 2) < 1 .or. &
         size(endogenous, 1) <= max(ar_order, exogenous_order) .or. &
         ar_order < 1 .or. exogenous_order < 1 .or. &
         (penalty /= bigtime_penalty_l1 .and. &
         penalty /= bigtime_penalty_hlag)) then
         out%info = 1
         return
      end if
      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      call build_varx_data(endogenous, exogenous, ar_order, exogenous_order, &
         response, ar_design, x_design)
      allocate(centered_response(size(response, 1), variables))
      allocate(centered_ar(size(ar_design, 1), size(ar_design, 2)))
      allocate(centered_x(size(x_design, 1), size(x_design, 2)))
      centered_response = response - spread( &
         sum(response, dim=1)/real(size(response, 1), dp), &
         1, size(response, 1))
      centered_ar = ar_design - spread( &
         sum(ar_design, dim=1)/real(size(ar_design, 1), dp), &
         1, size(ar_design, 1))
      centered_x = x_design - spread( &
         sum(x_design, dim=1)/real(size(x_design, 1), dp), &
         1, size(x_design, 1))
      allocate(phi_gradient(variables, size(ar_design, 2)))
      allocate(beta_gradient(variables, size(x_design, 2)))
      phi_gradient = matmul(transpose(centered_response), centered_ar)
      beta_gradient = matmul(transpose(centered_response), centered_x)
      phi_upper = max(1.0_dp, maxval(abs(phi_gradient)))
      beta_upper = max(1.0_dp, maxval(abs(beta_gradient)))
      phi_lower = 0.0_dp
      beta_lower = 0.0_dp
      do iteration = 1, 60
         phi_middle = 0.5_dp*(phi_lower + phi_upper)
         beta_middle = 0.5_dp*(beta_lower + beta_upper)
         phi_zero = .true.
         beta_zero = .true.
         do equation = 1, variables
            phi_zero = phi_zero .and. maxval(abs(apply_penalty_prox( &
               phi_gradient(equation, :), phi_middle, penalty, variables, &
               ar_order))) <= 100.0_dp*epsilon(1.0_dp)
            beta_zero = beta_zero .and. maxval(abs(apply_penalty_prox( &
               beta_gradient(equation, :), beta_middle, penalty, x_variables, &
               exogenous_order))) <= 100.0_dp*epsilon(1.0_dp)
         end do
         if (phi_zero) then
            phi_upper = phi_middle
         else
            phi_lower = phi_middle
         end if
         if (beta_zero) then
            beta_upper = beta_middle
         else
            beta_lower = beta_middle
         end if
      end do
      phi_upper = max(phi_upper*(1.0_dp + 1.0e-8_dp), tiny(1.0_dp))
      beta_upper = max(beta_upper*(1.0_dp + 1.0e-8_dp), tiny(1.0_dp))
      do iteration = 1, phi_size
         if (phi_size == 1) then
            exponent = 0.0_dp
         else
            exponent = real(iteration - 1, dp)/real(phi_size - 1, dp)
         end if
         out%lambda_phi(iteration) = phi_upper/phi_ratio**exponent
      end do
      do iteration = 1, beta_size
         if (beta_size == 1) then
            exponent = 0.0_dp
         else
            exponent = real(iteration - 1, dp)/real(beta_size - 1, dp)
         end if
         out%lambda_beta(iteration) = beta_upper/beta_ratio**exponent
      end do
      out%info = 0
   end function bigtime_varx_lambda_grid

   pure function bigtime_varma_innovation_proxy(series, phase1) result(proxy)
      ! Extend and demean Phase I residuals using bigtime's boundary convention.
      real(dp), intent(in) :: series(:, :)
      type(bigtime_var_fit_t), intent(in) :: phase1
      real(dp), allocatable :: proxy(:, :)
      integer :: observations, variables, leading

      observations = size(series, 1)
      variables = size(series, 2)
      allocate(proxy(observations, variables))
      proxy = 0.0_dp
      leading = phase1%lag_order
      if (.not. allocated(phase1%residuals)) return
      if (leading < 1 .or. size(phase1%residuals, 1) /= &
         observations - leading .or. size(phase1%residuals, 2) /= variables .or. &
         size(phase1%residuals, 1) < leading) return
      proxy(leading + 1:, :) = phase1%residuals
      proxy(:leading, :) = phase1%residuals(:leading, :)
      proxy = proxy - spread(sum(proxy, dim=1)/real(observations, dp), &
         1, observations)
   end function bigtime_varma_innovation_proxy

   pure function bigtime_sparse_varma(series, phase1_order, phase1_lambda, &
      phase1_penalty, ar_order, ma_order, lambda_phi, lambda_theta, &
      phase2_penalty, innovations, alpha, tolerance, phase1_max_iterations, &
      phase2_max_iterations) result(out)
      ! Estimate bigtime's two-stage sparse VARMA at fixed penalty values.
      real(dp), intent(in) :: series(:, :)
      integer, intent(in) :: phase1_order, phase1_penalty, ar_order, ma_order
      real(dp), intent(in) :: phase1_lambda, lambda_phi, lambda_theta
      integer, intent(in) :: phase2_penalty
      real(dp), intent(in), optional :: innovations(:, :), alpha, tolerance
      integer, intent(in), optional :: phase1_max_iterations
      integer, intent(in), optional :: phase2_max_iterations
      type(bigtime_varma_fit_t) :: out
      real(dp) :: selected_alpha, selected_tolerance
      integer :: selected_phase1_iterations, selected_phase2_iterations

      out%phase1_order = phase1_order
      out%ar_order = ar_order
      out%ma_order = ma_order
      if (size(series, 1) <= max(phase1_order, max(ar_order, ma_order)) .or. &
         size(series, 2) < 1 .or. phase1_order < 1 .or. ar_order < 1 .or. &
         ma_order < 1 .or. phase1_lambda < 0.0_dp .or. lambda_phi < 0.0_dp .or. &
         lambda_theta < 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(innovations)) then
         if (any(shape(innovations) /= shape(series))) then
            out%info = 2
            return
         end if
      else if (size(series, 1) < 2*phase1_order) then
         out%info = 1
         return
      end if
      selected_alpha = 0.0_dp
      if (present(alpha)) selected_alpha = alpha
      selected_tolerance = 1.0e-3_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_phase1_iterations = 10000
      if (present(phase1_max_iterations)) then
         selected_phase1_iterations = phase1_max_iterations
      end if
      selected_phase2_iterations = 100
      if (present(phase2_max_iterations)) then
         selected_phase2_iterations = phase2_max_iterations
      end if

      if (present(innovations)) then
         allocate(out%innovations(size(innovations, 1), size(innovations, 2)))
         out%innovations = innovations
         out%estimated_innovations = .false.
      else
         out%phase1 = bigtime_sparse_var(series, phase1_order, phase1_lambda, &
            phase1_penalty, selected_tolerance, selected_phase1_iterations)
         if (out%phase1%info /= 0 .and. out%phase1%info /= 4) then
            out%info = 100 + out%phase1%info
            return
         end if
         out%innovations = bigtime_varma_innovation_proxy(series, out%phase1)
         out%estimated_innovations = .true.
      end if

      out%phase2 = bigtime_sparse_varx(series, out%innovations, ar_order, &
         ma_order, lambda_phi, lambda_theta, phase2_penalty, selected_alpha, &
         selected_tolerance, selected_phase2_iterations)
      if (out%phase2%info /= 0 .and. out%phase2%info /= 4) then
         out%info = 200 + out%phase2%info
         return
      end if
      out%phi = out%phase2%phi
      out%theta = out%phase2%beta
      out%intercept = out%phase2%intercept
      out%fitted = out%phase2%fitted
      out%residuals = out%phase2%residuals
      if (out%estimated_innovations) then
         out%converged = out%phase1%converged .and. out%phase2%converged
      else
         out%converged = out%phase2%converged
      end if
      if (out%converged) then
         out%info = 0
      else
         out%info = 4
      end if
   end function bigtime_sparse_varma

   pure function bigtime_sparse_varma_path(series, phase1_order, &
      phase1_lambda, phase1_penalty, ar_order, ma_order, lambda_phi, &
      lambda_theta, phase2_penalty, innovations, alpha, tolerance, &
      phase1_max_iterations, phase2_max_iterations, warm_start) result(out)
      ! Estimate a Phase II VARMA path from one shared Phase I residual proxy.
      real(dp), intent(in) :: series(:, :)
      integer, intent(in) :: phase1_order, phase1_penalty, ar_order, ma_order
      real(dp), intent(in) :: phase1_lambda, lambda_phi(:), lambda_theta(:)
      integer, intent(in) :: phase2_penalty
      real(dp), intent(in), optional :: innovations(:, :), alpha, tolerance
      integer, intent(in), optional :: phase1_max_iterations
      integer, intent(in), optional :: phase2_max_iterations
      logical, intent(in), optional :: warm_start
      type(bigtime_varma_path_t) :: out
      real(dp) :: selected_alpha, selected_tolerance
      integer :: selected_phase1_iterations, selected_phase2_iterations
      logical :: selected_warm_start

      out%phase1_order = phase1_order
      out%ar_order = ar_order
      out%ma_order = ma_order
      if (size(lambda_phi) < 1 .or. size(lambda_phi) /= size(lambda_theta) .or. &
         any(lambda_phi < 0.0_dp) .or. any(lambda_theta < 0.0_dp) .or. &
         size(series, 1) <= max(phase1_order, max(ar_order, ma_order)) .or. &
         size(series, 2) < 1 .or. phase1_order < 1 .or. ar_order < 1 .or. &
         ma_order < 1 .or. phase1_lambda < 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(innovations)) then
         if (any(shape(innovations) /= shape(series))) then
            out%info = 2
            return
         end if
      else if (size(series, 1) < 2*phase1_order) then
         out%info = 1
         return
      end if
      selected_alpha = 0.0_dp
      if (present(alpha)) selected_alpha = alpha
      selected_tolerance = 1.0e-3_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_phase1_iterations = 10000
      if (present(phase1_max_iterations)) then
         selected_phase1_iterations = phase1_max_iterations
      end if
      selected_phase2_iterations = 100
      if (present(phase2_max_iterations)) then
         selected_phase2_iterations = phase2_max_iterations
      end if
      selected_warm_start = .true.
      if (present(warm_start)) selected_warm_start = warm_start

      if (present(innovations)) then
         allocate(out%innovations(size(innovations, 1), size(innovations, 2)))
         out%innovations = innovations
         out%estimated_innovations = .false.
      else
         out%phase1 = bigtime_sparse_var(series, phase1_order, phase1_lambda, &
            phase1_penalty, selected_tolerance, selected_phase1_iterations)
         if (out%phase1%info /= 0 .and. out%phase1%info /= 4) then
            out%info = 100 + out%phase1%info
            return
         end if
         out%innovations = bigtime_varma_innovation_proxy(series, out%phase1)
         out%estimated_innovations = .true.
      end if

      out%phase2 = bigtime_sparse_varx_path(series, out%innovations, ar_order, &
         ma_order, lambda_phi, lambda_theta, phase2_penalty, selected_alpha, &
         selected_tolerance, selected_phase2_iterations, selected_warm_start)
      if (out%phase2%info /= 0 .and. out%phase2%info /= 4) then
         out%info = 200 + out%phase2%info
         return
      end if
      if (out%estimated_innovations .and. .not. out%phase1%converged) then
         out%info = 4
      else
         out%info = out%phase2%info
      end if
   end function bigtime_sparse_varma_path

   pure function bigtime_var_cv(series, lag_order, lambdas, penalty, &
      training_fraction, tolerance, max_iterations) result(out)
      ! Select a sparse VAR penalty by expanding-window one-step forecasts.
      real(dp), intent(in) :: series(:, :), lambdas(:)
      integer, intent(in) :: lag_order, penalty
      real(dp), intent(in), optional :: training_fraction, tolerance
      integer, intent(in), optional :: max_iterations
      type(bigtime_cv_result_t) :: out
      type(bigtime_var_path_t) :: path
      real(dp) :: selected_fraction, selected_tolerance
      real(dp), allocatable :: prediction(:)
      integer :: selected_iterations, effective, first_train, folds
      integer :: fold, train_rows, origin, candidate, variables, lag

      selected_fraction = 0.9_dp
      if (present(training_fraction)) selected_fraction = training_fraction
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_iterations = 10000
      if (present(max_iterations)) selected_iterations = max_iterations
      variables = size(series, 2)
      effective = size(series, 1) - lag_order
      first_train = floor(selected_fraction*real(effective, dp))
      folds = effective - first_train
      if (variables < 1 .or. lag_order < 1 .or. size(lambdas) < 1 .or. &
         any(lambdas < 0.0_dp) .or. selected_fraction <= 0.0_dp .or. &
         selected_fraction >= 1.0_dp .or. first_train < 1 .or. folds < 1) then
         out%info = 1
         return
      end if

      allocate(out%squared_error(folds, size(lambdas)))
      allocate(out%mean_nonzero(size(lambdas)))
      out%mean_nonzero = 0.0_dp
      allocate(prediction(variables))
      do fold = 1, folds
         train_rows = first_train + fold - 1
         origin = lag_order + train_rows
         path = bigtime_sparse_var_path(series(:origin, :), lag_order, &
            lambdas, penalty, selected_tolerance, selected_iterations)
         if (path%info /= 0 .and. path%info /= 4) then
            out%info = 100*fold + path%info
            return
         end if
         do candidate = 1, size(lambdas)
            prediction = path%intercept(:, candidate)
            do lag = 1, lag_order
               prediction = prediction + matmul( &
                  path%phi(:, (lag - 1)*variables + 1:lag*variables, candidate), &
                  series(origin + 1 - lag, :))
            end do
            out%squared_error(fold, candidate) = sum( &
               (series(origin + 1, :) - prediction)**2)/real(variables, dp)
            out%mean_nonzero(candidate) = out%mean_nonzero(candidate) + &
               real(path%nonzero(candidate), dp)
         end do
      end do
      out%mean_nonzero = out%mean_nonzero/real(folds, dp)
      call summarize_cv(out)
   end function bigtime_var_cv

   pure function bigtime_varx_cv(endogenous, exogenous, ar_order, &
      exogenous_order, lambda_phi, lambda_beta, penalty, alpha, &
      training_fraction, tolerance, max_iterations) result(out)
      ! Select paired sparse VARX penalties by expanding-window forecasts.
      real(dp), intent(in) :: endogenous(:, :), exogenous(:, :)
      integer, intent(in) :: ar_order, exogenous_order, penalty
      real(dp), intent(in) :: lambda_phi(:), lambda_beta(:)
      real(dp), intent(in), optional :: alpha, training_fraction, tolerance
      integer, intent(in), optional :: max_iterations
      type(bigtime_cv_result_t) :: out
      type(bigtime_varx_path_t) :: path
      real(dp) :: selected_alpha, selected_fraction, selected_tolerance
      real(dp), allocatable :: prediction(:)
      integer :: selected_iterations, maximum_order, effective, first_train
      integer :: folds, fold, train_rows, origin, candidate, variables
      integer :: x_variables, lag

      selected_alpha = 0.0_dp
      if (present(alpha)) selected_alpha = alpha
      selected_fraction = 0.9_dp
      if (present(training_fraction)) selected_fraction = training_fraction
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_iterations = 100
      if (present(max_iterations)) selected_iterations = max_iterations
      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      maximum_order = max(ar_order, exogenous_order)
      effective = size(endogenous, 1) - maximum_order
      first_train = floor(selected_fraction*real(effective, dp))
      folds = effective - first_train
      if (size(endogenous, 1) /= size(exogenous, 1) .or. variables < 1 .or. &
         x_variables < 1 .or. ar_order < 1 .or. exogenous_order < 1 .or. &
         size(lambda_phi) < 1 .or. size(lambda_phi) /= size(lambda_beta) .or. &
         any(lambda_phi < 0.0_dp) .or. any(lambda_beta < 0.0_dp) .or. &
         selected_fraction <= 0.0_dp .or. selected_fraction >= 1.0_dp .or. &
         first_train < 1 .or. folds < 1) then
         out%info = 1
         return
      end if

      allocate(out%squared_error(folds, size(lambda_phi)))
      allocate(out%mean_nonzero(size(lambda_phi)))
      out%mean_nonzero = 0.0_dp
      allocate(prediction(variables))
      do fold = 1, folds
         train_rows = first_train + fold - 1
         origin = maximum_order + train_rows
         path = bigtime_sparse_varx_path(endogenous(:origin, :), &
            exogenous(:origin, :), ar_order, exogenous_order, lambda_phi, &
            lambda_beta, penalty, selected_alpha, selected_tolerance, &
            selected_iterations)
         if (path%info /= 0 .and. path%info /= 4) then
            out%info = 100*fold + path%info
            return
         end if
         do candidate = 1, size(lambda_phi)
            prediction = path%intercept(:, candidate)
            do lag = 1, ar_order
               prediction = prediction + matmul(path%phi(:, &
                  (lag - 1)*variables + 1:lag*variables, candidate), &
                  endogenous(origin + 1 - lag, :))
            end do
            do lag = 1, exogenous_order
               prediction = prediction + matmul(path%beta(:, &
                  (lag - 1)*x_variables + 1:lag*x_variables, candidate), &
                  exogenous(origin + 1 - lag, :))
            end do
            out%squared_error(fold, candidate) = sum( &
               (endogenous(origin + 1, :) - prediction)**2)/real(variables, dp)
            out%mean_nonzero(candidate) = out%mean_nonzero(candidate) + &
               real(path%nonzero(candidate), dp)
         end do
      end do
      out%mean_nonzero = out%mean_nonzero/real(folds, dp)
      call summarize_cv(out, lambda_phi, lambda_beta)
   end function bigtime_varx_cv

   pure function bigtime_varma_cv(series, phase1_order, phase1_lambda, &
      phase1_penalty, ar_order, ma_order, lambda_phi, lambda_theta, &
      phase2_penalty, innovations, alpha, training_fraction, tolerance, &
      phase1_max_iterations, phase2_max_iterations) result(out)
      ! Select Phase II VARMA penalties using one shared Phase I proxy.
      real(dp), intent(in) :: series(:, :)
      integer, intent(in) :: phase1_order, phase1_penalty, ar_order, ma_order
      real(dp), intent(in) :: phase1_lambda, lambda_phi(:), lambda_theta(:)
      integer, intent(in) :: phase2_penalty
      real(dp), intent(in), optional :: innovations(:, :)
      real(dp), intent(in), optional :: alpha, training_fraction, tolerance
      integer, intent(in), optional :: phase1_max_iterations
      integer, intent(in), optional :: phase2_max_iterations
      type(bigtime_varma_cv_t) :: out
      real(dp) :: selected_alpha, selected_fraction, selected_tolerance
      integer :: selected_phase1_iterations, selected_phase2_iterations

      selected_alpha = 0.0_dp
      if (present(alpha)) selected_alpha = alpha
      selected_fraction = 0.9_dp
      if (present(training_fraction)) selected_fraction = training_fraction
      selected_tolerance = 1.0e-3_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_phase1_iterations = 10000
      if (present(phase1_max_iterations)) then
         selected_phase1_iterations = phase1_max_iterations
      end if
      selected_phase2_iterations = 100
      if (present(phase2_max_iterations)) then
         selected_phase2_iterations = phase2_max_iterations
      end if
      if (present(innovations)) then
         if (any(shape(innovations) /= shape(series))) then
            out%info = 1
            return
         end if
         allocate(out%innovations(size(series, 1), size(series, 2)))
         out%innovations = innovations
      else
         if (size(series, 1) < 2*phase1_order) then
            out%info = 1
            return
         end if
         out%phase1 = bigtime_sparse_var(series, phase1_order, phase1_lambda, &
            phase1_penalty, selected_tolerance, selected_phase1_iterations)
         if (out%phase1%info /= 0 .and. out%phase1%info /= 4) then
            out%info = 100 + out%phase1%info
            return
         end if
         out%innovations = bigtime_varma_innovation_proxy(series, out%phase1)
         out%estimated_innovations = .true.
      end if
      out%phase2 = bigtime_varx_cv(series, out%innovations, ar_order, ma_order, &
         lambda_phi, lambda_theta, phase2_penalty, selected_alpha, &
         selected_fraction, selected_tolerance, selected_phase2_iterations)
      if (out%phase2%info /= 0) then
         out%info = 200 + out%phase2%info
      else
         out%info = 0
      end if
   end function bigtime_varma_cv

   pure function bigtime_var_path_ic(series, path) result(out)
      ! Compute AIC, BIC, and HQ for every sparse VAR path slice.
      real(dp), intent(in) :: series(:, :)
      type(bigtime_var_path_t), intent(in) :: path
      type(bigtime_ic_result_t) :: out
      real(dp), allocatable :: response(:, :), design(:, :), residuals(:, :)
      integer :: candidate, candidate_info

      if (.not. allocated(path%phi) .or. .not. allocated(path%intercept) .or. &
         path%lag_order < 1) then
         out%info = 1
         return
      end if
      call build_var_data(series, path%lag_order, response, design)
      allocate(out%value(size(path%phi, 3), 3))
      allocate(residuals(size(response, 1), size(response, 2)))
      do candidate = 1, size(path%phi, 3)
         residuals = response - matmul(design, &
            transpose(path%phi(:, :, candidate))) - &
            spread(path%intercept(:, candidate), 1, size(response, 1))
         call information_criteria(residuals, size(series, 1), &
            count(abs(path%phi(:, :, candidate)) > &
            100.0_dp*epsilon(1.0_dp)), out%value(candidate, :), candidate_info)
         if (candidate_info /= 0) then
            out%info = 100*candidate + candidate_info
            return
         end if
      end do
      call select_information_criteria(out)
   end function bigtime_var_path_ic

   pure function bigtime_varx_path_ic(endogenous, exogenous, path) result(out)
      ! Compute AIC, BIC, and HQ for every sparse VARX path slice.
      real(dp), intent(in) :: endogenous(:, :), exogenous(:, :)
      type(bigtime_varx_path_t), intent(in) :: path
      type(bigtime_ic_result_t) :: out
      real(dp), allocatable :: response(:, :), ar_design(:, :), x_design(:, :)
      real(dp), allocatable :: residuals(:, :)
      integer :: candidate, candidate_info, degrees

      if (.not. allocated(path%phi) .or. .not. allocated(path%beta) .or. &
         .not. allocated(path%intercept)) then
         out%info = 1
         return
      end if
      call build_varx_data(endogenous, exogenous, path%ar_order, &
         path%exogenous_order, response, ar_design, x_design)
      allocate(out%value(size(path%phi, 3), 3))
      allocate(residuals(size(response, 1), size(response, 2)))
      do candidate = 1, size(path%phi, 3)
         residuals = response - matmul(ar_design, &
            transpose(path%phi(:, :, candidate))) - matmul(x_design, &
            transpose(path%beta(:, :, candidate))) - &
            spread(path%intercept(:, candidate), 1, size(response, 1))
         degrees = count(abs(path%phi(:, :, candidate)) > &
            100.0_dp*epsilon(1.0_dp)) + count(abs(path%beta(:, :, candidate)) > &
            100.0_dp*epsilon(1.0_dp))
         call information_criteria(residuals, size(endogenous, 1), degrees, &
            out%value(candidate, :), candidate_info)
         if (candidate_info /= 0) then
            out%info = 100*candidate + candidate_info
            return
         end if
      end do
      call select_information_criteria(out)
   end function bigtime_varx_path_ic

   pure function bigtime_varma_path_ic(series, path) result(out)
      ! Compute information criteria for the Phase II sparse VARMA path.
      real(dp), intent(in) :: series(:, :)
      type(bigtime_varma_path_t), intent(in) :: path
      type(bigtime_ic_result_t) :: out

      if (.not. allocated(path%innovations)) then
         out%info = 1
         return
      end if
      out = bigtime_varx_path_ic(series, path%innovations, path%phase2)
   end function bigtime_varma_path_ic

   pure function bigtime_select_var_path(series, path, index) result(out)
      ! Materialize one fitted sparse VAR from a regularization path.
      real(dp), intent(in) :: series(:, :)
      type(bigtime_var_path_t), intent(in) :: path
      integer, intent(in) :: index
      type(bigtime_var_fit_t) :: out
      real(dp), allocatable :: response(:, :), design(:, :)

      if (.not. allocated(path%phi) .or. index < 1 .or. &
         index > size(path%phi, 3)) then
         out%info = 1
         return
      end if
      call build_var_data(series, path%lag_order, response, design)
      out%phi = path%phi(:, :, index)
      out%intercept = path%intercept(:, index)
      out%fitted = matmul(design, transpose(out%phi)) + &
         spread(out%intercept, 1, size(design, 1))
      out%residuals = response - out%fitted
      out%lag_order = path%lag_order
      if (allocated(path%lambda)) out%lambda = path%lambda(index)
      out%objective = path%objective(index)
      out%iterations = path%iterations(index)
      out%nonzero = path%nonzero(index)
      out%converged = path%converged(index)
      if (out%converged) then
         out%info = 0
      else
         out%info = 4
      end if
   end function bigtime_select_var_path

   pure function bigtime_select_varx_path(endogenous, exogenous, path, index) &
      result(out)
      ! Materialize one fitted sparse VARX from a paired regularization path.
      real(dp), intent(in) :: endogenous(:, :), exogenous(:, :)
      type(bigtime_varx_path_t), intent(in) :: path
      integer, intent(in) :: index
      type(bigtime_varx_fit_t) :: out
      real(dp), allocatable :: response(:, :), ar_design(:, :), x_design(:, :)

      if (.not. allocated(path%phi) .or. index < 1 .or. &
         index > size(path%phi, 3)) then
         out%info = 1
         return
      end if
      call build_varx_data(endogenous, exogenous, path%ar_order, &
         path%exogenous_order, response, ar_design, x_design)
      out%phi = path%phi(:, :, index)
      out%beta = path%beta(:, :, index)
      out%intercept = path%intercept(:, index)
      out%fitted = matmul(ar_design, transpose(out%phi)) + &
         matmul(x_design, transpose(out%beta)) + &
         spread(out%intercept, 1, size(response, 1))
      out%residuals = response - out%fitted
      out%ar_order = path%ar_order
      out%exogenous_order = path%exogenous_order
      if (allocated(path%lambda_phi)) then
         out%lambda_phi = path%lambda_phi(index)
         out%lambda_beta = path%lambda_beta(index)
      end if
      out%alpha = path%alpha
      out%objective = path%objective(index)
      out%iterations = path%iterations(index)
      out%nonzero = path%nonzero(index)
      out%converged = path%converged(index)
      if (out%converged) then
         out%info = 0
      else
         out%info = 4
      end if
   end function bigtime_select_varx_path

   pure function bigtime_var_forecast(fit, history, horizon) result(out)
      ! Recursively forecast a sparse VAR from its latest observations.
      type(bigtime_var_fit_t), intent(in) :: fit
      real(dp), intent(in) :: history(:, :)
      integer, intent(in) :: horizon
      type(bigtime_forecast_t) :: out
      real(dp), allocatable :: extended(:, :)
      integer :: observations, variables, step, lag, target

      observations = size(history, 1)
      variables = size(history, 2)
      if (.not. allocated(fit%phi) .or. .not. allocated(fit%intercept)) then
         out%info = 1
         return
      end if
      if (horizon < 1 .or. fit%lag_order < 1 .or. &
         observations < fit%lag_order .or. size(fit%phi, 1) /= variables .or. &
         size(fit%phi, 2) /= variables*fit%lag_order .or. &
         size(fit%intercept) /= variables) then
         out%info = 1
         return
      end if
      allocate(extended(observations + horizon, variables))
      extended(:observations, :) = history
      allocate(out%mean(horizon, variables))
      do step = 1, horizon
         target = observations + step
         extended(target, :) = fit%intercept
         do lag = 1, fit%lag_order
            extended(target, :) = extended(target, :) + matmul( &
               fit%phi(:, (lag - 1)*variables + 1:lag*variables), &
               extended(target - lag, :))
         end do
         out%mean(step, :) = extended(target, :)
      end do
      out%info = 0
   end function bigtime_var_forecast

   pure function bigtime_var_path_forecast(path, history, horizon) result(out)
      ! Recursively forecast every sparse VAR regularization-path slice.
      type(bigtime_var_path_t), intent(in) :: path
      real(dp), intent(in) :: history(:, :)
      integer, intent(in) :: horizon
      type(bigtime_path_forecast_t) :: out
      type(bigtime_var_fit_t) :: fit
      type(bigtime_forecast_t) :: forecast
      integer :: candidate

      if (.not. allocated(path%phi) .or. .not. allocated(path%intercept)) then
         out%info = 1
         return
      end if
      if (horizon < 1) then
         out%info = 1
         return
      end if
      allocate(out%mean(horizon, size(history, 2), size(path%phi, 3)))
      do candidate = 1, size(path%phi, 3)
         fit%phi = path%phi(:, :, candidate)
         fit%intercept = path%intercept(:, candidate)
         fit%lag_order = path%lag_order
         forecast = bigtime_var_forecast(fit, history, horizon)
         if (forecast%info /= 0) then
            out%info = 100*candidate + forecast%info
            return
         end if
         out%mean(:, :, candidate) = forecast%mean
      end do
      out%info = 0
   end function bigtime_var_path_forecast

   pure function bigtime_varx_forecast(fit, endogenous_history, &
      exogenous_values, horizon) result(out)
      ! Recursively forecast a sparse VARX using aligned exogenous values.
      type(bigtime_varx_fit_t), intent(in) :: fit
      real(dp), intent(in) :: endogenous_history(:, :), exogenous_values(:, :)
      integer, intent(in) :: horizon
      type(bigtime_forecast_t) :: out
      real(dp), allocatable :: extended(:, :)
      integer :: observations, variables, x_variables, step, lag, target

      observations = size(endogenous_history, 1)
      variables = size(endogenous_history, 2)
      x_variables = size(exogenous_values, 2)
      if (.not. allocated(fit%phi) .or. .not. allocated(fit%beta) .or. &
         .not. allocated(fit%intercept)) then
         out%info = 1
         return
      end if
      if (horizon < 1 .or. fit%ar_order < 1 .or. &
         fit%exogenous_order < 1 .or. &
         observations < fit%ar_order .or. &
         size(exogenous_values, 1) < observations + horizon - 1 .or. &
         size(fit%phi, 1) /= variables .or. &
         size(fit%phi, 2) /= variables*fit%ar_order .or. &
         size(fit%beta, 1) /= variables .or. &
         size(fit%beta, 2) /= x_variables*fit%exogenous_order) then
         out%info = 1
         return
      end if
      allocate(extended(observations + horizon, variables))
      extended(:observations, :) = endogenous_history
      allocate(out%mean(horizon, variables))
      do step = 1, horizon
         target = observations + step
         extended(target, :) = fit%intercept
         do lag = 1, fit%ar_order
            extended(target, :) = extended(target, :) + matmul( &
               fit%phi(:, (lag - 1)*variables + 1:lag*variables), &
               extended(target - lag, :))
         end do
         do lag = 1, fit%exogenous_order
            extended(target, :) = extended(target, :) + matmul( &
               fit%beta(:, (lag - 1)*x_variables + 1:lag*x_variables), &
               exogenous_values(target - lag, :))
         end do
         out%mean(step, :) = extended(target, :)
      end do
      out%info = 0
   end function bigtime_varx_forecast

   pure function bigtime_varx_path_forecast(path, endogenous_history, &
      exogenous_values, horizon) result(out)
      ! Recursively forecast every sparse VARX regularization-path slice.
      type(bigtime_varx_path_t), intent(in) :: path
      real(dp), intent(in) :: endogenous_history(:, :)
      real(dp), intent(in) :: exogenous_values(:, :)
      integer, intent(in) :: horizon
      type(bigtime_path_forecast_t) :: out
      type(bigtime_varx_fit_t) :: fit
      type(bigtime_forecast_t) :: forecast
      integer :: candidate

      if (.not. allocated(path%phi) .or. .not. allocated(path%beta) .or. &
         .not. allocated(path%intercept)) then
         out%info = 1
         return
      end if
      if (horizon < 1) then
         out%info = 1
         return
      end if
      allocate(out%mean(horizon, size(endogenous_history, 2), &
         size(path%phi, 3)))
      do candidate = 1, size(path%phi, 3)
         fit%phi = path%phi(:, :, candidate)
         fit%beta = path%beta(:, :, candidate)
         fit%intercept = path%intercept(:, candidate)
         fit%ar_order = path%ar_order
         fit%exogenous_order = path%exogenous_order
         forecast = bigtime_varx_forecast(fit, endogenous_history, &
            exogenous_values, horizon)
         if (forecast%info /= 0) then
            out%info = 100*candidate + forecast%info
            return
         end if
         out%mean(:, :, candidate) = forecast%mean
      end do
      out%info = 0
   end function bigtime_varx_path_forecast

   pure function bigtime_varma_forecast(fit, history, horizon) result(out)
      ! Recursively forecast sparse VARMA levels under zero future innovations.
      type(bigtime_varma_fit_t), intent(in) :: fit
      real(dp), intent(in) :: history(:, :)
      integer, intent(in) :: horizon
      type(bigtime_forecast_t) :: out
      real(dp), allocatable :: extended(:, :), innovation(:)
      integer :: observations, variables, step, lag, target, innovation_index

      observations = size(history, 1)
      variables = size(history, 2)
      if (.not. allocated(fit%phi) .or. .not. allocated(fit%theta) .or. &
         .not. allocated(fit%intercept) .or. &
         .not. allocated(fit%innovations)) then
         out%info = 1
         return
      end if
      if (horizon < 1 .or. fit%ar_order < 1 .or. fit%ma_order < 1 .or. &
         observations < max(fit%ar_order, fit%ma_order) .or. &
         any(shape(fit%innovations) /= shape(history)) .or. &
         size(fit%phi, 1) /= variables .or. &
         size(fit%phi, 2) /= variables*fit%ar_order .or. &
         size(fit%theta, 1) /= variables .or. &
         size(fit%theta, 2) /= variables*fit%ma_order) then
         out%info = 1
         return
      end if
      allocate(extended(observations + horizon, variables))
      allocate(innovation(variables))
      extended(:observations, :) = history
      allocate(out%mean(horizon, variables))
      do step = 1, horizon
         target = observations + step
         extended(target, :) = fit%intercept
         do lag = 1, fit%ar_order
            extended(target, :) = extended(target, :) + matmul( &
               fit%phi(:, (lag - 1)*variables + 1:lag*variables), &
               extended(target - lag, :))
         end do
         do lag = 1, fit%ma_order
            innovation_index = target - lag
            innovation = 0.0_dp
            if (innovation_index <= observations) then
               innovation = fit%innovations(innovation_index, :)
            end if
            extended(target, :) = extended(target, :) + matmul( &
               fit%theta(:, (lag - 1)*variables + 1:lag*variables), innovation)
         end do
         out%mean(step, :) = extended(target, :)
      end do
      out%info = 0
   end function bigtime_varma_forecast

   pure function bigtime_var_stability(phi, lag_order, margin) result(out)
      ! Diagnose VAR stability from the companion-matrix root moduli.
      real(dp), intent(in) :: phi(:, :)
      integer, intent(in) :: lag_order
      real(dp), intent(in), optional :: margin
      type(bigtime_stability_t) :: out
      real(dp), allocatable :: companion(:, :), perturbed(:, :)
      real(dp) :: selected_margin, perturbation
      integer :: variables, states, row, eigen_info

      variables = size(phi, 1)
      selected_margin = 0.0_dp
      if (present(margin)) selected_margin = margin
      if (variables < 1 .or. lag_order < 1 .or. &
         size(phi, 2) /= variables*lag_order .or. selected_margin < 0.0_dp .or. &
         selected_margin >= 1.0_dp) then
         out%info = 1
         return
      end if
      states = variables*lag_order
      allocate(companion(states, states))
      companion = 0.0_dp
      companion(:variables, :) = phi
      do row = variables + 1, states
         companion(row, row - variables) = 1.0_dp
      end do
      call general_eigenvalues(companion, out%roots, eigen_info)
      if (eigen_info /= 0) then
         allocate(perturbed(states, states))
         perturbed = companion
         perturbation = 1000.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(companion)))
         do row = 1, states
            perturbed(row, row) = perturbed(row, row) + &
               perturbation*real(row, dp)/real(states, dp)
         end do
         call general_eigenvalues(perturbed, out%roots, eigen_info)
         if (eigen_info /= 0) then
            out%info = 10 + eigen_info
            return
         end if
      end if
      out%maximum_modulus = maxval(abs(out%roots))
      out%stable = out%maximum_modulus < 1.0_dp - selected_margin
      out%info = 0
   end function bigtime_var_stability

   pure function bigtime_active_lags(primary, primary_predictors, &
      primary_order, secondary, secondary_predictors, secondary_order) &
      result(out)
      ! Report each response-predictor pair's largest nonzero lag.
      real(dp), intent(in) :: primary(:, :)
      integer, intent(in) :: primary_predictors, primary_order
      real(dp), intent(in), optional :: secondary(:, :)
      integer, intent(in), optional :: secondary_predictors, secondary_order
      type(bigtime_lag_order_t) :: out
      integer :: response, predictor, lag, responses
      real(dp), parameter :: zero_tolerance = 100.0_dp*epsilon(1.0_dp)

      responses = size(primary, 1)
      if (responses < 1 .or. primary_predictors < 1 .or. primary_order < 1 .or. &
         size(primary, 2) /= primary_predictors*primary_order) then
         out%info = 1
         return
      end if
      allocate(out%primary(responses, primary_predictors))
      out%primary = 0
      do response = 1, responses
         do predictor = 1, primary_predictors
            do lag = 1, primary_order
               if (abs(primary(response, predictor + &
                  (lag - 1)*primary_predictors)) > zero_tolerance) then
                  out%primary(response, predictor) = lag
               end if
            end do
         end do
      end do
      if (present(secondary)) then
         if (.not. present(secondary_predictors)) then
            out%info = 2
            return
         end if
         if (.not. present(secondary_order)) then
            out%info = 2
            return
         end if
         if (secondary_predictors < 1 .or. secondary_order < 1 .or. &
            size(secondary, 1) /= responses .or. &
            size(secondary, 2) /= secondary_predictors*secondary_order) then
            out%info = 2
            return
         end if
         allocate(out%secondary(responses, secondary_predictors))
         out%secondary = 0
         do response = 1, responses
            do predictor = 1, secondary_predictors
               do lag = 1, secondary_order
                  if (abs(secondary(response, predictor + &
                     (lag - 1)*secondary_predictors)) > zero_tolerance) then
                     out%secondary(response, predictor) = lag
                  end if
               end do
            end do
         end do
      end if
      out%info = 0
   end function bigtime_active_lags

   pure function bigtime_var_coefficients_from_draws(draws, lag_order, &
      maximum_modulus, sparsity, decay, zero_indices, trailing_zeros, &
      zero_self) result(out)
      ! Build and stabilize sparse VAR coefficients from supplied Gaussian draws.
      real(dp), intent(in) :: draws(:, :)
      integer, intent(in) :: lag_order, sparsity
      real(dp), intent(in) :: maximum_modulus, decay
      integer, intent(in), optional :: zero_indices(:)
      integer, intent(in), optional :: trailing_zeros(:, :)
      logical, intent(in), optional :: zero_self
      type(bigtime_coefficient_t) :: out
      type(bigtime_stability_t) :: stability
      logical :: allow_self_zeros
      integer :: variables, coefficient_count, lag, index, response, predictor
      integer :: zeros, maximum_iterations, iteration

      variables = size(draws, 1)
      coefficient_count = variables*lag_order
      if (variables < 1 .or. lag_order < 1 .or. &
         size(draws, 2) /= coefficient_count .or. maximum_modulus <= 0.0_dp .or. &
         maximum_modulus >= 1.0_dp .or. decay <= 0.0_dp .or. decay > 1.0_dp .or. &
         sparsity < bigtime_sparsity_dense .or. &
         sparsity > bigtime_sparsity_hlag) then
         out%info = 1
         return
      end if
      if (sparsity == bigtime_sparsity_l1 .and. .not. present(zero_indices)) then
         out%info = 2
         return
      end if
      if (sparsity == bigtime_sparsity_hlag .and. &
         .not. present(trailing_zeros)) then
         out%info = 2
         return
      end if

      allocate(out%phi(variables, coefficient_count))
      out%phi = draws
      do lag = 1, lag_order
         out%phi(:, (lag - 1)*variables + 1:lag*variables) = &
            decay**real(lag - 1, dp)* &
            out%phi(:, (lag - 1)*variables + 1:lag*variables)
      end do
      if (sparsity == bigtime_sparsity_l1) then
         if (any(zero_indices < 1) .or. &
            any(zero_indices > size(out%phi))) then
            out%info = 2
            return
         end if
         do index = 1, size(zero_indices)
            response = 1 + mod(zero_indices(index) - 1, variables)
            predictor = 1 + (zero_indices(index) - 1)/variables
            out%phi(response, predictor) = 0.0_dp
         end do
      else if (sparsity == bigtime_sparsity_hlag) then
         if (any(shape(trailing_zeros) /= [variables, variables]) .or. &
            any(trailing_zeros < 0) .or. any(trailing_zeros > lag_order)) then
            out%info = 2
            return
         end if
         allow_self_zeros = .true.
         if (present(zero_self)) allow_self_zeros = zero_self
         do response = 1, variables
            do predictor = 1, variables
               if (response == predictor .and. .not. allow_self_zeros) cycle
               zeros = trailing_zeros(response, predictor)
               if (zeros < 1) cycle
               do lag = lag_order - zeros + 1, lag_order
                  out%phi(response, predictor + (lag - 1)*variables) = 0.0_dp
               end do
            end do
         end do
      end if

      maximum_iterations = 100000
      do iteration = 0, maximum_iterations
         stability = bigtime_var_stability(out%phi, lag_order)
         if (stability%info /= 0) then
            out%info = 10 + stability%info
            return
         end if
         out%maximum_modulus = stability%maximum_modulus
         if (out%maximum_modulus <= maximum_modulus) exit
         out%phi = 0.99_dp*out%phi
      end do
      out%scaling_iterations = iteration
      if (iteration > maximum_iterations) then
         out%info = 3
         return
      end if
      out%companion = bigtime_companion_matrix(out%phi, lag_order)
      out%info = 0
   end function bigtime_var_coefficients_from_draws

   function bigtime_random_var_coefficients(variables, lag_order, &
      maximum_modulus, sparsity, decay, zero_count, zero_min, zero_max, &
      zero_self) result(out)
      ! Randomly generate and stabilize bigtime-style sparse VAR coefficients.
      integer, intent(in) :: variables, lag_order, sparsity
      real(dp), intent(in) :: maximum_modulus, decay
      integer, intent(in), optional :: zero_count, zero_min, zero_max
      logical, intent(in), optional :: zero_self
      type(bigtime_coefficient_t) :: out
      real(dp), allocatable :: draws(:, :)
      integer, allocatable :: indices(:), selected(:), trailing(:, :)
      integer :: total, requested, lower, upper, index, swap_index, temporary
      integer :: response, predictor
      logical :: allow_self_zeros

      if (variables < 1 .or. lag_order < 1) then
         out%info = 1
         return
      end if
      allocate(draws(variables, variables*lag_order))
      call random_standard_normal_matrix(draws)
      draws = 10.0_dp*draws
      if (sparsity == bigtime_sparsity_l1) then
         requested = variables*lag_order
         if (present(zero_count)) requested = zero_count
         total = size(draws)
         if (requested < 0 .or. requested > total) then
            out%info = 2
            return
         end if
         allocate(indices(total), selected(requested))
         indices = [(index, index=1, total)]
         do index = total, 2, -1
            swap_index = 1 + int(random_uniform()*real(index, dp))
            swap_index = min(index, swap_index)
            temporary = indices(index)
            indices(index) = indices(swap_index)
            indices(swap_index) = temporary
         end do
         if (requested > 0) selected = indices(:requested)
         out = bigtime_var_coefficients_from_draws(draws, lag_order, &
            maximum_modulus, sparsity, decay, zero_indices=selected)
      else if (sparsity == bigtime_sparsity_hlag) then
         lower = 0
         if (present(zero_min)) lower = zero_min
         upper = lag_order/2
         if (present(zero_max)) upper = zero_max
         allow_self_zeros = .true.
         if (present(zero_self)) allow_self_zeros = zero_self
         if (lower < 0 .or. upper < lower .or. upper > lag_order) then
            out%info = 2
            return
         end if
         allocate(trailing(variables, variables))
         do predictor = 1, variables
            do response = 1, variables
               trailing(response, predictor) = lower + int(random_uniform()* &
                  real(upper - lower + 1, dp))
               trailing(response, predictor) = min(upper, &
                  trailing(response, predictor))
            end do
         end do
         out = bigtime_var_coefficients_from_draws(draws, lag_order, &
            maximum_modulus, sparsity, decay, trailing_zeros=trailing, &
            zero_self=allow_self_zeros)
      else
         out = bigtime_var_coefficients_from_draws(draws, lag_order, &
            maximum_modulus, sparsity, decay)
      end if
   end function bigtime_random_var_coefficients

   pure function bigtime_var_simulate_from_innovations(phi, intercept, &
      innovations, initial_state, burnin) result(out)
      ! Simulate a VAR recursion from supplied innovations and companion state.
      real(dp), intent(in) :: phi(:, :), intercept(:), innovations(:, :)
      real(dp), intent(in) :: initial_state(:)
      integer, intent(in) :: burnin
      type(bigtime_simulation_t) :: out
      real(dp), allocatable :: companion(:, :), state(:), constant(:)
      integer :: variables, lag_order, states, total, periods, time

      variables = size(phi, 1)
      if (variables < 1 .or. mod(size(phi, 2), variables) /= 0) then
         out%info = 1
         return
      end if
      lag_order = size(phi, 2)/variables
      states = variables*lag_order
      total = size(innovations, 1)
      periods = total - burnin
      if (lag_order < 1 .or. size(intercept) /= variables .or. &
         size(innovations, 2) /= variables .or. size(initial_state) /= states .or. &
         burnin < 0 .or. periods < 1) then
         out%info = 1
         return
      end if
      companion = bigtime_companion_matrix(phi, lag_order)
      allocate(state(states), constant(states))
      state = initial_state
      constant = 0.0_dp
      constant(:variables) = intercept
      allocate(out%series(periods, variables))
      do time = 1, total
         state = constant + matmul(companion, state)
         state(:variables) = state(:variables) + innovations(time, :)
         if (time > burnin) out%series(time - burnin, :) = state(:variables)
      end do
      out%innovations = innovations
      out%initial_state = initial_state
      out%intercept = intercept
      out%phi = phi
      out%burnin = burnin
      out%info = 0
   end function bigtime_var_simulate_from_innovations

   function bigtime_var_simulate(phi, intercept, periods, burnin, &
      initial_state, innovation_scale) result(out)
      ! Simulate a Gaussian VAR using the shared random-number stream.
      real(dp), intent(in) :: phi(:, :), intercept(:)
      integer, intent(in) :: periods, burnin
      real(dp), intent(in), optional :: initial_state(:), innovation_scale
      type(bigtime_simulation_t) :: out
      real(dp), allocatable :: innovations(:, :), state(:)
      real(dp) :: scale
      integer :: variables, lag_order

      variables = size(phi, 1)
      if (variables < 1 .or. mod(size(phi, 2), variables) /= 0 .or. &
         periods < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      lag_order = size(phi, 2)/variables
      scale = 1.0_dp
      if (present(innovation_scale)) scale = innovation_scale
      if (scale < 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(innovations(periods + burnin, variables))
      call random_standard_normal_matrix(innovations)
      innovations = scale*innovations
      allocate(state(variables*lag_order))
      state = 0.0_dp
      if (present(initial_state)) then
         if (size(initial_state) /= size(state)) then
            out%info = 1
            return
         end if
         state = initial_state
      end if
      out = bigtime_var_simulate_from_innovations(phi, intercept, &
         innovations, state, burnin)
   end function bigtime_var_simulate

   pure function bigtime_companion_matrix(phi, lag_order) result(companion)
      ! Construct the block companion matrix for lag-major VAR coefficients.
      real(dp), intent(in) :: phi(:, :)
      integer, intent(in) :: lag_order
      real(dp) :: companion(size(phi, 1)*lag_order, &
         size(phi, 1)*lag_order)
      integer :: variables, row

      variables = size(phi, 1)
      companion = 0.0_dp
      companion(:variables, :) = phi
      do row = variables + 1, variables*lag_order
         companion(row, row - variables) = 1.0_dp
      end do
   end function bigtime_companion_matrix

   pure subroutine summarize_cv(out, lambda_phi, lambda_beta)
      ! Summarize fold errors and apply bigtime's one-standard-error rule.
      type(bigtime_cv_result_t), intent(inout) :: out
      real(dp), intent(in), optional :: lambda_phi(:), lambda_beta(:)
      real(dp) :: threshold, mean_value
      integer :: candidates, folds, candidate
      logical, allocatable :: eligible(:)

      folds = size(out%squared_error, 1)
      candidates = size(out%squared_error, 2)
      allocate(out%mean_squared_error(candidates))
      allocate(out%standard_error(candidates))
      allocate(eligible(candidates))
      do candidate = 1, candidates
         mean_value = sum(out%squared_error(:, candidate))/real(folds, dp)
         out%mean_squared_error(candidate) = mean_value
         if (folds > 1) then
            out%standard_error(candidate) = sqrt(sum( &
               (out%squared_error(:, candidate) - mean_value)**2)/ &
               real(folds - 1, dp))/sqrt(real(folds, dp))
         else
            out%standard_error(candidate) = 0.0_dp
         end if
      end do
      out%best_index = minloc(out%mean_squared_error, dim=1)
      threshold = out%mean_squared_error(out%best_index) + &
         out%standard_error(out%best_index)
      eligible = out%mean_squared_error < threshold
      if (present(lambda_phi) .and. present(lambda_beta)) then
         eligible = eligible .and. &
            lambda_phi >= lambda_phi(out%best_index) .and. &
            lambda_beta >= lambda_beta(out%best_index)
      end if
      out%one_se_index = 0
      do candidate = 1, candidates
         if (.not. eligible(candidate)) cycle
         if (out%one_se_index == 0) then
            out%one_se_index = candidate
         else if (out%mean_nonzero(candidate) < &
            out%mean_nonzero(out%one_se_index)) then
            out%one_se_index = candidate
         end if
      end do
      if (out%one_se_index == 0) out%one_se_index = out%best_index
      out%info = 0
   end subroutine summarize_cv

   pure subroutine information_criteria(residuals, observations, degrees, &
      values, info)
      ! Evaluate bigtime's covariance-determinant AIC, BIC, and HQ formulas.
      real(dp), intent(in) :: residuals(:, :)
      integer, intent(in) :: observations, degrees
      real(dp), intent(out) :: values(3)
      integer, intent(out) :: info
      real(dp), allocatable :: centered(:, :), covariance(:, :), inverse(:, :)
      real(dp) :: logdet, sample_size
      integer :: residual_count, variables

      residual_count = size(residuals, 1)
      variables = size(residuals, 2)
      if (residual_count < 2 .or. variables < 1 .or. observations < 2 .or. &
         degrees < 0) then
         info = 1
         values = huge(1.0_dp)
         return
      end if
      allocate(centered(residual_count, variables))
      centered = residuals - spread( &
         sum(residuals, dim=1)/real(residual_count, dp), 1, residual_count)
      allocate(covariance(variables, variables), inverse(variables, variables))
      covariance = matmul(transpose(centered), centered)/ &
         real(residual_count - 1, dp)
      call inverse_logdet(covariance, inverse, logdet, info, &
         100.0_dp*epsilon(1.0_dp))
      if (info /= 0) then
         values = huge(1.0_dp)
         return
      end if
      sample_size = real(observations, dp)
      values(bigtime_ic_aic) = logdet + 2.0_dp*real(degrees, dp)/sample_size
      values(bigtime_ic_bic) = logdet + log(sample_size)* &
         real(degrees, dp)/sample_size
      values(bigtime_ic_hq) = logdet + 2.0_dp*log(log(sample_size))* &
         real(degrees, dp)/sample_size
   end subroutine information_criteria

   pure elemental subroutine select_information_criteria(out)
      ! Record the minimizing path index for each information criterion.
      type(bigtime_ic_result_t), intent(inout) :: out
      integer :: criterion

      do criterion = 1, 3
         out%selected(criterion) = minloc(out%value(:, criterion), dim=1)
      end do
      out%info = 0
   end subroutine select_information_criteria

   pure subroutine build_var_data(series, lag_order, response, design)
      ! Form the response and lag-major design matrices used by bigtime.
      real(dp), intent(in) :: series(:, :)
      integer, intent(in) :: lag_order
      real(dp), allocatable, intent(out) :: response(:, :), design(:, :)
      integer :: observations, variables, effective, row, lag

      observations = size(series, 1)
      variables = size(series, 2)
      effective = observations - lag_order
      allocate(response(effective, variables))
      allocate(design(effective, variables*lag_order))
      response = series(lag_order + 1:, :)
      do row = 1, effective
         do lag = 1, lag_order
            design(row, (lag - 1)*variables + 1:lag*variables) = &
               series(lag_order + row - lag, :)
         end do
      end do
   end subroutine build_var_data

   pure subroutine build_varx_data(endogenous, exogenous, ar_order, &
      exogenous_order, response, ar_design, x_design)
      ! Align a VARX response with lag-major endogenous and exogenous designs.
      real(dp), intent(in) :: endogenous(:, :), exogenous(:, :)
      integer, intent(in) :: ar_order, exogenous_order
      real(dp), allocatable, intent(out) :: response(:, :)
      real(dp), allocatable, intent(out) :: ar_design(:, :), x_design(:, :)
      integer :: observations, variables, x_variables, maximum_order
      integer :: effective, row, lag

      observations = size(endogenous, 1)
      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      maximum_order = max(ar_order, exogenous_order)
      effective = observations - maximum_order
      allocate(response(effective, variables))
      allocate(ar_design(effective, variables*ar_order))
      allocate(x_design(effective, x_variables*exogenous_order))
      response = endogenous(maximum_order + 1:, :)
      do row = 1, effective
         do lag = 1, ar_order
            ar_design(row, (lag - 1)*variables + 1:lag*variables) = &
               endogenous(maximum_order + row - lag, :)
         end do
         do lag = 1, exogenous_order
            x_design(row, (lag - 1)*x_variables + 1:lag*x_variables) = &
               exogenous(maximum_order + row - lag, :)
         end do
      end do
   end subroutine build_varx_data

   pure subroutine sparse_var_row(response, design, initial, lambda, penalty, &
      variables, lag_order, step, tolerance, max_iterations, coefficient, &
      iterations, converged)
      ! Run one response row of bigtime's accelerated proximal iteration.
      real(dp), intent(in) :: response(:), design(:, :), initial(:)
      real(dp), intent(in) :: lambda, step, tolerance
      integer, intent(in) :: penalty, variables, lag_order, max_iterations
      real(dp), intent(out) :: coefficient(:)
      integer, intent(out) :: iterations
      logical, intent(out) :: converged
      real(dp) :: previous(size(initial)), older(size(initial))
      real(dp) :: extrapolated(size(initial)), gradient(size(initial))
      real(dp) :: acceleration

      coefficient = initial
      previous = initial
      older = initial
      converged = .false.
      do iterations = 1, max_iterations
         acceleration = real(iterations - 1, dp)/real(iterations + 2, dp)
         extrapolated = previous + acceleration*(previous - older)
         gradient = matmul(response - matmul(design, extrapolated), design)
         coefficient = apply_penalty_prox(extrapolated + step*gradient, &
            step*lambda, penalty, variables, lag_order)
         if (maxval(abs(coefficient - extrapolated)) <= tolerance) then
            converged = .true.
            exit
         end if
         older = previous
         previous = coefficient
      end do
      if (.not. converged) iterations = max_iterations
   end subroutine sparse_var_row

   pure subroutine sparse_varx_row(response, ar_design, x_design, initial_phi, &
      initial_beta, lambda_phi, lambda_beta, penalty, variables, ar_order, &
      x_variables, exogenous_order, step, alpha, tolerance, max_iterations, &
      phi, beta, iterations, converged)
      ! Run one joint endogenous-exogenous proximal iteration.
      real(dp), intent(in) :: response(:), ar_design(:, :), x_design(:, :)
      real(dp), intent(in) :: initial_phi(:), initial_beta(:)
      real(dp), intent(in) :: lambda_phi, lambda_beta, step, alpha, tolerance
      integer, intent(in) :: penalty, variables, ar_order, x_variables
      integer, intent(in) :: exogenous_order, max_iterations
      real(dp), intent(out) :: phi(:), beta(:)
      integer, intent(out) :: iterations
      logical, intent(out) :: converged
      real(dp) :: phi_previous(size(initial_phi)), phi_older(size(initial_phi))
      real(dp) :: beta_previous(size(initial_beta))
      real(dp) :: beta_older(size(initial_beta))
      real(dp) :: phi_extrapolated(size(initial_phi))
      real(dp) :: beta_extrapolated(size(initial_beta))
      real(dp) :: residual(size(response)), acceleration, threshold

      phi = initial_phi
      beta = initial_beta
      phi_previous = initial_phi
      phi_older = initial_phi
      beta_previous = initial_beta
      beta_older = initial_beta
      converged = .false.
      do iterations = 3, max_iterations - 1
         acceleration = real(iterations - 2, dp)/real(iterations + 1, dp)
         phi_extrapolated = phi_previous + &
            acceleration*(phi_previous - phi_older)
         beta_extrapolated = beta_previous + &
            acceleration*(beta_previous - beta_older)
         residual = response - matmul(ar_design, phi_extrapolated) - &
            matmul(x_design, beta_extrapolated)
         phi = apply_penalty_prox(phi_extrapolated + &
            step*matmul(residual, ar_design), step*lambda_phi, penalty, &
            variables, ar_order)/(1.0_dp + alpha)
         beta = apply_penalty_prox(beta_extrapolated + &
            step*matmul(residual, x_design), step*lambda_beta, penalty, &
            x_variables, exogenous_order)/(1.0_dp + alpha)
         threshold = max(maxval(abs(phi - phi_extrapolated)), &
            maxval(abs(beta - beta_extrapolated)))
         if (threshold <= tolerance) then
            converged = .true.
            exit
         end if
         phi_older = phi_previous
         phi_previous = phi
         beta_older = beta_previous
         beta_previous = beta
      end do
      if (.not. converged) iterations = max_iterations - 1
   end subroutine sparse_varx_row

   pure function apply_penalty_prox(values, threshold, penalty, variables, &
      lag_order) result(out)
      ! Dispatch to the lasso or hierarchical-lag proximal operator.
      real(dp), intent(in) :: values(:), threshold
      integer, intent(in) :: penalty, variables, lag_order
      real(dp) :: out(size(values))

      if (penalty == bigtime_penalty_hlag) then
         out = bigtime_hlag_prox(values, threshold, variables, lag_order)
      else
         out = bigtime_soft_threshold(values, threshold)
      end if
   end function apply_penalty_prox

   pure real(dp) function penalty_value(phi, penalty, variables, lag_order) &
      result(value)
      ! Evaluate the lasso or nested hierarchical-lag penalty.
      real(dp), intent(in) :: phi(:, :)
      integer, intent(in) :: penalty, variables, lag_order
      integer :: equation, predictor, first_lag, lag
      real(dp) :: group_norm

      if (penalty == bigtime_penalty_l1) then
         value = sum(abs(phi))
         return
      end if
      value = 0.0_dp
      do equation = 1, size(phi, 1)
         do predictor = 1, variables
            do first_lag = 1, lag_order
               group_norm = 0.0_dp
               do lag = first_lag, lag_order
                  group_norm = group_norm + &
                     phi(equation, predictor + (lag - 1)*variables)**2
               end do
               value = value + sqrt(group_norm)
            end do
         end do
      end do
   end function penalty_value

end module bigtime_mod
