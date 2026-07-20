! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R BigVAR package.
! Structured sparse-VAR algorithms translated from the R BigVAR package.
module bigvar_mod
   use kind_mod, only: dp
   use linalg_mod, only: symmetric_eigen, invert_matrix
   use linalg_mod, only: symmetric_pseudoinverse, inverse_logdet
   use stats_mod, only: normal_quantile
   use random_mod, only: random_multivariate_normal_matrix
   use time_series_var_utils_mod, only: build_var_data, build_varx_data
   use bigtime_mod, only: bigtime_var_fit_t
   use bigtime_mod, only: bigtime_forecast_t
   use bigtime_mod, only: bigtime_simulation_t
   use bigtime_mod, only: bigtime_stability_t, bigtime_var_forecast
   use bigtime_mod, only: bigtime_var_stability, bigtime_soft_threshold
   use bigtime_mod, only: bigtime_var_simulate_from_innovations
   use bigtime_mod, only: bigtime_companion_matrix
   implicit none
   private

   integer, parameter, public :: bigvar_structure_lag = 1
   integer, parameter, public :: bigvar_structure_own_other = 2
   integer, parameter, public :: bigvar_structure_sparse_lag = 3
   integer, parameter, public :: bigvar_structure_sparse_own_other = 4
   integer, parameter, public :: bigvar_structure_hlag_component = 5
   integer, parameter, public :: bigvar_structure_hlag_own_other = 6
   integer, parameter, public :: bigvar_structure_basic_en = 7
   integer, parameter, public :: bigvar_structure_tapered = 8
   integer, parameter, public :: bigvar_structure_mcp = 9
   integer, parameter, public :: bigvar_structure_scad = 10
   integer, parameter, public :: bigvar_structure_efx = 11
   integer, parameter, public :: bigvar_structure_bgr = 12
   integer, parameter, public :: bigvar_structure_hlag_element = 13
   integer, parameter, public :: bigvar_structure_basic = 14
   integer, parameter, public :: bigvar_loss_l1 = 1
   integer, parameter, public :: bigvar_loss_l2 = 2
   integer, parameter, public :: bigvar_loss_huber = 3
   integer, parameter, public :: bigvar_ic_aic = 1
   integer, parameter, public :: bigvar_ic_bic = 2

   type, public :: bigvar_fit_t
      ! Structured VAR coefficients and in-sample fit diagnostics.
      real(dp), allocatable :: phi(:, :)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: lambda_by_response(:)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: alpha = 0.0_dp
      real(dp) :: gamma = 3.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: lag_order = 0
      integer :: structure = 0
      integer :: iterations = 0
      integer :: nonzero = 0
      integer :: active_groups = 0
      integer :: forecast_horizon = 1
      integer :: info = 0
      logical :: converged = .false.
      logical :: relaxed = .false.
      logical :: direct = .false.
      logical :: minnesota = .false.
      logical :: include_intercept = .true.
      real(dp), allocatable :: shrinkage_target(:)
      real(dp) :: refit_fraction = 0.0_dp
   end type bigvar_fit_t

   type, public :: bigvar_path_t
      ! Warm-started structured VAR estimates over a penalty path.
      real(dp), allocatable :: phi(:, :, :)
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: objective(:)
      real(dp), allocatable :: lambda(:)
      real(dp) :: alpha = 0.0_dp
      real(dp) :: gamma = 3.0_dp
      integer, allocatable :: iterations(:)
      integer, allocatable :: nonzero(:)
      integer, allocatable :: active_groups(:)
      logical, allocatable :: converged(:)
      integer :: lag_order = 0
      integer :: forecast_horizon = 1
      integer :: structure = 0
      logical :: direct = .false.
      logical :: minnesota = .false.
      logical :: include_intercept = .true.
      real(dp), allocatable :: shrinkage_target(:)
      integer :: info = 0
   end type bigvar_path_t

   type, public :: bigvar_varx_fit_t
      ! Structured VARX coefficients and in-sample fit diagnostics.
      real(dp), allocatable :: phi(:, :)
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: alpha = 0.0_dp
      real(dp) :: gamma = 3.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: ar_order = 0
      integer :: exogenous_order = 0
      integer :: structure = 0
      integer :: iterations = 0
      integer :: nonzero = 0
      integer :: active_groups = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: relaxed = .false.
      logical :: contemporaneous = .false.
      logical :: transfer_function = .false.
      logical :: minnesota = .false.
      logical :: include_intercept = .true.
      real(dp), allocatable :: shrinkage_target(:)
      real(dp) :: refit_fraction = 0.0_dp
   end type bigvar_varx_fit_t

   type, public :: bigvar_varx_path_t
      ! Warm-started structured VARX estimates over a penalty path.
      real(dp), allocatable :: phi(:, :, :)
      real(dp), allocatable :: beta(:, :, :)
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: objective(:)
      real(dp), allocatable :: lambda(:)
      integer, allocatable :: iterations(:)
      integer, allocatable :: nonzero(:)
      integer, allocatable :: active_groups(:)
      logical, allocatable :: converged(:)
      real(dp) :: alpha = 0.0_dp
      real(dp) :: gamma = 3.0_dp
      integer :: ar_order = 0
      integer :: exogenous_order = 0
      integer :: structure = 0
      integer :: info = 0
      logical :: contemporaneous = .false.
      logical :: transfer_function = .false.
      logical :: minnesota = .false.
      logical :: include_intercept = .true.
      real(dp), allocatable :: shrinkage_target(:)
   end type bigvar_varx_path_t

   type, public :: bigvar_dual_path_t
      ! Structured VAR estimates over alpha-specific lambda paths.
      real(dp), allocatable :: phi(:, :, :, :)
      real(dp), allocatable :: intercept(:, :, :)
      real(dp), allocatable :: objective(:, :)
      real(dp), allocatable :: lambda(:, :)
      real(dp), allocatable :: alpha(:)
      integer, allocatable :: iterations(:, :)
      integer, allocatable :: nonzero(:, :)
      integer, allocatable :: active_groups(:, :)
      logical, allocatable :: converged(:, :)
      integer :: lag_order = 0
      integer :: forecast_horizon = 1
      integer :: structure = 0
      logical :: direct = .false.
      logical :: minnesota = .false.
      logical :: include_intercept = .true.
      real(dp), allocatable :: shrinkage_target(:)
      integer :: info = 0
   end type bigvar_dual_path_t

   type, public :: bigvar_varx_dual_path_t
      ! Structured VARX estimates over alpha-specific lambda paths.
      real(dp), allocatable :: phi(:, :, :, :)
      real(dp), allocatable :: beta(:, :, :, :)
      real(dp), allocatable :: intercept(:, :, :)
      real(dp), allocatable :: objective(:, :)
      real(dp), allocatable :: lambda(:, :)
      real(dp), allocatable :: alpha(:)
      integer, allocatable :: iterations(:, :)
      integer, allocatable :: nonzero(:, :)
      integer, allocatable :: active_groups(:, :)
      logical, allocatable :: converged(:, :)
      integer :: ar_order = 0
      integer :: exogenous_order = 0
      integer :: structure = 0
      logical :: contemporaneous = .false.
      logical :: transfer_function = .false.
      logical :: minnesota = .false.
      logical :: include_intercept = .true.
      real(dp), allocatable :: shrinkage_target(:)
      integer :: info = 0
   end type bigvar_varx_dual_path_t

   type, public :: bigvar_validation_t
      ! Rolling forecast losses, predictions, and selected penalty parameters.
      real(dp), allocatable :: loss(:, :)
      real(dp), allocatable :: forecasts(:, :, :)
      real(dp), allocatable :: mean_loss(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: lambda_grid(:, :)
      real(dp), allocatable :: mean_loss_surface(:, :)
      real(dp), allocatable :: standard_error_surface(:, :)
      integer, allocatable :: nonzero(:, :)
      logical, allocatable :: valid(:, :)
      real(dp) :: selection_standard_error = 0.0_dp
      real(dp) :: selected_lambda = 0.0_dp
      real(dp) :: selected_alpha = 0.0_dp
      integer :: best_index = 0
      integer :: one_se_index = 0
      integer :: selected_index = 0
      integer :: best_lambda_index = 0
      integer :: best_alpha_index = 0
      integer :: one_se_lambda_index = 0
      integer :: one_se_alpha_index = 0
      integer :: selected_lambda_index = 0
      integer :: selected_alpha_index = 0
      integer :: first_origin = 0
      integer :: last_origin = 0
      integer :: horizon = 0
      integer :: loss_type = bigvar_loss_l2
      integer :: info = 0
   end type bigvar_validation_t

   type, public :: bigvar_separate_path_t
      ! Warm-started VAR estimates over response-specific penalty paths.
      real(dp), allocatable :: phi(:, :, :)
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: objective(:)
      real(dp), allocatable :: lambda(:, :)
      integer, allocatable :: iterations(:)
      integer, allocatable :: nonzero(:)
      integer, allocatable :: active_groups(:)
      logical, allocatable :: converged(:)
      real(dp) :: alpha = 0.0_dp
      real(dp) :: gamma = 3.0_dp
      integer :: lag_order = 0
      integer :: forecast_horizon = 1
      integer :: structure = 0
      logical :: direct = .false.
      logical :: minnesota = .false.
      logical :: include_intercept = .true.
      real(dp), allocatable :: shrinkage_target(:)
      integer :: info = 0
   end type bigvar_separate_path_t

   type, public :: bigvar_separate_validation_t
      ! Per-response rolling losses and selected penalty parameters.
      real(dp), allocatable :: loss(:, :, :)
      real(dp), allocatable :: forecasts(:, :, :)
      real(dp), allocatable :: mean_loss(:, :)
      real(dp), allocatable :: standard_error(:, :)
      real(dp), allocatable :: lambda(:, :)
      integer, allocatable :: selected_index(:)
      integer, allocatable :: best_index(:)
      integer, allocatable :: one_se_index(:)
      real(dp), allocatable :: selected_lambda(:)
      logical, allocatable :: valid(:, :, :)
      integer :: first_origin = 0
      integer :: last_origin = 0
      integer :: horizon = 0
      integer :: loss_type = bigvar_loss_l2
      integer :: info = 0
   end type bigvar_separate_validation_t

   type, public :: bigvar_reselection_t
      ! Forecast evaluation with penalty parameters reselected at every origin.
      real(dp), allocatable :: loss(:, :)
      real(dp), allocatable :: forecasts(:, :)
      real(dp), allocatable :: selected_lambda(:, :)
      real(dp), allocatable :: selected_alpha(:)
      integer, allocatable :: selected_index(:, :)
      integer, allocatable :: selected_lambda_index(:)
      integer, allocatable :: selected_alpha_index(:)
      integer, allocatable :: nonzero(:)
      integer, allocatable :: fit_info(:)
      logical, allocatable :: converged(:)
      logical, allocatable :: valid(:)
      integer :: validation_first_origin = 0
      integer :: first_origin = 0
      integer :: last_origin = 0
      integer :: horizon = 0
      integer :: loss_type = bigvar_loss_l2
      integer :: info = 0
   end type bigvar_reselection_t

   type, public :: bigvar_interval_forecast_t
      ! Point forecasts, error covariances, and marginal normal intervals.
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: standard_error(:, :)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
      real(dp), allocatable :: innovation_covariance(:, :)
      real(dp) :: level = 0.95_dp
      integer :: info = 0
   end type bigvar_interval_forecast_t

   type, public :: bigvar_ls_varx_fit_t
      ! Least-squares VARX estimates and information-criterion diagnostics.
      real(dp), allocatable :: phi(:, :)
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: innovation_covariance(:, :)
      real(dp) :: criterion = huge(1.0_dp)
      integer :: ar_order = 0
      integer :: exogenous_order = 0
      integer :: forecast_horizon = 1
      integer :: information_criterion = 0
      integer :: info = 0
      logical :: contemporaneous = .false.
      logical :: direct = .false.
      logical :: transfer_function = .false.
      logical :: include_intercept = .true.
   end type bigvar_ls_varx_fit_t

   type, public :: bigvar_ic_search_t
      ! AIC or BIC surface and its selected least-squares VARX fit.
      type(bigvar_ls_varx_fit_t) :: fit
      real(dp), allocatable :: criterion(:, :)
      integer :: selected_ar_order = 0
      integer :: selected_exogenous_order = 0
      integer :: information_criterion = 0
      integer :: info = 0
   end type bigvar_ic_search_t

   type, public :: bigvar_ic_evaluation_t
      ! Rolling AIC or BIC benchmark forecasts and selected lag orders.
      real(dp), allocatable :: loss(:)
      real(dp), allocatable :: forecasts(:, :)
      integer, allocatable :: ar_order(:)
      integer, allocatable :: exogenous_order(:)
      logical, allocatable :: valid(:)
      integer :: first_origin = 0
      integer :: last_origin = 0
      integer :: horizon = 1
      integer :: information_criterion = 0
      integer :: loss_type = bigvar_loss_l2
      integer :: info = 0
   end type bigvar_ic_evaluation_t

   type, public :: bigvar_benchmark_t
      ! Benchmark forecasts, losses, and their sampling summaries.
      real(dp), allocatable :: loss(:)
      real(dp), allocatable :: forecasts(:, :)
      logical, allocatable :: valid(:)
      real(dp) :: mean_loss = huge(1.0_dp)
      real(dp) :: standard_error = huge(1.0_dp)
      integer :: first_origin = 0
      integer :: last_origin = 0
      integer :: horizon = 1
      integer :: loss_type = bigvar_loss_l2
      integer :: window_size = 0
      integer :: info = 0
   end type bigvar_benchmark_t

   type, public :: bigvar_simulation_t
      ! Simulated VAR observations and BigVAR recursion metadata.
      real(dp), allocatable :: series(:, :)
      real(dp), allocatable :: innovations(:, :)
      real(dp), allocatable :: initial_state(:)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: phi(:, :)
      real(dp), allocatable :: companion(:, :)
      real(dp), allocatable :: innovation_covariance(:, :)
      integer :: burnin = 0
      integer :: info = 0
   end type bigvar_simulation_t

   public :: bigvar_group_prox
   public :: bigvar_efx_prox
   public :: bigvar_mcp_update
   public :: bigvar_scad_update
   public :: bigvar_structured_var
   public :: bigvar_structured_var_path
   public :: bigvar_structured_var_dual_path
   public :: bigvar_structured_var_separate
   public :: bigvar_structured_var_separate_path
   public :: bigvar_lambda_grid
   public :: bigvar_lambda_alpha_grid
   public :: bigvar_separate_lambda_grid
   public :: bigvar_structured_varx
   public :: bigvar_structured_varx_path
   public :: bigvar_structured_varx_dual_path
   public :: bigvar_varx_lambda_grid
   public :: bigvar_varx_lambda_alpha_grid
   public :: bigvar_bgr
   public :: bigvar_bgr_path
   public :: bigvar_bgr_default_grid
   public :: bigvar_forecast
   public :: bigvar_direct_forecast
   public :: bigvar_varx_forecast
   public :: bigvar_innovation_covariance
   public :: bigvar_forecast_covariance
   public :: bigvar_var_forecast_interval
   public :: bigvar_varx_forecast_interval
   public :: bigvar_relaxed_var
   public :: bigvar_relaxed_varx
   public :: bigvar_forecast_loss
   public :: bigvar_var_validate
   public :: bigvar_var_validate_dual
   public :: bigvar_var_validate_separate
   public :: bigvar_varx_validate
   public :: bigvar_varx_validate_dual
   public :: bigvar_var_evaluate
   public :: bigvar_var_reselect
   public :: bigvar_var_reselect_dual
   public :: bigvar_var_reselect_separate
   public :: bigvar_varx_evaluate
   public :: bigvar_varx_reselect
   public :: bigvar_varx_reselect_dual
   public :: bigvar_stability
   public :: bigvar_least_squares_varx
   public :: bigvar_varx_ic_select
   public :: bigvar_ls_varx_forecast
   public :: bigvar_varx_ic_evaluate
   public :: bigvar_var_validate_loo
   public :: bigvar_var_validate_loo_dual
   public :: bigvar_varx_validate_loo
   public :: bigvar_varx_validate_loo_dual
   public :: bigvar_var_validate_separate_loo
   public :: bigvar_mean_benchmark
   public :: bigvar_random_walk_benchmark
   public :: bigvar_var_to_companion
   public :: bigvar_var_simulate_from_innovations
   public :: bigvar_var_simulate

contains

   pure elemental real(dp) function bigvar_mcp_update(score, lambda, gamma, &
      curvature) result(value)
      !! Apply BigVAR's scalar MCP coordinate update.
      real(dp), intent(in) :: score !! Score.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in) :: curvature !! Curvature.
      real(dp) :: magnitude, sign_value

      magnitude = abs(score)
      sign_value = sign(1.0_dp, score)
      if (magnitude <= lambda) then
         value = 0.0_dp
      else if (magnitude <= gamma*lambda*(1.0_dp + lambda)) then
         value = sign_value*(magnitude - lambda)/(curvature* &
            (1.0_dp + lambda - 1.0_dp/gamma))
      else
         value = score/(curvature*(1.0_dp + lambda))
      end if
   end function bigvar_mcp_update

   pure elemental real(dp) function bigvar_scad_update(score, lambda, gamma, &
      curvature) result(value)
      !! Apply BigVAR's scalar SCAD coordinate update.
      real(dp), intent(in) :: score !! Score.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in) :: curvature !! Curvature.
      real(dp) :: magnitude, sign_value

      magnitude = abs(score)
      sign_value = sign(1.0_dp, score)
      if (magnitude <= lambda) then
         value = 0.0_dp
      else if (magnitude <= lambda*(2.0_dp + lambda)) then
         value = sign_value*(magnitude - lambda)/(curvature* &
            (1.0_dp + lambda))
      else if (magnitude <= gamma*lambda*(1.0_dp + lambda)) then
         value = sign_value*(magnitude - gamma*lambda/(gamma - 1.0_dp))/ &
            (curvature*(1.0_dp - 1.0_dp/(gamma - 1.0_dp) + lambda))
      else
         value = score/(curvature*(1.0_dp + lambda))
      end if
   end function bigvar_scad_update

   pure function bigvar_group_prox(values, threshold, structure, variables, &
      lag_order, alpha) result(out)
      !! Apply a BigVAR group or sparse-group proximal operator.
      real(dp), intent(in) :: values(:, :) !! Input values.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp) :: out(size(values, 1), size(values, 2))
      real(dp) :: group_norm, selected_alpha, shrinkage, weight
      integer :: group_structure
      integer :: lag, suffix_lag, first_column, last_column
      integer :: predictor, response

      out = values
      if (threshold <= 0.0_dp) return
      if (variables < 1 .or. lag_order < 1) then
         out = 0.0_dp
         return
      end if
      if (any(shape(values) /= [variables, variables*lag_order])) then
         out = 0.0_dp
         return
      end if

      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
         if (selected_alpha < 0.0_dp .or. selected_alpha > 1.0_dp) then
            out = 0.0_dp
            return
         end if
      end if
      if (is_sparse_structure(structure)) then
         out = bigtime_soft_threshold(out, threshold*selected_alpha)
      end if
      group_structure = base_group_structure(structure)
      if (is_sparse_structure(structure) .and. selected_alpha >= 1.0_dp) return

      select case (group_structure)
      case (bigvar_structure_lag)
         if (structure == bigvar_structure_sparse_lag) then
            weight = real(variables, dp)
         else
            weight = sqrt(real(variables, dp))
         end if
         do lag = 1, lag_order
            first_column = (lag - 1)*variables + 1
            last_column = lag*variables
            group_norm = sqrt(sum(out(:, first_column:last_column)**2))
            call shrink_block(out(:, first_column:last_column), &
               threshold*(1.0_dp - selected_alpha)*weight, group_norm)
         end do
      case (bigvar_structure_own_other)
         weight = sqrt(real(variables, dp))
         do lag = 1, lag_order
            first_column = (lag - 1)*variables
            group_norm = 0.0_dp
            do response = 1, variables
               group_norm = group_norm + &
                  out(response, first_column + response)**2
            end do
            if (group_norm > 0.0_dp) then
               shrinkage = max(1.0_dp - threshold* &
                  (1.0_dp - selected_alpha)*weight/ &
                  sqrt(group_norm), 0.0_dp)
               do response = 1, variables
                  out(response, first_column + response) = shrinkage* &
                     out(response, first_column + response)
               end do
            end if

            if (variables > 1) then
               weight = sqrt(real(variables*(variables - 1), dp))
               group_norm = 0.0_dp
               do response = 1, variables
                  group_norm = group_norm + &
                     sum(out(response, first_column + 1: &
                     first_column + variables)**2) - &
                     out(response, first_column + response)**2
               end do
               if (group_norm > 0.0_dp) then
                  shrinkage = max(1.0_dp - threshold* &
                     (1.0_dp - selected_alpha)*weight/ &
                     sqrt(max(group_norm, 0.0_dp)), 0.0_dp)
                  do response = 1, variables
                     out(response, first_column + 1: &
                        first_column + response - 1) = shrinkage* &
                        out(response, first_column + 1: &
                        first_column + response - 1)
                     out(response, first_column + response + 1: &
                        first_column + variables) = shrinkage* &
                        out(response, first_column + response + 1: &
                        first_column + variables)
                  end do
               else
                  do response = 1, variables
                     out(response, first_column + 1: &
                        first_column + response - 1) = 0.0_dp
                     out(response, first_column + response + 1: &
                        first_column + variables) = 0.0_dp
                  end do
               end if
               weight = sqrt(real(variables, dp))
            end if
         end do
      case (bigvar_structure_hlag_component)
         weight = sqrt(real(variables, dp))
         do response = 1, variables
            do lag = lag_order, 1, -1
               first_column = (lag - 1)*variables + 1
               group_norm = sqrt(sum(out(response, first_column: &
                  variables*lag_order)**2))
               call shrink_vector(out(response, first_column: &
                  variables*lag_order), threshold*weight, group_norm)
            end do
         end do
      case (bigvar_structure_hlag_own_other)
         do response = 1, variables
            do lag = lag_order, 1, -1
               first_column = (lag - 1)*variables + 1
               group_norm = 0.0_dp
               do predictor = first_column, variables*lag_order
                  if (predictor /= first_column + response - 1) then
                     group_norm = group_norm + out(response, predictor)**2
                  end if
               end do
               call shrink_own_other_suffix(out(response, :), threshold* &
                  sqrt(real(max(variables - 1, 0), dp)), group_norm, &
                  first_column, response, variables*lag_order)
               group_norm = sqrt(sum(out(response, first_column: &
                  variables*lag_order)**2))
               call shrink_vector(out(response, first_column: &
                  variables*lag_order), threshold, group_norm)
            end do
         end do
      case (bigvar_structure_hlag_element)
         do response = 1, variables
            do predictor = 1, variables
               do lag = lag_order, 1, -1
                  group_norm = 0.0_dp
                  do suffix_lag = lag, lag_order
                     group_norm = group_norm + out(response, &
                        (suffix_lag - 1)*variables + predictor)**2
                  end do
                  call shrink_lag_suffix(out(response, :), threshold, &
                     group_norm, predictor, lag, lag_order, variables)
               end do
            end do
         end do
      case (bigvar_structure_basic)
         out = bigtime_soft_threshold(out, threshold)
      case (bigvar_structure_basic_en)
         out = bigtime_soft_threshold(out, threshold*selected_alpha)/ &
            (1.0_dp + threshold*(1.0_dp - selected_alpha))
      case (bigvar_structure_tapered)
         do lag = 1, lag_order
            first_column = (lag - 1)*variables + 1
            last_column = lag*variables
            out(:, first_column:last_column) = bigtime_soft_threshold( &
               out(:, first_column:last_column), threshold* &
               real(lag, dp)**selected_alpha)
         end do
      case default
         out = 0.0_dp
      end select
   end function bigvar_group_prox

   pure function bigvar_efx_prox(values, threshold, variables, ar_order, &
      x_variables, exogenous_order) result(out)
      !! Apply BigVAR's endogenous-first nested VARX proximal operator.
      real(dp), intent(in) :: values(:, :) !! Input values.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: x_variables !! X variables.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      real(dp) :: out(size(values, 1), size(values, 2))
      real(dp) :: group_norm, shrinkage
      integer :: response, lag, ar_first, ar_last, x_first, x_last
      integer :: ar_predictors

      out = values
      if (threshold <= 0.0_dp) return
      ar_predictors = variables*ar_order
      if (variables < 1 .or. x_variables < 1 .or. ar_order < 1 .or. &
         exogenous_order < 1 .or. exogenous_order > ar_order .or. &
         any(shape(values) /= [variables, ar_predictors + &
         x_variables*exogenous_order])) then
         out = 0.0_dp
         return
      end if
      do response = 1, variables
         do lag = 1, ar_order
            ar_first = (lag - 1)*variables + 1
            ar_last = lag*variables
            if (lag <= exogenous_order) then
               x_first = ar_predictors + (lag - 1)*x_variables + 1
               x_last = ar_predictors + lag*x_variables
               group_norm = sqrt(sum(out(response, x_first:x_last)**2))
               call shrink_vector(out(response, x_first:x_last), threshold, &
                  group_norm)
               group_norm = sqrt(sum(out(response, ar_first:ar_last)**2) + &
                  sum(out(response, x_first:x_last)**2))
               if (group_norm <= threshold) then
                  out(response, ar_first:ar_last) = 0.0_dp
                  out(response, x_first:x_last) = 0.0_dp
               else
                  shrinkage = 1.0_dp - threshold/group_norm
                  out(response, ar_first:ar_last) = shrinkage* &
                     out(response, ar_first:ar_last)
                  out(response, x_first:x_last) = shrinkage* &
                     out(response, x_first:x_last)
               end if
            else
               group_norm = sqrt(sum(out(response, ar_first:ar_last)**2))
               call shrink_vector(out(response, ar_first:ar_last), threshold, &
                  group_norm)
            end if
         end do
      end do
   end function bigvar_efx_prox

   pure function bigvar_structured_var(series, lag_order, lambda, structure, &
      tolerance, max_iterations, initial_phi, alpha, gamma, refit_fraction, &
      direct_horizon, minnesota_target, include_intercept) result(out)
      !! Estimate a structured sparse VAR with accelerated proximal gradients.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: initial_phi(:, :) !! Initial phi.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      integer, intent(in), optional :: direct_horizon !! Direct horizon.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_fit_t) :: out
      real(dp), allocatable :: response(:, :), original_response(:, :), design(:, :)
      real(dp), allocatable :: centered_response(:, :), centered_design(:, :)
      real(dp), allocatable :: response_mean(:), design_mean(:)
      real(dp), allocatable :: gram(:, :), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: coefficient(:, :), previous(:, :)
      real(dp), allocatable :: extrapolated(:, :), candidate(:, :)
      real(dp), allocatable :: target_matrix(:, :)
      real(dp) :: selected_alpha, selected_gamma, selected_tolerance, step
      real(dp) :: momentum, next_momentum, change
      integer :: selected_max_iterations, variables, predictors, observations
      integer :: selected_horizon
      integer :: eigen_info, iteration

      out%lag_order = lag_order
      selected_horizon = 1
      if (present(direct_horizon)) selected_horizon = direct_horizon
      out%forecast_horizon = selected_horizon
      out%direct = selected_horizon > 1
      out%lambda = lambda
      out%structure = structure
      if (present(include_intercept)) out%include_intercept = include_intercept
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      out%gamma = selected_gamma
      observations = size(series, 1)
      variables = size(series, 2)
      predictors = variables*lag_order
      if (variables < 1 .or. lag_order < 1 .or. &
         observations < lag_order + selected_horizon .or. lambda < 0.0_dp .or. &
         selected_horizon < 1) then
         out%info = 1
         return
      end if
      if (.not. is_supported_structure(structure)) then
         out%info = 2
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables) then
            out%info = 1
            return
         end if
         out%minnesota = .true.
         out%include_intercept = .false.
         out%shrinkage_target = minnesota_target
      end if
      if ((structure == bigvar_structure_mcp .and. selected_gamma <= 1.0_dp) &
         .or. (structure == bigvar_structure_scad .and. &
         selected_gamma <= 2.0_dp)) then
         out%info = 1
         return
      end if
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
         if (selected_alpha < 0.0_dp .or. selected_alpha > 1.0_dp) then
            out%info = 1
            return
         end if
      end if
      out%alpha = selected_alpha
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

      call build_var_direct_data(series, lag_order, selected_horizon, &
         response, design)
      original_response = response
      allocate(target_matrix(variables, predictors))
      target_matrix = 0.0_dp
      if (out%minnesota) then
         call set_minnesota_target(target_matrix, out%shrinkage_target)
         response = response - matmul(design, transpose(target_matrix))
      end if
      allocate(response_mean(variables), design_mean(predictors))
      response_mean = 0.0_dp
      design_mean = 0.0_dp
      if (out%include_intercept) then
         response_mean = sum(response, dim=1)/real(size(response, 1), dp)
         design_mean = sum(design, dim=1)/real(size(design, 1), dp)
      end if
      allocate(centered_response(size(response, 1), variables))
      allocate(centered_design(size(design, 1), predictors))
      centered_response = response - spread(response_mean, 1, size(response, 1))
      centered_design = design - spread(design_mean, 1, size(design, 1))

      allocate(coefficient(variables, predictors))
      allocate(previous(variables, predictors))
      allocate(extrapolated(variables, predictors))
      allocate(candidate(variables, predictors))
      coefficient = 0.0_dp
      if (present(initial_phi)) coefficient = initial_phi - target_matrix
      previous = coefficient
      momentum = 1.0_dp
      out%converged = .false.
      if (is_nonconvex_structure(structure)) then
         call nonconvex_coordinate_descent(centered_response, centered_design, &
            coefficient, lambda, selected_gamma, structure, &
            selected_tolerance, selected_max_iterations, out%iterations, &
            out%converged)
      else
         allocate(gram(predictors, predictors))
         gram = matmul(transpose(centered_design), centered_design)
         call symmetric_eigen(gram, eigenvalues, eigenvectors, eigen_info)
         if (eigen_info /= 0) then
            out%info = 10 + eigen_info
            return
         end if
         if (maxval(abs(eigenvalues)) <= tiny(1.0_dp)) then
            coefficient = 0.0_dp
            out%iterations = 0
            out%converged = .true.
         else
            step = 1.0_dp/maxval(eigenvalues)
            do iteration = 1, selected_max_iterations
               next_momentum = 0.5_dp*(1.0_dp + &
                  sqrt(1.0_dp + 4.0_dp*momentum*momentum))
               extrapolated = coefficient + ((momentum - 1.0_dp)/ &
                  next_momentum)*(coefficient - previous)
               candidate = extrapolated - step*matmul( &
                  matmul(extrapolated, transpose(centered_design)) - &
                  transpose(centered_response), centered_design)
               candidate = bigvar_group_prox(candidate, step*lambda, &
                  structure, variables, lag_order, selected_alpha)
               change = maxval(abs(candidate - coefficient)/ &
                  (1.0_dp + abs(coefficient)))
               previous = coefficient
               coefficient = candidate
               momentum = next_momentum
               if (change <= selected_tolerance) then
                  out%converged = .true.
                  exit
               end if
            end do
            out%iterations = min(iteration, selected_max_iterations)
         end if
      end if

      allocate(out%phi(variables, predictors))
      out%phi = coefficient + target_matrix
      allocate(out%intercept(variables))
      out%intercept = response_mean - matmul(coefficient, design_mean)
      allocate(out%fitted(size(response, 1), variables))
      out%fitted = matmul(design, transpose(out%phi)) + &
         spread(out%intercept, 1, size(response, 1))
      allocate(out%residuals(size(response, 1), variables))
      out%residuals = original_response - out%fitted
      if (is_nonconvex_structure(structure)) then
         out%objective = nonconvex_objective(out%residuals, coefficient, lambda, &
            selected_gamma, structure)
      else
         out%objective = 0.5_dp*sum(out%residuals**2) + lambda* &
            structured_penalty(coefficient, structure, variables, lag_order, &
            selected_alpha, selected_gamma)
      end if
      out%nonzero = count(abs(coefficient) > 100.0_dp*epsilon(1.0_dp))
      out%active_groups = count_active_groups(coefficient, structure, variables, &
         lag_order)
      if (out%converged) then
         out%info = 0
      else
         out%info = 4
      end if
      if (present(refit_fraction)) then
         out = bigvar_relaxed_var(series, out, refit_fraction)
      end if
   end function bigvar_structured_var

   pure function bigvar_structured_var_path(series, lag_order, lambdas, &
      structure, tolerance, max_iterations, warm_start, alpha, gamma, &
      refit_fraction, direct_horizon, minnesota_target, include_intercept) &
      result(out)
      !! Estimate a warm-started path of BigVAR structured models.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      logical, intent(in), optional :: warm_start !! Flag controlling warm start.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      integer, intent(in), optional :: direct_horizon !! Direct horizon.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_path_t) :: out
      type(bigvar_fit_t) :: fit, penalized_fit
      real(dp), allocatable :: starting_phi(:, :)
      real(dp) :: selected_alpha, selected_gamma, selected_tolerance
      integer :: selected_max_iterations, variables, predictors, index
      integer :: selected_horizon
      logical :: use_warm_start

      out%lag_order = lag_order
      selected_horizon = 1
      if (present(direct_horizon)) selected_horizon = direct_horizon
      out%forecast_horizon = selected_horizon
      out%direct = selected_horizon > 1
      out%structure = structure
      if (present(include_intercept)) out%include_intercept = include_intercept
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      out%gamma = selected_gamma
      variables = size(series, 2)
      predictors = variables*lag_order
      if (size(lambdas) < 1 .or. any(lambdas < 0.0_dp) .or. &
         variables < 1 .or. lag_order < 1) then
         out%info = 1
         return
      end if
      if (.not. is_supported_structure(structure)) then
         out%info = 2
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables) then
            out%info = 1
            return
         end if
         out%minnesota = .true.
         out%include_intercept = .false.
         out%shrinkage_target = minnesota_target
      end if
      if ((structure == bigvar_structure_mcp .and. selected_gamma <= 1.0_dp) &
         .or. (structure == bigvar_structure_scad .and. &
         selected_gamma <= 2.0_dp)) then
         out%info = 1
         return
      end if
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
         if (selected_alpha < 0.0_dp .or. selected_alpha > 1.0_dp) then
            out%info = 1
            return
         end if
      end if
      out%alpha = selected_alpha
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      use_warm_start = .true.
      if (present(warm_start)) use_warm_start = warm_start
      allocate(out%phi(variables, predictors, size(lambdas)))
      allocate(out%intercept(variables, size(lambdas)))
      allocate(out%objective(size(lambdas)))
      allocate(out%lambda(size(lambdas)))
      allocate(out%iterations(size(lambdas)))
      allocate(out%nonzero(size(lambdas)))
      allocate(out%active_groups(size(lambdas)))
      allocate(out%converged(size(lambdas)))
      out%lambda = lambdas
      allocate(starting_phi(variables, predictors))
      starting_phi = 0.0_dp
      do index = 1, size(lambdas)
         if (present(minnesota_target)) then
            penalized_fit = bigvar_structured_var(series, lag_order, &
               lambdas(index), structure, selected_tolerance, &
               selected_max_iterations, starting_phi, selected_alpha, &
               selected_gamma, direct_horizon=selected_horizon, &
               minnesota_target=minnesota_target, &
               include_intercept=include_intercept)
         else
            penalized_fit = bigvar_structured_var(series, lag_order, &
               lambdas(index), structure, selected_tolerance, &
               selected_max_iterations, starting_phi, selected_alpha, &
               selected_gamma, direct_horizon=selected_horizon, &
               include_intercept=include_intercept)
         end if
         if (penalized_fit%info /= 0 .and. penalized_fit%info /= 4) then
            out%info = 100*index + penalized_fit%info
            return
         end if
         fit = penalized_fit
         if (present(refit_fraction)) then
            fit = bigvar_relaxed_var(series, penalized_fit, refit_fraction)
         end if
         out%phi(:, :, index) = fit%phi
         out%intercept(:, index) = fit%intercept
         out%objective(index) = fit%objective
         out%iterations(index) = fit%iterations
         out%nonzero(index) = fit%nonzero
         out%active_groups(index) = fit%active_groups
         out%converged(index) = fit%converged
         if (use_warm_start) then
            starting_phi = penalized_fit%phi
         else
            starting_phi = 0.0_dp
         end if
      end do
      if (all(out%converged)) then
         out%info = 0
      else
         out%info = 4
      end if
   end function bigvar_structured_var_path

   pure function bigvar_structured_var_dual_path(series, lag_order, lambdas, &
      alphas, structure, tolerance, max_iterations, warm_start, gamma, &
      refit_fraction, direct_horizon, minnesota_target, include_intercept) &
      result(out)
      !! Estimate VAR paths with a distinct lambda grid for each alpha value.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: direct_horizon !! Direct horizon.
      logical, intent(in), optional :: warm_start !! Flag controlling warm start.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_dual_path_t) :: out
      type(bigvar_path_t) :: path
      integer :: alpha_index, variables, predictors

      variables = size(series, 2)
      predictors = variables*lag_order
      if (size(lambdas, 1) < 1 .or. size(lambdas, 2) /= size(alphas) .or. &
         size(alphas) < 1 .or. any(lambdas < 0.0_dp) .or. &
         any(alphas < 0.0_dp) .or. any(alphas > 1.0_dp) .or. &
         .not. uses_alpha(structure)) then
         out%info = 1
         return
      end if
      allocate(out%phi(variables, predictors, size(lambdas, 1), size(alphas)))
      allocate(out%intercept(variables, size(lambdas, 1), size(alphas)))
      allocate(out%objective(size(lambdas, 1), size(alphas)))
      allocate(out%iterations(size(lambdas, 1), size(alphas)))
      allocate(out%nonzero(size(lambdas, 1), size(alphas)))
      allocate(out%active_groups(size(lambdas, 1), size(alphas)))
      allocate(out%converged(size(lambdas, 1), size(alphas)))
      out%lambda = lambdas
      out%alpha = alphas
      out%lag_order = lag_order
      out%structure = structure
      if (present(include_intercept)) out%include_intercept = include_intercept
      out%forecast_horizon = 1
      if (present(direct_horizon)) out%forecast_horizon = direct_horizon
      out%direct = out%forecast_horizon > 1
      if (present(minnesota_target)) then
         out%minnesota = .true.
         out%include_intercept = .false.
         out%shrinkage_target = minnesota_target
      end if
      out%info = 0
      do alpha_index = 1, size(alphas)
         path = bigvar_structured_var_path(series, lag_order, &
            lambdas(:, alpha_index), structure, tolerance, max_iterations, &
            warm_start, alphas(alpha_index), gamma, refit_fraction, &
            direct_horizon, minnesota_target, include_intercept)
         if (path%info /= 0 .and. path%info /= 4) then
            out%info = path%info
            return
         end if
         out%phi(:, :, :, alpha_index) = path%phi
         out%intercept(:, :, alpha_index) = path%intercept
         out%objective(:, alpha_index) = path%objective
         out%iterations(:, alpha_index) = path%iterations
         out%nonzero(:, alpha_index) = path%nonzero
         out%active_groups(:, alpha_index) = path%active_groups
         out%converged(:, alpha_index) = path%converged
         if (path%info == 4) out%info = 4
      end do
   end function bigvar_structured_var_dual_path

   pure function bigvar_structured_var_separate(series, lag_order, lambdas, &
      structure, tolerance, max_iterations, initial_phi, alpha, gamma, &
      refit_fraction, direct_horizon, minnesota_target, include_intercept) &
      result(out)
      !! Estimate a VAR with a separate penalty for each response equation.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: initial_phi(:, :) !! Initial phi.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      integer, intent(in), optional :: direct_horizon !! Direct horizon.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_fit_t) :: out
      real(dp), allocatable :: response(:, :), original_response(:, :), design(:, :)
      real(dp), allocatable :: centered_response(:, :), centered_design(:, :)
      real(dp), allocatable :: response_mean(:), design_mean(:)
      real(dp), allocatable :: gram(:, :), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: coefficient(:, :), previous(:, :)
      real(dp), allocatable :: extrapolated(:, :), candidate(:, :)
      real(dp), allocatable :: target_matrix(:, :)
      real(dp) :: selected_alpha, selected_gamma, selected_tolerance
      real(dp) :: step, momentum, next_momentum, change
      integer :: variables, predictors, observations, selected_max_iterations
      integer :: selected_horizon
      integer :: eigen_info, iteration

      variables = size(series, 2)
      observations = size(series, 1)
      predictors = variables*lag_order
      out%lag_order = lag_order
      selected_horizon = 1
      if (present(direct_horizon)) selected_horizon = direct_horizon
      out%forecast_horizon = selected_horizon
      out%direct = selected_horizon > 1
      out%structure = structure
      if (present(include_intercept)) out%include_intercept = include_intercept
      out%lambda = 0.0_dp
      if (size(lambdas) > 0) out%lambda = maxval(lambdas)
      out%lambda_by_response = lambdas
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
      end if
      out%alpha = selected_alpha
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      out%gamma = selected_gamma
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      if (variables < 1 .or. lag_order < 1 .or. &
         observations < lag_order + selected_horizon .or. &
         selected_horizon < 1 .or. &
         size(lambdas) /= variables .or. any(lambdas < 0.0_dp) .or. &
         .not. supports_separate_lambdas(structure) .or. &
         selected_alpha < 0.0_dp .or. selected_alpha > 1.0_dp .or. &
         selected_tolerance <= 0.0_dp .or. selected_max_iterations < 1) then
         out%info = 1
         return
      end if
      if ((structure == bigvar_structure_mcp .and. selected_gamma <= 1.0_dp) &
         .or. (structure == bigvar_structure_scad .and. &
         selected_gamma <= 2.0_dp)) then
         out%info = 1
         return
      end if
      if (present(initial_phi)) then
         if (any(shape(initial_phi) /= [variables, predictors])) then
            out%info = 3
            return
         end if
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables) then
            out%info = 1
            return
         end if
         out%minnesota = .true.
         out%include_intercept = .false.
         out%shrinkage_target = minnesota_target
      end if

      call build_var_direct_data(series, lag_order, selected_horizon, &
         response, design)
      original_response = response
      allocate(target_matrix(variables, predictors))
      target_matrix = 0.0_dp
      if (out%minnesota) then
         call set_minnesota_target(target_matrix, out%shrinkage_target)
         response = response - matmul(design, transpose(target_matrix))
      end if
      allocate(response_mean(variables), design_mean(predictors))
      response_mean = 0.0_dp
      design_mean = 0.0_dp
      if (out%include_intercept) then
         response_mean = sum(response, dim=1)/real(size(response, 1), dp)
         design_mean = sum(design, dim=1)/real(size(design, 1), dp)
      end if
      allocate(centered_response(size(response, 1), variables))
      allocate(centered_design(size(design, 1), predictors))
      centered_response = response - spread(response_mean, 1, size(response, 1))
      centered_design = design - spread(design_mean, 1, size(design, 1))
      allocate(coefficient(variables, predictors))
      coefficient = 0.0_dp
      if (present(initial_phi)) coefficient = initial_phi - target_matrix
      out%converged = .false.
      if (is_nonconvex_structure(structure)) then
         call nonconvex_coordinate_descent_separate(centered_response, &
            centered_design, coefficient, lambdas, selected_gamma, structure, &
            selected_tolerance, selected_max_iterations, out%iterations, &
            out%converged)
      else
         allocate(gram(predictors, predictors))
         gram = matmul(transpose(centered_design), centered_design)
         call symmetric_eigen(gram, eigenvalues, eigenvectors, eigen_info)
         if (eigen_info /= 0) then
            out%info = 10 + eigen_info
            return
         end if
         if (maxval(abs(eigenvalues)) <= tiny(1.0_dp)) then
            coefficient = 0.0_dp
            out%converged = .true.
         else
            allocate(previous(variables, predictors))
            allocate(extrapolated(variables, predictors))
            allocate(candidate(variables, predictors))
            previous = coefficient
            momentum = 1.0_dp
            step = 1.0_dp/maxval(eigenvalues)
            do iteration = 1, selected_max_iterations
               next_momentum = 0.5_dp*(1.0_dp + &
                  sqrt(1.0_dp + 4.0_dp*momentum*momentum))
               extrapolated = coefficient + ((momentum - 1.0_dp)/ &
                  next_momentum)*(coefficient - previous)
               candidate = extrapolated - step*matmul( &
                  matmul(extrapolated, transpose(centered_design)) - &
                  transpose(centered_response), centered_design)
               candidate = separate_response_prox(candidate, step*lambdas, &
                  structure, variables, lag_order, selected_alpha)
               change = maxval(abs(candidate - coefficient)/ &
                  (1.0_dp + abs(coefficient)))
               previous = coefficient
               coefficient = candidate
               momentum = next_momentum
               if (change <= selected_tolerance) then
                  out%converged = .true.
                  exit
               end if
            end do
            out%iterations = min(iteration, selected_max_iterations)
         end if
      end if
      out%phi = coefficient + target_matrix
      out%intercept = response_mean - matmul(coefficient, design_mean)
      out%fitted = matmul(design, transpose(out%phi)) + &
         spread(out%intercept, 1, size(response, 1))
      out%residuals = original_response - out%fitted
      out%objective = separate_penalized_objective(out%residuals, coefficient, &
         lambdas, structure, variables, lag_order, selected_alpha, &
         selected_gamma)
      out%nonzero = count(abs(coefficient) > 100.0_dp*epsilon(1.0_dp))
      out%active_groups = count_active_groups(coefficient, structure, variables, &
         lag_order)
      if (out%converged) then
         out%info = 0
      else
         out%info = 4
      end if
      if (present(refit_fraction)) then
         out = bigvar_relaxed_var(series, out, refit_fraction)
      end if
   end function bigvar_structured_var_separate

   pure function bigvar_structured_var_separate_path(series, lag_order, &
      lambdas, structure, tolerance, max_iterations, warm_start, alpha, gamma, &
      refit_fraction, direct_horizon, minnesota_target, include_intercept) &
      result(out)
      !! Estimate a warm-started path of response-specific VAR penalties.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      integer, intent(in), optional :: direct_horizon !! Direct horizon.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      logical, intent(in), optional :: warm_start !! Flag controlling warm start.
      type(bigvar_separate_path_t) :: out
      type(bigvar_fit_t) :: fit, penalized_fit
      real(dp), allocatable :: starting_phi(:, :)
      real(dp) :: selected_alpha, selected_gamma, selected_tolerance
      integer :: variables, predictors, index, selected_max_iterations
      integer :: selected_horizon
      logical :: use_warm_start

      variables = size(series, 2)
      predictors = variables*lag_order
      out%lag_order = lag_order
      selected_horizon = 1
      if (present(direct_horizon)) selected_horizon = direct_horizon
      out%forecast_horizon = selected_horizon
      out%direct = selected_horizon > 1
      out%structure = structure
      if (present(include_intercept)) out%include_intercept = include_intercept
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
      end if
      out%alpha = selected_alpha
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      out%gamma = selected_gamma
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      use_warm_start = .true.
      if (present(warm_start)) use_warm_start = warm_start
      if (size(lambdas, 1) < 1 .or. size(lambdas, 2) /= variables .or. &
         any(lambdas < 0.0_dp) .or. &
         .not. supports_separate_lambdas(structure)) then
         out%info = 1
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables) then
            out%info = 1
            return
         end if
         out%minnesota = .true.
         out%include_intercept = .false.
         out%shrinkage_target = minnesota_target
      end if
      allocate(out%phi(variables, predictors, size(lambdas, 1)))
      allocate(out%intercept(variables, size(lambdas, 1)))
      allocate(out%objective(size(lambdas, 1)))
      allocate(out%lambda(size(lambdas, 1), variables))
      allocate(out%iterations(size(lambdas, 1)))
      allocate(out%nonzero(size(lambdas, 1)))
      allocate(out%active_groups(size(lambdas, 1)))
      allocate(out%converged(size(lambdas, 1)))
      out%lambda = lambdas
      allocate(starting_phi(variables, predictors))
      starting_phi = 0.0_dp
      do index = 1, size(lambdas, 1)
         if (present(minnesota_target)) then
            penalized_fit = bigvar_structured_var_separate(series, lag_order, &
               lambdas(index, :), structure, selected_tolerance, &
               selected_max_iterations, starting_phi, selected_alpha, &
               selected_gamma, direct_horizon=selected_horizon, &
               minnesota_target=minnesota_target, &
               include_intercept=include_intercept)
         else
            penalized_fit = bigvar_structured_var_separate(series, lag_order, &
               lambdas(index, :), structure, selected_tolerance, &
               selected_max_iterations, starting_phi, selected_alpha, &
               selected_gamma, direct_horizon=selected_horizon, &
               include_intercept=include_intercept)
         end if
         if (penalized_fit%info /= 0 .and. penalized_fit%info /= 4) then
            out%info = 100*index + penalized_fit%info
            return
         end if
         fit = penalized_fit
         if (present(refit_fraction)) then
            fit = bigvar_relaxed_var(series, penalized_fit, refit_fraction)
         end if
         out%phi(:, :, index) = fit%phi
         out%intercept(:, index) = fit%intercept
         out%objective(index) = fit%objective
         out%iterations(index) = fit%iterations
         out%nonzero(index) = fit%nonzero
         out%active_groups(index) = fit%active_groups
         out%converged(index) = fit%converged
         if (use_warm_start) then
            starting_phi = penalized_fit%phi
         else
            starting_phi = 0.0_dp
         end if
      end do
      if (all(out%converged)) then
         out%info = 0
      else
         out%info = 4
      end if
   end function bigvar_structured_var_separate_path

   pure function bigvar_separate_lambda_grid(series, lag_order, structure, &
      grid_ratio, grid_size, alpha, direct_horizon, minnesota_target, linear, &
      include_intercept) result(lambdas)
      !! Construct descending zero-model penalty grids for each response.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: grid_ratio !! Grid ratio.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: grid_size !! Grid size.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      integer, intent(in), optional :: direct_horizon !! Direct horizon.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: linear !! Flag controlling linear.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      real(dp), allocatable :: lambdas(:, :)
      real(dp), allocatable :: response(:, :)
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: centered_response(:, :), centered_design(:, :)
      real(dp), allocatable :: gradient(:, :), response_mean(:), design_mean(:)
      real(dp), allocatable :: row_gradient(:, :)
      real(dp), allocatable :: target_matrix(:, :)
      real(dp) :: maximum_lambda, selected_alpha
      integer :: variables, response_index, selected_horizon
      logical :: use_linear, use_intercept

      variables = size(series, 2)
      selected_horizon = 1
      if (present(direct_horizon)) selected_horizon = direct_horizon
      use_linear = .false.
      if (present(linear)) use_linear = linear
      use_intercept = .true.
      if (present(include_intercept)) use_intercept = include_intercept
      if (present(minnesota_target)) use_intercept = .false.
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
      end if
      if (variables < 1 .or. lag_order < 1 .or. &
         size(series, 1) < lag_order + selected_horizon .or. &
         selected_horizon < 1 .or. grid_size < 1 .or. &
         grid_ratio < 1.0_dp .or. .not. supports_separate_lambdas(structure) &
         .or. selected_alpha < 0.0_dp .or. selected_alpha > 1.0_dp) then
         allocate(lambdas(0, 0))
         return
      end if
      call build_var_direct_data(series, lag_order, selected_horizon, &
         response, design)
      allocate(response_mean(variables), design_mean(variables*lag_order))
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables) then
            allocate(lambdas(0, 0))
            return
         end if
         allocate(target_matrix(variables, variables*lag_order))
         target_matrix = 0.0_dp
         call set_minnesota_target(target_matrix, minnesota_target)
         response = response - matmul(design, transpose(target_matrix))
         response_mean = 0.0_dp
         design_mean = 0.0_dp
      else if (use_intercept) then
         response_mean = sum(response, dim=1)/real(size(response, 1), dp)
         design_mean = sum(design, dim=1)/real(size(design, 1), dp)
      end if
      centered_response = response - spread(response_mean, 1, size(response, 1))
      centered_design = design - spread(design_mean, 1, size(design, 1))
      gradient = matmul(transpose(centered_response), centered_design)
      allocate(lambdas(grid_size, variables))
      allocate(row_gradient(variables, variables*lag_order))
      do response_index = 1, variables
         row_gradient = 0.0_dp
         row_gradient(response_index, :) = gradient(response_index, :)
         if (is_nonconvex_structure(structure)) then
            maximum_lambda = maxval(abs(gradient(response_index, :)))/ &
               real(size(response, 1), dp)
         else
            maximum_lambda = zero_model_bound(row_gradient, structure, &
               variables, lag_order, selected_alpha)
         end if
         maximum_lambda = max(maximum_lambda*(1.0_dp + 1.0e-10_dp), &
            tiny(1.0_dp))
         call fill_lambda_grid(lambdas(:, response_index), maximum_lambda, &
            grid_ratio, use_linear)
      end do
   end function bigvar_separate_lambda_grid

   pure function bigvar_lambda_grid(series, lag_order, structure, grid_ratio, &
      grid_size, alpha, direct_horizon, minnesota_target, linear, &
      include_intercept) result(lambdas)
      !! Construct a descending grid whose first value gives the zero model.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: grid_ratio !! Grid ratio.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: grid_size !! Grid size.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      integer, intent(in), optional :: direct_horizon !! Direct horizon.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: linear !! Flag controlling linear.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      real(dp), allocatable :: lambdas(:)
      real(dp), allocatable :: response(:, :), design(:, :)
      real(dp), allocatable :: centered_response(:, :), centered_design(:, :)
      real(dp), allocatable :: gradient(:, :), response_mean(:), design_mean(:)
      real(dp), allocatable :: target_matrix(:, :)
      real(dp) :: maximum_lambda, selected_alpha
      integer :: variables, selected_horizon
      logical :: use_linear, use_intercept

      selected_horizon = 1
      if (present(direct_horizon)) selected_horizon = direct_horizon
      use_linear = .false.
      if (present(linear)) use_linear = linear
      use_intercept = .true.
      if (present(include_intercept)) use_intercept = include_intercept
      if (present(minnesota_target)) use_intercept = .false.
      if (size(series, 2) < 1 .or. &
         size(series, 1) < lag_order + selected_horizon .or. &
         lag_order < 1 .or. grid_size < 1 .or. grid_ratio < 1.0_dp .or. &
         selected_horizon < 1 .or. &
         .not. is_supported_structure(structure)) then
         allocate(lambdas(0))
         return
      end if
      variables = size(series, 2)
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables) then
            allocate(lambdas(0))
            return
         end if
      end if
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
         if (selected_alpha < 0.0_dp .or. selected_alpha > 1.0_dp) then
            allocate(lambdas(0))
            return
         end if
      end if
      call build_var_direct_data(series, lag_order, selected_horizon, &
         response, design)
      allocate(response_mean(variables), design_mean(variables*lag_order))
      response_mean = 0.0_dp
      design_mean = 0.0_dp
      if (present(minnesota_target)) then
         allocate(target_matrix(variables, variables*lag_order))
         target_matrix = 0.0_dp
         call set_minnesota_target(target_matrix, minnesota_target)
         response = response - matmul(design, transpose(target_matrix))
      else if (use_intercept) then
         response_mean = sum(response, dim=1)/real(size(response, 1), dp)
         design_mean = sum(design, dim=1)/real(size(design, 1), dp)
      end if
      allocate(centered_response(size(response, 1), variables))
      allocate(centered_design(size(design, 1), variables*lag_order))
      centered_response = response - spread(response_mean, 1, size(response, 1))
      centered_design = design - spread(design_mean, 1, size(design, 1))
      allocate(gradient(variables, variables*lag_order))
      gradient = matmul(transpose(centered_response), centered_design)
      if (is_nonconvex_structure(structure)) then
         maximum_lambda = maxval(abs(gradient))/real(size(response, 1), dp)
      else
         maximum_lambda = zero_model_bound(gradient, structure, variables, &
            lag_order, selected_alpha)
      end if
      maximum_lambda = max(maximum_lambda*(1.0_dp + 1.0e-10_dp), &
         tiny(1.0_dp))
      allocate(lambdas(grid_size))
      call fill_lambda_grid(lambdas, maximum_lambda, grid_ratio, use_linear)
   end function bigvar_lambda_grid

   pure function bigvar_lambda_alpha_grid(series, lag_order, structure, &
      grid_ratio, grid_size, alphas, direct_horizon, minnesota_target, linear, &
      include_intercept) result(lambdas)
      !! Construct an alpha-specific matrix of descending VAR lambda grids.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: grid_ratio !! Grid ratio.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: grid_size !! Grid size.
      integer, intent(in), optional :: direct_horizon !! Direct horizon.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: linear !! Flag controlling linear.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      real(dp), allocatable :: lambdas(:, :)
      real(dp), allocatable :: column(:)
      integer :: alpha_index

      if (size(alphas) < 1 .or. .not. uses_alpha(structure) .or. &
         any(alphas < 0.0_dp) .or. any(alphas > 1.0_dp)) then
         allocate(lambdas(0, 0))
         return
      end if
      allocate(lambdas(grid_size, size(alphas)))
      do alpha_index = 1, size(alphas)
         column = bigvar_lambda_grid(series, lag_order, structure, grid_ratio, &
            grid_size, alphas(alpha_index), direct_horizon, minnesota_target, &
            linear, include_intercept)
         if (size(column) /= grid_size) then
            deallocate(lambdas)
            allocate(lambdas(0, 0))
            return
         end if
         lambdas(:, alpha_index) = column
      end do
   end function bigvar_lambda_alpha_grid

   pure function bigvar_structured_varx(endogenous, exogenous, ar_order, &
      exogenous_order, lambda, structure, tolerance, max_iterations, &
      initial_phi, initial_beta, alpha, gamma, refit_fraction, &
      contemporaneous, minnesota_target, include_intercept) result(out)
      !! Estimate a BigVAR structured VARX with joint endogenous and exogenous updates.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: initial_phi(:, :) !! Initial phi.
      real(dp), intent(in), optional :: initial_beta(:, :) !! Initial beta.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_varx_fit_t) :: out
      real(dp), allocatable :: response(:, :), ar_design(:, :), x_design(:, :)
      real(dp), allocatable :: centered_response(:, :), combined(:, :)
      real(dp), allocatable :: response_mean(:), design_mean(:)
      real(dp), allocatable :: gram(:, :), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: coefficient(:, :), previous(:, :)
      real(dp), allocatable :: extrapolated(:, :), candidate(:, :)
      real(dp), allocatable :: target_matrix(:, :), original_response(:, :)
      real(dp) :: selected_alpha, selected_gamma, selected_tolerance
      real(dp) :: step, momentum, next_momentum, change
      integer :: selected_max_iterations, observations, variables, x_variables
      integer :: ar_predictors, x_predictors, eigen_info, iteration
      integer :: x_blocks
      logical :: use_contemporaneous

      out%ar_order = ar_order
      out%exogenous_order = exogenous_order
      out%transfer_function = ar_order == 0
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      out%contemporaneous = use_contemporaneous
      out%lambda = lambda
      out%structure = structure
      if (present(include_intercept)) out%include_intercept = include_intercept
      observations = size(endogenous, 1)
      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      ar_predictors = variables*ar_order
      x_blocks = exogenous_order
      if (use_contemporaneous) x_blocks = x_blocks + 1
      x_predictors = x_variables*x_blocks
      if (observations /= size(exogenous, 1) .or. variables < 1 .or. &
         x_variables < 1 .or. ar_order < 0 .or. exogenous_order < 0 .or. &
         (.not. use_contemporaneous .and. exogenous_order < 1) .or. &
         observations <= max(ar_order, exogenous_order) .or. lambda < 0.0_dp) then
         out%info = 1
         return
      end if
      if (ar_order == 0 .and. .not. supports_transfer_function(structure)) then
         out%info = 1
         return
      end if
      if (structure == bigvar_structure_efx .and. &
         (exogenous_order > ar_order .or. use_contemporaneous)) then
         out%info = 1
         return
      end if
      if (.not. is_varx_structure(structure)) then
         out%info = 2
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables .or. ar_order == 0) then
            out%info = 1
            return
         end if
         out%minnesota = .true.
         out%include_intercept = .false.
         out%shrinkage_target = minnesota_target
      end if
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
         if (selected_alpha < 0.0_dp .or. selected_alpha > 1.0_dp) then
            out%info = 1
            return
         end if
      end if
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      if ((structure == bigvar_structure_mcp .and. selected_gamma <= 1.0_dp) &
         .or. (structure == bigvar_structure_scad .and. &
         selected_gamma <= 2.0_dp)) then
         out%info = 1
         return
      end if
      out%alpha = selected_alpha
      out%gamma = selected_gamma
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
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      if (selected_tolerance <= 0.0_dp .or. selected_max_iterations < 1) then
         out%info = 1
         return
      end if

      call build_varx_data(endogenous, exogenous, ar_order, exogenous_order, &
         response, ar_design, x_design, use_contemporaneous)
      original_response = response
      allocate(combined(size(response, 1), ar_predictors + x_predictors))
      if (ar_predictors > 0) combined(:, :ar_predictors) = ar_design
      combined(:, ar_predictors + 1:) = x_design
      allocate(target_matrix(variables, size(combined, 2)))
      target_matrix = 0.0_dp
      if (out%minnesota) then
         call set_minnesota_target(target_matrix, out%shrinkage_target)
         response = response - matmul(combined, transpose(target_matrix))
      end if
      allocate(response_mean(variables), design_mean(size(combined, 2)))
      response_mean = 0.0_dp
      design_mean = 0.0_dp
      if (out%include_intercept) then
         response_mean = sum(response, dim=1)/real(size(response, 1), dp)
         design_mean = sum(combined, dim=1)/real(size(combined, 1), dp)
      end if
      allocate(centered_response(size(response, 1), variables))
      centered_response = response - spread(response_mean, 1, size(response, 1))
      combined = combined - spread(design_mean, 1, size(combined, 1))

      allocate(coefficient(variables, size(combined, 2)))
      coefficient = 0.0_dp
      if (present(initial_phi) .and. ar_predictors > 0) then
         coefficient(:, :ar_predictors) = initial_phi - &
            target_matrix(:, :ar_predictors)
      end if
      if (present(initial_beta)) coefficient(:, ar_predictors + 1:) = &
         initial_beta
      out%converged = .false.
      if (is_nonconvex_structure(structure)) then
         call nonconvex_coordinate_descent(centered_response, combined, &
            coefficient, lambda, selected_gamma, structure, &
            selected_tolerance, selected_max_iterations, out%iterations, &
            out%converged)
      else
         allocate(gram(size(combined, 2), size(combined, 2)))
         gram = matmul(transpose(combined), combined)
         call symmetric_eigen(gram, eigenvalues, eigenvectors, eigen_info)
         if (eigen_info /= 0) then
            out%info = 10 + eigen_info
            return
         end if
         if (maxval(abs(eigenvalues)) <= tiny(1.0_dp)) then
            coefficient = 0.0_dp
            out%converged = .true.
         else
            allocate(previous(variables, size(combined, 2)))
            allocate(extrapolated(variables, size(combined, 2)))
            allocate(candidate(variables, size(combined, 2)))
            previous = coefficient
            momentum = 1.0_dp
            step = 1.0_dp/maxval(eigenvalues)
            do iteration = 1, selected_max_iterations
               next_momentum = 0.5_dp*(1.0_dp + &
                  sqrt(1.0_dp + 4.0_dp*momentum*momentum))
               if (structure == bigvar_structure_efx) then
                  extrapolated = coefficient
               else
                  extrapolated = coefficient + ((momentum - 1.0_dp)/ &
                     next_momentum)*(coefficient - previous)
               end if
               candidate = extrapolated - step*matmul( &
                  matmul(extrapolated, transpose(combined)) - &
                  transpose(centered_response), combined)
               candidate = varx_prox(candidate, step*lambda, structure, &
                  variables, ar_order, x_variables, x_blocks, &
                  selected_alpha)
               change = maxval(abs(candidate - coefficient)/ &
                  (1.0_dp + abs(coefficient)))
               previous = coefficient
               coefficient = candidate
               momentum = next_momentum
               if (change <= selected_tolerance) then
                  out%converged = .true.
                  exit
               end if
            end do
            out%iterations = min(iteration, selected_max_iterations)
         end if
      end if

      allocate(out%phi(variables, ar_predictors))
      allocate(out%beta(variables, x_predictors))
      if (ar_predictors > 0) then
         out%phi = coefficient(:, :ar_predictors) + &
            target_matrix(:, :ar_predictors)
      end if
      out%beta = coefficient(:, ar_predictors + 1:)
      allocate(out%intercept(variables))
      out%intercept = response_mean - matmul(coefficient, design_mean)
      allocate(out%fitted(size(response, 1), variables))
      out%fitted = matmul(combined + spread(design_mean, 1, &
         size(combined, 1)), transpose(coefficient + target_matrix)) + &
         spread(out%intercept, 1, size(response, 1))
      allocate(out%residuals(size(response, 1), variables))
      out%residuals = original_response - out%fitted
      if (is_nonconvex_structure(structure)) then
         out%objective = nonconvex_objective(out%residuals, coefficient, &
            lambda, selected_gamma, structure)
      else
         out%objective = 0.5_dp*sum(out%residuals**2) + lambda* &
            varx_penalty(coefficient(:, :ar_predictors), out%beta, structure, &
            variables, ar_order, &
            x_variables, x_blocks, selected_alpha)
      end if
      out%nonzero = count(abs(coefficient) > 100.0_dp*epsilon(1.0_dp))
      if (structure == bigvar_structure_efx) then
         out%active_groups = count_active_efx_groups( &
            coefficient(:, :ar_predictors), out%beta, &
            variables, ar_order, x_variables, x_blocks)
      else
         out%active_groups = count_active_x_groups(out%beta)
         if (ar_predictors > 0) out%active_groups = out%active_groups + &
            count_active_groups(coefficient(:, :ar_predictors), structure, &
            variables, ar_order)
      end if
      if (out%converged) then
         out%info = 0
      else
         out%info = 4
      end if
      if (present(refit_fraction)) then
         out = bigvar_relaxed_varx(endogenous, exogenous, out, refit_fraction)
      end if
   end function bigvar_structured_varx

   pure function bigvar_structured_varx_path(endogenous, exogenous, ar_order, &
      exogenous_order, lambdas, structure, tolerance, max_iterations, &
      warm_start, alpha, gamma, refit_fraction, contemporaneous, &
      minnesota_target, include_intercept) result(out)
      !! Estimate a warm-started path of BigVAR structured VARX models.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      logical, intent(in), optional :: warm_start !! Flag controlling warm start.
      type(bigvar_varx_path_t) :: out
      type(bigvar_varx_fit_t) :: fit, penalized_fit
      real(dp), allocatable :: starting_phi(:, :), starting_beta(:, :)
      real(dp) :: selected_alpha, selected_gamma, selected_tolerance
      integer :: selected_max_iterations, variables, x_variables, index
      integer :: x_blocks
      logical :: use_warm_start, use_contemporaneous

      out%ar_order = ar_order
      out%exogenous_order = exogenous_order
      out%transfer_function = ar_order == 0
      out%structure = structure
      if (present(include_intercept)) out%include_intercept = include_intercept
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      out%contemporaneous = use_contemporaneous
      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      if (ar_order < 0 .or. &
         (ar_order == 0 .and. .not. supports_transfer_function(structure))) then
         out%info = 1
         return
      end if
      if (size(lambdas) < 1 .or. any(lambdas < 0.0_dp) .or. &
         .not. is_varx_structure(structure) .or. variables < 1 .or. &
         x_variables < 1 .or. ar_order < 0 .or. exogenous_order < 0 .or. &
         (.not. use_contemporaneous .and. exogenous_order < 1) .or. &
         (structure == bigvar_structure_efx .and. &
         (exogenous_order > ar_order .or. use_contemporaneous))) then
         out%info = 1
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables .or. ar_order == 0) then
            out%info = 1
            return
         end if
         out%minnesota = .true.
         out%include_intercept = .false.
         out%shrinkage_target = minnesota_target
      end if
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
      end if
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      use_warm_start = .true.
      if (present(warm_start)) use_warm_start = warm_start
      out%alpha = selected_alpha
      out%gamma = selected_gamma
      allocate(out%phi(variables, variables*ar_order, size(lambdas)))
      x_blocks = exogenous_order
      if (use_contemporaneous) x_blocks = x_blocks + 1
      allocate(out%beta(variables, x_variables*x_blocks, size(lambdas)))
      allocate(out%intercept(variables, size(lambdas)))
      allocate(out%objective(size(lambdas)), out%lambda(size(lambdas)))
      allocate(out%iterations(size(lambdas)), out%nonzero(size(lambdas)))
      allocate(out%active_groups(size(lambdas)), out%converged(size(lambdas)))
      out%lambda = lambdas
      allocate(starting_phi(variables, variables*ar_order))
      allocate(starting_beta(variables, x_variables*x_blocks))
      starting_phi = 0.0_dp
      starting_beta = 0.0_dp
      do index = 1, size(lambdas)
         if (present(minnesota_target)) then
            penalized_fit = bigvar_structured_varx(endogenous, exogenous, &
               ar_order, exogenous_order, lambdas(index), structure, &
               selected_tolerance, selected_max_iterations, starting_phi, &
               starting_beta, selected_alpha, selected_gamma, &
               contemporaneous=use_contemporaneous, &
               minnesota_target=minnesota_target, &
               include_intercept=include_intercept)
         else
            penalized_fit = bigvar_structured_varx(endogenous, exogenous, &
               ar_order, exogenous_order, lambdas(index), structure, &
               selected_tolerance, selected_max_iterations, starting_phi, &
               starting_beta, selected_alpha, selected_gamma, &
               contemporaneous=use_contemporaneous, &
               include_intercept=include_intercept)
         end if
         if (penalized_fit%info /= 0 .and. penalized_fit%info /= 4) then
            out%info = 100*index + penalized_fit%info
            return
         end if
         fit = penalized_fit
         if (present(refit_fraction)) then
            fit = bigvar_relaxed_varx(endogenous, exogenous, penalized_fit, &
               refit_fraction)
         end if
         out%phi(:, :, index) = fit%phi
         out%beta(:, :, index) = fit%beta
         out%intercept(:, index) = fit%intercept
         out%objective(index) = fit%objective
         out%iterations(index) = fit%iterations
         out%nonzero(index) = fit%nonzero
         out%active_groups(index) = fit%active_groups
         out%converged(index) = fit%converged
         if (use_warm_start) then
            starting_phi = penalized_fit%phi
            starting_beta = penalized_fit%beta
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
   end function bigvar_structured_varx_path

   pure function bigvar_structured_varx_dual_path(endogenous, exogenous, &
      ar_order, exogenous_order, lambdas, alphas, structure, tolerance, &
      max_iterations, warm_start, gamma, refit_fraction, contemporaneous, &
      minnesota_target, include_intercept) result(out)
      !! Estimate VARX paths with a distinct lambda grid for each alpha value.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      logical, intent(in), optional :: warm_start !! Flag controlling warm start.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_varx_dual_path_t) :: out
      type(bigvar_varx_path_t) :: path
      integer :: alpha_index, variables, x_variables, x_blocks

      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      if (ar_order < 0 .or. &
         (ar_order == 0 .and. .not. supports_transfer_function(structure))) then
         out%info = 1
         return
      end if
      x_blocks = exogenous_order
      if (present(contemporaneous)) then
         if (contemporaneous) x_blocks = x_blocks + 1
      end if
      if (size(lambdas, 1) < 1 .or. size(lambdas, 2) /= size(alphas) .or. &
         size(alphas) < 1 .or. any(lambdas < 0.0_dp) .or. &
         any(alphas < 0.0_dp) .or. any(alphas > 1.0_dp) .or. &
         .not. uses_alpha(structure)) then
         out%info = 1
         return
      end if
      if (present(minnesota_target)) then
         if (ar_order == 0 .or. size(minnesota_target) /= variables) then
            out%info = 1
            return
         end if
      end if
      allocate(out%phi(variables, variables*ar_order, size(lambdas, 1), &
         size(alphas)))
      allocate(out%beta(variables, x_variables*x_blocks, size(lambdas, 1), &
         size(alphas)))
      allocate(out%intercept(variables, size(lambdas, 1), size(alphas)))
      allocate(out%objective(size(lambdas, 1), size(alphas)))
      allocate(out%iterations(size(lambdas, 1), size(alphas)))
      allocate(out%nonzero(size(lambdas, 1), size(alphas)))
      allocate(out%active_groups(size(lambdas, 1), size(alphas)))
      allocate(out%converged(size(lambdas, 1), size(alphas)))
      out%lambda = lambdas
      out%alpha = alphas
      out%ar_order = ar_order
      out%exogenous_order = exogenous_order
      out%transfer_function = ar_order == 0
      out%structure = structure
      if (present(include_intercept)) out%include_intercept = include_intercept
      if (present(contemporaneous)) out%contemporaneous = contemporaneous
      if (present(minnesota_target)) then
         out%minnesota = .true.
         out%include_intercept = .false.
         out%shrinkage_target = minnesota_target
      end if
      out%info = 0
      do alpha_index = 1, size(alphas)
         path = bigvar_structured_varx_path(endogenous, exogenous, ar_order, &
            exogenous_order, lambdas(:, alpha_index), structure, tolerance, &
            max_iterations, warm_start, alphas(alpha_index), gamma, &
            refit_fraction, contemporaneous, minnesota_target, &
            include_intercept)
         if (path%info /= 0 .and. path%info /= 4) then
            out%info = path%info
            return
         end if
         out%phi(:, :, :, alpha_index) = path%phi
         out%beta(:, :, :, alpha_index) = path%beta
         out%intercept(:, :, alpha_index) = path%intercept
         out%objective(:, alpha_index) = path%objective
         out%iterations(:, alpha_index) = path%iterations
         out%nonzero(:, alpha_index) = path%nonzero
         out%active_groups(:, alpha_index) = path%active_groups
         out%converged(:, alpha_index) = path%converged
         if (path%info == 4) out%info = 4
      end do
   end function bigvar_structured_varx_dual_path

   pure function bigvar_varx_lambda_grid(endogenous, exogenous, ar_order, &
      exogenous_order, structure, grid_ratio, grid_size, alpha, &
      contemporaneous, minnesota_target, linear, include_intercept) &
      result(lambdas)
      !! Construct a descending BigVAR VARX grid from the joint zero-model bound.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: grid_ratio !! Grid ratio.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: grid_size !! Grid size.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: linear !! Flag controlling linear.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      real(dp), allocatable :: lambdas(:)
      real(dp), allocatable :: response(:, :), ar_design(:, :), x_design(:, :)
      real(dp), allocatable :: combined(:, :), gradient(:, :)
      real(dp), allocatable :: target_matrix(:, :)
      real(dp) :: selected_alpha, upper
      integer :: variables, x_variables, x_blocks
      logical :: use_contemporaneous, use_linear, use_intercept

      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      use_linear = .false.
      if (present(linear)) use_linear = linear
      use_intercept = .true.
      if (present(include_intercept)) use_intercept = include_intercept
      if (present(minnesota_target)) use_intercept = .false.
      x_blocks = exogenous_order
      if (use_contemporaneous) x_blocks = x_blocks + 1
      if (size(endogenous, 1) /= size(exogenous, 1) .or. variables < 1 .or. &
         x_variables < 1 .or. ar_order < 0 .or. exogenous_order < 0 .or. &
         (.not. use_contemporaneous .and. exogenous_order < 1) .or. &
         grid_ratio < 1.0_dp .or. grid_size < 1 .or. &
         .not. is_varx_structure(structure) .or. &
         (structure == bigvar_structure_efx .and. &
         (exogenous_order > ar_order .or. use_contemporaneous))) then
         allocate(lambdas(0))
         return
      end if
      if (ar_order == 0 .and. .not. supports_transfer_function(structure)) then
         allocate(lambdas(0))
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables .or. ar_order == 0) then
            allocate(lambdas(0))
            return
         end if
      end if
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
      end if
      call build_varx_data(endogenous, exogenous, ar_order, exogenous_order, &
         response, ar_design, x_design, use_contemporaneous)
      allocate(combined(size(response, 1), size(ar_design, 2) + &
         size(x_design, 2)))
      if (size(ar_design, 2) > 0) then
         combined(:, :size(ar_design, 2)) = ar_design
      end if
      combined(:, size(ar_design, 2) + 1:) = x_design
      if (present(minnesota_target)) then
         allocate(target_matrix(variables, size(combined, 2)))
         target_matrix = 0.0_dp
         call set_minnesota_target(target_matrix, minnesota_target)
         response = response - matmul(combined, transpose(target_matrix))
      else if (use_intercept) then
         response = response - spread(sum(response, dim=1)/ &
            real(size(response, 1), dp), 1, size(response, 1))
         combined = combined - spread(sum(combined, dim=1)/ &
            real(size(combined, 1), dp), 1, size(combined, 1))
      end if
      allocate(gradient(variables, size(combined, 2)))
      gradient = matmul(transpose(response), combined)
      if (is_nonconvex_structure(structure)) then
         upper = maxval(abs(gradient))/real(size(response, 1), dp)
      else
         upper = varx_zero_model_bound(gradient, structure, variables, &
            ar_order, x_variables, x_blocks, selected_alpha)
      end if
      upper = max(upper*(1.0_dp + 1.0e-10_dp), tiny(1.0_dp))
      allocate(lambdas(grid_size))
      call fill_lambda_grid(lambdas, upper, grid_ratio, use_linear)
   end function bigvar_varx_lambda_grid

   pure function bigvar_varx_lambda_alpha_grid(endogenous, exogenous, ar_order, &
      exogenous_order, structure, grid_ratio, grid_size, alphas, &
      contemporaneous, minnesota_target, linear, include_intercept) &
      result(lambdas)
      !! Construct an alpha-specific matrix of descending VARX lambda grids.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: grid_ratio !! Grid ratio.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: grid_size !! Grid size.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: linear !! Flag controlling linear.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      real(dp), allocatable :: lambdas(:, :)
      real(dp), allocatable :: column(:)
      integer :: alpha_index

      if (size(alphas) < 1 .or. .not. uses_alpha(structure) .or. &
         any(alphas < 0.0_dp) .or. any(alphas > 1.0_dp)) then
         allocate(lambdas(0, 0))
         return
      end if
      allocate(lambdas(grid_size, size(alphas)))
      do alpha_index = 1, size(alphas)
         column = bigvar_varx_lambda_grid(endogenous, exogenous, ar_order, &
            exogenous_order, structure, grid_ratio, grid_size, &
            alphas(alpha_index), contemporaneous, minnesota_target, linear, &
            include_intercept)
         if (size(column) /= grid_size) then
            deallocate(lambdas)
            allocate(lambdas(0, 0))
            return
         end if
         lambdas(:, alpha_index) = column
      end do
   end function bigvar_varx_lambda_alpha_grid

   pure function bigvar_bgr(series, lag_order, tau, random_walk) result(out)
      !! Estimate BigVAR's BGR Bayesian VAR from dummy observations.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: tau !! Tau.
      integer, intent(in) :: lag_order !! Model lag order.
      logical, intent(in), optional :: random_walk(:) !! Flag controlling random walk.
      type(bigvar_fit_t) :: out
      real(dp), allocatable :: response(:, :), lag_design(:, :), design(:, :)
      real(dp), allocatable :: scales(:), means(:), indicator(:)
      real(dp), allocatable :: prior_response(:, :), prior_design(:, :)
      real(dp), allocatable :: normal_matrix(:, :), inverse(:, :)
      real(dp), allocatable :: right_hand_side(:, :), coefficient(:, :)
      real(dp) :: mu, intercept_epsilon
      integer :: variables, predictors, prior_rows, status
      integer :: lag, variable, row

      out%lag_order = lag_order
      out%lambda = tau
      out%structure = bigvar_structure_bgr
      variables = size(series, 2)
      predictors = variables*lag_order
      if (variables < 1 .or. lag_order < 1 .or. &
         size(series, 1) <= 2*lag_order .or. tau <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(random_walk)) then
         if (size(random_walk) /= variables) then
            out%info = 1
            return
         end if
      end if
      call build_var_data(series, lag_order, response, lag_design)
      allocate(design(size(response, 1), predictors + 1))
      design(:, :predictors) = lag_design
      design(:, predictors + 1) = 1.0_dp
      scales = univariate_ar_scales(series, lag_order)
      means = sum(series, dim=1)/real(size(series, 1), dp)
      allocate(indicator(variables))
      indicator = 0.0_dp
      if (present(random_walk)) then
         where (random_walk) indicator = 1.0_dp
      end if

      prior_rows = predictors + 2*variables + 1
      allocate(prior_response(prior_rows, variables))
      allocate(prior_design(prior_rows, predictors + 1))
      prior_response = 0.0_dp
      prior_design = 0.0_dp
      do lag = 1, lag_order
         do variable = 1, variables
            row = (lag - 1)*variables + variable
            prior_design(row, row) = tau*real(lag, dp)*scales(variable)
            if (lag == 1) prior_response(row, variable) = &
               tau*scales(variable)*indicator(variable)
         end do
      end do
      mu = 0.1_dp*tau
      do variable = 1, variables
         row = predictors + variable
         prior_response(row, variable) = mu*means(variable)* &
            indicator(variable)
         do lag = 1, lag_order
            prior_design(row, (lag - 1)*variables + variable) = &
               mu*means(variable)*indicator(variable)
         end do
         row = predictors + variables + variable
         prior_response(row, variable) = scales(variable)
      end do
      intercept_epsilon = 1.0e-5_dp
      prior_design(prior_rows, predictors + 1) = intercept_epsilon

      allocate(normal_matrix(predictors + 1, predictors + 1))
      allocate(right_hand_side(predictors + 1, variables))
      normal_matrix = matmul(transpose(prior_design), prior_design) + &
         matmul(transpose(design), design)
      right_hand_side = matmul(transpose(prior_design), prior_response) + &
         matmul(transpose(design), response)
      call invert_matrix(normal_matrix, inverse, status)
      if (status /= 0) then
         if (allocated(inverse)) deallocate(inverse)
         allocate(inverse(predictors + 1, predictors + 1))
         call symmetric_pseudoinverse(normal_matrix, inverse, status)
      end if
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      coefficient = matmul(inverse, right_hand_side)
      allocate(out%phi(variables, predictors))
      allocate(out%intercept(variables))
      out%phi = transpose(coefficient(:predictors, :))
      out%intercept = coefficient(predictors + 1, :)
      allocate(out%fitted(size(response, 1), variables))
      out%fitted = matmul(design, coefficient)
      allocate(out%residuals(size(response, 1), variables))
      out%residuals = response - out%fitted
      out%objective = 0.5_dp*(sum(out%residuals**2) + &
         sum((prior_response - matmul(prior_design, coefficient))**2))
      out%nonzero = count(abs(out%phi) > 100.0_dp*epsilon(1.0_dp))
      out%active_groups = out%nonzero
      out%iterations = 1
      out%converged = .true.
      out%info = 0
   end function bigvar_bgr

   pure function bigvar_bgr_path(series, lag_order, tau, random_walk) &
      result(out)
      !! Estimate BGR fits over a supplied prior-tightness grid.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: tau(:) !! Tau.
      integer, intent(in) :: lag_order !! Model lag order.
      logical, intent(in), optional :: random_walk(:) !! Flag controlling random walk.
      type(bigvar_path_t) :: out
      type(bigvar_fit_t) :: fit
      integer :: index, variables, predictors

      variables = size(series, 2)
      predictors = variables*lag_order
      out%lag_order = lag_order
      out%structure = bigvar_structure_bgr
      if (size(tau) < 1 .or. any(tau <= 0.0_dp) .or. variables < 1 .or. &
         lag_order < 1) then
         out%info = 1
         return
      end if
      allocate(out%phi(variables, predictors, size(tau)))
      allocate(out%intercept(variables, size(tau)))
      allocate(out%objective(size(tau)), out%lambda(size(tau)))
      allocate(out%iterations(size(tau)), out%nonzero(size(tau)))
      allocate(out%active_groups(size(tau)), out%converged(size(tau)))
      out%lambda = tau
      do index = 1, size(tau)
         if (present(random_walk)) then
            fit = bigvar_bgr(series, lag_order, tau(index), random_walk)
         else
            fit = bigvar_bgr(series, lag_order, tau(index))
         end if
         if (fit%info /= 0) then
            out%info = 100*index + fit%info
            return
         end if
         out%phi(:, :, index) = fit%phi
         out%intercept(:, index) = fit%intercept
         out%objective(index) = fit%objective
         out%iterations(index) = fit%iterations
         out%nonzero(index) = fit%nonzero
         out%active_groups(index) = fit%active_groups
         out%converged(index) = fit%converged
      end do
      out%info = 0
   end function bigvar_bgr_path

   pure function bigvar_bgr_default_grid(variables, lag_order) result(tau)
      !! Construct BigVAR's package-default BGR tightness grid.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp), allocatable :: tau(:)
      integer :: index

      if (variables < 1 .or. lag_order < 1) then
         allocate(tau(0))
         return
      end if
      allocate(tau(161))
      do index = 1, 161
         tau(index) = (1.0_dp + 0.025_dp*real(index - 1, dp))* &
            sqrt(real(variables*lag_order, dp))
      end do
   end function bigvar_bgr_default_grid

   pure function bigvar_forecast(fit, history, horizon) result(out)
      !! Recursively forecast a BigVAR fit through the shared bigtime routine.
      type(bigvar_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: history(:, :) !! History.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bigtime_forecast_t) :: out
      type(bigtime_var_fit_t) :: adapter

      if (.not. allocated(fit%phi) .or. .not. allocated(fit%intercept)) then
         out%info = 1
         return
      end if
      if (fit%direct) then
         if (horizon /= fit%forecast_horizon) then
            out%info = 1
            return
         end if
         allocate(out%mean(horizon, size(fit%phi, 1)))
         out%mean = 0.0_dp
         out%mean(horizon, :) = bigvar_direct_forecast(fit, history)
         out%info = 0
         return
      end if
      adapter%phi = fit%phi
      adapter%intercept = fit%intercept
      adapter%lag_order = fit%lag_order
      adapter%info = fit%info
      out = bigtime_var_forecast(adapter, history, horizon)
   end function bigvar_forecast

   pure function bigvar_direct_forecast(fit, history) result(prediction)
      !! Forecast the fitted direct horizon from the latest observed VAR lags.
      type(bigvar_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: history(:, :) !! History.
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: lag_vector(:)
      integer :: variables, lag

      variables = size(history, 2)
      if (.not. allocated(fit%phi) .or. .not. allocated(fit%intercept) .or. &
         .not. fit%direct .or. variables /= size(fit%phi, 1) .or. &
         size(history, 1) < fit%lag_order) then
         allocate(prediction(0))
         return
      end if
      allocate(lag_vector(variables*fit%lag_order))
      do lag = 1, fit%lag_order
         lag_vector((lag - 1)*variables + 1:lag*variables) = &
            history(size(history, 1) + 1 - lag, :)
      end do
      prediction = fit%intercept + matmul(fit%phi, lag_vector)
   end function bigvar_direct_forecast

   pure function bigvar_innovation_covariance(residuals) result(covariance)
      !! Estimate the innovation covariance using BigVAR's MLE divisor.
      real(dp), intent(in) :: residuals(:, :) !! Model residuals.
      real(dp), allocatable :: covariance(:, :)
      integer :: observations, variables

      observations = size(residuals, 1)
      variables = size(residuals, 2)
      if (observations < 1 .or. variables < 1) then
         allocate(covariance(0, 0))
         return
      end if
      allocate(covariance(variables, variables))
      covariance = matmul(transpose(residuals), residuals)/ &
         real(observations, dp)
      covariance = 0.5_dp*(covariance + transpose(covariance))
   end function bigvar_innovation_covariance

   pure function bigvar_forecast_covariance(phi, lag_order, &
      innovation_covariance, horizon) result(covariance)
      !! Propagate VAR innovation covariance through moving-average coefficients.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: innovation_covariance(:, :) !! Innovation covariance matrix.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: psi(:, :, :), current(:, :)
      integer :: variables, step, lag

      variables = size(phi, 1)
      if (variables < 1 .or. lag_order < 1 .or. horizon < 1 .or. &
         size(phi, 2) /= variables*lag_order .or. &
         any(shape(innovation_covariance) /= [variables, variables])) then
         allocate(covariance(0, 0, 0))
         return
      end if
      allocate(covariance(variables, variables, horizon))
      allocate(psi(variables, variables, 0:horizon - 1))
      allocate(current(variables, variables))
      covariance = 0.0_dp
      psi = 0.0_dp
      do step = 1, variables
         psi(step, step, 0) = 1.0_dp
      end do
      do step = 1, horizon - 1
         do lag = 1, min(lag_order, step)
            psi(:, :, step) = psi(:, :, step) + matmul(phi(:, &
               (lag - 1)*variables + 1:lag*variables), &
               psi(:, :, step - lag))
         end do
      end do
      current = 0.0_dp
      do step = 1, horizon
         current = current + matmul(matmul(psi(:, :, step - 1), &
            innovation_covariance), transpose(psi(:, :, step - 1)))
         covariance(:, :, step) = 0.5_dp*(current + transpose(current))
      end do
   end function bigvar_forecast_covariance

   pure function bigvar_var_forecast_interval(fit, history, horizon, level) &
      result(out)
      !! Forecast a VAR with horizon-specific covariance and marginal intervals.
      type(bigvar_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: history(:, :) !! History.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in), optional :: level !! Model level or confidence level.
      type(bigvar_interval_forecast_t) :: out
      type(bigtime_forecast_t) :: point
      real(dp) :: selected_level, cutoff
      integer :: variables, step, variable

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      if (.not. allocated(fit%residuals) .or. horizon < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      if (fit%direct .and. horizon /= fit%forecast_horizon) then
         out%info = 1
         return
      end if
      point = bigvar_forecast(fit, history, horizon)
      if (point%info /= 0) then
         out%info = point%info
         return
      end if
      variables = size(fit%phi, 1)
      out%level = selected_level
      out%mean = point%mean
      out%innovation_covariance = &
         bigvar_innovation_covariance(fit%residuals)
      allocate(out%covariance(variables, variables, horizon))
      out%covariance = 0.0_dp
      if (fit%direct) then
         out%covariance(:, :, horizon) = out%innovation_covariance
      else
         out%covariance = bigvar_forecast_covariance(fit%phi, fit%lag_order, &
            out%innovation_covariance, horizon)
      end if
      allocate(out%standard_error(horizon, variables))
      allocate(out%lower(horizon, variables), out%upper(horizon, variables))
      out%standard_error = 0.0_dp
      do step = 1, horizon
         do variable = 1, variables
            out%standard_error(step, variable) = sqrt(max( &
               out%covariance(variable, variable, step), 0.0_dp))
         end do
      end do
      cutoff = normal_quantile(0.5_dp + 0.5_dp*selected_level)
      out%lower = out%mean - cutoff*out%standard_error
      out%upper = out%mean + cutoff*out%standard_error
      out%info = 0
   end function bigvar_var_forecast_interval

   pure function bigvar_varx_forecast_interval(fit, endogenous_history, &
      exogenous_values, horizon, level) result(out)
      !! Forecast VARX conditionally on exogenous values with normal intervals.
      type(bigvar_varx_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: endogenous_history(:, :) !! Endogenous history.
      real(dp), intent(in) :: exogenous_values(:, :) !! Exogenous values.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in), optional :: level !! Model level or confidence level.
      type(bigvar_interval_forecast_t) :: out
      type(bigtime_forecast_t) :: point
      real(dp) :: selected_level, cutoff
      integer :: variables, step, variable

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      if (.not. allocated(fit%residuals) .or. horizon < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      point = bigvar_varx_forecast(fit, endogenous_history, &
         exogenous_values, horizon)
      if (point%info /= 0) then
         out%info = point%info
         return
      end if
      variables = size(fit%phi, 1)
      out%level = selected_level
      out%mean = point%mean
      out%innovation_covariance = &
         bigvar_innovation_covariance(fit%residuals)
      if (fit%ar_order == 0) then
         allocate(out%covariance(variables, variables, horizon))
         do step = 1, horizon
            out%covariance(:, :, step) = out%innovation_covariance
         end do
      else
         out%covariance = bigvar_forecast_covariance(fit%phi, fit%ar_order, &
            out%innovation_covariance, horizon)
      end if
      allocate(out%standard_error(horizon, variables))
      allocate(out%lower(horizon, variables), out%upper(horizon, variables))
      do step = 1, horizon
         do variable = 1, variables
            out%standard_error(step, variable) = sqrt(max( &
               out%covariance(variable, variable, step), 0.0_dp))
         end do
      end do
      cutoff = normal_quantile(0.5_dp + 0.5_dp*selected_level)
      out%lower = out%mean - cutoff*out%standard_error
      out%upper = out%mean + cutoff*out%standard_error
      out%info = 0
   end function bigvar_varx_forecast_interval

   pure function bigvar_varx_forecast(fit, endogenous_history, &
      exogenous_values, horizon) result(out)
      !! Forecast a VARX fit, including optional contemporaneous predictors.
      type(bigvar_varx_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: endogenous_history(:, :) !! Endogenous history.
      real(dp), intent(in) :: exogenous_values(:, :) !! Exogenous values.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bigtime_forecast_t) :: out
      real(dp), allocatable :: working(:, :)
      integer :: observations, variables, x_variables
      integer :: step, lag, target, block

      if (.not. allocated(fit%phi) .or. .not. allocated(fit%beta) .or. &
         .not. allocated(fit%intercept)) then
         out%info = 1
         return
      end if
      observations = size(endogenous_history, 1)
      variables = size(endogenous_history, 2)
      x_variables = size(exogenous_values, 2)
      if (horizon < 1 .or. fit%ar_order < 0 .or. &
         fit%exogenous_order < 0 .or. observations < fit%ar_order .or. &
         size(fit%phi, 1) /= variables .or. &
         size(fit%phi, 2) /= variables*fit%ar_order .or. &
         size(fit%beta, 1) /= variables .or. &
         size(fit%beta, 2) /= x_variables*(fit%exogenous_order + &
         merge(1, 0, fit%contemporaneous)) .or. &
         size(fit%intercept) /= variables) then
         out%info = 1
         return
      end if
      if (fit%contemporaneous) then
         if (size(exogenous_values, 1) < observations + horizon) then
            out%info = 1
            return
         end if
      else
         if (fit%exogenous_order < 1 .or. &
            size(exogenous_values, 1) < observations + horizon - 1) then
            out%info = 1
            return
         end if
      end if
      allocate(working(observations + horizon, variables))
      working(:observations, :) = endogenous_history
      allocate(out%mean(horizon, variables))
      do step = 1, horizon
         target = observations + step
         working(target, :) = fit%intercept
         do lag = 1, fit%ar_order
            working(target, :) = working(target, :) + matmul( &
               fit%phi(:, (lag - 1)*variables + 1:lag*variables), &
               working(target - lag, :))
         end do
         do lag = 0, fit%exogenous_order
            if (.not. fit%contemporaneous .and. lag == 0) cycle
            block = lag
            if (fit%contemporaneous) block = block + 1
            working(target, :) = working(target, :) + matmul( &
               fit%beta(:, (block - 1)*x_variables + 1:block*x_variables), &
               exogenous_values(target - lag, :))
         end do
         out%mean(step, :) = working(target, :)
      end do
      out%info = 0
   end function bigvar_varx_forecast

   pure function bigvar_relaxed_var(series, fit, refit_fraction) result(out)
      !! Refit selected VAR coefficients by least squares and blend the result.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      type(bigvar_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      type(bigvar_fit_t) :: out
      real(dp), allocatable :: response(:, :), original_response(:, :)
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: centered_response(:, :), centered_design(:, :)
      real(dp), allocatable :: response_mean(:), design_mean(:)
      real(dp), allocatable :: refitted(:, :), target_matrix(:, :)
      real(dp), allocatable :: deviation(:, :)
      real(dp) :: fraction

      out = fit
      fraction = 1.0_dp
      if (present(refit_fraction)) fraction = refit_fraction
      if (.not. allocated(fit%phi) .or. .not. allocated(fit%intercept) .or. &
         fit%lag_order < 1 .or. size(series, 2) /= size(fit%phi, 1) .or. &
         size(series, 1) <= fit%lag_order .or. fraction < 0.0_dp .or. &
         fraction > 1.0_dp) then
         out%info = 1
         return
      end if
      call build_var_direct_data(series, fit%lag_order, fit%forecast_horizon, &
         response, design)
      original_response = response
      allocate(target_matrix(size(fit%phi, 1), size(fit%phi, 2)))
      target_matrix = 0.0_dp
      if (fit%minnesota) then
         call set_minnesota_target(target_matrix, fit%shrinkage_target)
         response = response - matmul(design, transpose(target_matrix))
      end if
      allocate(response_mean(size(response, 2)), design_mean(size(design, 2)))
      response_mean = 0.0_dp
      design_mean = 0.0_dp
      if (fit%include_intercept) then
         response_mean = sum(response, dim=1)/real(size(response, 1), dp)
         design_mean = sum(design, dim=1)/real(size(design, 1), dp)
      end if
      centered_response = response - spread(response_mean, 1, size(response, 1))
      centered_design = design - spread(design_mean, 1, size(design, 1))
      deviation = fit%phi - target_matrix
      refitted = relaxed_coefficients(centered_response, centered_design, &
         deviation)
      deviation = (1.0_dp - fraction)*deviation + fraction*refitted
      out%phi = deviation + target_matrix
      out%intercept = response_mean - matmul(deviation, design_mean)
      out%fitted = matmul(design, transpose(out%phi)) + &
         spread(out%intercept, 1, size(response, 1))
      out%residuals = original_response - out%fitted
      if (allocated(out%lambda_by_response)) then
         out%objective = separate_penalized_objective(out%residuals, deviation, &
            out%lambda_by_response, out%structure, size(out%phi, 1), &
            out%lag_order, out%alpha, out%gamma)
      else if (is_nonconvex_structure(out%structure)) then
         out%objective = nonconvex_objective(out%residuals, deviation, &
            out%lambda, out%gamma, out%structure)
      else
         out%objective = 0.5_dp*sum(out%residuals**2) + out%lambda* &
            structured_penalty(deviation, out%structure, size(out%phi, 1), &
            out%lag_order, out%alpha, out%gamma)
      end if
      out%nonzero = count(abs(deviation) > 100.0_dp*epsilon(1.0_dp))
      out%active_groups = count_active_groups(deviation, out%structure, &
         size(out%phi, 1), out%lag_order)
      out%relaxed = .true.
      out%refit_fraction = fraction
   end function bigvar_relaxed_var

   pure function bigvar_relaxed_varx(endogenous, exogenous, fit, &
      refit_fraction) result(out)
      !! Refit selected VARX coefficients by least squares and blend the result.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      type(bigvar_varx_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      type(bigvar_varx_fit_t) :: out
      real(dp), allocatable :: response(:, :), original_response(:, :)
      real(dp), allocatable :: ar_design(:, :), x_design(:, :)
      real(dp), allocatable :: design(:, :), centered_response(:, :)
      real(dp), allocatable :: centered_design(:, :), response_mean(:)
      real(dp), allocatable :: design_mean(:), coefficient(:, :), refitted(:, :)
      real(dp), allocatable :: target_matrix(:, :)
      real(dp) :: fraction
      integer :: ar_predictors, x_blocks, x_variables

      out = fit
      fraction = 1.0_dp
      if (present(refit_fraction)) fraction = refit_fraction
      x_variables = size(exogenous, 2)
      x_blocks = fit%exogenous_order
      if (fit%contemporaneous) x_blocks = x_blocks + 1
      if (.not. allocated(fit%phi) .or. .not. allocated(fit%beta) .or. &
         .not. allocated(fit%intercept) .or. &
         size(endogenous, 1) /= size(exogenous, 1) .or. &
         size(endogenous, 2) /= size(fit%phi, 1) .or. &
         size(fit%phi, 2) /= size(endogenous, 2)*fit%ar_order .or. &
         size(fit%beta, 1) /= size(endogenous, 2) .or. &
         size(fit%beta, 2) /= x_variables*x_blocks .or. &
         fit%ar_order < 0 .or. fit%exogenous_order < 0 .or. &
         (.not. fit%contemporaneous .and. fit%exogenous_order < 1) .or. &
         size(endogenous, 1) <= max(fit%ar_order, fit%exogenous_order) .or. &
         fraction < 0.0_dp .or. fraction > 1.0_dp) then
         out%info = 1
         return
      end if
      call build_varx_data(endogenous, exogenous, fit%ar_order, &
         fit%exogenous_order, response, ar_design, x_design, &
         fit%contemporaneous)
      original_response = response
      ar_predictors = size(ar_design, 2)
      allocate(design(size(response, 1), ar_predictors + size(x_design, 2)))
      if (ar_predictors > 0) design(:, :ar_predictors) = ar_design
      design(:, ar_predictors + 1:) = x_design
      allocate(target_matrix(size(fit%phi, 1), size(design, 2)))
      target_matrix = 0.0_dp
      if (fit%minnesota) then
         call set_minnesota_target(target_matrix, fit%shrinkage_target)
         response = response - matmul(design, transpose(target_matrix))
      end if
      allocate(response_mean(size(response, 2)), design_mean(size(design, 2)))
      response_mean = 0.0_dp
      design_mean = 0.0_dp
      if (fit%include_intercept) then
         response_mean = sum(response, dim=1)/real(size(response, 1), dp)
         design_mean = sum(design, dim=1)/real(size(design, 1), dp)
      end if
      centered_response = response - spread(response_mean, 1, size(response, 1))
      centered_design = design - spread(design_mean, 1, size(design, 1))
      allocate(coefficient(size(fit%phi, 1), size(design, 2)))
      if (ar_predictors > 0) coefficient(:, :ar_predictors) = fit%phi - &
         target_matrix(:, :ar_predictors)
      coefficient(:, ar_predictors + 1:) = fit%beta
      refitted = relaxed_coefficients(centered_response, centered_design, &
         coefficient)
      coefficient = (1.0_dp - fraction)*coefficient + fraction*refitted
      if (ar_predictors > 0) out%phi = coefficient(:, :ar_predictors) + &
         target_matrix(:, :ar_predictors)
      out%beta = coefficient(:, ar_predictors + 1:)
      out%intercept = response_mean - matmul(coefficient, design_mean)
      out%fitted = matmul(design, transpose(coefficient + target_matrix)) + &
         spread(out%intercept, 1, size(response, 1))
      out%residuals = original_response - out%fitted
      if (is_nonconvex_structure(out%structure)) then
         out%objective = nonconvex_objective(out%residuals, coefficient, &
            out%lambda, out%gamma, out%structure)
      else
         out%objective = 0.5_dp*sum(out%residuals**2) + out%lambda* &
            varx_penalty(coefficient(:, :ar_predictors), out%beta, &
            out%structure, &
            size(out%phi, 1), out%ar_order, size(exogenous, 2), &
            x_blocks, out%alpha)
      end if
      out%nonzero = count(abs(coefficient) > 100.0_dp*epsilon(1.0_dp))
      if (out%structure == bigvar_structure_efx) then
         out%active_groups = count_active_efx_groups( &
            coefficient(:, :ar_predictors), out%beta, &
            size(out%phi, 1), out%ar_order, size(exogenous, 2), &
            x_blocks)
      else
      out%active_groups = count_active_x_groups(out%beta)
      if (ar_predictors > 0) out%active_groups = out%active_groups + &
         count_active_groups(coefficient(:, :ar_predictors), out%structure, &
         size(out%phi, 1), out%ar_order)
      end if
      out%relaxed = .true.
      out%refit_fraction = fraction
   end function bigvar_relaxed_varx

   pure function bigvar_forecast_loss(error, loss_type, delta) result(value)
      !! Sum the selected BigVAR forecast loss over response components.
      real(dp), intent(in) :: error(:) !! Error.
      integer, intent(in) :: loss_type !! Loss type.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp) :: value
      real(dp) :: selected_delta

      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      select case (loss_type)
      case (bigvar_loss_l1)
         value = sum(abs(error))
      case (bigvar_loss_l2)
         value = sum(error**2)
      case (bigvar_loss_huber)
         value = sum(merge(0.5_dp*error**2, &
            selected_delta*(abs(error) - 0.5_dp*selected_delta), &
            abs(error) < selected_delta))
      case default
         value = huge(1.0_dp)
      end select
   end function bigvar_forecast_loss

   pure elemental real(dp) function huber_component(error, delta) result(value)
      !! Return one component of the Huber forecast loss.
      real(dp), intent(in) :: error !! Error.
      real(dp), intent(in) :: delta !! Model increment or differencing parameter.

      if (abs(error) <= delta) then
         value = 0.5_dp*error**2
      else
         value = delta*(abs(error) - 0.5_dp*delta)
      end if
   end function huber_component

   pure function bigvar_var_validate(series, lag_order, lambdas, structure, &
      first_origin, last_origin, horizon, loss_type, delta, &
      one_standard_error, window_size, alphas, tolerance, max_iterations, &
      gamma, random_walk, refit_fraction, recursive, minnesota_target, &
      include_intercept) &
      result(out)
      !! Select BigVAR penalties by rolling multi-step forecast validation.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alphas(:) !! Alphas.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: random_walk(:) !! Flag controlling random walk.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      type(bigvar_fit_t) :: fit
      type(bigtime_forecast_t) :: forecast
      real(dp), allocatable :: alpha_grid(:)
      real(dp) :: selected_delta, selected_tolerance, selected_gamma
      logical :: use_one_se, use_recursive
      integer :: selected_loss, selected_window, selected_max_iterations
      integer :: origins, candidates, alpha_count, lambda_count
      integer :: origin_index, candidate, alpha_index, lambda_index
      integer :: origin, training_start

      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_window = 0
      if (present(window_size)) selected_window = window_size
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      use_recursive = .true.
      if (present(recursive)) use_recursive = recursive
      lambda_count = size(lambdas)
      alpha_count = 1
      if (present(alphas)) alpha_count = size(alphas)
      allocate(alpha_grid(alpha_count))
      if (present(alphas)) then
         alpha_grid = alphas
      else
         alpha_grid = 1.0_dp/real(size(series, 2) + 1, dp)
      end if
      origins = last_origin - first_origin + 1
      candidates = lambda_count*alpha_count
      if (size(series, 2) < 1 .or. lag_order < 1 .or. lambda_count < 1 .or. &
         alpha_count < 1 .or. origins < 1 .or. horizon < 1 .or. &
         first_origin <= lag_order .or. last_origin + horizon > &
         size(series, 1) .or. selected_window < 0 .or. &
         selected_loss < bigvar_loss_l1 .or. &
         selected_loss > bigvar_loss_huber .or. selected_delta <= 0.0_dp .or. &
         any(lambdas < 0.0_dp) .or. any(alpha_grid < 0.0_dp) .or. &
         any(alpha_grid > 1.0_dp)) then
         out%info = 1
         return
      end if
      if (structure == bigvar_structure_bgr .and. alpha_count > 1) then
         out%info = 1
         return
      end if
      if (structure == bigvar_structure_bgr .and. present(include_intercept)) then
         if (.not. include_intercept) then
            out%info = 1
            return
         end if
      end if
      if (structure == bigvar_structure_bgr .and. horizon > 1 .and. &
         .not. use_recursive) then
         out%info = 1
         return
      end if
      if (present(random_walk)) then
         if (size(random_walk) /= size(series, 2)) then
            out%info = 1
            return
         end if
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= size(series, 2) .or. &
            structure == bigvar_structure_bgr) then
            out%info = 1
            return
         end if
      end if

      call initialize_validation(out, origins, candidates, size(series, 2), &
         first_origin, last_origin, horizon, selected_loss)
      call fill_candidate_grid(out, lambdas, alpha_grid)
      do origin_index = 1, origins
         origin = first_origin + origin_index - 1
         training_start = 1
         if (selected_window > 0) then
            if (horizon > 1 .and. .not. use_recursive) then
               training_start = max(1, origin - selected_window - horizon + 1)
            else
               training_start = max(1, origin - selected_window + 1)
            end if
         end if
         if (origin - training_start + 1 <= lag_order) cycle
         do alpha_index = 1, alpha_count
            do lambda_index = 1, lambda_count
               candidate = lambda_index + lambda_count*(alpha_index - 1)
               if (structure == bigvar_structure_bgr) then
                  if (present(random_walk)) then
                     fit = bigvar_bgr(series(training_start:origin, :), &
                        lag_order, lambdas(lambda_index), random_walk)
                  else
                     fit = bigvar_bgr(series(training_start:origin, :), &
                        lag_order, lambdas(lambda_index))
                  end if
               else
                  if (present(refit_fraction)) then
                     fit = bigvar_structured_var( &
                        series(training_start:origin, :), lag_order, &
                        lambdas(lambda_index), structure, selected_tolerance, &
                        selected_max_iterations, alpha=alpha_grid(alpha_index), &
                        gamma=selected_gamma, &
                        refit_fraction=refit_fraction, &
                        direct_horizon=merge(1, horizon, use_recursive), &
                        minnesota_target=minnesota_target, &
                        include_intercept=include_intercept)
                  else
                     if (present(minnesota_target)) then
                        fit = bigvar_structured_var( &
                           series(training_start:origin, :), lag_order, &
                           lambdas(lambda_index), structure, &
                           selected_tolerance, selected_max_iterations, &
                           alpha=alpha_grid(alpha_index), gamma=selected_gamma, &
                           direct_horizon=merge(1, horizon, use_recursive), &
                           minnesota_target=minnesota_target, &
                           include_intercept=include_intercept)
                     else
                        fit = bigvar_structured_var( &
                           series(training_start:origin, :), lag_order, &
                           lambdas(lambda_index), structure, &
                           selected_tolerance, selected_max_iterations, &
                           alpha=alpha_grid(alpha_index), gamma=selected_gamma, &
                           direct_horizon=merge(1, horizon, use_recursive), &
                           include_intercept=include_intercept)
                     end if
                  end if
               end if
               if (fit%info /= 0 .and. fit%info /= 4) cycle
               forecast = bigvar_forecast(fit, &
                  series(training_start:origin, :), horizon)
               if (forecast%info /= 0) cycle
               out%forecasts(origin_index, candidate, :) = &
                  forecast%mean(horizon, :)
               out%loss(origin_index, candidate) = bigvar_forecast_loss( &
                  series(origin + horizon, :) - forecast%mean(horizon, :), &
                  selected_loss, selected_delta)
               out%nonzero(origin_index, candidate) = fit%nonzero
               out%valid(origin_index, candidate) = .true.
               if (fit%info == 4) out%info = 4
            end do
         end do
      end do
      call summarize_validation(out, use_one_se)
   end function bigvar_var_validate

   pure function bigvar_varx_validate(endogenous, exogenous, ar_order, &
      exogenous_order, lambdas, structure, first_origin, last_origin, &
      horizon, loss_type, delta, one_standard_error, window_size, alphas, &
      tolerance, max_iterations, gamma, refit_fraction, contemporaneous, &
      minnesota_target, include_intercept) result(out)
      !! Select BigVAR VARX penalties by rolling multi-step validation.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alphas(:) !! Alphas.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      type(bigvar_varx_fit_t) :: fit
      type(bigtime_forecast_t) :: forecast
      real(dp), allocatable :: alpha_grid(:)
      real(dp) :: selected_delta, selected_tolerance, selected_gamma
      logical :: use_one_se, use_contemporaneous
      integer :: selected_loss, selected_window, selected_max_iterations
      integer :: origins, candidates, alpha_count, lambda_count
      integer :: origin_index, candidate, alpha_index, lambda_index
      integer :: origin, training_start, required_order

      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_window = 0
      if (present(window_size)) selected_window = window_size
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      lambda_count = size(lambdas)
      alpha_count = 1
      if (present(alphas)) alpha_count = size(alphas)
      allocate(alpha_grid(alpha_count))
      if (present(alphas)) then
         alpha_grid = alphas
      else
         alpha_grid = 1.0_dp/real(size(endogenous, 2) + 1, dp)
      end if
      origins = last_origin - first_origin + 1
      candidates = lambda_count*alpha_count
      required_order = max(ar_order, exogenous_order)
      if (size(endogenous, 1) /= size(exogenous, 1) .or. &
         size(endogenous, 2) < 1 .or. size(exogenous, 2) < 1 .or. &
         ar_order < 0 .or. exogenous_order < 0 .or. &
         (.not. use_contemporaneous .and. exogenous_order < 1) .or. &
         lambda_count < 1 .or. &
         alpha_count < 1 .or. origins < 1 .or. horizon < 1 .or. &
         first_origin <= required_order .or. last_origin + horizon > &
         size(endogenous, 1) .or. selected_window < 0 .or. &
         selected_loss < bigvar_loss_l1 .or. &
         selected_loss > bigvar_loss_huber .or. selected_delta <= 0.0_dp .or. &
         any(lambdas < 0.0_dp) .or. any(alpha_grid < 0.0_dp) .or. &
         any(alpha_grid > 1.0_dp) .or. &
         (structure == bigvar_structure_efx .and. &
         (exogenous_order > ar_order .or. use_contemporaneous))) then
         out%info = 1
         return
      end if
      if (ar_order == 0 .and. .not. supports_transfer_function(structure)) then
         out%info = 1
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= size(endogenous, 2) .or. &
            ar_order == 0) then
            out%info = 1
            return
         end if
      end if

      call initialize_validation(out, origins, candidates, &
         size(endogenous, 2), first_origin, last_origin, horizon, selected_loss)
      call fill_candidate_grid(out, lambdas, alpha_grid)
      do origin_index = 1, origins
         origin = first_origin + origin_index - 1
         training_start = 1
         if (selected_window > 0) then
            training_start = max(1, origin - selected_window + 1)
         end if
         if (origin - training_start + 1 <= required_order) cycle
         do alpha_index = 1, alpha_count
            do lambda_index = 1, lambda_count
               candidate = lambda_index + lambda_count*(alpha_index - 1)
               if (present(refit_fraction)) then
                  fit = bigvar_structured_varx( &
                     endogenous(training_start:origin, :), &
                     exogenous(training_start:origin, :), ar_order, &
                     exogenous_order, lambdas(lambda_index), structure, &
                     selected_tolerance, selected_max_iterations, &
                     alpha=alpha_grid(alpha_index), gamma=selected_gamma, &
                     refit_fraction=refit_fraction, &
                     contemporaneous=use_contemporaneous, &
                     minnesota_target=minnesota_target, &
                     include_intercept=include_intercept)
               else
                  if (present(minnesota_target)) then
                     fit = bigvar_structured_varx( &
                        endogenous(training_start:origin, :), &
                        exogenous(training_start:origin, :), ar_order, &
                        exogenous_order, lambdas(lambda_index), structure, &
                        selected_tolerance, selected_max_iterations, &
                        alpha=alpha_grid(alpha_index), gamma=selected_gamma, &
                        contemporaneous=use_contemporaneous, &
                        minnesota_target=minnesota_target, &
                        include_intercept=include_intercept)
                  else
                     fit = bigvar_structured_varx( &
                        endogenous(training_start:origin, :), &
                        exogenous(training_start:origin, :), ar_order, &
                        exogenous_order, lambdas(lambda_index), structure, &
                        selected_tolerance, selected_max_iterations, &
                        alpha=alpha_grid(alpha_index), gamma=selected_gamma, &
                        contemporaneous=use_contemporaneous, &
                        include_intercept=include_intercept)
                  end if
               end if
               if (fit%info /= 0 .and. fit%info /= 4) cycle
               if (use_contemporaneous) then
                  forecast = bigvar_varx_forecast(fit, &
                     endogenous(training_start:origin, :), &
                     exogenous(training_start:origin + horizon, :), horizon)
               else
                  forecast = bigvar_varx_forecast(fit, &
                     endogenous(training_start:origin, :), &
                     exogenous(training_start:origin + horizon - 1, :), &
                     horizon)
               end if
               if (forecast%info /= 0) cycle
               out%forecasts(origin_index, candidate, :) = &
                  forecast%mean(horizon, :)
               out%loss(origin_index, candidate) = bigvar_forecast_loss( &
                  endogenous(origin + horizon, :) - &
                  forecast%mean(horizon, :), selected_loss, selected_delta)
               out%nonzero(origin_index, candidate) = fit%nonzero
               out%valid(origin_index, candidate) = .true.
               if (fit%info == 4) out%info = 4
            end do
         end do
      end do
      call summarize_validation(out, use_one_se)
   end function bigvar_varx_validate

   pure function bigvar_var_validate_separate(series, lag_order, lambdas, &
      structure, first_origin, last_origin, horizon, loss_type, delta, &
      one_standard_error, window_size, alpha, tolerance, max_iterations, &
      gamma, refit_fraction, recursive, minnesota_target, include_intercept) &
      result(out)
      !! Select a separate rolling-validation penalty for each VAR response.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_separate_validation_t) :: out
      type(bigvar_fit_t) :: fit
      type(bigtime_forecast_t) :: forecast
      real(dp) :: selected_delta, selected_alpha
      real(dp) :: selected_tolerance, selected_gamma, error_value
      integer :: selected_loss, selected_window, selected_max_iterations
      integer :: origins, candidates, variables, origin_index, candidate
      integer :: response, origin, training_start
      logical :: use_one_se, use_recursive

      variables = size(series, 2)
      origins = last_origin - first_origin + 1
      candidates = size(lambdas, 1)
      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_window = 0
      if (present(window_size)) selected_window = window_size
      selected_alpha = 0.0_dp
      if (uses_alpha(structure)) then
         selected_alpha = 1.0_dp/real(variables + 1, dp)
         if (present(alpha)) selected_alpha = alpha
      end if
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      use_recursive = .true.
      if (present(recursive)) use_recursive = recursive
      if (variables < 1 .or. lag_order < 1 .or. origins < 1 .or. &
         candidates < 1 .or. size(lambdas, 2) /= variables .or. &
         any(lambdas < 0.0_dp) .or. .not. supports_separate_lambdas(structure) &
         .or. first_origin <= lag_order .or. last_origin + horizon > &
         size(series, 1) .or. horizon < 1 .or. selected_window < 0 .or. &
         selected_loss < bigvar_loss_l1 .or. &
         selected_loss > bigvar_loss_huber .or. selected_delta <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables) then
            out%info = 1
            return
         end if
      end if
      allocate(out%loss(origins, candidates, variables))
      allocate(out%forecasts(origins, candidates, variables))
      allocate(out%mean_loss(candidates, variables))
      allocate(out%standard_error(candidates, variables))
      allocate(out%lambda(candidates, variables))
      allocate(out%selected_index(variables), out%best_index(variables))
      allocate(out%one_se_index(variables), out%selected_lambda(variables))
      allocate(out%valid(origins, candidates, variables))
      out%loss = huge(1.0_dp)
      out%forecasts = 0.0_dp
      out%mean_loss = huge(1.0_dp)
      out%standard_error = huge(1.0_dp)
      out%lambda = lambdas
      out%selected_index = 0
      out%best_index = 0
      out%one_se_index = 0
      out%selected_lambda = 0.0_dp
      out%valid = .false.
      out%first_origin = first_origin
      out%last_origin = last_origin
      out%horizon = horizon
      out%loss_type = selected_loss
      do origin_index = 1, origins
         origin = first_origin + origin_index - 1
         training_start = 1
         if (selected_window > 0) then
            if (horizon > 1 .and. .not. use_recursive) then
               training_start = max(1, origin - selected_window - horizon + 1)
            else
               training_start = max(1, origin - selected_window + 1)
            end if
         end if
         if (origin - training_start + 1 <= lag_order) cycle
         do candidate = 1, candidates
            if (present(refit_fraction)) then
               fit = bigvar_structured_var_separate( &
                  series(training_start:origin, :), lag_order, &
                  lambdas(candidate, :), structure, selected_tolerance, &
                  selected_max_iterations, alpha=selected_alpha, &
                  gamma=selected_gamma, refit_fraction=refit_fraction, &
                  direct_horizon=merge(1, horizon, use_recursive), &
                  minnesota_target=minnesota_target, &
                  include_intercept=include_intercept)
            else
               if (present(minnesota_target)) then
                  fit = bigvar_structured_var_separate( &
                     series(training_start:origin, :), lag_order, &
                     lambdas(candidate, :), structure, selected_tolerance, &
                     selected_max_iterations, alpha=selected_alpha, &
                     gamma=selected_gamma, &
                     direct_horizon=merge(1, horizon, use_recursive), &
                     minnesota_target=minnesota_target, &
                     include_intercept=include_intercept)
               else
                  fit = bigvar_structured_var_separate( &
                     series(training_start:origin, :), lag_order, &
                     lambdas(candidate, :), structure, selected_tolerance, &
                     selected_max_iterations, alpha=selected_alpha, &
                     gamma=selected_gamma, &
                     direct_horizon=merge(1, horizon, use_recursive), &
                     include_intercept=include_intercept)
               end if
            end if
            if (fit%info /= 0 .and. fit%info /= 4) cycle
            forecast = bigvar_forecast(fit, &
               series(training_start:origin, :), horizon)
            if (forecast%info /= 0) cycle
            out%forecasts(origin_index, candidate, :) = &
               forecast%mean(horizon, :)
            do response = 1, variables
               error_value = series(origin + horizon, response) - &
                  forecast%mean(horizon, response)
               out%loss(origin_index, candidate, response) = &
                  bigvar_forecast_loss([error_value], selected_loss, &
                  selected_delta)
               out%valid(origin_index, candidate, response) = .true.
            end do
            if (fit%info == 4) out%info = 4
         end do
      end do
      call summarize_separate_validation(out, use_one_se)
   end function bigvar_var_validate_separate

   pure function bigvar_var_evaluate(series, lag_order, lambda, structure, &
      first_origin, last_origin, horizon, loss_type, delta, window_size, &
      alpha, tolerance, max_iterations, gamma, random_walk, refit_fraction, &
      recursive, minnesota_target, include_intercept) result(out)
      !! Refit a selected BigVAR model over a rolling evaluation period.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: random_walk(:) !! Flag controlling random walk.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      real(dp) :: selected_delta, selected_alpha
      real(dp) :: selected_tolerance, selected_gamma
      integer :: selected_loss, selected_window, selected_max_iterations
      logical :: use_recursive

      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_window = 0
      if (present(window_size)) selected_window = window_size
      selected_alpha = 1.0_dp/real(size(series, 2) + 1, dp)
      if (present(alpha)) selected_alpha = alpha
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      use_recursive = .true.
      if (present(recursive)) use_recursive = recursive
      if (present(random_walk)) then
         out = bigvar_var_validate(series, lag_order, [lambda], structure, &
            first_origin, last_origin, horizon, selected_loss, selected_delta, &
            .false., selected_window, [selected_alpha], selected_tolerance, &
            selected_max_iterations, selected_gamma, random_walk, &
            recursive=use_recursive, minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
      else if (present(refit_fraction)) then
         out = bigvar_var_validate(series, lag_order, [lambda], structure, &
            first_origin, last_origin, horizon, selected_loss, selected_delta, &
            .false., selected_window, [selected_alpha], selected_tolerance, &
            selected_max_iterations, selected_gamma, &
            refit_fraction=refit_fraction, recursive=use_recursive, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
      else
         out = bigvar_var_validate(series, lag_order, [lambda], structure, &
            first_origin, last_origin, horizon, selected_loss, selected_delta, &
            .false., selected_window, [selected_alpha], selected_tolerance, &
            selected_max_iterations, selected_gamma, recursive=use_recursive, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
      end if
   end function bigvar_var_evaluate

   pure function bigvar_varx_evaluate(endogenous, exogenous, ar_order, &
      exogenous_order, lambda, structure, first_origin, last_origin, horizon, &
      loss_type, delta, window_size, alpha, tolerance, max_iterations, gamma, &
      refit_fraction, contemporaneous, minnesota_target, include_intercept) &
      result(out)
      !! Refit a selected BigVAR VARX model over a rolling evaluation period.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      real(dp) :: selected_delta, selected_alpha
      real(dp) :: selected_tolerance, selected_gamma
      integer :: selected_loss, selected_window, selected_max_iterations
      logical :: use_contemporaneous

      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_window = 0
      if (present(window_size)) selected_window = window_size
      selected_alpha = 1.0_dp/real(size(endogenous, 2) + 1, dp)
      if (present(alpha)) selected_alpha = alpha
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      if (present(refit_fraction)) then
         out = bigvar_varx_validate(endogenous, exogenous, ar_order, &
            exogenous_order, [lambda], structure, first_origin, last_origin, &
            horizon, selected_loss, selected_delta, .false., selected_window, &
            [selected_alpha], selected_tolerance, selected_max_iterations, &
            selected_gamma, refit_fraction, use_contemporaneous, &
            minnesota_target, include_intercept)
      else
         out = bigvar_varx_validate(endogenous, exogenous, ar_order, &
            exogenous_order, [lambda], structure, first_origin, last_origin, &
            horizon, selected_loss, selected_delta, .false., selected_window, &
            [selected_alpha], selected_tolerance, selected_max_iterations, &
            selected_gamma, contemporaneous=use_contemporaneous, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
      end if
   end function bigvar_varx_evaluate

   pure function bigvar_var_reselect(series, lag_order, lambdas, structure, &
      validation_first_origin, first_origin, last_origin, horizon, loss_type, &
      delta, one_standard_error, selection_window, window_size, alphas, &
      tolerance, max_iterations, gamma, random_walk, refit_fraction, recursive, &
      minnesota_target, include_intercept) result(out)
      !! Reselect VAR penalties cumulatively before every evaluation forecast.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: validation_first_origin !! Validation first origin.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: selection_window !! Selection window.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alphas(:) !! Alphas.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: random_walk(:) !! Flag controlling random walk.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_reselection_t) :: out
      type(bigvar_validation_t) :: selection, evaluation
      integer :: origin, index, selection_first, selection_last

      call initialize_reselection(out, first_origin, last_origin, horizon, &
         validation_first_origin, size(series, 2), 1, loss_type)
      if (out%info /= 0) return
      do index = 1, size(out%valid)
         origin = first_origin + index - 1
         selection_last = origin - horizon
         selection_first = validation_first_origin
         if (present(selection_window)) then
            if (selection_window > 0) then
               selection_first = max(selection_first, &
                  selection_last - selection_window + 1)
            end if
         end if
         if (selection_last < selection_first) cycle
         selection = bigvar_var_validate(series, lag_order, lambdas, structure, &
            selection_first, selection_last, horizon, loss_type=loss_type, &
            delta=delta, one_standard_error=one_standard_error, &
            window_size=window_size, alphas=alphas, tolerance=tolerance, &
            max_iterations=max_iterations, gamma=gamma, random_walk=random_walk, &
            refit_fraction=refit_fraction, recursive=recursive, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (selection%selected_index < 1) cycle
         evaluation = bigvar_var_evaluate(series, lag_order, &
            selection%selected_lambda, structure, origin, origin, horizon, &
            loss_type=loss_type, delta=delta, window_size=window_size, &
            alpha=selection%selected_alpha, tolerance=tolerance, &
            max_iterations=max_iterations, gamma=gamma, random_walk=random_walk, &
            refit_fraction=refit_fraction, recursive=recursive, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         call store_joint_reselection(out, index, selection, evaluation)
      end do
      call finalize_reselection(out)
   end function bigvar_var_reselect

   pure function bigvar_var_reselect_dual(series, lag_order, lambdas, alphas, &
      structure, validation_first_origin, first_origin, last_origin, horizon, &
      loss_type, delta, one_standard_error, selection_window, window_size, &
      tolerance, max_iterations, gamma, refit_fraction, recursive, &
      minnesota_target, include_intercept) result(out)
      !! Reselect an alpha-specific VAR lambda surface at every origin.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: validation_first_origin !! Validation first origin.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: selection_window !! Selection window.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_reselection_t) :: out
      type(bigvar_validation_t) :: selection, evaluation
      integer :: origin, index, selection_first, selection_last

      call initialize_reselection(out, first_origin, last_origin, horizon, &
         validation_first_origin, size(series, 2), 1, loss_type)
      if (out%info /= 0 .or. .not. valid_dual_grid(lambdas, alphas, structure)) &
         then
         out%info = 1
         return
      end if
      do index = 1, size(out%valid)
         origin = first_origin + index - 1
         selection_last = origin - horizon
         selection_first = validation_first_origin
         if (present(selection_window)) then
            if (selection_window > 0) then
               selection_first = max(selection_first, &
                  selection_last - selection_window + 1)
            end if
         end if
         if (selection_last < selection_first) cycle
         selection = bigvar_var_validate_dual(series, lag_order, lambdas, &
            alphas, structure, selection_first, selection_last, horizon, &
            loss_type=loss_type, delta=delta, &
            one_standard_error=one_standard_error, window_size=window_size, &
            tolerance=tolerance, max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, recursive=recursive, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (selection%selected_index < 1) cycle
         evaluation = bigvar_var_evaluate(series, lag_order, &
            selection%selected_lambda, structure, origin, origin, horizon, &
            loss_type=loss_type, delta=delta, window_size=window_size, &
            alpha=selection%selected_alpha, tolerance=tolerance, &
            max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, recursive=recursive, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         call store_joint_reselection(out, index, selection, evaluation)
      end do
      call finalize_reselection(out)
   end function bigvar_var_reselect_dual

   pure function bigvar_varx_reselect(endogenous, exogenous, ar_order, &
      exogenous_order, lambdas, structure, validation_first_origin, &
      first_origin, last_origin, horizon, loss_type, delta, one_standard_error, &
      selection_window, window_size, alphas, tolerance, max_iterations, gamma, &
      refit_fraction, contemporaneous, minnesota_target, include_intercept) &
      result(out)
      !! Reselect VARX penalties cumulatively before every evaluation forecast.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: validation_first_origin !! Validation first origin.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: selection_window !! Selection window.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alphas(:) !! Alphas.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_reselection_t) :: out
      type(bigvar_validation_t) :: selection, evaluation
      integer :: origin, index, selection_first, selection_last

      call initialize_reselection(out, first_origin, last_origin, horizon, &
         validation_first_origin, size(endogenous, 2), 1, loss_type)
      if (out%info /= 0) return
      do index = 1, size(out%valid)
         origin = first_origin + index - 1
         selection_last = origin - horizon
         selection_first = validation_first_origin
         if (present(selection_window)) then
            if (selection_window > 0) then
               selection_first = max(selection_first, &
                  selection_last - selection_window + 1)
            end if
         end if
         if (selection_last < selection_first) cycle
         selection = bigvar_varx_validate(endogenous, exogenous, ar_order, &
            exogenous_order, lambdas, structure, selection_first, &
            selection_last, horizon, loss_type=loss_type, delta=delta, &
            one_standard_error=one_standard_error, window_size=window_size, &
            alphas=alphas, tolerance=tolerance, &
            max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, contemporaneous=contemporaneous, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (selection%selected_index < 1) cycle
         evaluation = bigvar_varx_evaluate(endogenous, exogenous, ar_order, &
            exogenous_order, selection%selected_lambda, structure, origin, &
            origin, horizon, loss_type=loss_type, delta=delta, &
            window_size=window_size, alpha=selection%selected_alpha, &
            tolerance=tolerance, max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, contemporaneous=contemporaneous, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         call store_joint_reselection(out, index, selection, evaluation)
      end do
      call finalize_reselection(out)
   end function bigvar_varx_reselect

   pure function bigvar_varx_reselect_dual(endogenous, exogenous, ar_order, &
      exogenous_order, lambdas, alphas, structure, validation_first_origin, &
      first_origin, last_origin, horizon, loss_type, delta, one_standard_error, &
      selection_window, window_size, tolerance, max_iterations, gamma, &
      refit_fraction, contemporaneous, minnesota_target, include_intercept) &
      result(out)
      !! Reselect an alpha-specific VARX lambda surface at every origin.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: validation_first_origin !! Validation first origin.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: selection_window !! Selection window.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_reselection_t) :: out
      type(bigvar_validation_t) :: selection, evaluation
      integer :: origin, index, selection_first, selection_last

      call initialize_reselection(out, first_origin, last_origin, horizon, &
         validation_first_origin, size(endogenous, 2), 1, loss_type)
      if (out%info /= 0 .or. .not. valid_dual_grid(lambdas, alphas, structure)) &
         then
         out%info = 1
         return
      end if
      do index = 1, size(out%valid)
         origin = first_origin + index - 1
         selection_last = origin - horizon
         selection_first = validation_first_origin
         if (present(selection_window)) then
            if (selection_window > 0) then
               selection_first = max(selection_first, &
                  selection_last - selection_window + 1)
            end if
         end if
         if (selection_last < selection_first) cycle
         selection = bigvar_varx_validate_dual(endogenous, exogenous, ar_order, &
            exogenous_order, lambdas, alphas, structure, selection_first, &
            selection_last, horizon, loss_type=loss_type, delta=delta, &
            one_standard_error=one_standard_error, window_size=window_size, &
            tolerance=tolerance, max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, contemporaneous=contemporaneous, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (selection%selected_index < 1) cycle
         evaluation = bigvar_varx_evaluate(endogenous, exogenous, ar_order, &
            exogenous_order, selection%selected_lambda, structure, origin, &
            origin, horizon, loss_type=loss_type, delta=delta, &
            window_size=window_size, alpha=selection%selected_alpha, &
            tolerance=tolerance, max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, contemporaneous=contemporaneous, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         call store_joint_reselection(out, index, selection, evaluation)
      end do
      call finalize_reselection(out)
   end function bigvar_varx_reselect_dual

   pure function bigvar_var_reselect_separate(series, lag_order, lambdas, &
      structure, validation_first_origin, first_origin, last_origin, horizon, &
      loss_type, delta, one_standard_error, selection_window, window_size, &
      alpha, tolerance, max_iterations, gamma, refit_fraction, recursive, &
      minnesota_target, include_intercept) result(out)
      !! Reselect response-specific VAR penalties before every forecast.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: validation_first_origin !! Validation first origin.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: selection_window !! Selection window.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_reselection_t) :: out
      type(bigvar_separate_validation_t) :: selection
      type(bigvar_fit_t) :: fit
      type(bigtime_forecast_t) :: forecast
      real(dp) :: selected_alpha, selected_delta
      integer :: origin, index, selection_first, selection_last, training_start
      integer :: selected_loss, selected_window, fit_horizon
      logical :: use_recursive

      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_window = 0
      if (present(window_size)) selected_window = window_size
      selected_alpha = 1.0_dp/real(size(series, 2) + 1, dp)
      if (present(alpha)) selected_alpha = alpha
      use_recursive = .true.
      if (present(recursive)) use_recursive = recursive
      fit_horizon = merge(1, horizon, use_recursive)
      call initialize_reselection(out, first_origin, last_origin, horizon, &
         validation_first_origin, size(series, 2), size(series, 2), loss_type)
      if (out%info /= 0) return
      out%selected_alpha = selected_alpha
      do index = 1, size(out%valid)
         origin = first_origin + index - 1
         selection_last = origin - horizon
         selection_first = validation_first_origin
         if (present(selection_window)) then
            if (selection_window > 0) then
               selection_first = max(selection_first, &
                  selection_last - selection_window + 1)
            end if
         end if
         if (selection_last < selection_first) cycle
         selection = bigvar_var_validate_separate(series, lag_order, lambdas, &
            structure, selection_first, selection_last, horizon, &
            loss_type=selected_loss, delta=selected_delta, &
            one_standard_error=one_standard_error, window_size=selected_window, &
            alpha=selected_alpha, tolerance=tolerance, &
            max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, recursive=use_recursive, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (.not. allocated(selection%selected_lambda)) cycle
         training_start = 1
         if (selected_window > 0) then
            if (horizon > 1 .and. .not. use_recursive) then
               training_start = max(1, origin - selected_window - horizon + 1)
            else
               training_start = max(1, origin - selected_window + 1)
            end if
         end if
         fit = bigvar_structured_var_separate(series(training_start:origin, :), &
            lag_order, selection%selected_lambda, structure, &
            tolerance=tolerance, max_iterations=max_iterations, &
            alpha=selected_alpha, gamma=gamma, refit_fraction=refit_fraction, &
            direct_horizon=fit_horizon, minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (fit%info /= 0 .and. fit%info /= 4) cycle
         forecast = bigvar_forecast(fit, series(training_start:origin, :), &
            horizon)
         if (forecast%info /= 0) cycle
         out%forecasts(index, :) = forecast%mean(horizon, :)
         out%loss(index, :) = abs(series(origin + horizon, :) - &
            out%forecasts(index, :))
         if (selected_loss == bigvar_loss_l2) out%loss(index, :) = &
            (series(origin + horizon, :) - out%forecasts(index, :))**2
         if (selected_loss == bigvar_loss_huber) out%loss(index, :) = &
            huber_component(series(origin + horizon, :) - &
            out%forecasts(index, :), selected_delta)
         out%selected_lambda(index, :) = selection%selected_lambda
         out%selected_index(index, :) = selection%selected_index
         out%nonzero(index) = fit%nonzero
         out%fit_info(index) = fit%info
         out%converged(index) = fit%converged
         out%valid(index) = .true.
      end do
      call finalize_reselection(out)
   end function bigvar_var_reselect_separate

   pure function bigvar_least_squares_varx(endogenous, exogenous, ar_order, &
      exogenous_order, forecast_horizon, information_criterion, &
      contemporaneous, include_intercept) result(out)
      !! Fit a VARX model by stable least squares with an optional intercept.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in), optional :: forecast_horizon !! Forecast horizon.
      integer, intent(in), optional :: information_criterion !! Information criterion.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_ls_varx_fit_t) :: out
      real(dp), allocatable :: response(:, :), predictors(:, :), design(:, :)
      real(dp), allocatable :: normal_matrix(:, :), inverse(:, :)
      real(dp), allocatable :: coefficient(:, :), covariance_inverse(:, :)
      real(dp) :: logdet, sample_size, penalty
      integer :: variables, x_variables, predictor_count, design_count, x_blocks
      integer :: cases, status, criterion_type, horizon
      logical :: use_contemporaneous, use_intercept

      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      horizon = 1
      if (present(forecast_horizon)) horizon = forecast_horizon
      criterion_type = 0
      if (present(information_criterion)) criterion_type = &
         information_criterion
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      use_intercept = .true.
      if (present(include_intercept)) use_intercept = include_intercept
      out%ar_order = ar_order
      out%exogenous_order = exogenous_order
      out%transfer_function = ar_order == 0
      out%forecast_horizon = horizon
      out%direct = horizon > 1
      out%information_criterion = criterion_type
      out%contemporaneous = use_contemporaneous
      out%include_intercept = use_intercept
      x_blocks = exogenous_order
      if (use_contemporaneous) x_blocks = x_blocks + 1
      if (size(endogenous, 1) /= size(exogenous, 1) .or. variables < 1 .or. &
         ar_order < 0 .or. exogenous_order < 0 .or. horizon < 1 .or. &
         (x_blocks > 0 .and. x_variables < 1) .or. &
         (criterion_type /= 0 .and. criterion_type /= bigvar_ic_aic .and. &
         criterion_type /= bigvar_ic_bic)) then
         out%info = 1
         return
      end if
      call build_ls_varx_data(endogenous, exogenous, ar_order, &
         exogenous_order, horizon, use_contemporaneous, response, predictors)
      cases = size(response, 1)
      predictor_count = size(predictors, 2)
      if (cases < 2 .or. predictor_count >= cases) then
         out%info = 1
         return
      end if
      design_count = predictor_count + merge(1, 0, use_intercept)
      allocate(design(cases, design_count))
      if (predictor_count > 0) design(:, :predictor_count) = predictors
      if (use_intercept) design(:, design_count) = 1.0_dp
      allocate(coefficient(design_count, variables))
      coefficient = 0.0_dp
      if (design_count > 0) then
         normal_matrix = matmul(transpose(design), design)
         allocate(inverse(design_count, design_count))
         call symmetric_pseudoinverse(normal_matrix, inverse, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         coefficient = matmul(inverse, matmul(transpose(design), response))
      end if
      allocate(out%phi(variables, variables*ar_order))
      allocate(out%beta(variables, x_variables*x_blocks))
      allocate(out%intercept(variables))
      if (ar_order > 0) out%phi = transpose( &
         coefficient(:variables*ar_order, :))
      if (x_blocks > 0) out%beta = transpose(coefficient( &
         variables*ar_order + 1:predictor_count, :))
      out%intercept = 0.0_dp
      if (use_intercept) out%intercept = coefficient(design_count, :)
      allocate(out%fitted(cases, variables), out%residuals(cases, variables))
      out%fitted = 0.0_dp
      if (design_count > 0) out%fitted = matmul(design, coefficient)
      out%residuals = response - out%fitted
      allocate(out%innovation_covariance(variables, variables))
      out%innovation_covariance = matmul(transpose(out%residuals), &
         out%residuals)/real(cases, dp)
      if (criterion_type /= 0) then
         allocate(covariance_inverse(variables, variables))
         call inverse_logdet(out%innovation_covariance, covariance_inverse, &
            logdet, status, 100.0_dp*epsilon(1.0_dp))
         if (status == 0) then
            sample_size = real(cases, dp)
            penalty = real(variables*predictor_count + &
               merge(variables, 0, use_intercept), dp)
            if (criterion_type == bigvar_ic_aic) then
               out%criterion = logdet + 2.0_dp*penalty/sample_size
            else
               out%criterion = logdet + log(sample_size)*penalty/sample_size
            end if
         end if
      end if
      out%info = 0
   end function bigvar_least_squares_varx

   pure function bigvar_varx_ic_select(endogenous, exogenous, max_ar_order, &
      max_exogenous_order, information_criterion, forecast_horizon, &
      contemporaneous, include_intercept) result(out)
      !! Select least-squares VARX lag orders by AIC or BIC.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      integer, intent(in) :: max_ar_order !! Maximum autoregressive order.
      integer, intent(in) :: max_exogenous_order !! Maximum exogenous order.
      integer, intent(in) :: information_criterion !! Information criterion.
      integer, intent(in), optional :: forecast_horizon !! Forecast horizon.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_ic_search_t) :: out
      type(bigvar_ls_varx_fit_t) :: candidate
      real(dp) :: best_value
      integer :: ar_order, exogenous_order, horizon
      logical :: use_contemporaneous

      horizon = 1
      if (present(forecast_horizon)) horizon = forecast_horizon
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      out%information_criterion = information_criterion
      if (max_ar_order < 0 .or. max_exogenous_order < 0 .or. horizon < 1 .or. &
         (information_criterion /= bigvar_ic_aic .and. &
         information_criterion /= bigvar_ic_bic)) then
         out%info = 1
         return
      end if
      allocate(out%criterion(0:max_ar_order, 0:max_exogenous_order))
      out%criterion = huge(1.0_dp)
      best_value = huge(1.0_dp)
      do ar_order = 0, max_ar_order
         do exogenous_order = 0, max_exogenous_order
            candidate = bigvar_least_squares_varx(endogenous, exogenous, &
               ar_order, exogenous_order, horizon, information_criterion, &
               use_contemporaneous, include_intercept)
            if (candidate%info /= 0) cycle
            out%criterion(ar_order, exogenous_order) = candidate%criterion
            if (candidate%criterion < best_value) then
               best_value = candidate%criterion
               out%selected_ar_order = ar_order
               out%selected_exogenous_order = exogenous_order
               out%fit = candidate
            end if
         end do
      end do
      if (best_value == huge(1.0_dp)) then
         out%info = 2
      else
         out%info = 0
      end if
   end function bigvar_varx_ic_select

   pure function bigvar_ls_varx_forecast(fit, endogenous_history, &
      exogenous_values, horizon) result(out)
      !! Forecast a least-squares VARX fit recursively or at its direct horizon.
      type(bigvar_ls_varx_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: endogenous_history(:, :) !! Endogenous history.
      real(dp), intent(in) :: exogenous_values(:, :) !! Exogenous values.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bigtime_forecast_t) :: out
      real(dp), allocatable :: working(:, :)
      integer :: observations, variables, x_variables, step, lag, target
      integer :: block, exogenous_index

      if (.not. allocated(fit%phi) .or. .not. allocated(fit%beta) .or. &
         .not. allocated(fit%intercept) .or. horizon < 1) then
         out%info = 1
         return
      end if
      observations = size(endogenous_history, 1)
      variables = size(endogenous_history, 2)
      x_variables = size(exogenous_values, 2)
      if (size(fit%phi, 1) /= variables .or. &
         size(fit%phi, 2) /= variables*fit%ar_order .or. &
         size(fit%beta, 1) /= variables .or. &
         observations < fit%ar_order .or. &
         (fit%direct .and. horizon /= fit%forecast_horizon)) then
         out%info = 1
         return
      end if
      if (fit%contemporaneous) then
         if (size(fit%beta, 2) /= &
            x_variables*(fit%exogenous_order + 1) .or. &
            size(exogenous_values, 1) < observations + horizon) then
            out%info = 1
            return
         end if
      else if (size(fit%beta, 2) /= x_variables*fit%exogenous_order .or. &
         (fit%exogenous_order > 0 .and. size(exogenous_values, 1) < &
         observations + merge(0, horizon - 1, fit%direct))) then
         out%info = 1
         return
      end if
      allocate(out%mean(horizon, variables))
      out%mean = 0.0_dp
      allocate(working(observations + horizon, variables))
      working(:observations, :) = endogenous_history
      if (fit%direct) then
         target = observations + horizon
         working(target, :) = fit%intercept
         do lag = 1, fit%ar_order
            working(target, :) = working(target, :) + matmul( &
               fit%phi(:, (lag - 1)*variables + 1:lag*variables), &
               endogenous_history(observations + 1 - lag, :))
         end do
         do lag = 0, fit%exogenous_order
            if (.not. fit%contemporaneous .and. lag == 0) cycle
            block = lag
            if (fit%contemporaneous) block = block + 1
            if (fit%contemporaneous) then
               exogenous_index = target - lag
            else
               exogenous_index = observations + 1 - lag
            end if
            working(target, :) = working(target, :) + matmul(fit%beta(:, &
               (block - 1)*x_variables + 1:block*x_variables), &
               exogenous_values(exogenous_index, :))
         end do
         out%mean(horizon, :) = working(target, :)
      else
         do step = 1, horizon
            target = observations + step
            working(target, :) = fit%intercept
            do lag = 1, fit%ar_order
               working(target, :) = working(target, :) + matmul(fit%phi(:, &
                  (lag - 1)*variables + 1:lag*variables), &
                  working(target - lag, :))
            end do
            do lag = 0, fit%exogenous_order
               if (.not. fit%contemporaneous .and. lag == 0) cycle
               block = lag
               if (fit%contemporaneous) block = block + 1
               working(target, :) = working(target, :) + matmul(fit%beta(:, &
                  (block - 1)*x_variables + 1:block*x_variables), &
                  exogenous_values(target - lag, :))
            end do
            out%mean(step, :) = working(target, :)
         end do
      end if
      out%info = 0
   end function bigvar_ls_varx_forecast

   pure function bigvar_varx_ic_evaluate(endogenous, exogenous, &
      max_ar_order, max_exogenous_order, information_criterion, first_origin, &
      last_origin, horizon, iterated, loss_type, delta, contemporaneous, &
      include_intercept) result(out)
      !! Evaluate rolling forecasts from AIC- or BIC-selected VARX models.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      integer, intent(in) :: max_ar_order !! Maximum autoregressive order.
      integer, intent(in) :: max_exogenous_order !! Maximum exogenous order.
      integer, intent(in) :: information_criterion !! Information criterion.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      logical, intent(in), optional :: iterated !! Flag controlling iterated.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      integer, intent(in), optional :: loss_type !! Loss type.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      type(bigvar_ic_evaluation_t) :: out
      type(bigvar_ic_search_t) :: search
      type(bigtime_forecast_t) :: forecast
      real(dp) :: selected_delta
      integer :: origins, index, origin, fit_horizon, selected_loss
      logical :: use_iterated, use_contemporaneous

      use_iterated = .false.
      if (present(iterated)) use_iterated = iterated
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      origins = last_origin - first_origin + 1
      out%first_origin = first_origin
      out%last_origin = last_origin
      out%horizon = horizon
      out%information_criterion = information_criterion
      out%loss_type = selected_loss
      if (size(endogenous, 1) /= size(exogenous, 1) .or. origins < 1 .or. &
         horizon < 1 .or. first_origin < 2 .or. &
         last_origin + horizon > size(endogenous, 1) .or. &
         selected_loss < bigvar_loss_l1 .or. &
         selected_loss > bigvar_loss_huber .or. selected_delta <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%loss(origins), out%forecasts(origins, &
         size(endogenous, 2)))
      allocate(out%ar_order(origins), out%exogenous_order(origins))
      allocate(out%valid(origins))
      out%loss = huge(1.0_dp)
      out%forecasts = 0.0_dp
      out%ar_order = 0
      out%exogenous_order = 0
      out%valid = .false.
      fit_horizon = horizon
      if (use_iterated) fit_horizon = 1
      do index = 1, origins
         origin = first_origin + index - 1
         search = bigvar_varx_ic_select(endogenous(:origin, :), &
            exogenous(:origin, :), max_ar_order, max_exogenous_order, &
            information_criterion, fit_horizon, use_contemporaneous, &
            include_intercept)
         if (search%info /= 0) cycle
         forecast = bigvar_ls_varx_forecast(search%fit, &
            endogenous(:origin, :), exogenous, horizon)
         if (forecast%info /= 0) cycle
         out%forecasts(index, :) = forecast%mean(horizon, :)
         out%loss(index) = bigvar_forecast_loss( &
            endogenous(origin + horizon, :) - out%forecasts(index, :), &
            selected_loss, selected_delta)
         out%ar_order(index) = search%selected_ar_order
         out%exogenous_order(index) = search%selected_exogenous_order
         out%valid(index) = .true.
      end do
      if (all(out%valid)) then
         out%info = 0
      else
         out%info = 2
      end if
   end function bigvar_varx_ic_evaluate

   pure function bigvar_mean_benchmark(series, first_origin, last_origin, &
      horizon, loss_type, delta, window_size) result(out)
      !! Evaluate expanding- or fixed-window unconditional-mean forecasts.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: window_size !! Window size.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      type(bigvar_benchmark_t) :: out
      real(dp) :: selected_delta
      integer :: selected_loss, selected_window, origins
      integer :: index, origin, training_start

      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_window = 0
      if (present(window_size)) selected_window = window_size
      origins = last_origin - first_origin + 1
      out%first_origin = first_origin
      out%last_origin = last_origin
      out%horizon = horizon
      out%loss_type = selected_loss
      out%window_size = selected_window
      if (size(series, 2) < 1 .or. origins < 1 .or. horizon < 1 .or. &
         first_origin < 1 .or. last_origin + horizon > size(series, 1) .or. &
         selected_window < 0 .or. selected_loss < bigvar_loss_l1 .or. &
         selected_loss > bigvar_loss_huber .or. selected_delta <= 0.0_dp) then
         out%info = 1
         return
      end if
      call initialize_benchmark(out, origins, size(series, 2))
      do index = 1, origins
         origin = first_origin + index - 1
         training_start = 1
         if (selected_window > 0) then
            training_start = max(1, origin - selected_window + 1)
         end if
         out%forecasts(index, :) = sum(series(training_start:origin, :), &
            dim=1)/real(origin - training_start + 1, dp)
         out%loss(index) = bigvar_forecast_loss( &
            series(origin + horizon, :) - out%forecasts(index, :), &
            selected_loss, selected_delta)
         out%valid(index) = .true.
      end do
      call summarize_benchmark(out)
   end function bigvar_mean_benchmark

   pure function bigvar_random_walk_benchmark(series, first_origin, &
      last_origin, horizon, loss_type, delta) result(out)
      !! Evaluate no-drift random-walk forecasts over a range of origins.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      type(bigvar_benchmark_t) :: out
      real(dp) :: selected_delta
      integer :: selected_loss, origins, index, origin

      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      origins = last_origin - first_origin + 1
      out%first_origin = first_origin
      out%last_origin = last_origin
      out%horizon = horizon
      out%loss_type = selected_loss
      if (size(series, 2) < 1 .or. origins < 1 .or. horizon < 1 .or. &
         first_origin < 1 .or. last_origin + horizon > size(series, 1) .or. &
         selected_loss < bigvar_loss_l1 .or. &
         selected_loss > bigvar_loss_huber .or. selected_delta <= 0.0_dp) then
         out%info = 1
         return
      end if
      call initialize_benchmark(out, origins, size(series, 2))
      do index = 1, origins
         origin = first_origin + index - 1
         out%forecasts(index, :) = series(origin, :)
         out%loss(index) = bigvar_forecast_loss( &
            series(origin + horizon, :) - out%forecasts(index, :), &
            selected_loss, selected_delta)
         out%valid(index) = .true.
      end do
      call summarize_benchmark(out)
   end function bigvar_random_walk_benchmark

   pure function bigvar_var_validate_loo(series, lag_order, lambdas, &
      structure, loss_type, delta, one_standard_error, alphas, tolerance, &
      max_iterations, gamma, refit_fraction, horizon, recursive, &
      first_observation, last_observation, minnesota_target, &
      include_intercept) result(out)
      !! Select structured VAR penalties by leave-one-out cross-validation.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: first_observation !! First observation.
      integer, intent(in), optional :: last_observation !! Last observation.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alphas(:) !! Alphas.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      type(bigvar_path_t) :: path
      real(dp), allocatable :: reduced(:, :), alpha_grid(:), prediction(:)
      real(dp) :: selected_delta, selected_tolerance, selected_gamma
      integer :: selected_loss, selected_max_iterations, selected_horizon
      integer :: first_case, last_case, cases, alpha_count, lambda_count
      integer :: fold, omitted, alpha_index, lambda_index, candidate
      integer :: fit_horizon
      logical :: use_one_se, use_recursive

      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      selected_horizon = 1
      if (present(horizon)) selected_horizon = horizon
      use_recursive = .true.
      if (present(recursive)) use_recursive = recursive
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      first_case = lag_order + selected_horizon
      if (present(first_observation)) first_case = first_observation
      last_case = size(series, 1)
      if (present(last_observation)) last_case = last_observation
      lambda_count = size(lambdas)
      alpha_count = 1
      if (present(alphas)) alpha_count = size(alphas)
      allocate(alpha_grid(alpha_count))
      if (present(alphas)) then
         alpha_grid = alphas
      else
         alpha_grid = 1.0_dp/real(size(series, 2) + 1, dp)
      end if
      cases = last_case - first_case + 1
      if (size(series, 2) < 1 .or. lag_order < 1 .or. lambda_count < 1 .or. &
         alpha_count < 1 .or. cases < 1 .or. selected_horizon < 1 .or. &
         first_case < lag_order + selected_horizon .or. &
         last_case > size(series, 1) .or. size(series, 1) - 1 <= lag_order .or. &
         selected_loss < bigvar_loss_l1 .or. &
         selected_loss > bigvar_loss_huber .or. selected_delta <= 0.0_dp .or. &
         any(lambdas < 0.0_dp) .or. any(alpha_grid < 0.0_dp) .or. &
         any(alpha_grid > 1.0_dp)) then
         out%info = 1
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= size(series, 2) .or. &
            structure == bigvar_structure_bgr) then
            out%info = 1
            return
         end if
      end if
      call initialize_validation(out, cases, lambda_count*alpha_count, &
         size(series, 2), first_case, last_case, selected_horizon, &
         selected_loss)
      call fill_candidate_grid(out, lambdas, alpha_grid)
      fit_horizon = selected_horizon
      if (use_recursive) fit_horizon = 1
      do fold = 1, cases
         omitted = first_case + fold - 1
         reduced = matrix_without_row(series, omitted)
         do alpha_index = 1, alpha_count
            if (present(refit_fraction)) then
               path = bigvar_structured_var_path(reduced, lag_order, lambdas, &
                  structure, selected_tolerance, selected_max_iterations, &
                  alpha=alpha_grid(alpha_index), gamma=selected_gamma, &
                  refit_fraction=refit_fraction, direct_horizon=fit_horizon, &
                  minnesota_target=minnesota_target, &
                  include_intercept=include_intercept)
            else
               if (present(minnesota_target)) then
                  path = bigvar_structured_var_path(reduced, lag_order, &
                     lambdas, structure, selected_tolerance, &
                     selected_max_iterations, alpha=alpha_grid(alpha_index), &
                     gamma=selected_gamma, direct_horizon=fit_horizon, &
                     minnesota_target=minnesota_target, &
                     include_intercept=include_intercept)
               else
                  path = bigvar_structured_var_path(reduced, lag_order, &
                     lambdas, structure, selected_tolerance, &
                     selected_max_iterations, alpha=alpha_grid(alpha_index), &
                     gamma=selected_gamma, direct_horizon=fit_horizon, &
                     include_intercept=include_intercept)
               end if
            end if
            if (path%info /= 0 .and. path%info /= 4) cycle
            do lambda_index = 1, lambda_count
               candidate = lambda_index + lambda_count*(alpha_index - 1)
               prediction = predict_var_coefficients( &
                  path%phi(:, :, lambda_index), &
                  path%intercept(:, lambda_index), lag_order, &
                  series(:omitted - selected_horizon, :), selected_horizon, &
                  use_recursive)
               out%forecasts(fold, candidate, :) = prediction
               out%loss(fold, candidate) = bigvar_forecast_loss( &
                  series(omitted, :) - prediction, selected_loss, &
                  selected_delta)
               out%nonzero(fold, candidate) = path%nonzero(lambda_index)
               out%valid(fold, candidate) = .true.
            end do
            if (path%info == 4) out%info = 4
         end do
      end do
      call summarize_validation(out, use_one_se)
   end function bigvar_var_validate_loo

   pure function bigvar_varx_validate_loo(endogenous, exogenous, ar_order, &
      exogenous_order, lambdas, structure, loss_type, delta, &
      one_standard_error, alphas, tolerance, max_iterations, gamma, &
      refit_fraction, horizon, contemporaneous, first_observation, &
      last_observation, minnesota_target, include_intercept) result(out)
      !! Select structured VARX penalties by leave-one-out cross-validation.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: first_observation !! First observation.
      integer, intent(in), optional :: last_observation !! Last observation.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alphas(:) !! Alphas.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      type(bigvar_varx_path_t) :: path
      real(dp), allocatable :: reduced_y(:, :), reduced_x(:, :)
      real(dp), allocatable :: alpha_grid(:), prediction(:)
      real(dp) :: selected_delta, selected_tolerance, selected_gamma
      integer :: selected_loss, selected_max_iterations, selected_horizon
      integer :: first_case, last_case, cases, alpha_count, lambda_count
      integer :: fold, omitted, alpha_index, lambda_index, candidate
      integer :: required_order
      logical :: use_one_se, use_contemporaneous

      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      selected_horizon = 1
      if (present(horizon)) selected_horizon = horizon
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      required_order = max(ar_order, exogenous_order)
      first_case = required_order + selected_horizon
      if (present(first_observation)) first_case = first_observation
      last_case = size(endogenous, 1)
      if (present(last_observation)) last_case = last_observation
      lambda_count = size(lambdas)
      alpha_count = 1
      if (present(alphas)) alpha_count = size(alphas)
      allocate(alpha_grid(alpha_count))
      if (present(alphas)) then
         alpha_grid = alphas
      else
         alpha_grid = 1.0_dp/real(size(endogenous, 2) + 1, dp)
      end if
      cases = last_case - first_case + 1
      if (size(endogenous, 1) /= size(exogenous, 1) .or. &
         size(endogenous, 2) < 1 .or. size(exogenous, 2) < 1 .or. &
         ar_order < 0 .or. exogenous_order < 0 .or. &
         (.not. use_contemporaneous .and. exogenous_order < 1) .or. &
         lambda_count < 1 .or. alpha_count < 1 .or. cases < 1 .or. &
         selected_horizon < 1 .or. &
         first_case < required_order + selected_horizon .or. &
         last_case > size(endogenous, 1) .or. &
         size(endogenous, 1) - 1 <= required_order .or. &
         selected_loss < bigvar_loss_l1 .or. &
         selected_loss > bigvar_loss_huber .or. selected_delta <= 0.0_dp .or. &
         any(lambdas < 0.0_dp) .or. any(alpha_grid < 0.0_dp) .or. &
         any(alpha_grid > 1.0_dp) .or. &
         (structure == bigvar_structure_efx .and. use_contemporaneous)) then
         out%info = 1
         return
      end if
      if (ar_order == 0 .and. .not. supports_transfer_function(structure)) then
         out%info = 1
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= size(endogenous, 2) .or. &
            ar_order == 0) then
            out%info = 1
            return
         end if
      end if
      call initialize_validation(out, cases, lambda_count*alpha_count, &
         size(endogenous, 2), first_case, last_case, selected_horizon, &
         selected_loss)
      call fill_candidate_grid(out, lambdas, alpha_grid)
      do fold = 1, cases
         omitted = first_case + fold - 1
         reduced_y = matrix_without_row(endogenous, omitted)
         reduced_x = matrix_without_row(exogenous, omitted)
         do alpha_index = 1, alpha_count
            if (present(refit_fraction)) then
               path = bigvar_structured_varx_path(reduced_y, reduced_x, &
                  ar_order, exogenous_order, lambdas, structure, &
                  selected_tolerance, selected_max_iterations, &
                  alpha=alpha_grid(alpha_index), gamma=selected_gamma, &
                  refit_fraction=refit_fraction, &
                  contemporaneous=use_contemporaneous, &
                  minnesota_target=minnesota_target, &
                  include_intercept=include_intercept)
            else
               if (present(minnesota_target)) then
                  path = bigvar_structured_varx_path(reduced_y, reduced_x, &
                     ar_order, exogenous_order, lambdas, structure, &
                     selected_tolerance, selected_max_iterations, &
                     alpha=alpha_grid(alpha_index), gamma=selected_gamma, &
                     contemporaneous=use_contemporaneous, &
                     minnesota_target=minnesota_target, &
                     include_intercept=include_intercept)
               else
                  path = bigvar_structured_varx_path(reduced_y, reduced_x, &
                     ar_order, exogenous_order, lambdas, structure, &
                     selected_tolerance, selected_max_iterations, &
                     alpha=alpha_grid(alpha_index), gamma=selected_gamma, &
                     contemporaneous=use_contemporaneous, &
                     include_intercept=include_intercept)
               end if
            end if
            if (path%info /= 0 .and. path%info /= 4) cycle
            do lambda_index = 1, lambda_count
               candidate = lambda_index + lambda_count*(alpha_index - 1)
               prediction = predict_varx_coefficients( &
                  path%phi(:, :, lambda_index), &
                  path%beta(:, :, lambda_index), &
                  path%intercept(:, lambda_index), ar_order, &
                  exogenous_order, use_contemporaneous, &
                  endogenous(:omitted - selected_horizon, :), exogenous, &
                  selected_horizon)
               out%forecasts(fold, candidate, :) = prediction
               out%loss(fold, candidate) = bigvar_forecast_loss( &
                  endogenous(omitted, :) - prediction, selected_loss, &
                  selected_delta)
               out%nonzero(fold, candidate) = path%nonzero(lambda_index)
               out%valid(fold, candidate) = .true.
            end do
            if (path%info == 4) out%info = 4
         end do
      end do
      call summarize_validation(out, use_one_se)
   end function bigvar_varx_validate_loo

   pure function bigvar_var_validate_separate_loo(series, lag_order, &
      lambdas, structure, loss_type, delta, one_standard_error, alpha, &
      tolerance, max_iterations, gamma, refit_fraction, horizon, recursive, &
      first_observation, last_observation, minnesota_target, &
      include_intercept) result(out)
      !! Select response-specific VAR penalties by leave-one-out validation.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: first_observation !! First observation.
      integer, intent(in), optional :: last_observation !! Last observation.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_separate_validation_t) :: out
      type(bigvar_separate_path_t) :: path
      real(dp), allocatable :: reduced(:, :), prediction(:)
      real(dp) :: selected_delta, selected_alpha, selected_tolerance
      real(dp) :: selected_gamma
      integer :: selected_loss, selected_max_iterations, selected_horizon
      integer :: first_case, last_case, cases, candidates, variables
      integer :: fold, omitted, candidate, response, fit_horizon
      logical :: use_one_se, use_recursive

      variables = size(series, 2)
      candidates = size(lambdas, 1)
      selected_loss = bigvar_loss_l2
      if (present(loss_type)) selected_loss = loss_type
      selected_delta = 2.5_dp
      if (present(delta)) selected_delta = delta
      selected_alpha = 1.0_dp/real(variables + 1, dp)
      if (present(alpha)) selected_alpha = alpha
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_max_iterations = 10000
      if (present(max_iterations)) selected_max_iterations = max_iterations
      selected_gamma = 3.0_dp
      if (present(gamma)) selected_gamma = gamma
      selected_horizon = 1
      if (present(horizon)) selected_horizon = horizon
      use_recursive = .true.
      if (present(recursive)) use_recursive = recursive
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      first_case = lag_order + selected_horizon
      if (present(first_observation)) first_case = first_observation
      last_case = size(series, 1)
      if (present(last_observation)) last_case = last_observation
      cases = last_case - first_case + 1
      if (variables < 1 .or. lag_order < 1 .or. candidates < 1 .or. &
         size(lambdas, 2) /= variables .or. any(lambdas < 0.0_dp) .or. &
         cases < 1 .or. first_case < lag_order + selected_horizon .or. &
         last_case > size(series, 1) .or. selected_horizon < 1 .or. &
         selected_alpha < 0.0_dp .or. selected_alpha > 1.0_dp) then
         out%info = 1
         return
      end if
      if (present(minnesota_target)) then
         if (size(minnesota_target) /= variables) then
            out%info = 1
            return
         end if
      end if
      allocate(out%loss(cases, candidates, variables))
      allocate(out%forecasts(cases, candidates, variables))
      allocate(out%mean_loss(candidates, variables))
      allocate(out%standard_error(candidates, variables))
      allocate(out%lambda(candidates, variables))
      allocate(out%selected_index(variables), out%best_index(variables))
      allocate(out%one_se_index(variables), out%selected_lambda(variables))
      allocate(out%valid(cases, candidates, variables))
      out%loss = huge(1.0_dp)
      out%forecasts = 0.0_dp
      out%mean_loss = huge(1.0_dp)
      out%standard_error = huge(1.0_dp)
      out%lambda = lambdas
      out%selected_index = 0
      out%best_index = 0
      out%one_se_index = 0
      out%selected_lambda = 0.0_dp
      out%valid = .false.
      out%first_origin = first_case
      out%last_origin = last_case
      out%horizon = selected_horizon
      out%loss_type = selected_loss
      fit_horizon = selected_horizon
      if (use_recursive) fit_horizon = 1
      do fold = 1, cases
         omitted = first_case + fold - 1
         reduced = matrix_without_row(series, omitted)
         if (present(refit_fraction)) then
            path = bigvar_structured_var_separate_path(reduced, lag_order, &
               lambdas, structure, selected_tolerance, &
               selected_max_iterations, alpha=selected_alpha, &
               gamma=selected_gamma, refit_fraction=refit_fraction, &
               direct_horizon=fit_horizon, &
               minnesota_target=minnesota_target, &
               include_intercept=include_intercept)
         else
            if (present(minnesota_target)) then
               path = bigvar_structured_var_separate_path(reduced, lag_order, &
                  lambdas, structure, selected_tolerance, &
                  selected_max_iterations, alpha=selected_alpha, &
                  gamma=selected_gamma, direct_horizon=fit_horizon, &
                  minnesota_target=minnesota_target, &
                  include_intercept=include_intercept)
            else
               path = bigvar_structured_var_separate_path(reduced, lag_order, &
                  lambdas, structure, selected_tolerance, &
                  selected_max_iterations, alpha=selected_alpha, &
                  gamma=selected_gamma, direct_horizon=fit_horizon, &
                  include_intercept=include_intercept)
            end if
         end if
         if (path%info /= 0 .and. path%info /= 4) cycle
         do candidate = 1, candidates
            prediction = predict_var_coefficients(path%phi(:, :, candidate), &
               path%intercept(:, candidate), lag_order, &
               series(:omitted - selected_horizon, :), selected_horizon, &
               use_recursive)
            out%forecasts(fold, candidate, :) = prediction
            do response = 1, variables
               out%loss(fold, candidate, response) = bigvar_forecast_loss( &
                  [series(omitted, response) - prediction(response)], &
                  selected_loss, selected_delta)
               out%valid(fold, candidate, response) = .true.
            end do
         end do
         if (path%info == 4) out%info = 4
      end do
      call summarize_separate_validation(out, use_one_se)
   end function bigvar_var_validate_separate_loo

   pure function bigvar_var_validate_dual(series, lag_order, lambdas, alphas, &
      structure, first_origin, last_origin, horizon, loss_type, delta, &
      one_standard_error, window_size, tolerance, max_iterations, gamma, &
      refit_fraction, recursive, minnesota_target, include_intercept) &
      result(out)
      !! Validate alpha-specific VAR lambda paths over rolling forecast origins.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      type(bigvar_validation_t) :: partial
      logical :: use_one_se
      integer :: alpha_index, first_candidate, last_candidate

      if (.not. valid_dual_grid(lambdas, alphas, structure)) then
         out%info = 1
         return
      end if
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      do alpha_index = 1, size(alphas)
         partial = bigvar_var_validate(series, lag_order, &
            lambdas(:, alpha_index), structure, first_origin, last_origin, &
            horizon, loss_type=loss_type, delta=delta, &
            one_standard_error=.false., window_size=window_size, &
            alphas=[alphas(alpha_index)], tolerance=tolerance, &
            max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, recursive=recursive, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (.not. allocated(partial%loss)) then
            out%info = partial%info
            return
         end if
         if (alpha_index == 1) then
            call initialize_validation(out, size(partial%loss, 1), &
               size(lambdas), size(series, 2), first_origin, last_origin, &
               horizon, partial%loss_type)
            call configure_dual_grid(out, lambdas, alphas)
         end if
         first_candidate = 1 + size(lambdas, 1)*(alpha_index - 1)
         last_candidate = size(lambdas, 1)*alpha_index
         call copy_validation_candidates(out, partial, first_candidate, &
            last_candidate)
         if (partial%info == 4) out%info = 4
      end do
      call summarize_dual_validation(out, use_one_se)
   end function bigvar_var_validate_dual

   pure function bigvar_varx_validate_dual(endogenous, exogenous, ar_order, &
      exogenous_order, lambdas, alphas, structure, first_origin, last_origin, &
      horizon, loss_type, delta, one_standard_error, window_size, tolerance, &
      max_iterations, gamma, refit_fraction, contemporaneous, &
      minnesota_target, include_intercept) result(out)
      !! Validate alpha-specific VARX lambda paths over rolling forecast origins.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: window_size !! Window size.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      type(bigvar_validation_t) :: partial
      logical :: use_one_se
      integer :: alpha_index, first_candidate, last_candidate

      if (.not. valid_dual_grid(lambdas, alphas, structure)) then
         out%info = 1
         return
      end if
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      do alpha_index = 1, size(alphas)
         partial = bigvar_varx_validate(endogenous, exogenous, ar_order, &
            exogenous_order, lambdas(:, alpha_index), structure, first_origin, &
            last_origin, horizon, loss_type=loss_type, delta=delta, &
            one_standard_error=.false., window_size=window_size, &
            alphas=[alphas(alpha_index)], tolerance=tolerance, &
            max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, contemporaneous=contemporaneous, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (.not. allocated(partial%loss)) then
            out%info = partial%info
            return
         end if
         if (alpha_index == 1) then
            call initialize_validation(out, size(partial%loss, 1), &
               size(lambdas), size(endogenous, 2), first_origin, last_origin, &
               horizon, partial%loss_type)
            call configure_dual_grid(out, lambdas, alphas)
         end if
         first_candidate = 1 + size(lambdas, 1)*(alpha_index - 1)
         last_candidate = size(lambdas, 1)*alpha_index
         call copy_validation_candidates(out, partial, first_candidate, &
            last_candidate)
         if (partial%info == 4) out%info = 4
      end do
      call summarize_dual_validation(out, use_one_se)
   end function bigvar_varx_validate_dual

   pure function bigvar_var_validate_loo_dual(series, lag_order, lambdas, &
      alphas, structure, loss_type, delta, one_standard_error, tolerance, &
      max_iterations, gamma, refit_fraction, horizon, recursive, &
      first_observation, last_observation, minnesota_target, &
      include_intercept) result(out)
      !! Validate alpha-specific VAR lambda paths by leave-one-out selection.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: first_observation !! First observation.
      integer, intent(in), optional :: last_observation !! Last observation.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: recursive !! Flag controlling recursive.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      type(bigvar_validation_t) :: partial
      logical :: use_one_se
      integer :: alpha_index, first_candidate, last_candidate

      if (.not. valid_dual_grid(lambdas, alphas, structure)) then
         out%info = 1
         return
      end if
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      do alpha_index = 1, size(alphas)
         partial = bigvar_var_validate_loo(series, lag_order, &
            lambdas(:, alpha_index), structure, loss_type=loss_type, &
            delta=delta, one_standard_error=.false., &
            alphas=[alphas(alpha_index)], tolerance=tolerance, &
            max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, horizon=horizon, recursive=recursive, &
            first_observation=first_observation, &
            last_observation=last_observation, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (.not. allocated(partial%loss)) then
            out%info = partial%info
            return
         end if
         if (alpha_index == 1) then
            call initialize_validation(out, size(partial%loss, 1), &
               size(lambdas), size(series, 2), partial%first_origin, &
               partial%last_origin, partial%horizon, partial%loss_type)
            call configure_dual_grid(out, lambdas, alphas)
         end if
         first_candidate = 1 + size(lambdas, 1)*(alpha_index - 1)
         last_candidate = size(lambdas, 1)*alpha_index
         call copy_validation_candidates(out, partial, first_candidate, &
            last_candidate)
         if (partial%info == 4) out%info = 4
      end do
      call summarize_dual_validation(out, use_one_se)
   end function bigvar_var_validate_loo_dual

   pure function bigvar_varx_validate_loo_dual(endogenous, exogenous, ar_order, &
      exogenous_order, lambdas, alphas, structure, loss_type, delta, &
      one_standard_error, tolerance, max_iterations, gamma, refit_fraction, &
      horizon, contemporaneous, first_observation, last_observation, &
      minnesota_target, include_intercept) result(out)
      !! Validate alpha-specific VARX paths by leave-one-out selection.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: first_observation !! First observation.
      integer, intent(in), optional :: last_observation !! Last observation.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(in), optional :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in), optional :: refit_fraction !! Refit fraction.
      logical, intent(in), optional :: one_standard_error !! Flag controlling one standard error.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), intent(in), optional :: minnesota_target(:) !! Minnesota target.
      logical, intent(in), optional :: include_intercept !! Whether to include an intercept.
      type(bigvar_validation_t) :: out
      type(bigvar_validation_t) :: partial
      logical :: use_one_se
      integer :: alpha_index, first_candidate, last_candidate

      if (.not. valid_dual_grid(lambdas, alphas, structure)) then
         out%info = 1
         return
      end if
      use_one_se = .false.
      if (present(one_standard_error)) use_one_se = one_standard_error
      do alpha_index = 1, size(alphas)
         partial = bigvar_varx_validate_loo(endogenous, exogenous, ar_order, &
            exogenous_order, lambdas(:, alpha_index), structure, &
            loss_type=loss_type, delta=delta, one_standard_error=.false., &
            alphas=[alphas(alpha_index)], tolerance=tolerance, &
            max_iterations=max_iterations, gamma=gamma, &
            refit_fraction=refit_fraction, horizon=horizon, &
            contemporaneous=contemporaneous, &
            first_observation=first_observation, &
            last_observation=last_observation, &
            minnesota_target=minnesota_target, &
            include_intercept=include_intercept)
         if (.not. allocated(partial%loss)) then
            out%info = partial%info
            return
         end if
         if (alpha_index == 1) then
            call initialize_validation(out, size(partial%loss, 1), &
               size(lambdas), size(endogenous, 2), partial%first_origin, &
               partial%last_origin, partial%horizon, partial%loss_type)
            call configure_dual_grid(out, lambdas, alphas)
         end if
         first_candidate = 1 + size(lambdas, 1)*(alpha_index - 1)
         last_candidate = size(lambdas, 1)*alpha_index
         call copy_validation_candidates(out, partial, first_candidate, &
            last_candidate)
         if (partial%info == 4) out%info = 4
      end do
      call summarize_dual_validation(out, use_one_se)
   end function bigvar_varx_validate_loo_dual

   pure function bigvar_stability(fit, margin) result(out)
      !! Diagnose BigVAR companion roots through the shared bigtime routine.
      type(bigvar_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: margin !! Margin.
      type(bigtime_stability_t) :: out
      real(dp) :: selected_margin

      if (.not. allocated(fit%phi)) then
         out%info = 1
         return
      end if
      selected_margin = 0.0_dp
      if (present(margin)) selected_margin = margin
      out = bigtime_var_stability(fit%phi, fit%lag_order, selected_margin)
   end function bigvar_stability

   pure function bigvar_var_to_companion(phi, lag_order) result(companion)
      !! Convert lag-major VAR coefficients to multiple companion form.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp), allocatable :: companion(:, :)
      integer :: variables

      variables = size(phi, 1)
      if (variables < 1 .or. lag_order < 1 .or. &
         size(phi, 2) /= variables*lag_order) then
         allocate(companion(0, 0))
         return
      end if
      companion = bigtime_companion_matrix(phi, lag_order)
   end function bigvar_var_to_companion

   pure function bigvar_var_simulate_from_innovations(phi, innovations, &
      burnin, intercept, initial_state, innovation_covariance) result(out)
      !! Simulate a BigVAR recursion from supplied innovations.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: innovations(:, :) !! Model innovations.
      integer, intent(in), optional :: burnin !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: intercept(:) !! Model intercept.
      real(dp), intent(in), optional :: initial_state(:) !! Initial state vector.
      real(dp), intent(in), optional :: innovation_covariance(:, :) !! Innovation covariance matrix.
      type(bigvar_simulation_t) :: out
      type(bigtime_simulation_t) :: shared
      real(dp), allocatable :: selected_intercept(:), selected_state(:)
      integer :: variables, lag_order, discard

      variables = size(phi, 1)
      if (variables < 1 .or. mod(size(phi, 2), variables) /= 0) then
         out%info = 1
         return
      end if
      lag_order = size(phi, 2)/variables
      discard = 0
      if (present(burnin)) discard = burnin
      allocate(selected_intercept(variables))
      selected_intercept = 0.0_dp
      if (present(intercept)) then
         if (size(intercept) /= variables) then
            out%info = 1
            return
         end if
         selected_intercept = intercept
      end if
      allocate(selected_state(variables*lag_order))
      selected_state = 0.0_dp
      if (present(initial_state)) then
         if (size(initial_state) /= size(selected_state)) then
            out%info = 1
            return
         end if
         selected_state = initial_state
      end if
      if (present(innovation_covariance)) then
         if (any(shape(innovation_covariance) /= [variables, variables])) then
            out%info = 1
            return
         end if
      end if
      shared = bigtime_var_simulate_from_innovations(phi, selected_intercept, &
         innovations, selected_state, discard)
      out%info = shared%info
      if (shared%info /= 0) return
      out%series = shared%series
      out%innovations = shared%innovations
      out%initial_state = shared%initial_state
      out%intercept = shared%intercept
      out%phi = shared%phi
      out%companion = bigtime_companion_matrix(phi, lag_order)
      if (present(innovation_covariance)) then
         out%innovation_covariance = innovation_covariance
      end if
      out%burnin = discard
   end function bigvar_var_simulate_from_innovations

   function bigvar_var_simulate(phi, innovation_covariance, periods, burnin, &
      intercept, initial_state) result(out)
      !! Simulate a stationary Gaussian VAR with BigVAR's default burn-in.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: innovation_covariance(:, :) !! Innovation covariance matrix.
      integer, intent(in) :: periods !! Periods.
      integer, intent(in), optional :: burnin !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: intercept(:) !! Model intercept.
      real(dp), intent(in), optional :: initial_state(:) !! Initial state vector.
      type(bigvar_simulation_t) :: out
      type(bigtime_stability_t) :: stability
      real(dp), allocatable :: innovations(:, :), zero(:)
      integer :: variables, lag_order, discard, status

      variables = size(phi, 1)
      if (variables < 1 .or. mod(size(phi, 2), variables) /= 0 .or. &
         periods < 1 .or. &
         any(shape(innovation_covariance) /= [variables, variables])) then
         out%info = 1
         return
      end if
      lag_order = size(phi, 2)/variables
      discard = 500
      if (present(burnin)) discard = burnin
      if (lag_order < 1 .or. discard < 0) then
         out%info = 1
         return
      end if
      stability = bigtime_var_stability(phi, lag_order)
      if (stability%info /= 0 .or. .not. stability%stable) then
         out%info = 2
         return
      end if
      allocate(innovations(periods + discard, variables), zero(variables))
      zero = 0.0_dp
      call random_multivariate_normal_matrix(zero, innovation_covariance, &
         innovations, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out = bigvar_var_simulate_from_innovations(phi, innovations, discard, &
         intercept, initial_state, innovation_covariance)
   end function bigvar_var_simulate

   pure subroutine initialize_validation(out, origins, candidates, variables, &
      first_origin, last_origin, horizon, loss_type)
      !! Allocate and initialize a rolling validation result.
      type(bigvar_validation_t), intent(out) :: out !! Procedure result.
      integer, intent(in) :: origins !! Origins.
      integer, intent(in) :: candidates !! Candidates.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: loss_type !! Loss type.

      allocate(out%loss(origins, candidates))
      allocate(out%forecasts(origins, candidates, variables))
      allocate(out%mean_loss(candidates), out%standard_error(candidates))
      allocate(out%lambda(candidates), out%alpha(candidates))
      allocate(out%nonzero(origins, candidates))
      allocate(out%valid(origins, candidates))
      out%loss = huge(1.0_dp)
      out%forecasts = 0.0_dp
      out%mean_loss = huge(1.0_dp)
      out%standard_error = huge(1.0_dp)
      out%lambda = 0.0_dp
      out%alpha = 0.0_dp
      out%nonzero = 0
      out%valid = .false.
      out%first_origin = first_origin
      out%last_origin = last_origin
      out%horizon = horizon
      out%loss_type = loss_type
   end subroutine initialize_validation

   pure subroutine initialize_reselection(out, first_origin, last_origin, &
      horizon, validation_first_origin, variables, selections, loss_type)
      !! Allocate a rolling-reselection result and initialize its metadata.
      type(bigvar_reselection_t), intent(out) :: out !! Procedure result.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: last_origin !! Last origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: validation_first_origin !! Validation first origin.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: selections !! Selections.
      integer, intent(in), optional :: loss_type !! Loss type.
      integer :: origins

      origins = last_origin - first_origin + 1
      if (origins < 1 .or. horizon < 1 .or. validation_first_origin < 1 .or. &
         variables < 1 .or. selections < 1) then
         out%info = 1
         return
      end if
      allocate(out%loss(origins, selections))
      allocate(out%forecasts(origins, variables))
      allocate(out%selected_lambda(origins, selections))
      allocate(out%selected_alpha(origins))
      allocate(out%selected_index(origins, selections))
      allocate(out%selected_lambda_index(origins))
      allocate(out%selected_alpha_index(origins))
      allocate(out%nonzero(origins), out%fit_info(origins))
      allocate(out%converged(origins), out%valid(origins))
      out%loss = huge(1.0_dp)
      out%forecasts = 0.0_dp
      out%selected_lambda = 0.0_dp
      out%selected_alpha = 0.0_dp
      out%selected_index = 0
      out%selected_lambda_index = 0
      out%selected_alpha_index = 0
      out%nonzero = 0
      out%fit_info = 0
      out%converged = .false.
      out%valid = .false.
      out%validation_first_origin = validation_first_origin
      out%first_origin = first_origin
      out%last_origin = last_origin
      out%horizon = horizon
      out%loss_type = bigvar_loss_l2
      if (present(loss_type)) out%loss_type = loss_type
      if (out%loss_type < bigvar_loss_l1 .or. &
         out%loss_type > bigvar_loss_huber) out%info = 1
   end subroutine initialize_reselection

   pure subroutine store_joint_reselection(out, index, selection, evaluation)
      !! Store one jointly selected penalty and its evaluation forecast.
      type(bigvar_reselection_t), intent(inout) :: out !! Procedure result, updated in place.
      integer, intent(in) :: index !! Element or observation index.
      type(bigvar_validation_t), intent(in) :: selection !! Selection.
      type(bigvar_validation_t), intent(in) :: evaluation !! Evaluation.

      if (.not. allocated(evaluation%valid)) return
      if (.not. evaluation%valid(1, 1)) return
      out%loss(index, 1) = evaluation%loss(1, 1)
      out%forecasts(index, :) = evaluation%forecasts(1, 1, :)
      out%selected_lambda(index, 1) = selection%selected_lambda
      out%selected_alpha(index) = selection%selected_alpha
      out%selected_index(index, 1) = selection%selected_index
      out%selected_lambda_index(index) = selection%selected_lambda_index
      out%selected_alpha_index(index) = selection%selected_alpha_index
      out%nonzero(index) = evaluation%nonzero(1, 1)
      out%fit_info(index) = evaluation%info
      out%converged(index) = evaluation%info == 0
      out%valid(index) = .true.
   end subroutine store_joint_reselection

   pure subroutine finalize_reselection(out)
      !! Set the aggregate status after all rolling-reselection origins.
      type(bigvar_reselection_t), intent(inout) :: out !! Procedure result, updated in place.

      if (.not. any(out%valid)) then
         out%info = 2
      else if (all(out%valid) .and. all(out%converged)) then
         out%info = 0
      else
         out%info = 4
      end if
   end subroutine finalize_reselection

   pure subroutine initialize_benchmark(out, origins, variables)
      !! Allocate benchmark forecasts and losses while retaining metadata.
      type(bigvar_benchmark_t), intent(inout) :: out !! Procedure result, updated in place.
      integer, intent(in) :: origins !! Origins.
      integer, intent(in) :: variables !! Number or indices of variables.

      allocate(out%loss(origins), out%forecasts(origins, variables))
      allocate(out%valid(origins))
      out%loss = huge(1.0_dp)
      out%forecasts = 0.0_dp
      out%valid = .false.
   end subroutine initialize_benchmark

   pure subroutine summarize_benchmark(out)
      !! Compute the mean loss and its standard error over valid origins.
      type(bigvar_benchmark_t), intent(inout) :: out !! Procedure result, updated in place.
      real(dp) :: sum_squares
      integer :: valid_count, index

      valid_count = count(out%valid)
      if (valid_count < 1) then
         out%info = 2
         return
      end if
      out%mean_loss = 0.0_dp
      do index = 1, size(out%loss)
         if (out%valid(index)) out%mean_loss = out%mean_loss + out%loss(index)
      end do
      out%mean_loss = out%mean_loss/real(valid_count, dp)
      out%standard_error = 0.0_dp
      if (valid_count > 1) then
         sum_squares = 0.0_dp
         do index = 1, size(out%loss)
            if (out%valid(index)) then
               sum_squares = sum_squares + &
                  (out%loss(index) - out%mean_loss)**2
            end if
         end do
         out%standard_error = sqrt(sum_squares/real(valid_count - 1, dp))/ &
            sqrt(real(valid_count, dp))
      end if
      out%info = 0
   end subroutine summarize_benchmark

   pure subroutine fill_candidate_grid(out, lambdas, alphas)
      !! Store candidates with lambda varying fastest within each alpha.
      type(bigvar_validation_t), intent(inout) :: out !! Procedure result, updated in place.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer :: alpha_index, lambda_index, candidate

      do alpha_index = 1, size(alphas)
         do lambda_index = 1, size(lambdas)
            candidate = lambda_index + size(lambdas)*(alpha_index - 1)
            out%lambda(candidate) = lambdas(lambda_index)
            out%alpha(candidate) = alphas(alpha_index)
         end do
      end do
   end subroutine fill_candidate_grid

   pure subroutine configure_dual_grid(out, lambdas, alphas)
      !! Store an alpha-specific lambda matrix and its flattened candidates.
      type(bigvar_validation_t), intent(inout) :: out !! Procedure result, updated in place.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer :: alpha_index, lambda_index, candidate

      allocate(out%lambda_grid(size(lambdas, 1), size(lambdas, 2)))
      allocate(out%mean_loss_surface(size(lambdas, 1), size(lambdas, 2)))
      allocate(out%standard_error_surface(size(lambdas, 1), size(lambdas, 2)))
      out%lambda_grid = lambdas
      out%mean_loss_surface = huge(1.0_dp)
      out%standard_error_surface = huge(1.0_dp)
      do alpha_index = 1, size(alphas)
         do lambda_index = 1, size(lambdas, 1)
            candidate = lambda_index + size(lambdas, 1)*(alpha_index - 1)
            out%lambda(candidate) = lambdas(lambda_index, alpha_index)
            out%alpha(candidate) = alphas(alpha_index)
         end do
      end do
   end subroutine configure_dual_grid

   pure subroutine copy_validation_candidates(out, partial, first_candidate, &
      last_candidate)
      !! Copy one alpha column of validation results into a combined surface.
      type(bigvar_validation_t), intent(inout) :: out !! Procedure result, updated in place.
      type(bigvar_validation_t), intent(in) :: partial !! Partial.
      integer, intent(in) :: first_candidate !! First candidate.
      integer, intent(in) :: last_candidate !! Last candidate.

      out%loss(:, first_candidate:last_candidate) = partial%loss
      out%forecasts(:, first_candidate:last_candidate, :) = partial%forecasts
      out%nonzero(:, first_candidate:last_candidate) = partial%nonzero
      out%valid(:, first_candidate:last_candidate) = partial%valid
   end subroutine copy_validation_candidates

   pure subroutine summarize_dual_validation(out, use_one_se)
      !! Summarize a flattened dual grid and recover its surface indices.
      type(bigvar_validation_t), intent(inout) :: out !! Procedure result, updated in place.
      logical, intent(in) :: use_one_se !! Whether to use the one se.
      integer :: lambda_count

      call summarize_validation(out, use_one_se)
      if (out%selected_index < 1) return
      lambda_count = size(out%lambda_grid, 1)
      out%mean_loss_surface = reshape(out%mean_loss, &
         shape(out%mean_loss_surface))
      out%standard_error_surface = reshape(out%standard_error, &
         shape(out%standard_error_surface))
      call split_candidate_index(out%best_index, lambda_count, &
         out%best_lambda_index, out%best_alpha_index)
      call split_candidate_index(out%one_se_index, lambda_count, &
         out%one_se_lambda_index, out%one_se_alpha_index)
      call split_candidate_index(out%selected_index, lambda_count, &
         out%selected_lambda_index, out%selected_alpha_index)
   end subroutine summarize_dual_validation

   pure subroutine split_candidate_index(candidate, lambda_count, lambda_index, &
      alpha_index)
      !! Convert a lambda-fast candidate index into two surface indices.
      integer, intent(in) :: candidate !! Candidate.
      integer, intent(in) :: lambda_count !! Number of lambda.
      integer, intent(out) :: lambda_index !! Index of lambda.
      integer, intent(out) :: alpha_index !! Index of alpha.

      lambda_index = modulo(candidate - 1, lambda_count) + 1
      alpha_index = (candidate - 1)/lambda_count + 1
   end subroutine split_candidate_index

   pure logical function valid_dual_grid(lambdas, alphas, structure) &
      result(value)
      !! Check dimensions and parameter ranges for alpha-specific lambda grids.
      real(dp), intent(in) :: lambdas(:, :) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alphas(:) !! Alphas.
      integer, intent(in) :: structure !! Model-structure specification.

      value = size(lambdas, 1) > 0 .and. &
         size(lambdas, 2) == size(alphas) .and. size(alphas) > 0 .and. &
         all(lambdas >= 0.0_dp) .and. all(alphas >= 0.0_dp) .and. &
         all(alphas <= 1.0_dp) .and. uses_alpha(structure)
   end function valid_dual_grid

   pure subroutine fill_lambda_grid(values, maximum_lambda, grid_ratio, linear)
      !! Fill a descending lambda grid with linear or geometric spacing.
      real(dp), intent(out) :: values(:) !! Input values.
      real(dp), intent(in) :: maximum_lambda !! Maximum lambda.
      real(dp), intent(in) :: grid_ratio !! Grid ratio.
      logical, intent(in) :: linear !! Flag controlling linear.
      real(dp) :: fraction, minimum_lambda
      integer :: index

      minimum_lambda = maximum_lambda/grid_ratio
      if (size(values) == 1) then
         values(1) = maximum_lambda
         return
      end if
      do index = 1, size(values)
         fraction = real(index - 1, dp)/real(size(values) - 1, dp)
         if (linear) then
            values(index) = maximum_lambda + &
               fraction*(minimum_lambda - maximum_lambda)
         else
            values(index) = maximum_lambda/grid_ratio**fraction
         end if
      end do
   end subroutine fill_lambda_grid

   pure subroutine summarize_validation(out, use_one_se)
      !! Compute validation summaries and BigVAR's one-standard-error choice.
      type(bigvar_validation_t), intent(inout) :: out !! Procedure result, updated in place.
      logical, intent(in) :: use_one_se !! Whether to use the one se.
      real(dp) :: mean_value, sum_squares, pooled_mean, threshold
      real(dp) :: minimum_loss, tie_tolerance
      integer :: candidate, observation, valid_count, pooled_count

      pooled_mean = 0.0_dp
      pooled_count = count(out%valid)
      if (pooled_count < 1) then
         out%info = 2
         return
      end if
      do candidate = 1, size(out%mean_loss)
         valid_count = count(out%valid(:, candidate))
         if (valid_count < 1) cycle
         mean_value = 0.0_dp
         do observation = 1, size(out%loss, 1)
            if (out%valid(observation, candidate)) then
               mean_value = mean_value + out%loss(observation, candidate)
               pooled_mean = pooled_mean + out%loss(observation, candidate)
            end if
         end do
         mean_value = mean_value/real(valid_count, dp)
         out%mean_loss(candidate) = mean_value
         if (valid_count > 1) then
            sum_squares = 0.0_dp
            do observation = 1, size(out%loss, 1)
               if (out%valid(observation, candidate)) then
                  sum_squares = sum_squares + &
                     (out%loss(observation, candidate) - mean_value)**2
               end if
            end do
            out%standard_error(candidate) = sqrt(sum_squares/ &
               real(valid_count - 1, dp))/sqrt(real(valid_count, dp))
         else
            out%standard_error(candidate) = 0.0_dp
         end if
      end do
      pooled_mean = pooled_mean/real(pooled_count, dp)
      if (pooled_count > 1) then
         sum_squares = 0.0_dp
         do candidate = 1, size(out%mean_loss)
            do observation = 1, size(out%loss, 1)
               if (out%valid(observation, candidate)) then
                  sum_squares = sum_squares + &
                     (out%loss(observation, candidate) - pooled_mean)**2
               end if
            end do
         end do
         out%selection_standard_error = sqrt(sum_squares/ &
            real(pooled_count - 1, dp))/sqrt(real(size(out%loss, 1), dp))
      end if

      minimum_loss = minval(out%mean_loss)
      tie_tolerance = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, minimum_loss)
      do candidate = 1, size(out%mean_loss)
         if (abs(out%mean_loss(candidate) - minimum_loss) <= tie_tolerance) &
            out%best_index = candidate
      end do
      threshold = minimum_loss + out%selection_standard_error
      do candidate = 1, size(out%mean_loss)
         if (out%mean_loss(candidate) < threshold) then
            out%one_se_index = candidate
            exit
         end if
      end do
      if (out%one_se_index == 0) out%one_se_index = out%best_index
      out%selected_index = out%best_index
      if (use_one_se) out%selected_index = out%one_se_index
      out%selected_lambda = out%lambda(out%selected_index)
      out%selected_alpha = out%alpha(out%selected_index)
   end subroutine summarize_validation

   pure subroutine summarize_separate_validation(out, use_one_se)
      !! Summarize and select each response's rolling penalty independently.
      type(bigvar_separate_validation_t), intent(inout) :: out !! Procedure result, updated in place.
      logical, intent(in) :: use_one_se !! Whether to use the one se.
      real(dp) :: mean_value, sum_squares, minimum_loss, threshold
      real(dp) :: tie_tolerance
      integer :: response, candidate, observation, valid_count

      if (count(out%valid) < 1) then
         out%info = 2
         return
      end if
      do response = 1, size(out%loss, 3)
         do candidate = 1, size(out%loss, 2)
            valid_count = count(out%valid(:, candidate, response))
            if (valid_count < 1) cycle
            mean_value = 0.0_dp
            do observation = 1, size(out%loss, 1)
               if (out%valid(observation, candidate, response)) then
                  mean_value = mean_value + &
                     out%loss(observation, candidate, response)
               end if
            end do
            mean_value = mean_value/real(valid_count, dp)
            out%mean_loss(candidate, response) = mean_value
            if (valid_count > 1) then
               sum_squares = 0.0_dp
               do observation = 1, size(out%loss, 1)
                  if (out%valid(observation, candidate, response)) then
                     sum_squares = sum_squares + &
                        (out%loss(observation, candidate, response) - &
                        mean_value)**2
                  end if
               end do
               out%standard_error(candidate, response) = sqrt(sum_squares/ &
                  real(valid_count - 1, dp))/sqrt(real(valid_count, dp))
            else
               out%standard_error(candidate, response) = 0.0_dp
            end if
         end do
         minimum_loss = minval(out%mean_loss(:, response))
         tie_tolerance = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, minimum_loss)
         do candidate = 1, size(out%loss, 2)
            if (abs(out%mean_loss(candidate, response) - minimum_loss) <= &
               tie_tolerance) out%best_index(response) = candidate
         end do
         threshold = minimum_loss + &
            out%standard_error(out%best_index(response), response)
         do candidate = 1, size(out%loss, 2)
            if (out%mean_loss(candidate, response) < threshold) then
               out%one_se_index(response) = candidate
               exit
            end if
         end do
         if (out%one_se_index(response) == 0) then
            out%one_se_index(response) = out%best_index(response)
         end if
         out%selected_index(response) = out%best_index(response)
         if (use_one_se) then
            out%selected_index(response) = out%one_se_index(response)
         end if
         out%selected_lambda(response) = out%lambda( &
            out%selected_index(response), response)
      end do
   end subroutine summarize_separate_validation

   pure function varx_prox(values, threshold, structure, variables, ar_order, &
      x_variables, exogenous_order, alpha) result(out)
      !! Apply BigVAR's joint endogenous and exogenous VARX proximal operator.
      real(dp), intent(in) :: values(:, :) !! Input values.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: x_variables !! X variables.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      real(dp) :: out(size(values, 1), size(values, 2))
      real(dp) :: group_norm, group_threshold
      integer :: ar_predictors, column

      ar_predictors = variables*ar_order
      out = values
      select case (structure)
      case (bigvar_structure_lag, bigvar_structure_own_other, &
         bigvar_structure_sparse_lag, bigvar_structure_sparse_own_other)
         out(:, :ar_predictors) = bigvar_group_prox(values(:, :ar_predictors), &
            threshold, structure, variables, ar_order, alpha)
         if (is_sparse_structure(structure)) then
            out(:, ar_predictors + 1:) = bigtime_soft_threshold( &
               values(:, ar_predictors + 1:), threshold*alpha)
            group_threshold = threshold*(1.0_dp - alpha)* &
               sqrt(real(variables, dp))
         else
            group_threshold = threshold*sqrt(real(variables, dp))
         end if
         do column = ar_predictors + 1, ar_predictors + &
            x_variables*exogenous_order
            group_norm = sqrt(sum(out(:, column)**2))
            call shrink_vector(out(:, column), group_threshold, group_norm)
         end do
      case (bigvar_structure_basic)
         out = bigtime_soft_threshold(values, threshold)
      case (bigvar_structure_basic_en)
         out = bigtime_soft_threshold(values, threshold*alpha)/ &
            (1.0_dp + threshold*(1.0_dp - alpha))
      case (bigvar_structure_efx)
         out = bigvar_efx_prox(values, threshold, variables, ar_order, &
            x_variables, exogenous_order)
      case default
         out = 0.0_dp
      end select
   end function varx_prox

   pure real(dp) function varx_penalty(phi, beta, structure, variables, &
      ar_order, x_variables, exogenous_order, alpha) result(value)
      !! Evaluate BigVAR's joint VARX structured penalty.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: beta(:, :) !! Regression or model coefficients.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: x_variables !! X variables.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp) :: x_weight
      integer :: column

      if (structure == bigvar_structure_basic) then
         value = sum(abs(phi)) + sum(abs(beta))
         return
      end if
      if (structure == bigvar_structure_basic_en) then
         value = alpha*(sum(abs(phi)) + sum(abs(beta))) + &
            0.5_dp*(1.0_dp - alpha)*(sum(phi**2) + sum(beta**2))
         return
      end if
      if (structure == bigvar_structure_efx) then
         value = efx_penalty(phi, beta, variables, ar_order, x_variables, &
            exogenous_order)
         return
      end if
      value = structured_penalty(phi, structure, variables, ar_order, alpha, &
         3.0_dp)
      x_weight = sqrt(real(variables, dp))
      do column = 1, size(beta, 2)
         value = value + (1.0_dp - alpha)*x_weight* &
            sqrt(sum(beta(:, column)**2))
      end do
      if (is_sparse_structure(structure)) value = value + &
         alpha*sum(abs(beta))
   end function varx_penalty

   pure real(dp) function efx_penalty(phi, beta, variables, ar_order, &
      x_variables, exogenous_order) result(value)
      !! Evaluate BigVAR's endogenous-first nested VARX penalty.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: beta(:, :) !! Regression or model coefficients.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: x_variables !! X variables.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      real(dp) :: joint_norm
      integer :: response, lag, ar_first, ar_last, x_first, x_last

      value = 0.0_dp
      do response = 1, variables
         do lag = 1, ar_order
            ar_first = (lag - 1)*variables + 1
            ar_last = lag*variables
            if (lag <= exogenous_order) then
               x_first = (lag - 1)*x_variables + 1
               x_last = lag*x_variables
               joint_norm = sqrt(sum(phi(response, ar_first:ar_last)**2) + &
                  sum(beta(response, x_first:x_last)**2))
               value = value + joint_norm + &
                  sqrt(sum(beta(response, x_first:x_last)**2))
            else
               value = value + &
                  sqrt(sum(phi(response, ar_first:ar_last)**2))
            end if
         end do
      end do
   end function efx_penalty

   pure integer function count_active_x_groups(beta) result(value)
      !! Count active response-vector groups in an exogenous coefficient block.
      real(dp), intent(in) :: beta(:, :) !! Regression or model coefficients.
      real(dp), parameter :: threshold = 100.0_dp*epsilon(1.0_dp)
      integer :: column

      value = 0
      do column = 1, size(beta, 2)
         if (maxval(abs(beta(:, column))) > threshold) value = value + 1
      end do
   end function count_active_x_groups

   pure integer function count_active_efx_groups(phi, beta, variables, &
      ar_order, x_variables, exogenous_order) result(value)
      !! Count active nested groups in an EFX coefficient estimate.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: beta(:, :) !! Regression or model coefficients.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: x_variables !! X variables.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      real(dp), parameter :: threshold = 100.0_dp*epsilon(1.0_dp)
      integer :: response, lag, ar_first, ar_last, x_first, x_last

      value = 0
      do response = 1, variables
         do lag = 1, ar_order
            ar_first = (lag - 1)*variables + 1
            ar_last = lag*variables
            if (lag <= exogenous_order) then
               x_first = (lag - 1)*x_variables + 1
               x_last = lag*x_variables
               if (maxval(abs(beta(response, x_first:x_last))) > &
                  threshold) value = value + 1
               if (max(maxval(abs(phi(response, ar_first:ar_last))), &
                  maxval(abs(beta(response, x_first:x_last)))) > &
                  threshold) value = value + 1
            else if (maxval(abs(phi(response, ar_first:ar_last))) > &
               threshold) then
               value = value + 1
            end if
         end do
      end do
   end function count_active_efx_groups

   pure real(dp) function varx_zero_model_bound(gradient, structure, variables, &
      ar_order, x_variables, exogenous_order, alpha) result(value)
      !! Find the smallest joint VARX penalty producing an all-zero update.
      real(dp), intent(in) :: gradient(:, :) !! Gradient.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: x_variables !! X variables.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      real(dp) :: lower, upper, midpoint
      integer :: iteration

      lower = 0.0_dp
      upper = max(1.0_dp, maxval(abs(gradient)))
      do while (.not. varx_prox_is_zero(gradient, upper, structure, variables, &
         ar_order, x_variables, exogenous_order, alpha))
         upper = 2.0_dp*upper
      end do
      do iteration = 1, 80
         midpoint = 0.5_dp*(lower + upper)
         if (varx_prox_is_zero(gradient, midpoint, structure, variables, &
            ar_order, x_variables, exogenous_order, alpha)) then
            upper = midpoint
         else
            lower = midpoint
         end if
      end do
      value = upper
   end function varx_zero_model_bound

   pure logical function varx_prox_is_zero(values, threshold, structure, &
      variables, ar_order, x_variables, exogenous_order, alpha) result(value)
      !! Test whether a joint VARX proximal update is identically zero.
      real(dp), intent(in) :: values(:, :) !! Input values.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: x_variables !! X variables.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      real(dp) :: transformed(size(values, 1), size(values, 2))

      transformed = varx_prox(values, threshold, structure, variables, &
         ar_order, x_variables, exogenous_order, alpha)
      value = maxval(abs(transformed)) <= &
         10.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(values)))
   end function varx_prox_is_zero

   pure subroutine shrink_block(block, threshold, block_norm)
      !! Shrink one nonoverlapping coefficient block toward zero.
      real(dp), intent(inout) :: block(:, :) !! Block, updated in place.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      real(dp), intent(in) :: block_norm !! Block norm.
      real(dp) :: shrinkage

      if (block_norm <= threshold) then
         block = 0.0_dp
      else
         shrinkage = 1.0_dp - threshold/block_norm
         block = shrinkage*block
      end if
   end subroutine shrink_block

   pure function relaxed_coefficients(response, design, coefficient) &
      result(out)
      !! Refit each response on its selected coefficient support.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: design(:, :) !! Design.
      real(dp), intent(in) :: coefficient(:, :) !! Coefficient.
      real(dp) :: out(size(coefficient, 1), size(coefficient, 2))
      real(dp), allocatable :: selected_design(:, :), normal_matrix(:, :)
      real(dp), allocatable :: inverse(:, :), right_hand_side(:), estimate(:)
      integer, allocatable :: support(:)
      integer :: equation, active, status

      out = coefficient
      do equation = 1, size(coefficient, 1)
         support = pack([(active, active=1, size(coefficient, 2))], &
            abs(coefficient(equation, :)) > 1.0e-8_dp)
         if (size(support) < 1) cycle
         allocate(selected_design(size(design, 1), size(support)))
         selected_design = design(:, support)
         allocate(normal_matrix(size(support), size(support)))
         normal_matrix = matmul(transpose(selected_design), selected_design)
         allocate(inverse(size(support), size(support)))
         call symmetric_pseudoinverse(normal_matrix, inverse, status)
         if (status == 0) then
            allocate(right_hand_side(size(support)), estimate(size(support)))
            right_hand_side = matmul(transpose(selected_design), &
               response(:, equation))
            estimate = matmul(inverse, right_hand_side)
            out(equation, :) = 0.0_dp
            out(equation, support) = estimate
         end if
         deallocate(selected_design, normal_matrix, inverse)
         if (allocated(right_hand_side)) deallocate(right_hand_side)
         if (allocated(estimate)) deallocate(estimate)
         if (allocated(support)) deallocate(support)
      end do
   end function relaxed_coefficients

   pure subroutine nonconvex_coordinate_descent(response, design, coefficient, &
      lambda, gamma, structure, tolerance, max_iterations, iterations, &
      converged)
      !! Fit MCP or SCAD coefficients by cyclic coordinate descent.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: design(:, :) !! Design.
      real(dp), intent(inout) :: coefficient(:, :) !! Coefficient, updated in place.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(out) :: iterations !! Number of algorithm iterations.
      logical, intent(out) :: converged !! Flag controlling converged.
      real(dp), allocatable :: residual(:), curvature(:)
      real(dp) :: score, old_value, new_value, shift, maximum_shift
      real(dp) :: response_scale
      integer :: equation, iteration, predictor, observations
      logical :: row_converged

      observations = size(response, 1)
      allocate(curvature(size(design, 2)))
      curvature = sum(design**2, dim=1)/real(observations, dp)
      allocate(residual(observations))
      iterations = 0
      converged = .true.
      do equation = 1, size(response, 2)
         residual = response(:, equation) - &
            matmul(design, coefficient(equation, :))
         response_scale = sqrt(sum(response(:, equation)**2)/ &
            real(observations, dp))
         row_converged = .false.
         do iteration = 1, max_iterations
            maximum_shift = 0.0_dp
            do predictor = 1, size(design, 2)
               old_value = coefficient(equation, predictor)
               if (curvature(predictor) <= tiny(1.0_dp)) then
                  new_value = 0.0_dp
               else
                  score = dot_product(residual, design(:, predictor))/ &
                     real(observations, dp) + old_value*curvature(predictor)
                  if (structure == bigvar_structure_mcp) then
                     new_value = bigvar_mcp_update(score, lambda, gamma, &
                        curvature(predictor))
                  else
                     new_value = bigvar_scad_update(score, lambda, gamma, &
                        curvature(predictor))
                  end if
               end if
               shift = new_value - old_value
               if (shift /= 0.0_dp) then
                  residual = residual - shift*design(:, predictor)
                  coefficient(equation, predictor) = new_value
                  maximum_shift = max(maximum_shift, abs(shift))
               end if
            end do
            if (maximum_shift <= tolerance*max(response_scale, 1.0_dp)) then
               row_converged = .true.
               exit
            end if
         end do
         iterations = max(iterations, min(iteration, max_iterations))
         converged = converged .and. row_converged
      end do
   end subroutine nonconvex_coordinate_descent

   pure subroutine nonconvex_coordinate_descent_separate(response, design, &
      coefficient, lambdas, gamma, structure, tolerance, max_iterations, &
      iterations, converged)
      !! Fit MCP or SCAD with a distinct coordinate penalty per response.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: design(:, :) !! Design.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(inout) :: coefficient(:, :) !! Coefficient, updated in place.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(out) :: iterations !! Number of algorithm iterations.
      logical, intent(out) :: converged !! Flag controlling converged.
      real(dp), allocatable :: residual(:), curvature(:)
      real(dp) :: score, old_value, new_value, shift, maximum_shift
      real(dp) :: response_scale
      integer :: equation, iteration, predictor, observations
      logical :: row_converged

      observations = size(response, 1)
      allocate(curvature(size(design, 2)))
      curvature = sum(design**2, dim=1)/real(observations, dp)
      allocate(residual(observations))
      iterations = 0
      converged = .true.
      do equation = 1, size(response, 2)
         residual = response(:, equation) - &
            matmul(design, coefficient(equation, :))
         response_scale = sqrt(sum(response(:, equation)**2)/ &
            real(observations, dp))
         row_converged = .false.
         do iteration = 1, max_iterations
            maximum_shift = 0.0_dp
            do predictor = 1, size(design, 2)
               old_value = coefficient(equation, predictor)
               if (curvature(predictor) <= tiny(1.0_dp)) then
                  new_value = 0.0_dp
               else
                  score = dot_product(residual, design(:, predictor))/ &
                     real(observations, dp) + old_value*curvature(predictor)
                  if (structure == bigvar_structure_mcp) then
                     new_value = bigvar_mcp_update(score, lambdas(equation), &
                        gamma, curvature(predictor))
                  else
                     new_value = bigvar_scad_update(score, lambdas(equation), &
                        gamma, curvature(predictor))
                  end if
               end if
               shift = new_value - old_value
               if (shift /= 0.0_dp) then
                  residual = residual - shift*design(:, predictor)
                  coefficient(equation, predictor) = new_value
                  maximum_shift = max(maximum_shift, abs(shift))
               end if
            end do
            if (maximum_shift <= tolerance*max(response_scale, 1.0_dp)) then
               row_converged = .true.
               exit
            end if
         end do
         iterations = max(iterations, min(iteration, max_iterations))
         converged = converged .and. row_converged
      end do
   end subroutine nonconvex_coordinate_descent_separate

   pure function separate_response_prox(values, thresholds, structure, &
      variables, lag_order, alpha) result(out)
      !! Apply row-specific thresholds to a response-separable VAR penalty.
      real(dp), intent(in) :: values(:, :) !! Input values.
      real(dp), intent(in) :: thresholds(:) !! Thresholds.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp) :: out(size(values, 1), size(values, 2))
      real(dp) :: work(size(values, 1), size(values, 2))
      integer :: response

      out = 0.0_dp
      do response = 1, variables
         work = 0.0_dp
         work(response, :) = values(response, :)
         work = bigvar_group_prox(work, thresholds(response), structure, &
            variables, lag_order, alpha)
         out(response, :) = work(response, :)
      end do
   end function separate_response_prox

   pure real(dp) function separate_penalized_objective(residuals, phi, &
      lambdas, structure, variables, lag_order, alpha, gamma) result(value)
      !! Evaluate an objective with a distinct penalty for each response.
      real(dp), intent(in) :: residuals(:, :) !! Model residuals.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: lambdas(:) !! Candidate penalty or shrinkage parameters.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp) :: work(size(phi, 1), size(phi, 2))
      integer :: response

      value = 0.5_dp*sum(residuals**2)
      do response = 1, variables
         work = 0.0_dp
         work(response, :) = phi(response, :)
         if (structure == bigvar_structure_mcp) then
            value = value + sum(mcp_penalty_value(phi(response, :), &
               lambdas(response), gamma)) + &
               0.5_dp*lambdas(response)*sum(phi(response, :)**2)
         else if (structure == bigvar_structure_scad) then
            value = value + sum(scad_penalty_value(phi(response, :), &
               lambdas(response), gamma)) + &
               0.5_dp*lambdas(response)*sum(phi(response, :)**2)
         else
            value = value + lambdas(response)*structured_penalty(work, &
               structure, variables, lag_order, alpha, gamma)
         end if
      end do
   end function separate_penalized_objective

   pure real(dp) function nonconvex_objective(residuals, phi, lambda, gamma, &
      structure) result(value)
      !! Evaluate a diagnostic MCP or SCAD penalized mean-square criterion.
      real(dp), intent(in) :: residuals(:, :) !! Model residuals.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      integer, intent(in) :: structure !! Model-structure specification.

      value = 0.5_dp*sum(residuals**2)/real(size(residuals, 1), dp) + &
         0.5_dp*lambda*sum(phi**2)
      if (structure == bigvar_structure_mcp) then
         value = value + sum(mcp_penalty_value(phi, lambda, gamma))
      else
         value = value + sum(scad_penalty_value(phi, lambda, gamma))
      end if
   end function nonconvex_objective

   pure elemental real(dp) function mcp_penalty_value(coefficient, lambda, &
      gamma) result(value)
      !! Evaluate the scalar minimax concave penalty.
      real(dp), intent(in) :: coefficient !! Coefficient.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      real(dp) :: magnitude

      magnitude = abs(coefficient)
      if (magnitude <= gamma*lambda) then
         value = lambda*magnitude - magnitude**2/(2.0_dp*gamma)
      else
         value = 0.5_dp*gamma*lambda**2
      end if
   end function mcp_penalty_value

   pure elemental real(dp) function scad_penalty_value(coefficient, lambda, &
      gamma) result(value)
      !! Evaluate the scalar smoothly clipped absolute-deviation penalty.
      real(dp), intent(in) :: coefficient !! Coefficient.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      real(dp) :: magnitude

      magnitude = abs(coefficient)
      if (magnitude <= lambda) then
         value = lambda*magnitude
      else if (magnitude <= gamma*lambda) then
         value = (-magnitude**2 + 2.0_dp*gamma*lambda*magnitude - &
            lambda**2)/(2.0_dp*(gamma - 1.0_dp))
      else
         value = 0.5_dp*(gamma + 1.0_dp)*lambda**2
      end if
   end function scad_penalty_value

   pure subroutine shrink_vector(vector, threshold, vector_norm)
      !! Shrink one vector group toward zero.
      real(dp), intent(inout) :: vector(:) !! Vector, updated in place.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      real(dp), intent(in) :: vector_norm !! Vector norm.
      real(dp) :: shrinkage

      if (vector_norm <= threshold) then
         vector = 0.0_dp
      else
         shrinkage = 1.0_dp - threshold/vector_norm
         vector = shrinkage*vector
      end if
   end subroutine shrink_vector

   pure subroutine shrink_own_other_suffix(vector, threshold, squared_norm, &
      first_column, response, last_column)
      !! Shrink an Own/Other suffix while preserving its current own lag.
      real(dp), intent(inout) :: vector(:) !! Vector, updated in place.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      real(dp), intent(in) :: squared_norm !! Squared norm.
      integer, intent(in) :: first_column !! First column.
      integer, intent(in) :: response !! Response observations.
      integer, intent(in) :: last_column !! Last column.
      real(dp) :: shrinkage
      integer :: own_column

      own_column = first_column + response - 1
      if (squared_norm <= threshold**2) then
         vector(first_column:own_column - 1) = 0.0_dp
         vector(own_column + 1:last_column) = 0.0_dp
      else
         shrinkage = 1.0_dp - threshold/sqrt(squared_norm)
         vector(first_column:own_column - 1) = shrinkage* &
            vector(first_column:own_column - 1)
         vector(own_column + 1:last_column) = shrinkage* &
            vector(own_column + 1:last_column)
      end if
   end subroutine shrink_own_other_suffix

   pure subroutine shrink_lag_suffix(vector, threshold, squared_norm, &
      predictor, first_lag, last_lag, variables)
      !! Shrink one predictor's coefficient suffix across successive lags.
      real(dp), intent(inout) :: vector(:) !! Vector, updated in place.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      real(dp), intent(in) :: squared_norm !! Squared norm.
      integer, intent(in) :: predictor !! Predictor.
      integer, intent(in) :: first_lag !! First lag.
      integer, intent(in) :: last_lag !! Last lag.
      integer, intent(in) :: variables !! Number or indices of variables.
      real(dp) :: shrinkage
      integer :: lag, column

      if (squared_norm <= threshold**2) then
         shrinkage = 0.0_dp
      else
         shrinkage = 1.0_dp - threshold/sqrt(squared_norm)
      end if
      do lag = first_lag, last_lag
         column = (lag - 1)*variables + predictor
         vector(column) = shrinkage*vector(column)
      end do
   end subroutine shrink_lag_suffix

   pure real(dp) function structured_penalty(phi, structure, variables, &
      lag_order, alpha, gamma) &
      result(value)
      !! Evaluate a BigVAR group or sparse-group penalty.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: gamma !! Model coefficient or scale parameter.
      real(dp) :: diagonal_norm, other_norm
      real(dp) :: group_weight
      integer :: group_structure, lag, suffix_lag, first_column
      integer :: predictor, response

      value = alpha*sum(abs(phi))
      group_structure = base_group_structure(structure)
      select case (group_structure)
      case (bigvar_structure_lag)
         if (structure == bigvar_structure_sparse_lag) then
            group_weight = real(variables, dp)
         else
            group_weight = sqrt(real(variables, dp))
         end if
         do lag = 1, lag_order
            first_column = (lag - 1)*variables + 1
            value = value + (1.0_dp - alpha)*group_weight*sqrt( &
               sum(phi(:, first_column:first_column + variables - 1)**2))
         end do
      case (bigvar_structure_own_other)
         do lag = 1, lag_order
            first_column = (lag - 1)*variables
            diagonal_norm = 0.0_dp
            other_norm = 0.0_dp
            do response = 1, variables
               diagonal_norm = diagonal_norm + &
                  phi(response, first_column + response)**2
               other_norm = other_norm + &
                  sum(phi(response, first_column + 1: &
                  first_column + variables)**2) - &
                  phi(response, first_column + response)**2
            end do
            value = value + (1.0_dp - alpha)* &
               sqrt(real(variables, dp))*sqrt(diagonal_norm)
            if (variables > 1) value = value + (1.0_dp - alpha)* &
               sqrt(real(variables*(variables - 1), dp))* &
               sqrt(max(other_norm, 0.0_dp))
         end do
      case (bigvar_structure_hlag_component)
         do response = 1, variables
            do lag = 1, lag_order
               first_column = (lag - 1)*variables + 1
               value = value + sqrt(real(variables, dp))*sqrt(sum( &
                  phi(response, first_column:variables*lag_order)**2))
            end do
         end do
      case (bigvar_structure_hlag_own_other)
         do response = 1, variables
            do lag = 1, lag_order
               first_column = (lag - 1)*variables + 1
               value = value + sqrt(sum(phi(response, &
                  first_column:variables*lag_order)**2))
               other_norm = sum(phi(response, first_column: &
                  variables*lag_order)**2) - &
                  phi(response, first_column + response - 1)**2
               if (variables > 1) value = value + &
                  sqrt(real(variables - 1, dp))*sqrt(max(other_norm, 0.0_dp))
            end do
         end do
      case (bigvar_structure_hlag_element)
         do response = 1, variables
            do predictor = 1, variables
               do lag = 1, lag_order
                  diagonal_norm = 0.0_dp
                  do suffix_lag = lag, lag_order
                     diagonal_norm = diagonal_norm + phi(response, &
                        (suffix_lag - 1)*variables + predictor)**2
                  end do
                  value = value + sqrt(diagonal_norm)
               end do
            end do
         end do
      case (bigvar_structure_basic)
         value = sum(abs(phi))
      case (bigvar_structure_basic_en)
         value = alpha*sum(abs(phi)) + &
            0.5_dp*(1.0_dp - alpha)*sum(phi**2)
      case (bigvar_structure_tapered)
         value = 0.0_dp
         do lag = 1, lag_order
            first_column = (lag - 1)*variables + 1
            value = value + real(lag, dp)**alpha*sum(abs(phi(:, &
               first_column:first_column + variables - 1)))
         end do
      case (bigvar_structure_mcp)
         value = sum(mcp_penalty_value(phi, 1.0_dp, gamma))
      case (bigvar_structure_scad)
         value = sum(scad_penalty_value(phi, 1.0_dp, gamma))
      end select
   end function structured_penalty

   pure integer function count_active_groups(phi, structure, variables, &
      lag_order) result(value)
      !! Count nonzero groups in a structured coefficient matrix.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp), parameter :: threshold = 100.0_dp*epsilon(1.0_dp)
      real(dp) :: other_norm
      integer :: lag, first_column, predictor, response, suffix_lag

      value = 0
      select case (base_group_structure(structure))
      case (bigvar_structure_lag)
         do lag = 1, lag_order
            first_column = (lag - 1)*variables + 1
            if (maxval(abs(phi(:, first_column: &
               first_column + variables - 1))) > threshold) value = value + 1
         end do
      case (bigvar_structure_own_other)
         do lag = 1, lag_order
            first_column = (lag - 1)*variables
            if (any([(abs(phi(response, first_column + response)) > &
               threshold, response=1, variables)])) value = value + 1
            other_norm = 0.0_dp
            do response = 1, variables
               other_norm = other_norm + sum(abs(phi(response, &
                  first_column + 1:first_column + variables)))
               other_norm = other_norm - &
                  abs(phi(response, first_column + response))
            end do
            if (other_norm > threshold) value = value + 1
         end do
      case (bigvar_structure_hlag_component)
         do response = 1, variables
            do lag = 1, lag_order
               first_column = (lag - 1)*variables + 1
               if (maxval(abs(phi(response, first_column: &
                  variables*lag_order))) > threshold) value = value + 1
            end do
         end do
      case (bigvar_structure_hlag_own_other)
         do response = 1, variables
            do lag = 1, lag_order
               first_column = (lag - 1)*variables + 1
               if (maxval(abs(phi(response, first_column: &
                  variables*lag_order))) > threshold) value = value + 1
               other_norm = sum(abs(phi(response, first_column: &
                  variables*lag_order))) - &
                  abs(phi(response, first_column + response - 1))
               if (other_norm > threshold) value = value + 1
            end do
         end do
      case (bigvar_structure_hlag_element)
         do response = 1, variables
            do predictor = 1, variables
               do lag = 1, lag_order
                  if (any([(abs(phi(response, &
                     (suffix_lag - 1)*variables + predictor)) > threshold, &
                     suffix_lag=lag, lag_order)])) value = value + 1
               end do
            end do
         end do
      case (bigvar_structure_basic, bigvar_structure_basic_en, &
         bigvar_structure_tapered, &
         bigvar_structure_mcp, bigvar_structure_scad)
         value = count(abs(phi) > threshold)
      end select
   end function count_active_groups

   pure real(dp) function zero_model_bound(gradient, structure, variables, &
      lag_order, alpha) result(value)
      !! Find the smallest penalty that satisfies every zero-group condition.
      real(dp), intent(in) :: gradient(:, :) !! Gradient.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp) :: lower, upper, midpoint
      integer :: iteration

      if (maxval(abs(gradient)) <= tiny(1.0_dp)) then
         value = 0.0_dp
         return
      end if
      lower = 0.0_dp
      upper = maxval(abs(gradient))
      do while (.not. prox_is_zero(gradient, upper, structure, variables, &
         lag_order, alpha))
         upper = 2.0_dp*upper
      end do
      do iteration = 1, 80
         midpoint = 0.5_dp*(lower + upper)
         if (prox_is_zero(gradient, midpoint, structure, variables, &
            lag_order, alpha)) then
            upper = midpoint
         else
            lower = midpoint
         end if
      end do
      value = upper
   end function zero_model_bound

   pure logical function prox_is_zero(values, threshold, structure, variables, &
      lag_order, alpha) result(value)
      !! Test whether a structured proximal update is identically zero.
      real(dp), intent(in) :: values(:, :) !! Input values.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      integer, intent(in) :: structure !! Model-structure specification.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp) :: transformed(size(values, 1), size(values, 2))

      transformed = bigvar_group_prox(values, threshold, structure, variables, &
         lag_order, alpha)
      value = maxval(abs(transformed)) <= &
         10.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(values)))
   end function prox_is_zero

   pure elemental logical function is_sparse_structure(structure) result(value)
      !! Identify structures that combine elementwise and group penalties.
      integer, intent(in) :: structure !! Model-structure specification.

      value = structure == bigvar_structure_sparse_lag .or. &
         structure == bigvar_structure_sparse_own_other
   end function is_sparse_structure

   pure elemental logical function is_supported_structure(structure) &
      result(value)
      !! Identify the BigVAR structures implemented by this module.
      integer, intent(in) :: structure !! Model-structure specification.

      value = (structure >= bigvar_structure_lag .and. &
         structure <= bigvar_structure_scad) .or. &
         structure == bigvar_structure_hlag_element .or. &
         structure == bigvar_structure_basic
   end function is_supported_structure

   pure elemental logical function supports_separate_lambdas(structure) &
      result(value)
      !! Identify BigVAR penalties that separate across response equations.
      integer, intent(in) :: structure !! Model-structure specification.

      value = structure == bigvar_structure_basic .or. &
         structure == bigvar_structure_basic_en .or. &
         structure == bigvar_structure_hlag_component .or. &
         structure == bigvar_structure_hlag_own_other .or. &
         structure == bigvar_structure_hlag_element .or. &
         structure == bigvar_structure_mcp .or. &
         structure == bigvar_structure_scad
   end function supports_separate_lambdas

   pure elemental logical function is_nonconvex_structure(structure) &
      result(value)
      !! Identify structures estimated by non-convex coordinate descent.
      integer, intent(in) :: structure !! Model-structure specification.

      value = structure == bigvar_structure_mcp .or. &
         structure == bigvar_structure_scad
   end function is_nonconvex_structure

   pure elemental logical function is_varx_structure(structure) result(value)
      !! Identify BigVAR structures that support exogenous predictors.
      integer, intent(in) :: structure !! Model-structure specification.

      value = structure == bigvar_structure_lag .or. &
         structure == bigvar_structure_own_other .or. &
         structure == bigvar_structure_sparse_lag .or. &
         structure == bigvar_structure_sparse_own_other .or. &
         structure == bigvar_structure_basic .or. &
         structure == bigvar_structure_basic_en .or. &
         structure == bigvar_structure_mcp .or. &
         structure == bigvar_structure_scad .or. &
         structure == bigvar_structure_efx
   end function is_varx_structure

   pure elemental logical function supports_transfer_function(structure) &
      result(value)
      !! Identify penalties supported by BigVAR transfer-function models.
      integer, intent(in) :: structure !! Model-structure specification.

      value = structure == bigvar_structure_basic .or. &
         structure == bigvar_structure_basic_en .or. &
         structure == bigvar_structure_mcp .or. &
         structure == bigvar_structure_scad
   end function supports_transfer_function

   pure elemental logical function uses_alpha(structure) result(value)
      !! Identify structures with an additional mixing or taper parameter.
      integer, intent(in) :: structure !! Model-structure specification.

      value = is_sparse_structure(structure) .or. &
         structure == bigvar_structure_basic_en .or. &
         structure == bigvar_structure_tapered
   end function uses_alpha

   pure elemental integer function base_group_structure(structure) &
      result(value)
      !! Map a sparse-group structure to its underlying group layout.
      integer, intent(in) :: structure !! Model-structure specification.

      select case (structure)
      case (bigvar_structure_sparse_lag)
         value = bigvar_structure_lag
      case (bigvar_structure_sparse_own_other)
         value = bigvar_structure_own_other
      case default
         value = structure
      end select
   end function base_group_structure

   pure function univariate_ar_scales(series, lag_order) result(scales)
      !! Estimate BGR prior scales from separate univariate AR regressions.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp), allocatable :: scales(:)
      real(dp), allocatable :: response(:, :), lag_design(:, :), design(:, :)
      real(dp), allocatable :: normal_matrix(:, :), inverse(:, :)
      real(dp), allocatable :: right(:), coefficient(:), residual(:)
      integer :: variable, status

      allocate(scales(size(series, 2)))
      do variable = 1, size(series, 2)
         call build_var_data(series(:, variable:variable), lag_order, &
            response, lag_design)
         allocate(design(size(response, 1), lag_order + 1))
         design(:, :lag_order) = lag_design
         design(:, lag_order + 1) = 1.0_dp
         allocate(normal_matrix(lag_order + 1, lag_order + 1))
         normal_matrix = matmul(transpose(design), design)
         call invert_matrix(normal_matrix, inverse, status)
         if (status /= 0) then
            if (allocated(inverse)) deallocate(inverse)
            allocate(inverse(lag_order + 1, lag_order + 1))
            call symmetric_pseudoinverse(normal_matrix, inverse, status)
         end if
         if (status == 0) then
            allocate(right(lag_order + 1), coefficient(lag_order + 1))
            right = matmul(transpose(design), response(:, 1))
            coefficient = matmul(inverse, right)
            allocate(residual(size(response, 1)))
            residual = response(:, 1) - matmul(design, coefficient)
            scales(variable) = sqrt(sum(residual**2)/ &
               real(size(residual), dp))
         else
            scales(variable) = sqrt(sum((series(:, variable) - &
               sum(series(:, variable))/real(size(series, 1), dp))**2)/ &
               real(size(series, 1), dp))
         end if
         scales(variable) = max(scales(variable), sqrt(epsilon(1.0_dp)))
         if (allocated(response)) deallocate(response)
         if (allocated(lag_design)) deallocate(lag_design)
         if (allocated(design)) deallocate(design)
         if (allocated(normal_matrix)) deallocate(normal_matrix)
         if (allocated(inverse)) deallocate(inverse)
         if (allocated(right)) deallocate(right)
         if (allocated(coefficient)) deallocate(coefficient)
         if (allocated(residual)) deallocate(residual)
      end do
   end function univariate_ar_scales

   pure function matrix_without_row(values, omitted) result(out)
      !! Return a matrix with one row removed while preserving row order.
      real(dp), intent(in) :: values(:, :) !! Input values.
      integer, intent(in) :: omitted !! Omitted.
      real(dp), allocatable :: out(:, :)

      allocate(out(size(values, 1) - 1, size(values, 2)))
      if (omitted > 1) out(:omitted - 1, :) = values(:omitted - 1, :)
      if (omitted < size(values, 1)) then
         out(omitted:, :) = values(omitted + 1:, :)
      end if
   end function matrix_without_row

   pure subroutine set_minnesota_target(matrix, target)
      !! Place per-series targets on the endogenous first-lag diagonal.
      real(dp), intent(inout) :: matrix(:, :) !! Input matrix, updated in place.
      real(dp), intent(in) :: target(:) !! Target.
      integer :: variable

      do variable = 1, size(target)
         matrix(variable, variable) = target(variable)
      end do
   end subroutine set_minnesota_target

   pure function predict_var_coefficients(phi, intercept, lag_order, history, &
      horizon, recursive) result(prediction)
      !! Forecast VAR coefficients recursively or at one direct horizon.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: intercept(:) !! Model intercept.
      real(dp), intent(in) :: history(:, :) !! History.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      logical, intent(in) :: recursive !! Flag controlling recursive.
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: working(:, :)
      integer :: observations, variables, step, lag, target

      observations = size(history, 1)
      variables = size(history, 2)
      allocate(prediction(variables))
      allocate(working(observations + horizon, variables))
      working(:observations, :) = history
      if (.not. recursive) then
         prediction = intercept
         do lag = 1, lag_order
            prediction = prediction + matmul(phi(:, &
               (lag - 1)*variables + 1:lag*variables), &
               history(observations + 1 - lag, :))
         end do
         return
      end if
      do step = 1, horizon
         target = observations + step
         working(target, :) = intercept
         do lag = 1, lag_order
            working(target, :) = working(target, :) + matmul(phi(:, &
               (lag - 1)*variables + 1:lag*variables), &
               working(target - lag, :))
         end do
      end do
      prediction = working(observations + horizon, :)
   end function predict_var_coefficients

   pure function predict_varx_coefficients(phi, beta, intercept, ar_order, &
      exogenous_order, contemporaneous, history, exogenous, horizon) &
      result(prediction)
      !! Forecast VARX coefficients using the supplied exogenous path.
      real(dp), intent(in) :: phi(:, :) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: beta(:, :) !! Regression or model coefficients.
      real(dp), intent(in) :: intercept(:) !! Model intercept.
      real(dp), intent(in) :: history(:, :) !! History.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      logical, intent(in) :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: working(:, :)
      integer :: observations, variables, x_variables
      integer :: step, lag, block, target

      observations = size(history, 1)
      variables = size(history, 2)
      x_variables = size(exogenous, 2)
      allocate(working(observations + horizon, variables))
      working(:observations, :) = history
      do step = 1, horizon
         target = observations + step
         working(target, :) = intercept
         do lag = 1, ar_order
            working(target, :) = working(target, :) + matmul(phi(:, &
               (lag - 1)*variables + 1:lag*variables), &
               working(target - lag, :))
         end do
         do lag = 0, exogenous_order
            if (.not. contemporaneous .and. lag == 0) cycle
            block = lag
            if (contemporaneous) block = block + 1
            working(target, :) = working(target, :) + matmul(beta(:, &
               (block - 1)*x_variables + 1:block*x_variables), &
               exogenous(target - lag, :))
         end do
      end do
      allocate(prediction(variables))
      prediction = working(observations + horizon, :)
   end function predict_varx_coefficients

   pure subroutine build_var_direct_data(series, lag_order, horizon, response, &
      design)
      !! Align lagged predictors with responses at a specified direct horizon.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), allocatable, intent(out) :: response(:, :) !! Response observations.
      real(dp), allocatable, intent(out) :: design(:, :) !! Design.
      integer :: observations, variables, cases, lag

      observations = size(series, 1)
      variables = size(series, 2)
      cases = observations - lag_order - horizon + 1
      allocate(response(cases, variables))
      allocate(design(cases, variables*lag_order))
      response = series(lag_order + horizon:observations, :)
      do lag = 1, lag_order
         design(:, (lag - 1)*variables + 1:lag*variables) = series( &
            lag_order + 1 - lag:observations - horizon + 1 - lag, :)
      end do
   end subroutine build_var_direct_data

   pure subroutine build_ls_varx_data(endogenous, exogenous, ar_order, &
      exogenous_order, horizon, contemporaneous, response, predictors)
      !! Align least-squares VARX predictors with direct-horizon responses.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      logical, intent(in) :: contemporaneous !! Flag controlling contemporaneous.
      real(dp), allocatable, intent(out) :: response(:, :) !! Response observations.
      real(dp), allocatable, intent(out) :: predictors(:, :) !! Predictor matrix.
      integer :: observations, variables, x_variables, leading, cases
      integer :: lag, block, x_blocks, ar_predictors

      observations = size(endogenous, 1)
      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      leading = max(ar_order, exogenous_order)
      cases = observations - leading - horizon + 1
      x_blocks = exogenous_order
      if (contemporaneous) x_blocks = x_blocks + 1
      ar_predictors = variables*ar_order
      allocate(response(max(cases, 0), variables))
      allocate(predictors(max(cases, 0), &
         ar_predictors + x_variables*x_blocks))
      if (cases < 1) return
      response = endogenous(leading + horizon:observations, :)
      do lag = 1, ar_order
         predictors(:, (lag - 1)*variables + 1:lag*variables) = endogenous( &
            leading + 1 - lag:observations - horizon + 1 - lag, :)
      end do
      if (contemporaneous) then
         predictors(:, ar_predictors + 1:ar_predictors + x_variables) = &
            exogenous(leading + horizon:observations, :)
      end if
      do lag = 1, exogenous_order
         block = lag
         if (contemporaneous) block = block + 1
         if (contemporaneous) then
            predictors(:, ar_predictors + (block - 1)*x_variables + 1: &
               ar_predictors + block*x_variables) = exogenous( &
               leading + horizon - lag:observations - lag, :)
         else
            predictors(:, ar_predictors + (block - 1)*x_variables + 1: &
               ar_predictors + block*x_variables) = exogenous( &
               leading + 1 - lag:observations - horizon + 1 - lag, :)
         end if
      end do
   end subroutine build_ls_varx_data

end module bigvar_mod
