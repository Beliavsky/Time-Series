! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R FCVAR package.
! Fractionally cointegrated VAR algorithms translated from the GPL-3 FCVAR package.
module fcvar_mod
   use kind_mod, only: dp
   use fourier_mod, only: fft_transform
   use linalg_mod, only: cholesky_lower, identity_matrix, &
      invert_matrix, symmetric_eigen, symmetrize, general_eigenvalues, kronecker_product
   use linalg_mod, only: symmetric_pseudoinverse
   use optimization_mod, only: optimization_result_t, &
      bfgs_minimize_fd, nelder_mead_minimize
   use random_mod, only: random_standard_normal_matrix, random_uniform
   use time_series_diagnostics_mod, only: multivariate_white_noise_test_t, &
      multivariate_white_noise_test
   use special_functions_mod, only: regularized_gamma_q
   use stats_mod, only: sort
   use fracdist_mod, only: fracdist_probability_t, fracdist_p_value
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
      ieee_value, ieee_positive_inf
   implicit none
   private

   type, public :: fcvar_options_t
      ! Deterministic terms and the number of conditioned initial observations.
      integer :: initial_values = 0
      logical :: level_parameter = .false.
      logical :: prefer_high_b_local_max = .false.
      logical :: restricted_constant = .false.
      logical :: unrestricted_constant = .false.
   end type fcvar_options_t

   type, public :: fcvar_transform_t
      ! Fractionally transformed response, levels, lag, and deterministic arrays.
      real(dp), allocatable :: z0(:, :), z1(:, :), z2(:, :), z3(:, :)
      integer :: info = 0
   end type fcvar_transform_t

   type, public :: fcvar_fit_t
      ! Concentrated FCVAR estimates at fixed fractional orders.
      real(dp), allocatable :: alpha(:, :), beta(:, :), rho(:)
      real(dp), allocatable :: pi(:, :), covariance(:, :), gamma(:, :)
      real(dp), allocatable :: level(:), unrestricted_constant(:), residuals(:, :)
      real(dp) :: d = 0.0_dp
      real(dp) :: b = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: lag_order = 0
      integer :: rank = 0
      integer :: observations = 0
      integer :: info = 0
   end type fcvar_fit_t

   type, public :: fcvar_grid_t
      ! Bounded fractional-order likelihood grid and its maximizing point.
      real(dp), allocatable :: d(:), b(:), log_likelihood(:, :)
      real(dp), allocatable :: level(:, :, :)
      real(dp), allocatable :: best_level(:)
      real(dp), allocatable :: local_d(:), local_b(:), local_log_likelihood(:)
      integer, allocatable :: local_d_index(:), local_b_index(:)
      real(dp) :: best_d = 0.0_dp
      real(dp) :: best_b = 0.0_dp
      real(dp) :: best_log_likelihood = -huge(1.0_dp)
      real(dp) :: global_d = 0.0_dp
      real(dp) :: global_b = 0.0_dp
      real(dp) :: global_log_likelihood = -huge(1.0_dp)
      integer :: info = 0
      integer :: best_d_index = 0
      integer :: best_b_index = 0
      integer :: global_d_index = 0
      integer :: global_b_index = 0
      logical :: used_local_maximum = .false.
   end type fcvar_grid_t

   type, public :: fcvar_estimation_t
      ! Grid-initialized fractional-order optimization and concentrated fit.
      type(fcvar_fit_t) :: fit
      type(fcvar_grid_t) :: grid
      real(dp), allocatable :: optimizer_coordinates(:)
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: converged = .false.
   end type fcvar_estimation_t

   type, public :: fcvar_path_t
      ! Recursive FCVAR values and the innovations applied to them.
      real(dp), allocatable :: series(:, :), innovations(:, :)
      integer :: info = 0
   end type fcvar_path_t

   type, public :: fcvar_lr_test_t
      ! Nested-model likelihood-ratio statistic, degrees of freedom, and tail.
      real(dp) :: unrestricted_log_likelihood = 0.0_dp
      real(dp) :: restricted_log_likelihood = 0.0_dp
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: degrees_of_freedom = 0
      integer :: info = 0
   end type fcvar_lr_test_t

   type, public :: fcvar_rank_tests_t
      ! Rankwise estimates, full-rank LR statistics, and information criteria.
      real(dp), allocatable :: d(:), b(:), log_likelihood(:), statistic(:)
      real(dp), allocatable :: p_value(:), aic(:), bic(:)
      integer, allocatable :: free_parameters(:)
      logical, allocatable :: p_value_available(:)
      integer :: aic_rank = 0
      integer :: bic_rank = 0
      integer :: lag_order = 0
      integer :: observations = 0
      integer :: info = 0
   end type fcvar_rank_tests_t

   type, public :: fcvar_lag_selection_t
      ! Candidate lag estimates, tests, criteria, and residual diagnostics.
      real(dp), allocatable :: d(:), b(:), log_likelihood(:)
      real(dp), allocatable :: statistic(:), p_value(:), aic(:), bic(:)
      real(dp), allocatable :: q_p_value(:, :), lm_p_value(:, :)
      real(dp), allocatable :: multivariate_q_p_value(:)
      integer, allocatable :: free_parameters(:)
      integer :: aic_lag = 0
      integer :: bic_lag = 0
      integer :: maximum_lag = 0
      integer :: rank = 0
      integer :: diagnostic_lag = 0
      integer :: observations = 0
      integer :: info = 0
   end type fcvar_lag_selection_t

   type, public :: fcvar_bootstrap_rank_t
      ! Wild-bootstrap LR distribution for a nested FCVAR rank comparison.
      type(fcvar_fit_t) :: null_fit, alternative_fit
      real(dp), allocatable :: statistic(:)
      real(dp) :: observed_statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: null_rank = 0
      integer :: alternative_rank = 0
      integer :: replications = 0
      integer :: exceedances = 0
      integer :: info = 0
   end type fcvar_bootstrap_rank_t

   type, public :: fcvar_standard_errors_t
      ! Numerical observed-information covariance and coefficient errors.
      real(dp), allocatable :: parameters(:), hessian(:, :), covariance(:, :)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: alpha(:, :), beta(:, :), rho(:), gamma(:, :)
      real(dp), allocatable :: level(:), unrestricted_constant(:)
      real(dp) :: d = 0.0_dp
      real(dp) :: b = 0.0_dp
      real(dp) :: increment = 0.0_dp
      integer :: info = 0
      logical :: equal_orders = .true.
   end type fcvar_standard_errors_t

   type, public :: fcvar_roots_t
      ! FCVAR companion eigenvalues and inverse characteristic roots.
      complex(dp), allocatable :: eigenvalues(:), roots(:)
      real(dp) :: b = 0.0_dp
      real(dp) :: minimum_root_modulus = 0.0_dp
      real(dp) :: maximum_inverse_root_modulus = 0.0_dp
      integer :: info = 0
      logical :: outside_unit_circle = .false.
   end type fcvar_roots_t

   type, public :: fcvar_restrictions_t
      ! Linear restrictions on packed, order, alpha, and beta-star parameters.
      real(dp), allocatable :: parameter_matrix(:, :), parameter_value(:)
      real(dp), allocatable :: order_matrix(:, :), order_value(:)
      real(dp), allocatable :: order_inequality_matrix(:, :)
      real(dp), allocatable :: order_inequality_value(:)
      real(dp), allocatable :: alpha_matrix(:, :), alpha_value(:)
      real(dp), allocatable :: beta_matrix(:, :), beta_value(:)
      real(dp), allocatable :: level_matrix(:, :), level_value(:)
   end type fcvar_restrictions_t

   type, public :: fcvar_restricted_estimation_t
      ! Unrestricted and linearly restricted fits with their LR test.
      type(fcvar_fit_t) :: unrestricted_fit, restricted_fit
      type(fcvar_lr_test_t) :: test
      real(dp), allocatable :: constraint_matrix(:, :), constraint_value(:)
      integer :: restriction_rank = 0
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: converged = .false.
   end type fcvar_restricted_estimation_t

   type, public :: fcvar_order_estimation_t
      ! FCVAR fit under affine equality and inequality restrictions on d and b.
      type(fcvar_fit_t) :: fit
      real(dp) :: starting_d = 0.0_dp
      real(dp) :: starting_b = 0.0_dp
      integer :: equality_rank = 0
      integer :: free_dimensions = 0
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: converged = .false.
   end type fcvar_order_estimation_t

   type, public :: fcvar_switching_estimation_t
      ! Restricted fixed-order fit from alternating alpha-beta GLS updates.
      type(fcvar_fit_t) :: fit
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type fcvar_switching_estimation_t

   type, public :: fcvar_bootstrap_hypothesis_t
      ! Wild-bootstrap distribution for a linearly restricted FCVAR LR test.
      type(fcvar_fit_t) :: unrestricted_fit, restricted_fit
      real(dp), allocatable :: statistic(:)
      real(dp) :: observed_statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: restriction_rank = 0
      integer :: replications = 0
      integer :: exceedances = 0
      integer :: info = 0
   end type fcvar_bootstrap_hypothesis_t

   public :: fcvar_fractional_difference, fcvar_fractional_lags
   public :: fcvar_transform, fcvar_estimate_fixed, fcvar_log_likelihood_fixed
   public :: fcvar_likelihood_grid, fcvar_estimate
   public :: fcvar_find_local_maxima
   public :: fcvar_estimate_equal_orders
   public :: fcvar_estimate_constrained_orders
   public :: fcvar_forecast, fcvar_simulate_from_innovations, fcvar_simulate
   public :: fcvar_bootstrap_from_signs, fcvar_bootstrap
   public :: fcvar_likelihood_ratio, fcvar_free_parameter_count
   public :: fcvar_rank_tests
   public :: fcvar_lag_select
   public :: fcvar_bootstrap_rank_from_signs, fcvar_bootstrap_rank
   public :: fcvar_pack_parameters, fcvar_unpack_parameters
   public :: fcvar_full_log_likelihood, fcvar_standard_errors
   public :: fcvar_characteristic_roots
   public :: fcvar_restricted_estimate
   public :: fcvar_estimate_restricted_fixed
   public :: fcvar_bootstrap_hypothesis_from_signs
   public :: fcvar_bootstrap_hypothesis

contains

   pure function fcvar_fractional_difference(x, d) result(differenced)
      !! Apply FCVAR's level-preserving fractional difference by FFT convolution.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), allocatable :: differenced(:, :)
      complex(dp), allocatable :: padded(:), weights(:), transformed(:)
      integer :: n, variables, nfft, lag, variable

      n = size(x, 1)
      variables = size(x, 2)
      if (n < 1 .or. variables < 1 .or. .not. ieee_is_finite(d) .or. &
         .not. all(ieee_is_finite(x))) then
         allocate(differenced(0, 0))
         return
      end if
      nfft = 1
      do while (nfft < 2*n - 1)
         nfft = 2*nfft
      end do
      allocate(padded(nfft), weights(nfft), differenced(n, variables))
      weights = cmplx(0.0_dp, 0.0_dp, dp)
      weights(1) = cmplx(1.0_dp, 0.0_dp, dp)
      do lag = 1, n - 1
         weights(lag + 1) = weights(lag)*real(lag - d - 1.0_dp, dp)/ &
            real(lag, dp)
      end do
      weights = fft_transform(weights)
      do variable = 1, variables
         padded = cmplx(0.0_dp, 0.0_dp, dp)
         padded(:n) = cmplx(x(:, variable), 0.0_dp, dp)
         transformed = fft_transform(fft_transform(padded)*weights, inverse=.true.)
         differenced(:, variable) = real(transformed(:n), dp)
      end do
   end function fcvar_fractional_difference

   pure function fcvar_fractional_lags(x, b, order) result(lagged)
      !! Form powers of L_b = 1 - (1-L)^b through the requested order.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: b !! B.
      integer, intent(in) :: order !! Model or polynomial order.
      real(dp), allocatable :: lagged(:, :)
      real(dp), allocatable :: previous(:, :)
      integer :: n, variables, lag

      n = size(x, 1)
      variables = size(x, 2)
      if (n < 1 .or. variables < 1 .or. order < 0 .or. &
         .not. ieee_is_finite(b) .or. .not. all(ieee_is_finite(x))) then
         allocate(lagged(0, 0))
         return
      end if
      allocate(lagged(n, variables*order))
      if (order == 0) return
      previous = x - fcvar_fractional_difference(x, b)
      lagged(:, :variables) = previous
      do lag = 2, order
         previous = previous - fcvar_fractional_difference(previous, b)
         lagged(:, 1 + (lag - 1)*variables:lag*variables) = previous
      end do
   end function fcvar_fractional_lags

   pure function fcvar_transform(x, lag_order, d, b, options) result(out)
      !! Construct FCVAR's Z0, Z1, Z2, and Z3 regression arrays.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: b !! B.
      integer, intent(in) :: lag_order !! Model lag order.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      type(fcvar_transform_t) :: out
      type(fcvar_options_t) :: selected
      real(dp), allocatable :: levels(:, :), work(:, :)
      integer :: n, variables, initial, rows

      selected = fcvar_options_t()
      if (present(options)) selected = options
      n = size(x, 1)
      variables = size(x, 2)
      initial = selected%initial_values
      rows = n - initial
      if (n < 1 .or. variables < 1 .or. lag_order < 0 .or. initial < 0 .or. &
         rows < 1 .or. .not. ieee_is_finite(d) .or. &
         .not. ieee_is_finite(b) .or. .not. all(ieee_is_finite(x))) then
         out%info = 1
         return
      end if
      out%z0 = fcvar_fractional_difference(x, d)
      allocate(levels(n, variables + merge(1, 0, selected%restricted_constant)))
      levels(:, :variables) = x
      if (selected%restricted_constant) levels(:, variables + 1) = 1.0_dp
      work = fcvar_fractional_lags(levels, b, 1)
      out%z1 = fcvar_fractional_difference(work, d - b)
      if (lag_order > 0) then
         work = fcvar_fractional_lags(x, b, lag_order)
         out%z2 = fcvar_fractional_difference(work, d)
      else
         allocate(out%z2(n, 0))
      end if
      out%z0 = out%z0(initial + 1:n, :)
      out%z1 = out%z1(initial + 1:n, :)
      out%z2 = out%z2(initial + 1:n, :)
      if (selected%unrestricted_constant) then
         allocate(out%z3(rows, 1))
         out%z3 = 1.0_dp
      else
         allocate(out%z3(rows, 0))
      end if
   end function fcvar_transform

   pure function fcvar_estimate_fixed(x, lag_order, rank, d, b, options, &
      level, level_initial) result(out)
      !! Concentrate the fixed-order FCVAR likelihood, profiling a level if enabled.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: b !! B.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      real(dp), intent(in), optional :: level(:) !! Model level or confidence level.
      real(dp), intent(in), optional :: level_initial(:) !! Level initial.
      type(fcvar_fit_t) :: out
      type(fcvar_options_t) :: selected, profile_options
      type(optimization_result_t) :: optimized
      real(dp), allocatable :: initial(:), shifted(:, :), selected_level(:)
      integer :: variables

      selected = fcvar_options_t()
      if (present(options)) selected = options
      profile_options = selected
      profile_options%level_parameter = .false.
      variables = size(x, 2)
      if (.not. selected%level_parameter) then
         out = fcvar_estimate_fixed_core(x, lag_order, rank, d, b, &
            profile_options)
         allocate(out%level(0))
         return
      end if
      if (variables < 1) then
         out%info = 1
         return
      end if
      if (present(level)) then
         if (size(level) /= variables .or. .not. all(ieee_is_finite(level))) then
            out%info = 1
            return
         end if
         selected_level = level
      else
         initial = x(1, :)
         if (present(level_initial)) then
            if (size(level_initial) /= variables .or. &
               .not. all(ieee_is_finite(level_initial))) then
               out%info = 1
               return
            end if
            initial = level_initial
         end if
         optimized = nelder_mead_minimize(level_objective, initial, 500, &
            1.0e-8_dp)
         if (.not. allocated(optimized%parameters) .or. &
            .not. ieee_is_finite(optimized%objective) .or. &
            optimized%objective >= 0.2_dp*huge(1.0_dp)) then
            out%info = 2
            return
         end if
         selected_level = optimized%parameters
      end if
      shifted = x - spread(selected_level, 1, size(x, 1))
      out = fcvar_estimate_fixed_core(shifted, lag_order, rank, d, b, &
         profile_options)
      if (out%info == 0) out%level = selected_level

   contains

      pure function level_objective(candidate) result(value)
         !! Return negative fixed-order likelihood for one level vector.
         real(dp), intent(in) :: candidate(:) !! Candidate.
         real(dp) :: value
         type(fcvar_fit_t) :: candidate_fit
         real(dp), allocatable :: candidate_shifted(:, :)

         candidate_shifted = x - spread(candidate, 1, size(x, 1))
         candidate_fit = fcvar_estimate_fixed_core(candidate_shifted, &
            lag_order, rank, d, b, profile_options)
         if (candidate_fit%info /= 0) then
            value = 0.25_dp*huge(1.0_dp)
         else
            value = -candidate_fit%log_likelihood
         end if
      end function level_objective

   end function fcvar_estimate_fixed

   pure function fcvar_estimate_fixed_core(x, lag_order, rank, d, b, options) &
      result(out)
      !! Concentrate the Gaussian FCVAR likelihood for already level-shifted data.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: b !! B.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      type(fcvar_fit_t) :: out
      type(fcvar_options_t) :: selected
      type(fcvar_transform_t) :: transformed
      real(dp), allocatable :: z0hat(:, :), z1hat(:, :), z2hat(:, :)
      real(dp), allocatable :: r0(:, :), r1(:, :), target(:, :)
      real(dp), allocatable :: s00(:, :), s01(:, :), s10(:, :), s11(:, :)
      real(dp), allocatable :: s00_inverse(:, :), s11_inverse(:, :)
      real(dp), allocatable :: cross_inverse(:, :), beta_star(:, :)
      real(dp), allocatable :: middle(:, :), middle_inverse(:, :)
      real(dp), allocatable :: lower(:, :), lower_inverse(:, :)
      real(dp), allocatable :: eigenmatrix(:, :), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: normalization(:, :), normalization_inverse(:, :)
      real(dp), allocatable :: pi_star(:, :)
      real(dp) :: sample, log_determinant, pi_value
      integer :: n, variables, p1, rows, status

      selected = fcvar_options_t()
      if (present(options)) selected = options
      n = size(x, 1)
      variables = size(x, 2)
      rows = n - selected%initial_values
      if (n < 2 .or. variables < 1 .or. lag_order < 0 .or. rank < 0 .or. &
         rank > variables .or. rows <= variables*lag_order) then
         out%info = 1
         return
      end if
      transformed = fcvar_transform(x, lag_order, d, b, selected)
      if (transformed%info /= 0) then
         out%info = 10 + transformed%info
         return
      end if
      z0hat = transformed%z0
      z1hat = transformed%z1
      z2hat = transformed%z2
      if (selected%unrestricted_constant) then
         z0hat = z0hat - spread(sum(z0hat, dim=1)/real(rows, dp), 1, rows)
         z1hat = z1hat - spread(sum(z1hat, dim=1)/real(rows, dp), 1, rows)
         if (lag_order > 0) then
            z2hat = z2hat - spread(sum(z2hat, dim=1)/real(rows, dp), 1, rows)
         end if
      end if
      if (lag_order > 0) then
         call invert_matrix(matmul(transpose(z2hat), z2hat), cross_inverse, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         r0 = z0hat - matmul(z2hat, matmul(cross_inverse, &
            matmul(transpose(z2hat), z0hat)))
         r1 = z1hat - matmul(z2hat, matmul(cross_inverse, &
            matmul(transpose(z2hat), z1hat)))
      else
         r0 = z0hat
         r1 = z1hat
      end if
      sample = real(rows, dp)
      s00 = matmul(transpose(r0), r0)/sample
      s01 = matmul(transpose(r0), r1)/sample
      s10 = transpose(s01)
      s11 = matmul(transpose(r1), r1)/sample
      p1 = size(s11, 1)
      allocate(out%alpha(variables, rank), out%beta(variables, rank))
      allocate(out%rho(merge(rank, 0, selected%restricted_constant)))
      allocate(out%pi(variables, variables))
      if (rank == 0) then
         out%alpha = 0.0_dp
         out%beta = 0.0_dp
         out%pi = 0.0_dp
         out%covariance = s00
         allocate(beta_star(p1, 0), pi_star(variables, p1))
         pi_star = 0.0_dp
      else if (rank == variables) then
         call invert_matrix(s11, s11_inverse, status)
         if (status /= 0) then
            out%info = 30 + status
            return
         end if
         pi_star = matmul(s01, s11_inverse)
         beta_star = transpose(pi_star)
         out%alpha = identity_matrix(variables)
         out%beta = beta_star(:variables, :)
         out%pi = pi_star(:, :variables)
         out%covariance = symmetrize(s00 - matmul(s01, matmul(s11_inverse, s10)))
      else
         call invert_matrix(s00, s00_inverse, status)
         if (status /= 0) then
            out%info = 40 + status
            return
         end if
         call cholesky_lower(s11, lower, status)
         if (status /= 0) then
            out%info = 50 + status
            return
         end if
         call invert_matrix(lower, lower_inverse, status)
         if (status /= 0) then
            out%info = 60 + status
            return
         end if
         eigenmatrix = matmul(lower_inverse, matmul(s10, &
            matmul(s00_inverse, matmul(s01, transpose(lower_inverse)))))
         call symmetric_eigen(symmetrize(eigenmatrix), eigenvalues, eigenvectors, status)
         if (status /= 0) then
            out%info = 70 + status
            return
         end if
         beta_star = matmul(transpose(lower_inverse), eigenvectors(:, :rank))
         middle = matmul(transpose(beta_star), matmul(s11, beta_star))
         call invert_matrix(middle, middle_inverse, status)
         if (status /= 0) then
            out%info = 80 + status
            return
         end if
         out%alpha = matmul(s01, matmul(beta_star, middle_inverse))
         out%covariance = s00 - matmul(out%alpha, &
            matmul(middle, transpose(out%alpha)))
         out%beta = beta_star(:variables, :)
         normalization = out%beta(:rank, :)
         call invert_matrix(normalization, normalization_inverse, status)
         if (status /= 0) then
            out%info = 90 + status
            return
         end if
         beta_star = matmul(beta_star, normalization_inverse)
         out%beta = beta_star(:variables, :)
         out%alpha = matmul(out%alpha, transpose(normalization))
         pi_star = matmul(out%alpha, transpose(beta_star))
         out%pi = pi_star(:, :variables)
         out%covariance = symmetrize(out%covariance)
      end if
      if (selected%restricted_constant) out%rho = beta_star(p1, :)
      if (rank == 0) then
         pi_star = 0.0_dp
      else
         pi_star = matmul(out%alpha, transpose(beta_star))
      end if
      allocate(out%gamma(variables, variables*lag_order))
      if (lag_order > 0) then
         target = z0hat - matmul(z1hat, transpose(pi_star))
         out%gamma = transpose(matmul(cross_inverse, &
            matmul(transpose(z2hat), target)))
      end if
      allocate(out%unrestricted_constant(merge(variables, 0, &
         selected%unrestricted_constant)))
      if (selected%unrestricted_constant) then
         target = transformed%z0 - matmul(transformed%z1, transpose(pi_star))
         if (lag_order > 0) then
            target = target - matmul(transformed%z2, transpose(out%gamma))
         end if
         out%unrestricted_constant = sum(target, dim=1)/sample
      end if
      out%residuals = transformed%z0 - matmul(transformed%z1, transpose(pi_star))
      if (lag_order > 0) then
         out%residuals = out%residuals - &
            matmul(transformed%z2, transpose(out%gamma))
      end if
      if (selected%unrestricted_constant) then
         out%residuals = out%residuals - &
            spread(out%unrestricted_constant, 1, rows)
      end if
      call covariance_log_determinant(out%covariance, log_determinant, status)
      if (status /= 0) then
         out%info = 100 + status
         return
      end if
      pi_value = acos(-1.0_dp)
      out%log_likelihood = -0.5_dp*sample*real(variables, dp)* &
         (log(2.0_dp*pi_value) + 1.0_dp) - 0.5_dp*sample*log_determinant
      out%d = d
      out%b = b
      out%lag_order = lag_order
      out%rank = rank
      out%observations = rows
   end function fcvar_estimate_fixed_core

   pure function fcvar_estimate_restricted_fixed(x, lag_order, rank, d, b, &
      restrictions, options, max_iterations, tolerance) result(out)
      !! Estimate restricted alpha and beta by the FCVAR switching algorithm.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: b !! B.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_restrictions_t), intent(in) :: restrictions !! Restrictions.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(fcvar_switching_estimation_t) :: out
      type(fcvar_options_t) :: selected
      type(fcvar_transform_t) :: transformed
      type(fcvar_fit_t) :: initial_fit, coefficients
      real(dp), allocatable :: z0hat(:, :), z1hat(:, :), z2hat(:, :)
      real(dp), allocatable :: r0(:, :), r1(:, :), cross_inverse(:, :)
      real(dp), allocatable :: s00(:, :), s01(:, :), s10(:, :), s11(:, :)
      real(dp), allocatable :: s11_inverse(:, :), omega(:, :), omega_inverse(:, :)
      real(dp), allocatable :: alpha(:, :), beta(:, :), beta_star(:, :)
      real(dp), allocatable :: old_alpha(:, :), old_beta(:, :), candidate_alpha(:, :)
      real(dp), allocatable :: candidate_beta(:, :), candidate_omega(:, :)
      real(dp), allocatable :: middle(:, :), middle_inverse(:, :)
      real(dp), allocatable :: alpha_cross(:, :), alpha_cross_inverse(:, :)
      real(dp), allocatable :: unrestricted(:), adjusted(:), weight(:, :)
      real(dp), allocatable :: alpha_matrix(:, :), alpha_value(:)
      real(dp), allocatable :: beta_matrix(:, :), beta_value(:), target(:, :)
      real(dp), allocatable :: supplied_beta_matrix(:, :), supplied_beta_value(:)
      real(dp), allocatable :: old_pi(:, :), current_pi(:, :)
      real(dp) :: selected_tolerance, old_criterion, criterion, candidate_criterion
      real(dp) :: likelihood_change, coefficient_change, sample
      real(dp), parameter :: lambda(5) = [1.0_dp, 1.2_dp, 2.0_dp, 4.0_dp, 8.0_dp]
      integer :: variables, beta_rows, rows, limit, iteration, status, candidate

      selected = fcvar_options_t()
      if (present(options)) selected = options
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-8_dp
      if (present(tolerance)) selected_tolerance = tolerance
      variables = size(x, 2)
      beta_rows = variables + merge(1, 0, selected%restricted_constant)
      if (rank < 1 .or. rank > variables .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. selected%level_parameter .or. &
         allocated(restrictions%parameter_matrix) .or. &
         allocated(restrictions%order_matrix) .or. &
         allocated(restrictions%level_matrix) .or. &
         .not. all_constraints_finite(restrictions)) then
         out%info = 1
         return
      end if
      if ((allocated(restrictions%alpha_matrix) .neqv. &
         allocated(restrictions%alpha_value)) .or. &
         (allocated(restrictions%beta_matrix) .neqv. &
         allocated(restrictions%beta_value))) then
         out%info = 2
         return
      end if
      if (allocated(restrictions%alpha_matrix)) then
         if (size(restrictions%alpha_matrix, 2) /= variables*rank .or. &
            size(restrictions%alpha_matrix, 1) /= &
            size(restrictions%alpha_value)) then
            out%info = 3
            return
         end if
         alpha_matrix = restrictions%alpha_matrix
         alpha_value = restrictions%alpha_value
      else
         allocate(alpha_matrix(0, variables*rank), alpha_value(0))
      end if
      if (allocated(restrictions%beta_matrix)) then
         if (size(restrictions%beta_matrix, 2) /= beta_rows*rank .or. &
            size(restrictions%beta_matrix, 1) /= &
            size(restrictions%beta_value)) then
            out%info = 4
            return
         end if
         beta_matrix = restrictions%beta_matrix
         beta_value = restrictions%beta_value
      else
         allocate(beta_matrix(0, beta_rows*rank), beta_value(0))
      end if
      supplied_beta_matrix = beta_matrix
      supplied_beta_value = beta_value
      deallocate(beta_matrix, beta_value)
      allocate(beta_matrix(size(supplied_beta_matrix, 1) + rank*rank, &
         beta_rows*rank))
      allocate(beta_value(size(supplied_beta_value) + rank*rank))
      beta_matrix = 0.0_dp
      beta_value = 0.0_dp
      if (size(supplied_beta_matrix, 1) > 0) then
         beta_matrix(:size(supplied_beta_matrix, 1), :) = supplied_beta_matrix
         beta_value(:size(supplied_beta_value)) = supplied_beta_value
      end if
      do candidate = 1, rank
         do iteration = 1, rank
            rows = size(supplied_beta_matrix, 1) + &
               (candidate - 1)*rank + iteration
            beta_matrix(rows, (candidate - 1)*beta_rows + iteration) = 1.0_dp
            beta_value(rows) = merge(1.0_dp, 0.0_dp, iteration == candidate)
         end do
      end do
      if (size(alpha_matrix, 1) + size(beta_matrix, 1) < 1) then
         out%info = 5
         return
      end if
      initial_fit = fcvar_estimate_fixed(x, lag_order, rank, d, b, selected)
      if (initial_fit%info /= 0) then
         out%info = 10 + initial_fit%info
         return
      end if
      transformed = fcvar_transform(x, lag_order, d, b, selected)
      if (transformed%info /= 0) then
         out%info = 20 + transformed%info
         return
      end if
      rows = size(transformed%z0, 1)
      sample = real(rows, dp)
      z0hat = transformed%z0
      z1hat = transformed%z1
      z2hat = transformed%z2
      if (selected%unrestricted_constant) then
         z0hat = z0hat - spread(sum(z0hat, dim=1)/sample, 1, rows)
         z1hat = z1hat - spread(sum(z1hat, dim=1)/sample, 1, rows)
         if (lag_order > 0) z2hat = z2hat - &
            spread(sum(z2hat, dim=1)/sample, 1, rows)
      end if
      if (lag_order > 0) then
         call invert_matrix(matmul(transpose(z2hat), z2hat), cross_inverse, status)
         if (status /= 0) then
            out%info = 30 + status
            return
         end if
         r0 = z0hat - matmul(z2hat, matmul(cross_inverse, &
            matmul(transpose(z2hat), z0hat)))
         r1 = z1hat - matmul(z2hat, matmul(cross_inverse, &
            matmul(transpose(z2hat), z1hat)))
      else
         r0 = z0hat
         r1 = z1hat
      end if
      s00 = matmul(transpose(r0), r0)/sample
      s01 = matmul(transpose(r0), r1)/sample
      s10 = transpose(s01)
      s11 = matmul(transpose(r1), r1)/sample
      call invert_matrix(s11, s11_inverse, status)
      if (status /= 0) then
         out%info = 40 + status
         return
      end if
      alpha = initial_fit%alpha
      allocate(beta_star(beta_rows, rank))
      beta_star(:variables, :) = initial_fit%beta
      if (selected%restricted_constant) beta_star(beta_rows, :) = initial_fit%rho
      beta = beta_star
      omega = initial_fit%covariance
      criterion = -huge(1.0_dp)
      do iteration = 1, limit
         old_alpha = alpha
         old_beta = beta
         old_pi = matmul(alpha, transpose(beta))
         old_criterion = criterion
         call invert_matrix(omega, omega_inverse, status)
         if (status /= 0) exit
         middle = matmul(transpose(beta), matmul(s11, beta))
         call invert_matrix(middle, middle_inverse, status)
         if (status /= 0) exit
         alpha = matmul(s01, matmul(beta, middle_inverse))
         weight = kronecker_product(middle, omega_inverse)
         unrestricted = reshape(alpha, [variables*rank])
         call constrained_quadratic_projection(unrestricted, weight, &
            alpha_matrix, alpha_value, adjusted, status)
         if (status /= 0) exit
         alpha = reshape(adjusted, [variables, rank])
         omega = switching_covariance(s00, s01, s11, alpha, beta)
         call invert_matrix(omega, omega_inverse, status)
         if (status /= 0) exit
         alpha_cross = matmul(transpose(alpha), matmul(omega_inverse, alpha))
         call invert_matrix(alpha_cross, alpha_cross_inverse, status)
         if (status /= 0) exit
         beta = matmul(s11_inverse, matmul(s10, matmul(omega_inverse, &
            matmul(alpha, alpha_cross_inverse))))
         weight = kronecker_product(alpha_cross, s11)
         unrestricted = reshape(beta, [beta_rows*rank])
         call constrained_quadratic_projection(unrestricted, weight, &
            beta_matrix, beta_value, adjusted, status)
         if (status /= 0) exit
         beta = reshape(adjusted, [beta_rows, rank])
         omega = switching_covariance(s00, s01, s11, alpha, beta)
         call covariance_log_determinant(omega, criterion, status)
         if (status /= 0) exit
         criterion = -criterion
         likelihood_change = huge(1.0_dp)
         if (iteration > 1) likelihood_change = &
            (criterion - old_criterion)/(1.0_dp + abs(old_criterion))
         if (iteration > 1 .and. likelihood_change < 0.01_dp) then
            candidate_alpha = alpha
            candidate_beta = beta
            do candidate = 2, size(lambda)
               candidate_alpha = old_alpha + lambda(candidate)*(alpha - old_alpha)
               candidate_beta = old_beta + lambda(candidate)*(beta - old_beta)
               candidate_omega = switching_covariance(s00, s01, s11, &
                  candidate_alpha, candidate_beta)
               call covariance_log_determinant(candidate_omega, &
                  candidate_criterion, status)
               if (status /= 0) cycle
               candidate_criterion = -candidate_criterion
               if (candidate_criterion > criterion) then
                  alpha = candidate_alpha
                  beta = candidate_beta
                  omega = candidate_omega
                  criterion = candidate_criterion
               end if
            end do
            status = 0
         end if
         current_pi = matmul(alpha, transpose(beta))
         coefficient_change = maxval(abs(current_pi - old_pi)/ &
            (1.0_dp + abs(old_pi)))
         likelihood_change = (criterion - old_criterion)/ &
            (1.0_dp + abs(old_criterion))
         if (iteration > 1 .and. abs(likelihood_change) <= selected_tolerance .and. &
            coefficient_change <= sqrt(selected_tolerance)) then
            out%converged = .true.
            exit
         end if
      end do
      out%iterations = min(iteration, limit)
      if (status /= 0) then
         out%info = 50 + status
         return
      end if
      coefficients = initial_fit
      coefficients%alpha = alpha
      coefficients%beta = beta(:variables, :)
      if (selected%restricted_constant) coefficients%rho = beta(beta_rows, :)
      coefficients%pi = matmul(alpha, transpose(beta(:variables, :)))
      target = z0hat - matmul(z1hat, &
         transpose(matmul(alpha, transpose(beta))))
      if (lag_order > 0) coefficients%gamma = transpose(matmul(cross_inverse, &
         matmul(transpose(z2hat), target)))
      if (selected%unrestricted_constant) then
         target = transformed%z0 - matmul(transformed%z1, &
            transpose(matmul(alpha, transpose(beta))))
         if (lag_order > 0) target = target - &
            matmul(transformed%z2, transpose(coefficients%gamma))
         coefficients%unrestricted_constant = sum(target, dim=1)/sample
      end if
      out%fit = refresh_fcvar_fit(x, coefficients, selected)
      if (out%fit%info /= 0) out%info = 60 + out%fit%info
   end function fcvar_estimate_restricted_fixed

   pure function fcvar_log_likelihood_fixed(x, lag_order, rank, d, b, options) &
      result(log_likelihood)
      !! Return the concentrated Gaussian log likelihood at fixed d and b.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: b !! B.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      real(dp) :: log_likelihood
      type(fcvar_fit_t) :: fit

      fit = fcvar_estimate_fixed(x, lag_order, rank, d, b, options)
      if (fit%info == 0) then
         log_likelihood = fit%log_likelihood
      else
         log_likelihood = -huge(1.0_dp)
      end if
   end function fcvar_log_likelihood_fixed

   pure function fcvar_likelihood_grid(x, lag_order, rank, d_grid, b_grid, &
      options, constrain_b) result(out)
      !! Evaluate a grid and select either its global or highest-b local maximum.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: d_grid(:) !! D grid.
      real(dp), intent(in) :: b_grid(:) !! B grid.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      type(fcvar_grid_t) :: out
      type(fcvar_options_t) :: selected
      type(fcvar_fit_t) :: fit
      logical :: constrained
      logical, allocatable :: feasible(:, :)
      real(dp), allocatable :: level_start(:)
      real(dp) :: value
      integer :: d_index, b_index, local_index, selected_local

      selected = fcvar_options_t()
      if (present(options)) selected = options
      constrained = .true.
      if (present(constrain_b)) constrained = constrain_b
      if (size(d_grid) < 1 .or. size(b_grid) < 1 .or. &
         .not. all(ieee_is_finite(d_grid)) .or. &
         .not. all(ieee_is_finite(b_grid))) then
         out%info = 1
         return
      end if
      out%d = d_grid
      out%b = b_grid
      allocate(out%log_likelihood(size(d_grid), size(b_grid)), &
         feasible(size(d_grid), size(b_grid)))
      out%log_likelihood = -huge(1.0_dp)
      feasible = .false.
      if (selected%level_parameter) then
         allocate(out%level(size(x, 2), size(d_grid), size(b_grid)))
         out%level = 0.0_dp
         level_start = x(1, :)
      else
         allocate(out%level(0, size(d_grid), size(b_grid)))
      end if
      do b_index = 1, size(b_grid)
         do d_index = 1, size(d_grid)
            if (constrained .and. b_grid(b_index) > d_grid(d_index)) cycle
            if (selected%level_parameter) then
               fit = fcvar_estimate_fixed(x, lag_order, rank, d_grid(d_index), &
                  b_grid(b_index), selected, level_initial=level_start)
            else
               fit = fcvar_estimate_fixed(x, lag_order, rank, d_grid(d_index), &
                  b_grid(b_index), selected)
            end if
            if (fit%info /= 0) cycle
            value = fit%log_likelihood
            out%log_likelihood(d_index, b_index) = value
            feasible(d_index, b_index) = .true.
            if (selected%level_parameter) then
               out%level(:, d_index, b_index) = fit%level
               level_start = fit%level
            end if
            if (value > out%global_log_likelihood .or. &
               (value == out%global_log_likelihood .and. &
               b_grid(b_index) > out%global_b)) then
               out%global_log_likelihood = value
               out%global_d = d_grid(d_index)
               out%global_b = b_grid(b_index)
               out%global_d_index = d_index
               out%global_b_index = b_index
            end if
         end do
      end do
      if (out%global_log_likelihood <= -0.5_dp*huge(1.0_dp)) then
         out%info = 2
         return
      end if
      call fcvar_find_local_maxima(out%log_likelihood, feasible, &
         out%local_d_index, out%local_b_index)
      allocate(out%local_d(size(out%local_d_index)), &
         out%local_b(size(out%local_d_index)), &
         out%local_log_likelihood(size(out%local_d_index)))
      do local_index = 1, size(out%local_d_index)
         out%local_d(local_index) = d_grid(out%local_d_index(local_index))
         out%local_b(local_index) = b_grid(out%local_b_index(local_index))
         out%local_log_likelihood(local_index) = out%log_likelihood( &
            out%local_d_index(local_index), out%local_b_index(local_index))
      end do
      out%best_d = out%global_d
      out%best_b = out%global_b
      out%best_log_likelihood = out%global_log_likelihood
      out%best_d_index = out%global_d_index
      out%best_b_index = out%global_b_index
      if (selected%prefer_high_b_local_max .and. size(out%local_d) > 0) then
         selected_local = 1
         do local_index = 2, size(out%local_d)
            if (out%local_b(local_index) > out%local_b(selected_local) .or. &
               (out%local_b(local_index) == out%local_b(selected_local) .and. &
               out%local_log_likelihood(local_index) > &
               out%local_log_likelihood(selected_local)) .or. &
               (out%local_b(local_index) == out%local_b(selected_local) .and. &
               out%local_log_likelihood(local_index) == &
               out%local_log_likelihood(selected_local) .and. &
               out%local_d(local_index) > out%local_d(selected_local))) then
               selected_local = local_index
            end if
         end do
         out%best_d = out%local_d(selected_local)
         out%best_b = out%local_b(selected_local)
         out%best_log_likelihood = out%local_log_likelihood(selected_local)
         out%best_d_index = out%local_d_index(selected_local)
         out%best_b_index = out%local_b_index(selected_local)
         out%used_local_maximum = .true.
      end if
      if (selected%level_parameter) then
         out%best_level = out%level(:, out%best_d_index, out%best_b_index)
      else
         allocate(out%best_level(0))
      end if
   end function fcvar_likelihood_grid

   pure subroutine fcvar_find_local_maxima(values, feasible, row_index, &
      column_index)
      !! Locate strict maxima against every feasible neighboring grid cell.
      real(dp), intent(in) :: values(:, :) !! Input values.
      logical, intent(in) :: feasible(:, :) !! Flag controlling feasible.
      integer, allocatable, intent(out) :: row_index(:) !! Index of row.
      integer, allocatable, intent(out) :: column_index(:) !! Index of column.
      logical, allocatable :: is_maximum(:, :)
      integer :: row, column, neighbor_row, neighbor_column, count_maximum

      if (any(shape(values) /= shape(feasible))) then
         allocate(row_index(0), column_index(0))
         return
      end if
      allocate(is_maximum, mold=feasible)
      is_maximum = feasible
      do column = 1, size(values, 2)
         do row = 1, size(values, 1)
            if (.not. feasible(row, column)) cycle
            do neighbor_column = max(1, column - 1), min(size(values, 2), &
               column + 1)
               do neighbor_row = max(1, row - 1), min(size(values, 1), row + 1)
                  if (neighbor_row == row .and. neighbor_column == column) cycle
                  if (feasible(neighbor_row, neighbor_column) .and. &
                     values(row, column) <= &
                     values(neighbor_row, neighbor_column)) then
                     is_maximum(row, column) = .false.
                     exit
                  end if
               end do
               if (.not. is_maximum(row, column)) exit
            end do
         end do
      end do
      count_maximum = count(is_maximum)
      allocate(row_index(count_maximum), column_index(count_maximum))
      count_maximum = 0
      do column = 1, size(values, 2)
         do row = 1, size(values, 1)
            if (.not. is_maximum(row, column)) cycle
            count_maximum = count_maximum + 1
            row_index(count_maximum) = row
            column_index(count_maximum) = column
         end do
      end do
   end subroutine fcvar_find_local_maxima

   pure function fcvar_estimate(x, lag_order, rank, d_bounds, b_bounds, &
      options, constrain_b, grid_points, max_iterations) result(out)
      !! Estimate d and b by grid initialization and bounded BFGS refinement.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: d_bounds(2) !! D bounds.
      real(dp), intent(in) :: b_bounds(2) !! B bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      integer, intent(in), optional :: grid_points !! Grid points.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_estimation_t) :: out
      type(fcvar_options_t) :: selected
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: d_grid(:), b_grid(:), initial(:)
      real(dp) :: d_value, b_value, refined_likelihood, best_likelihood
      integer :: points, limit, sweep
      logical :: constrained

      selected = fcvar_options_t()
      if (present(options)) selected = options
      constrained = .true.
      if (present(constrain_b)) constrained = constrain_b
      points = 11
      if (present(grid_points)) points = grid_points
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      if (points < 2 .or. limit < 1 .or. &
         d_bounds(1) >= d_bounds(2) .or. b_bounds(1) >= b_bounds(2) .or. &
         .not. all(ieee_is_finite(d_bounds)) .or. &
         .not. all(ieee_is_finite(b_bounds))) then
         out%info = 1
         return
      end if
      d_grid = regular_grid(d_bounds, points)
      b_grid = regular_grid(b_bounds, points)
      out%grid = fcvar_likelihood_grid(x, lag_order, rank, d_grid, b_grid, &
         options, constrained)
      if (out%grid%info /= 0) then
         out%info = 10 + out%grid%info
         return
      end if
      allocate(initial(2))
      call encode_orders(out%grid%best_d, out%grid%best_b, d_bounds, b_bounds, &
         constrained, initial)
      optimization = bfgs_minimize_fd(objective, initial, limit, 1.0e-7_dp, &
         1.0e-5_dp)
      out%optimizer_coordinates = optimization%parameters
      out%iterations = optimization%iterations
      out%optimizer_info = optimization%info
      out%converged = optimization%converged
      call decode_orders(optimization%parameters, d_bounds, b_bounds, &
         constrained, d_value, b_value)
      refined_likelihood = fcvar_log_likelihood_fixed(x, lag_order, rank, &
         d_value, b_value, options)
      if (refined_likelihood <= out%grid%best_log_likelihood) then
         d_value = out%grid%best_d
         b_value = out%grid%best_b
         refined_likelihood = out%grid%best_log_likelihood
      end if
      if (.not. selected%prefer_high_b_local_max) then
         do sweep = 1, 4
            call maximize_b_order(x, lag_order, rank, d_value, b_bounds, &
               selected, constrained, b_value, best_likelihood)
            refined_likelihood = best_likelihood
            call maximize_d_order(x, lag_order, rank, b_value, d_bounds, &
               selected, constrained, d_value, best_likelihood)
            refined_likelihood = best_likelihood
         end do
      end if
      if (selected%level_parameter) then
         out%fit = fcvar_estimate_fixed(x, lag_order, rank, d_value, b_value, &
            selected, level_initial=out%grid%best_level)
      else
         out%fit = fcvar_estimate_fixed(x, lag_order, rank, d_value, b_value, &
            selected)
      end if
      if (out%fit%info /= 0) out%info = 20 + out%fit%info

   contains

      pure function objective(coordinates) result(value)
         !! Return the negative concentrated likelihood in bounded coordinates.
         real(dp), intent(in) :: coordinates(:) !! Coordinates.
         real(dp) :: value, current_d, current_b, likelihood

         call decode_orders(coordinates, d_bounds, b_bounds, constrained, &
            current_d, current_b)
         likelihood = fcvar_log_likelihood_fixed(x, lag_order, rank, current_d, &
            current_b, options)
         if (likelihood <= -0.5_dp*huge(1.0_dp)) then
            value = 0.25_dp*huge(1.0_dp)
         else
            value = -likelihood
         end if
      end function objective

   end function fcvar_estimate

   pure function fcvar_estimate_equal_orders(x, lag_order, rank, bounds, &
      options, grid_points) result(out)
      !! Estimate the default FCVAR restriction d equals b on a bounded interval.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      integer, intent(in), optional :: grid_points !! Grid points.
      type(fcvar_estimation_t) :: out
      real(dp), allocatable :: grid(:)
      real(dp) :: order_value, likelihood
      integer :: points, index

      points = 21
      if (present(grid_points)) points = grid_points
      if (points < 2 .or. bounds(1) >= bounds(2) .or. &
         .not. all(ieee_is_finite(bounds))) then
         out%info = 1
         return
      end if
      grid = regular_grid(bounds, points)
      out%grid%d = grid
      out%grid%b = grid
      allocate(out%grid%log_likelihood(points, points))
      out%grid%log_likelihood = -huge(1.0_dp)
      do index = 1, points
         likelihood = fcvar_log_likelihood_fixed(x, lag_order, rank, &
            grid(index), grid(index), options)
         out%grid%log_likelihood(index, index) = likelihood
         if (likelihood > out%grid%best_log_likelihood) then
            out%grid%best_log_likelihood = likelihood
            out%grid%best_d = grid(index)
            out%grid%best_b = grid(index)
         end if
      end do
      if (out%grid%best_log_likelihood <= -0.5_dp*huge(1.0_dp)) then
         out%info = 2
         out%grid%info = 2
         return
      end if
      call maximize_equal_order(x, lag_order, rank, bounds, options, &
         order_value, likelihood)
      if (likelihood < out%grid%best_log_likelihood) then
         order_value = out%grid%best_d
      end if
      out%fit = fcvar_estimate_fixed(x, lag_order, rank, order_value, &
         order_value, options)
      allocate(out%optimizer_coordinates(0))
      out%iterations = 60
      out%converged = out%fit%info == 0
      if (out%fit%info /= 0) out%info = 10 + out%fit%info
   end function fcvar_estimate_equal_orders

   pure function fcvar_estimate_constrained_orders(x, lag_order, rank, bounds, &
      restrictions, options, constrain_b, grid_points, max_iterations) result(out)
      !! Estimate fractional orders under affine equalities and inequalities.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2, 2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_restrictions_t), intent(in) :: restrictions !! Restrictions.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      integer, intent(in), optional :: grid_points !! Grid points.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_order_estimation_t) :: out
      type(fcvar_options_t) :: selected
      type(optimization_result_t) :: optimized
      real(dp), allocatable :: equality(:, :), equality_value(:)
      real(dp), allocatable :: inequality(:, :), inequality_value(:)
      real(dp), allocatable :: particular(:), basis(:, :), initial(:)
      real(dp), allocatable :: reduced_matrix(:, :), reduced_value(:)
      real(dp), allocatable :: d_grid(:), b_grid(:), candidate(:)
      real(dp) :: lower, upper, value, best_value, tolerance
      real(dp) :: line_point(2), line_direction(2), denominator
      real(dp) :: edge_lower, edge_upper, edge_value
      integer :: equality_count, inequality_count, status, points, limit
      integer :: d_index, b_index, constraint, free_dimensions
      logical :: constrained

      selected = fcvar_options_t()
      if (present(options)) selected = options
      constrained = .true.
      if (present(constrain_b)) constrained = constrain_b
      points = 11
      if (present(grid_points)) points = grid_points
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      if (points < 2 .or. limit < 1 .or. lag_order < 0 .or. rank < 0 .or. &
         bounds(1, 1) >= bounds(2, 1) .or. &
         bounds(1, 2) >= bounds(2, 2) .or. &
         .not. all(ieee_is_finite(bounds)) .or. &
         (allocated(restrictions%order_matrix) .neqv. &
         allocated(restrictions%order_value)) .or. &
         (allocated(restrictions%order_inequality_matrix) .neqv. &
         allocated(restrictions%order_inequality_value))) then
         out%info = 1
         return
      end if
      equality_count = 0
      if (allocated(restrictions%order_matrix)) then
         equality_count = size(restrictions%order_matrix, 1)
         if (size(restrictions%order_matrix, 2) /= 2 .or. &
            size(restrictions%order_value) /= equality_count .or. &
            .not. all(ieee_is_finite(restrictions%order_matrix)) .or. &
            .not. all(ieee_is_finite(restrictions%order_value))) then
            out%info = 2
            return
         end if
         equality = restrictions%order_matrix
         equality_value = restrictions%order_value
      else
         allocate(equality(0, 2), equality_value(0))
      end if
      inequality_count = 4 + merge(1, 0, constrained)
      if (allocated(restrictions%order_inequality_matrix)) then
         if (size(restrictions%order_inequality_matrix, 2) /= 2 .or. &
            size(restrictions%order_inequality_value) /= &
            size(restrictions%order_inequality_matrix, 1) .or. &
            .not. all(ieee_is_finite( &
            restrictions%order_inequality_matrix)) .or. &
            .not. all(ieee_is_finite( &
            restrictions%order_inequality_value))) then
            out%info = 3
            return
         end if
         inequality_count = inequality_count + &
            size(restrictions%order_inequality_matrix, 1)
      end if
      allocate(inequality(inequality_count, 2), &
         inequality_value(inequality_count))
      inequality = 0.0_dp
      inequality(1, 1) = 1.0_dp
      inequality_value(1) = bounds(1, 1)
      inequality(2, 1) = -1.0_dp
      inequality_value(2) = -bounds(2, 1)
      inequality(3, 2) = 1.0_dp
      inequality_value(3) = bounds(1, 2)
      inequality(4, 2) = -1.0_dp
      inequality_value(4) = -bounds(2, 2)
      status = 4
      if (constrained) then
         status = 5
         inequality(5, :) = [1.0_dp, -1.0_dp]
         inequality_value(5) = 0.0_dp
      end if
      if (allocated(restrictions%order_inequality_matrix)) then
         inequality(status + 1:, :) = restrictions%order_inequality_matrix
         inequality_value(status + 1:) = &
            restrictions%order_inequality_value
      end if
      initial = 0.5_dp*(bounds(1, :) + bounds(2, :))
      if (equality_count > 0) then
         call constraint_geometry(equality, equality_value, initial, particular, &
            basis, out%equality_rank, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
      else
         particular = [0.0_dp, 0.0_dp]
         basis = identity_matrix(2)
         out%equality_rank = 0
      end if
      free_dimensions = size(basis, 2)
      out%free_dimensions = free_dimensions
      reduced_matrix = matmul(inequality, basis)
      reduced_value = inequality_value - matmul(inequality, particular)
      tolerance = 1.0e-10_dp
      if (free_dimensions == 0) then
         if (.not. orders_feasible(particular, inequality, inequality_value, &
            tolerance)) then
            out%info = 20
            return
         end if
         out%starting_d = particular(1)
         out%starting_b = particular(2)
         out%fit = fcvar_estimate_fixed(x, lag_order, rank, particular(1), &
            particular(2), selected)
         if (out%fit%info /= 0) out%info = 30 + out%fit%info
         out%converged = out%info == 0
         return
      end if
      if (free_dimensions == 1) then
         call feasible_line_interval(reduced_matrix(:, 1), reduced_value, &
            lower, upper, status)
         if (status /= 0 .or. lower > upper) then
            out%info = 20 + max(1, status)
            return
         end if
         call maximize_fcvar_order_line(x, lag_order, rank, selected, lower, &
            upper, particular, basis(:, 1), candidate, best_value)
         out%starting_d = candidate(1)
         out%starting_b = candidate(2)
         out%fit = fcvar_estimate_fixed(x, lag_order, rank, candidate(1), &
            candidate(2), selected)
         if (out%fit%info /= 0) out%info = 30 + out%fit%info
         out%iterations = 60
         out%converged = out%info == 0
         return
      end if
      d_grid = regular_grid(bounds(:, 1), points)
      b_grid = regular_grid(bounds(:, 2), points)
      best_value = -huge(1.0_dp)
      allocate(candidate(2))
      do b_index = 1, points
         do d_index = 1, points
            initial = [d_grid(d_index), b_grid(b_index)]
            if (.not. orders_feasible(initial, inequality, inequality_value, &
               tolerance)) cycle
            value = fcvar_log_likelihood_fixed(x, lag_order, rank, initial(1), &
               initial(2), selected)
            if (value > best_value) then
               best_value = value
               candidate = initial
            end if
         end do
      end do
      if (best_value <= -0.5_dp*huge(1.0_dp)) then
         out%info = 21
         return
      end if
      out%starting_d = candidate(1)
      out%starting_b = candidate(2)
      optimized = nelder_mead_minimize(order_objective, candidate, limit, &
         1.0e-8_dp, 0.05_dp)
      out%iterations = optimized%iterations
      out%optimizer_info = optimized%info
      if (allocated(optimized%parameters) .and. &
         orders_feasible(optimized%parameters, inequality, inequality_value, &
         1.0e-8_dp) .and. -optimized%objective > best_value) then
         candidate = optimized%parameters
         best_value = -optimized%objective
      end if
      do constraint = 1, inequality_count
         denominator = dot_product(inequality(constraint, :), &
            inequality(constraint, :))
         if (denominator <= tiny(1.0_dp)) cycle
         line_point = inequality_value(constraint)* &
            inequality(constraint, :)/denominator
         line_direction = [-inequality(constraint, 2), &
            inequality(constraint, 1)]/sqrt(denominator)
         call boundary_line_interval(line_point, line_direction, inequality, &
            inequality_value, edge_lower, edge_upper, status)
         if (status /= 0 .or. edge_lower > edge_upper) cycle
         call maximize_fcvar_order_line(x, lag_order, rank, selected, &
            edge_lower, edge_upper, line_point, line_direction, initial, &
            edge_value)
         if (edge_value > best_value) then
            candidate = initial
            best_value = edge_value
         end if
      end do
      out%fit = fcvar_estimate_fixed(x, lag_order, rank, candidate(1), &
         candidate(2), selected)
      if (out%fit%info /= 0) out%info = 30 + out%fit%info
      out%converged = out%info == 0 .and. &
         (optimized%converged .or. best_value > -0.5_dp*huge(1.0_dp))

   contains

      pure function order_objective(orders) result(objective)
         !! Return penalized negative likelihood for a candidate order pair.
         real(dp), intent(in) :: orders(:) !! Orders.
         real(dp) :: objective

         if (size(orders) /= 2 .or. .not. orders_feasible(orders, inequality, &
            inequality_value, 1.0e-10_dp)) then
            objective = 0.25_dp*huge(1.0_dp)
         else
            objective = -fcvar_log_likelihood_fixed(x, lag_order, rank, &
               orders(1), orders(2), selected)
         end if
      end function order_objective

   end function fcvar_estimate_constrained_orders

   pure function fcvar_simulate_from_innovations(history, model, innovations, &
      options, level) result(out)
      !! Generate an FCVAR path from supplied innovations and starting history.
      real(dp), intent(in) :: history(:, :) !! History.
      real(dp), intent(in) :: innovations(:, :) !! Model innovations.
      type(fcvar_fit_t), intent(in) :: model !! Model specification.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      real(dp), intent(in), optional :: level(:) !! Model level or confidence level.
      type(fcvar_path_t) :: out
      type(fcvar_options_t) :: selected
      real(dp), allocatable :: work(:, :), y(:, :), z(:, :), term(:, :)
      real(dp), allocatable :: ones(:, :), level_value(:), constant_effect(:)
      integer :: initial, periods, variables, time, current

      selected = fcvar_options_t()
      if (present(options)) selected = options
      initial = size(history, 1)
      variables = size(history, 2)
      periods = size(innovations, 1)
      allocate(level_value(variables))
      level_value = 0.0_dp
      if (allocated(model%level)) then
         if (size(model%level) == variables) then
            level_value = model%level
         else if (size(model%level) /= 0) then
            out%info = 1
            return
         end if
      end if
      if (present(level)) then
         if (size(level) /= variables) then
            out%info = 1
            return
         end if
         level_value = level
      end if
      if (initial < 1 .or. variables < 1 .or. periods < 1 .or. &
         size(innovations, 2) /= variables .or. &
         any(shape(model%pi) /= [variables, variables]) .or. &
         size(model%alpha, 1) /= variables .or. &
         size(model%beta, 1) /= variables .or. &
         size(model%gamma, 1) /= variables .or. &
         size(model%gamma, 2) /= variables*model%lag_order .or. &
         .not. all(ieee_is_finite(history)) .or. &
         .not. all(ieee_is_finite(innovations)) .or. &
         .not. all(ieee_is_finite(level_value))) then
         out%info = 1
         return
      end if
      if (selected%restricted_constant .and. &
         size(model%rho) /= model%rank) then
         out%info = 2
         return
      end if
      if (selected%unrestricted_constant .and. &
         size(model%unrestricted_constant) /= variables) then
         out%info = 3
         return
      end if
      allocate(work(initial + periods, variables))
      work(:initial, :) = history
      do time = 1, periods
         current = initial + time
         work(current, :) = 0.0_dp
         allocate(y(current, variables))
         y = work(:current, :) - spread(level_value, 1, current)
         z = fcvar_fractional_lags(y, model%d, 1)
         if (model%rank > 0) then
            term = fcvar_fractional_difference( &
               fcvar_fractional_lags(y, model%b, 1), model%d - model%b)
            z = z + matmul(term, transpose(model%pi))
            if (selected%restricted_constant) then
               allocate(ones(current, 1))
               ones = 1.0_dp
               term = fcvar_fractional_difference( &
                  fcvar_fractional_lags(ones, model%b, 1), &
                  model%d - model%b)
               constant_effect = matmul(model%alpha, model%rho)
               z = z + spread(term(:, 1), 2, variables)* &
                  spread(constant_effect, 1, current)
               deallocate(ones, constant_effect)
            end if
         end if
         if (selected%unrestricted_constant) then
            z = z + spread(model%unrestricted_constant, 1, current)
         end if
         if (model%lag_order > 0) then
            term = fcvar_fractional_difference( &
               fcvar_fractional_lags(y, model%b, model%lag_order), model%d)
            z = z + matmul(term, transpose(model%gamma))
         end if
         z = z + spread(level_value, 1, current)
         work(current, :) = z(current, :) + innovations(time, :)
         deallocate(y, z)
         if (allocated(term)) deallocate(term)
      end do
      allocate(out%series(periods, variables), out%innovations(periods, variables))
      out%series = work(initial + 1:, :)
      out%innovations = innovations
   end function fcvar_simulate_from_innovations

   pure function fcvar_forecast(history, model, periods, options, level) &
      result(out)
      !! Recursively forecast an FCVAR model with zero future innovations.
      real(dp), intent(in) :: history(:, :) !! History.
      type(fcvar_fit_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: periods !! Periods.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      real(dp), intent(in), optional :: level(:) !! Model level or confidence level.
      type(fcvar_path_t) :: out
      real(dp), allocatable :: innovations(:, :)

      if (periods < 1) then
         out%info = 1
         return
      end if
      allocate(innovations(periods, size(history, 2)))
      innovations = 0.0_dp
      out = fcvar_simulate_from_innovations(history, model, innovations, &
         options, level)
   end function fcvar_forecast

   function fcvar_simulate(history, model, periods, options, level) result(out)
      !! Simulate FCVAR with standard-normal disturbances as in FCVARsim.
      real(dp), intent(in) :: history(:, :) !! History.
      type(fcvar_fit_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: periods !! Periods.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      real(dp), intent(in), optional :: level(:) !! Model level or confidence level.
      type(fcvar_path_t) :: out
      real(dp), allocatable :: draws(:, :), innovations(:, :)

      if (periods < 1 .or. size(history, 2) < 1) then
         out%info = 1
         return
      end if
      allocate(draws(size(history, 2), periods))
      call random_standard_normal_matrix(draws)
      innovations = transpose(draws)
      out = fcvar_simulate_from_innovations(history, model, innovations, &
         options, level)
   end function fcvar_simulate

   pure function fcvar_bootstrap_from_signs(history, model, residuals, signs, &
      options, level) result(out)
      !! Generate FCVAR's centered-residual wild bootstrap for supplied signs.
      real(dp), intent(in) :: history(:, :) !! History.
      real(dp), intent(in) :: residuals(:, :) !! Model residuals.
      real(dp), intent(in) :: signs(:) !! Signs.
      type(fcvar_fit_t), intent(in) :: model !! Model specification.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      real(dp), intent(in), optional :: level(:) !! Model level or confidence level.
      type(fcvar_path_t) :: out
      real(dp), allocatable :: centered(:, :), innovations(:, :)
      integer :: periods, variables

      periods = size(signs)
      variables = size(history, 2)
      if (periods < 1 .or. size(residuals, 1) < periods .or. &
         size(residuals, 2) /= variables .or. &
         .not. all(ieee_is_finite(residuals)) .or. &
         .not. all(ieee_is_finite(signs))) then
         out%info = 1
         return
      end if
      centered = residuals - spread(sum(residuals, dim=1)/ &
         real(size(residuals, 1), dp), 1, size(residuals, 1))
      allocate(innovations(periods, variables))
      innovations = centered(:periods, :)*spread(signs, 2, variables)
      out = fcvar_simulate_from_innovations(history, model, innovations, &
         options, level)
   end function fcvar_bootstrap_from_signs

   function fcvar_bootstrap(history, model, residuals, periods, options, &
      level) result(out)
      !! Draw Rademacher signs and generate FCVAR's wild bootstrap path.
      real(dp), intent(in) :: history(:, :) !! History.
      real(dp), intent(in) :: residuals(:, :) !! Model residuals.
      type(fcvar_fit_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: periods !! Periods.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      real(dp), intent(in), optional :: level(:) !! Model level or confidence level.
      type(fcvar_path_t) :: out
      real(dp), allocatable :: signs(:)
      integer :: time

      if (periods < 1) then
         out%info = 1
         return
      end if
      allocate(signs(periods))
      do time = 1, periods
         if (random_uniform() > 0.5_dp) then
            signs(time) = 1.0_dp
         else
            signs(time) = -1.0_dp
         end if
      end do
      out = fcvar_bootstrap_from_signs(history, model, residuals, signs, &
         options, level)
   end function fcvar_bootstrap

   pure function fcvar_likelihood_ratio(unrestricted_log_likelihood, &
      restricted_log_likelihood, unrestricted_parameters, &
      restricted_parameters) result(out)
      !! Test ordinary nested FCVAR restrictions by a chi-square LR statistic.
      real(dp), intent(in) :: unrestricted_log_likelihood !! Unrestricted log likelihood.
      real(dp), intent(in) :: restricted_log_likelihood !! Restricted log likelihood.
      integer, intent(in) :: unrestricted_parameters !! Unrestricted parameters.
      integer, intent(in) :: restricted_parameters !! Restricted parameters.
      type(fcvar_lr_test_t) :: out

      out%unrestricted_log_likelihood = unrestricted_log_likelihood
      out%restricted_log_likelihood = restricted_log_likelihood
      if (.not. ieee_is_finite(unrestricted_log_likelihood) .or. &
         .not. ieee_is_finite(restricted_log_likelihood) .or. &
         unrestricted_log_likelihood < restricted_log_likelihood) then
         out%info = 1
         return
      end if
      out%degrees_of_freedom = unrestricted_parameters - restricted_parameters
      if (out%degrees_of_freedom <= 0) then
         out%info = 2
         return
      end if
      out%statistic = 2.0_dp*(unrestricted_log_likelihood - &
         restricted_log_likelihood)
      out%p_value = regularized_gamma_q( &
         0.5_dp*real(out%degrees_of_freedom, dp), 0.5_dp*out%statistic)
   end function fcvar_likelihood_ratio

   pure function fcvar_free_parameter_count(variables, lag_order, rank, &
      options, equal_orders) result(count)
      !! Count unrestricted FCVAR parameters after rank normalization.
      integer, intent(in) :: variables !! Number or indices of variables.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      integer :: count
      type(fcvar_options_t) :: selected
      logical :: equal

      selected = fcvar_options_t()
      if (present(options)) selected = options
      equal = .false.
      if (present(equal_orders)) equal = equal_orders
      if (variables < 1 .or. lag_order < 0 .or. rank < 0 .or. &
         rank > variables) then
         count = -1
         return
      end if
      count = merge(1, 2, equal) + 2*variables*rank + &
         variables*variables*lag_order - rank*rank
      if (selected%level_parameter) count = count + variables
      if (selected%restricted_constant) count = count + rank
      if (selected%unrestricted_constant) count = count + variables
   end function fcvar_free_parameter_count

   pure function fcvar_rank_tests(x, lag_order, bounds, options, equal_orders, &
      constrain_b, grid_points, max_iterations) result(out)
      !! Estimate every rank and compare each with the full-rank FCVAR model.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2, 2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      integer, intent(in), optional :: grid_points !! Grid points.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_rank_tests_t) :: out
      type(fcvar_estimation_t) :: estimation
      type(fcvar_options_t) :: selected
      type(fracdist_probability_t) :: rank_probability
      logical :: equal, constrained, probability_eligible
      real(dp) :: sample
      integer :: variables, rank, points, limit, constant_indicator

      selected = fcvar_options_t()
      if (present(options)) selected = options
      variables = size(x, 2)
      equal = .true.
      if (present(equal_orders)) equal = equal_orders
      constrained = .true.
      if (present(constrain_b)) constrained = constrain_b
      points = 11
      if (present(grid_points)) points = grid_points
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      if (size(x, 1) < 2 .or. variables < 1 .or. lag_order < 0 .or. &
         points < 2 .or. limit < 1) then
         out%info = 1
         return
      end if
      allocate(out%d(0:variables), out%b(0:variables))
      allocate(out%log_likelihood(0:variables), out%statistic(0:variables))
      allocate(out%p_value(0:variables), out%aic(0:variables), &
         out%bic(0:variables))
      allocate(out%free_parameters(0:variables))
      allocate(out%p_value_available(0:variables))
      out%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
      out%p_value_available = .false.
      do rank = 0, variables
         if (equal) then
            estimation = fcvar_estimate_equal_orders(x, lag_order, rank, &
               [max(bounds(1, 1), bounds(1, 2)), &
               min(bounds(2, 1), bounds(2, 2))], options, points)
         else
            estimation = fcvar_estimate(x, lag_order, rank, bounds(:, 1), &
               bounds(:, 2), options, constrained, points, limit)
         end if
         if (estimation%info /= 0) then
            out%info = 10 + rank
            return
         end if
         out%d(rank) = estimation%fit%d
         out%b(rank) = estimation%fit%b
         out%log_likelihood(rank) = estimation%fit%log_likelihood
         out%free_parameters(rank) = fcvar_free_parameter_count(variables, &
            lag_order, rank, options, equal)
      end do
      sample = real(size(x, 1), dp)
      if (present(options)) sample = sample - real(options%initial_values, dp)
      do rank = 0, variables
         out%statistic(rank) = 2.0_dp*(out%log_likelihood(variables) - &
            out%log_likelihood(rank))
         out%aic(rank) = -2.0_dp*out%log_likelihood(rank) + &
            2.0_dp*real(out%free_parameters(rank), dp)
         out%bic(rank) = -2.0_dp*out%log_likelihood(rank) + &
            log(sample)*real(out%free_parameters(rank), dp)
      end do
      probability_eligible = (.not. selected%restricted_constant .and. &
         .not. selected%unrestricted_constant .and. &
         .not. selected%level_parameter) .or. &
         (equal .and. .not. selected%unrestricted_constant .and. &
         (selected%restricted_constant .neqv. selected%level_parameter))
      constant_indicator = merge(1, 0, selected%restricted_constant .or. &
         selected%level_parameter)
      if (probability_eligible) then
         do rank = 0, variables - 1
            if (variables - rank > 12 .or. out%b(rank) <= 0.0_dp .or. &
               out%b(rank) >= 2.0_dp) cycle
            rank_probability = fracdist_p_value(variables - rank, &
               constant_indicator, out%b(rank), out%statistic(rank))
            if (rank_probability%info /= 0) cycle
            out%p_value(rank) = rank_probability%p_value
            out%p_value_available(rank) = .true.
         end do
      end if
      out%aic_rank = minloc(out%aic, dim=1) - 1
      out%bic_rank = minloc(out%bic, dim=1) - 1
      out%lag_order = lag_order
      out%observations = int(sample)
   end function fcvar_rank_tests

   pure function fcvar_lag_select(x, maximum_lag, rank, bounds, &
      diagnostic_lag, options, equal_orders, constrain_b, grid_points, &
      max_iterations) result(out)
      !! Select the FCVAR lag order by sequential LR tests, AIC, and BIC.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2, 2) !! Bounds.
      integer, intent(in) :: maximum_lag !! Maximum lag.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      integer, intent(in) :: diagnostic_lag !! Diagnostic lag.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      integer, intent(in), optional :: grid_points !! Grid points.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_lag_selection_t) :: out
      type(fcvar_estimation_t) :: estimation
      type(multivariate_white_noise_test_t) :: diagnostics
      type(fcvar_lr_test_t) :: test
      logical :: equal, constrained
      real(dp) :: sample
      integer :: variables, lag_order, points, limit

      variables = size(x, 2)
      equal = .true.
      if (present(equal_orders)) equal = equal_orders
      constrained = .true.
      if (present(constrain_b)) constrained = constrain_b
      points = 11
      if (present(grid_points)) points = grid_points
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      sample = real(size(x, 1), dp)
      if (present(options)) sample = sample - real(options%initial_values, dp)
      if (size(x, 1) < 2 .or. variables < 1 .or. maximum_lag < 0 .or. &
         rank < 0 .or. rank > variables .or. diagnostic_lag < 1 .or. &
         diagnostic_lag >= int(sample) .or. points < 2 .or. limit < 1) then
         out%info = 1
         return
      end if
      allocate(out%d(0:maximum_lag), out%b(0:maximum_lag))
      allocate(out%log_likelihood(0:maximum_lag))
      allocate(out%statistic(0:maximum_lag), out%p_value(0:maximum_lag))
      allocate(out%aic(0:maximum_lag), out%bic(0:maximum_lag))
      allocate(out%free_parameters(0:maximum_lag))
      allocate(out%q_p_value(0:maximum_lag, variables))
      allocate(out%lm_p_value(0:maximum_lag, variables))
      allocate(out%multivariate_q_p_value(0:maximum_lag))
      out%statistic = 0.0_dp
      out%p_value = 1.0_dp
      do lag_order = 0, maximum_lag
         if (equal) then
            estimation = fcvar_estimate_equal_orders(x, lag_order, rank, &
               [max(bounds(1, 1), bounds(1, 2)), &
               min(bounds(2, 1), bounds(2, 2))], options, points)
         else
            estimation = fcvar_estimate(x, lag_order, rank, bounds(:, 1), &
               bounds(:, 2), options, constrained, points, limit)
         end if
         if (estimation%info /= 0) then
            out%info = 10 + lag_order
            return
         end if
         out%d(lag_order) = estimation%fit%d
         out%b(lag_order) = estimation%fit%b
         out%log_likelihood(lag_order) = estimation%fit%log_likelihood
         out%free_parameters(lag_order) = fcvar_free_parameter_count( &
            variables, lag_order, rank, options, equal)
         out%aic(lag_order) = -2.0_dp*out%log_likelihood(lag_order) + &
            2.0_dp*real(out%free_parameters(lag_order), dp)
         out%bic(lag_order) = -2.0_dp*out%log_likelihood(lag_order) + &
            log(sample)*real(out%free_parameters(lag_order), dp)
         diagnostics = multivariate_white_noise_test( &
            estimation%fit%residuals, diagnostic_lag)
         if (diagnostics%info /= 0) then
            out%info = 100 + 10*lag_order + diagnostics%info
            return
         end if
         out%q_p_value(lag_order, :) = diagnostics%q_p_value
         out%lm_p_value(lag_order, :) = diagnostics%lm_p_value
         out%multivariate_q_p_value(lag_order) = &
            diagnostics%multivariate_q_p_value
         if (lag_order > 0) then
            test = fcvar_likelihood_ratio(out%log_likelihood(lag_order), &
               out%log_likelihood(lag_order - 1), &
               out%free_parameters(lag_order), &
               out%free_parameters(lag_order - 1))
            if (test%info /= 0) then
               out%info = 200 + lag_order
               return
            end if
            out%statistic(lag_order) = test%statistic
            out%p_value(lag_order) = test%p_value
         end if
      end do
      out%aic_lag = minloc(out%aic, dim=1) - 1
      out%bic_lag = minloc(out%bic, dim=1) - 1
      out%maximum_lag = maximum_lag
      out%rank = rank
      out%diagnostic_lag = diagnostic_lag
      out%observations = int(sample)
   end function fcvar_lag_select

   pure function fcvar_bootstrap_rank_from_signs(x, lag_order, null_rank, &
      alternative_rank, bounds, signs, options, equal_orders, constrain_b, &
      grid_points, max_iterations) result(out)
      !! Bootstrap a nested FCVAR rank LR statistic using supplied wild signs.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2, 2) !! Bounds.
      real(dp), intent(in) :: signs(:, :) !! Signs.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: null_rank !! Null rank.
      integer, intent(in) :: alternative_rank !! Alternative rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      integer, intent(in), optional :: grid_points !! Grid points.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_bootstrap_rank_t) :: out
      type(fcvar_options_t) :: selected
      type(fcvar_estimation_t) :: null_estimation, alternative_estimation
      type(fcvar_path_t) :: path
      real(dp), allocatable :: bootstrap_sample(:, :)
      logical :: equal, constrained
      real(dp) :: tolerance
      integer :: variables, periods, replications, replication, points, limit

      selected = fcvar_options_t()
      if (present(options)) selected = options
      variables = size(x, 2)
      periods = size(x, 1) - selected%initial_values
      replications = size(signs, 2)
      equal = .true.
      if (present(equal_orders)) equal = equal_orders
      constrained = .true.
      if (present(constrain_b)) constrained = constrain_b
      points = 11
      if (present(grid_points)) points = grid_points
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      out%null_rank = null_rank
      out%alternative_rank = alternative_rank
      out%replications = replications
      if (size(x, 1) < 2 .or. variables < 1 .or. lag_order < 0 .or. &
         lag_order + 1 > size(x, 1) .or. null_rank < 0 .or. &
         alternative_rank <= null_rank .or. alternative_rank > variables .or. &
         periods < 1 .or. replications < 1 .or. size(signs, 1) /= periods .or. &
         points < 2 .or. limit < 1 .or. .not. all(ieee_is_finite(x)) .or. &
         .not. all(ieee_is_finite(signs)) .or. &
         any(abs(abs(signs) - 1.0_dp) > 10.0_dp*epsilon(1.0_dp))) then
         out%info = 1
         return
      end if
      null_estimation = estimate_rank_model(x, lag_order, null_rank, bounds, &
         selected, equal, constrained, points, limit)
      if (null_estimation%info /= 0) then
         out%info = 10 + null_estimation%info
         return
      end if
      alternative_estimation = estimate_rank_model(x, lag_order, &
         alternative_rank, bounds, selected, equal, constrained, points, limit)
      if (alternative_estimation%info /= 0) then
         out%info = 20 + alternative_estimation%info
         return
      end if
      out%null_fit = null_estimation%fit
      out%alternative_fit = alternative_estimation%fit
      tolerance = 1.0e3_dp*epsilon(1.0_dp)*max(1.0_dp, &
         abs(out%null_fit%log_likelihood), &
         abs(out%alternative_fit%log_likelihood))
      if (out%alternative_fit%log_likelihood < &
         out%null_fit%log_likelihood - tolerance) then
         out%info = 2
         return
      end if
      out%observed_statistic = max(0.0_dp, 2.0_dp*( &
         out%alternative_fit%log_likelihood - out%null_fit%log_likelihood))
      allocate(out%statistic(replications))
      do replication = 1, replications
         path = fcvar_bootstrap_from_signs(x(:lag_order + 1, :), &
            out%null_fit, out%null_fit%residuals, signs(:, replication), selected)
         if (path%info /= 0) then
            out%info = 100 + replication
            return
         end if
         allocate(bootstrap_sample(lag_order + 1 + periods, variables))
         bootstrap_sample(:lag_order + 1, :) = x(:lag_order + 1, :)
         bootstrap_sample(lag_order + 2:, :) = path%series
         alternative_estimation = estimate_rank_model(bootstrap_sample, &
            lag_order, alternative_rank, bounds, selected, equal, constrained, &
            points, limit)
         if (alternative_estimation%info /= 0) then
            out%info = 200 + replication
            return
         end if
         null_estimation = estimate_rank_model(bootstrap_sample, lag_order, &
            null_rank, bounds, selected, equal, constrained, points, limit)
         if (null_estimation%info /= 0) then
            out%info = 300 + replication
            return
         end if
         tolerance = 1.0e3_dp*epsilon(1.0_dp)*max(1.0_dp, &
            abs(null_estimation%fit%log_likelihood), &
            abs(alternative_estimation%fit%log_likelihood))
         if (alternative_estimation%fit%log_likelihood < &
            null_estimation%fit%log_likelihood - tolerance) then
            out%info = 400 + replication
            return
         end if
         out%statistic(replication) = max(0.0_dp, 2.0_dp*( &
            alternative_estimation%fit%log_likelihood - &
            null_estimation%fit%log_likelihood))
         deallocate(bootstrap_sample)
      end do
      call sort(out%statistic)
      out%exceedances = count(out%statistic > out%observed_statistic)
      out%p_value = real(out%exceedances, dp)/real(replications, dp)
   end function fcvar_bootstrap_rank_from_signs

   function fcvar_bootstrap_rank(x, lag_order, null_rank, alternative_rank, &
      bounds, replications, options, equal_orders, constrain_b, grid_points, &
      max_iterations) result(out)
      !! Draw shared-RNG wild signs and bootstrap a nested FCVAR rank test.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2, 2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: null_rank !! Null rank.
      integer, intent(in) :: alternative_rank !! Alternative rank.
      integer, intent(in) :: replications !! Replications.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      integer, intent(in), optional :: grid_points !! Grid points.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_bootstrap_rank_t) :: out
      type(fcvar_options_t) :: selected
      real(dp), allocatable :: signs(:, :)
      integer :: periods, replication, time

      selected = fcvar_options_t()
      if (present(options)) selected = options
      periods = size(x, 1) - selected%initial_values
      if (periods < 1 .or. replications < 1) then
         out%info = 1
         return
      end if
      allocate(signs(periods, replications))
      do replication = 1, replications
         do time = 1, periods
            signs(time, replication) = merge(1.0_dp, -1.0_dp, &
               random_uniform() > 0.5_dp)
         end do
      end do
      out = fcvar_bootstrap_rank_from_signs(x, lag_order, null_rank, &
         alternative_rank, bounds, signs, selected, equal_orders, constrain_b, &
         grid_points, max_iterations)
   end function fcvar_bootstrap_rank

   pure function estimate_rank_model(x, lag_order, rank, bounds, options, &
      equal_orders, constrain_b, grid_points, max_iterations) result(out)
      !! Estimate one rank using the selected fractional-order restriction.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2, 2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      integer, intent(in) :: grid_points !! Grid points.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_options_t), intent(in) :: options !! Algorithm options.
      logical, intent(in) :: equal_orders !! Flag controlling equal orders.
      logical, intent(in) :: constrain_b !! Flag controlling constrain b.
      type(fcvar_estimation_t) :: out

      if (equal_orders) then
         out = fcvar_estimate_equal_orders(x, lag_order, rank, &
            [max(bounds(1, 1), bounds(1, 2)), &
            min(bounds(2, 1), bounds(2, 2))], options, grid_points)
      else
         out = fcvar_estimate(x, lag_order, rank, bounds(:, 1), bounds(:, 2), &
            options, constrain_b, grid_points, max_iterations)
      end if
   end function estimate_rank_model

   pure function fcvar_pack_parameters(fit, equal_orders) result(parameters)
      !! Pack identified FCVAR mean parameters into the Hessian vector order.
      type(fcvar_fit_t), intent(in) :: fit !! Previously fitted model.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      real(dp), allocatable :: parameters(:)
      logical :: equal
      integer :: count, index, block, variables

      equal = .true.
      if (present(equal_orders)) equal = equal_orders
      if (.not. valid_fit_layout(fit)) then
         allocate(parameters(0))
         return
      end if
      variables = size(fit%alpha, 1)
      count = inference_parameter_count(fit, equal)
      allocate(parameters(count))
      index = 1
      parameters(index) = fit%d
      index = index + 1
      if (.not. equal) then
         parameters(index) = fit%b
         index = index + 1
      end if
      block = size(fit%level)
      if (block > 0) then
         parameters(index:index + block - 1) = fit%level
         index = index + block
      end if
      block = size(fit%alpha)
      if (block > 0) then
         parameters(index:index + block - 1) = reshape(fit%alpha, [block])
         index = index + block
      end if
      block = max(0, variables - fit%rank)*fit%rank
      if (block > 0) then
         parameters(index:index + block - 1) = &
            reshape(fit%beta(fit%rank + 1:, :), [block])
         index = index + block
      end if
      block = size(fit%rho)
      if (block > 0) then
         parameters(index:index + block - 1) = fit%rho
         index = index + block
      end if
      block = size(fit%gamma)
      if (block > 0) then
         parameters(index:index + block - 1) = reshape(fit%gamma, [block])
         index = index + block
      end if
      block = size(fit%unrestricted_constant)
      if (block > 0) parameters(index:index + block - 1) = &
         fit%unrestricted_constant
   end function fcvar_pack_parameters

   pure function fcvar_unpack_parameters(parameters, template, equal_orders) &
      result(out)
      !! Restore identified FCVAR coefficient fields from a packed vector.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      type(fcvar_fit_t), intent(in) :: template !! Template.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      type(fcvar_fit_t) :: out
      logical :: equal
      integer :: count, index, block, variables

      equal = .true.
      if (present(equal_orders)) equal = equal_orders
      out = template
      out%info = 0
      if (.not. valid_fit_layout(template) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      count = inference_parameter_count(template, equal)
      if (size(parameters) /= count) then
         out%info = 1
         return
      end if
      variables = size(out%alpha, 1)
      index = 1
      out%d = parameters(index)
      index = index + 1
      if (equal) then
         out%b = out%d
      else
         out%b = parameters(index)
         index = index + 1
      end if
      block = size(out%level)
      if (block > 0) then
         out%level = parameters(index:index + block - 1)
         index = index + block
      end if
      block = size(out%alpha)
      if (block > 0) then
         out%alpha = reshape(parameters(index:index + block - 1), &
            shape(out%alpha))
         index = index + block
      end if
      block = max(0, variables - out%rank)*out%rank
      if (block > 0) then
         out%beta(out%rank + 1:, :) = reshape( &
            parameters(index:index + block - 1), &
            [variables - out%rank, out%rank])
         index = index + block
      end if
      block = size(out%rho)
      if (block > 0) then
         out%rho = parameters(index:index + block - 1)
         index = index + block
      end if
      block = size(out%gamma)
      if (block > 0) then
         out%gamma = reshape(parameters(index:index + block - 1), &
            shape(out%gamma))
         index = index + block
      end if
      block = size(out%unrestricted_constant)
      if (block > 0) out%unrestricted_constant = &
         parameters(index:index + block - 1)
      if (out%rank > 0) then
         out%pi = matmul(out%alpha, transpose(out%beta))
      else
         out%pi = 0.0_dp
      end if
   end function fcvar_unpack_parameters

   pure function fcvar_full_log_likelihood(x, fit, options) result(value)
      !! Evaluate the covariance-concentrated likelihood at all mean parameters.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      type(fcvar_fit_t), intent(in) :: fit !! Previously fitted model.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      real(dp) :: value
      type(fcvar_options_t) :: selected
      type(fcvar_transform_t) :: transformed
      real(dp), allocatable :: beta_star(:, :), residuals(:, :), covariance(:, :)
      real(dp), allocatable :: shifted(:, :)
      real(dp) :: log_determinant, sample, pi_value
      integer :: variables, status

      selected = fcvar_options_t()
      if (present(options)) selected = options
      value = -huge(1.0_dp)
      if (.not. valid_fit_layout(fit) .or. .not. all(ieee_is_finite(x))) return
      variables = size(x, 2)
      if (variables /= size(fit%alpha, 1) .or. &
         selected%level_parameter .neqv. (size(fit%level) > 0) .or. &
         selected%restricted_constant .neqv. (size(fit%rho) > 0) .or. &
         selected%unrestricted_constant .neqv. &
         (size(fit%unrestricted_constant) > 0)) return
      shifted = x
      if (selected%level_parameter) shifted = shifted - &
         spread(fit%level, 1, size(x, 1))
      transformed = fcvar_transform(shifted, fit%lag_order, fit%d, fit%b, &
         selected)
      if (transformed%info /= 0) return
      residuals = transformed%z0
      if (fit%rank > 0) then
         allocate(beta_star(variables + merge(1, 0, &
            selected%restricted_constant), fit%rank))
         beta_star(:variables, :) = fit%beta
         if (selected%restricted_constant) beta_star(variables + 1, :) = fit%rho
         residuals = residuals - matmul(transformed%z1, &
            matmul(beta_star, transpose(fit%alpha)))
      end if
      if (fit%lag_order > 0) residuals = residuals - &
         matmul(transformed%z2, transpose(fit%gamma))
      if (selected%unrestricted_constant) residuals = residuals - &
         spread(fit%unrestricted_constant, 1, size(residuals, 1))
      sample = real(size(residuals, 1), dp)
      covariance = matmul(transpose(residuals), residuals)/sample
      call covariance_log_determinant(covariance, log_determinant, status)
      if (status /= 0) return
      pi_value = acos(-1.0_dp)
      value = -0.5_dp*sample*real(variables, dp)* &
         (log(2.0_dp*pi_value) + 1.0_dp) - 0.5_dp*sample*log_determinant
   end function fcvar_full_log_likelihood

   pure function refresh_fcvar_fit(x, coefficients, options) result(out)
      !! Recompute residuals, covariance, and likelihood for supplied coefficients.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      type(fcvar_fit_t), intent(in) :: coefficients !! Model coefficients.
      type(fcvar_options_t), intent(in) :: options !! Algorithm options.
      type(fcvar_fit_t) :: out
      type(fcvar_transform_t) :: transformed
      real(dp), allocatable :: beta_star(:, :)
      real(dp) :: log_determinant, sample, pi_value
      integer :: variables, status

      out = coefficients
      out%info = 0
      if (.not. valid_fit_layout(coefficients) .or. &
         .not. all(ieee_is_finite(x))) then
         out%info = 1
         return
      end if
      variables = size(x, 2)
      if (variables /= size(out%alpha, 1) .or. &
         options%level_parameter .neqv. (size(out%level) > 0) .or. &
         options%restricted_constant .neqv. (size(out%rho) > 0) .or. &
         options%unrestricted_constant .neqv. &
         (size(out%unrestricted_constant) > 0)) then
         out%info = 1
         return
      end if
      if (options%level_parameter) then
         transformed = fcvar_transform(x - spread(out%level, 1, size(x, 1)), &
            out%lag_order, out%d, out%b, options)
      else
         transformed = fcvar_transform(x, out%lag_order, out%d, out%b, options)
      end if
      if (transformed%info /= 0) then
         out%info = 10 + transformed%info
         return
      end if
      out%residuals = transformed%z0
      if (out%rank > 0) then
         allocate(beta_star(variables + merge(1, 0, &
            options%restricted_constant), out%rank))
         beta_star(:variables, :) = out%beta
         if (options%restricted_constant) &
            beta_star(variables + 1, :) = out%rho
         out%residuals = out%residuals - matmul(transformed%z1, &
            matmul(beta_star, transpose(out%alpha)))
      end if
      if (out%lag_order > 0) out%residuals = out%residuals - &
         matmul(transformed%z2, transpose(out%gamma))
      if (options%unrestricted_constant) out%residuals = out%residuals - &
         spread(out%unrestricted_constant, 1, size(out%residuals, 1))
      sample = real(size(out%residuals, 1), dp)
      out%covariance = matmul(transpose(out%residuals), out%residuals)/sample
      call covariance_log_determinant(out%covariance, log_determinant, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      pi_value = acos(-1.0_dp)
      out%log_likelihood = -0.5_dp*sample*real(variables, dp)* &
         (log(2.0_dp*pi_value) + 1.0_dp) - 0.5_dp*sample*log_determinant
      out%observations = size(out%residuals, 1)
   end function refresh_fcvar_fit

   pure function fcvar_standard_errors(x, fit, options, equal_orders, &
      increment) result(out)
      !! Estimate coefficient covariance from a centered numerical Hessian.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      type(fcvar_fit_t), intent(in) :: fit !! Previously fitted model.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      real(dp), intent(in), optional :: increment !! Increment.
      type(fcvar_standard_errors_t) :: out
      real(dp), allocatable :: step(:), first(:), second(:), inverse(:, :)
      real(dp) :: center, first_value, second_value, third_value, fourth_value
      logical :: equal
      integer :: count, i, j, status

      equal = .true.
      if (present(equal_orders)) equal = equal_orders
      out%equal_orders = equal
      out%increment = 1.0e-4_dp
      if (present(increment)) out%increment = increment
      out%parameters = fcvar_pack_parameters(fit, equal)
      count = size(out%parameters)
      if (count < 1 .or. out%increment <= 0.0_dp .or. &
         .not. ieee_is_finite(out%increment) .or. &
         size(x, 2) /= size(fit%alpha, 1)) then
         out%info = 1
         return
      end if
      center = objective(out%parameters)
      if (center <= -0.5_dp*huge(1.0_dp)) then
         out%info = 2
         return
      end if
      allocate(out%hessian(count, count), step(count), first(count), second(count))
      step = out%increment*max(1.0_dp, abs(out%parameters))
      out%hessian = 0.0_dp
      do i = 1, count
         first = out%parameters
         second = out%parameters
         first(i) = first(i) + step(i)
         second(i) = second(i) - step(i)
         first_value = objective(first)
         second_value = objective(second)
         if (min(first_value, second_value) <= -0.5_dp*huge(1.0_dp)) then
            out%info = 2
            return
         end if
         out%hessian(i, i) = (first_value - 2.0_dp*center + second_value)/ &
            step(i)**2
         do j = 1, i - 1
            first = out%parameters
            second = out%parameters
            first(i) = first(i) + step(i)
            first(j) = first(j) + step(j)
            second(i) = second(i) - step(i)
            second(j) = second(j) + step(j)
            first_value = objective(first)
            second_value = objective(second)
            first = out%parameters
            second = out%parameters
            first(i) = first(i) + step(i)
            first(j) = first(j) - step(j)
            second(i) = second(i) - step(i)
            second(j) = second(j) - step(j)
            third_value = objective(first)
            fourth_value = objective(second)
            if (min(first_value, second_value, third_value, fourth_value) <= &
               -0.5_dp*huge(1.0_dp)) then
               out%info = 2
               return
            end if
            out%hessian(i, j) = (first_value - second_value - third_value + &
               fourth_value)/(4.0_dp*step(i)*step(j))
            out%hessian(j, i) = out%hessian(i, j)
         end do
      end do
      call invert_matrix(symmetrize(-out%hessian), inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      out%covariance = symmetrize(inverse)
      allocate(out%standard_error(count))
      do i = 1, count
         if (out%covariance(i, i) > 0.0_dp .and. &
            ieee_is_finite(out%covariance(i, i))) then
            out%standard_error(i) = sqrt(out%covariance(i, i))
         else
            out%standard_error(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            out%info = 4
         end if
      end do
      call map_standard_errors(out, fit)

   contains

      pure function objective(parameters) result(log_likelihood)
         !! Evaluate the full likelihood at one packed parameter vector.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: log_likelihood
         type(fcvar_fit_t) :: adjusted

         adjusted = fcvar_unpack_parameters(parameters, fit, equal)
         if (adjusted%info /= 0) then
            log_likelihood = -huge(1.0_dp)
         else
            log_likelihood = fcvar_full_log_likelihood(x, adjusted, options)
         end if
      end function objective

   end function fcvar_standard_errors

   pure function fcvar_characteristic_roots(fit) result(out)
      !! Compute inverse companion eigenvalues as FCVAR characteristic roots.
      type(fcvar_fit_t), intent(in) :: fit !! Previously fitted model.
      type(fcvar_roots_t) :: out
      real(dp), allocatable :: companion(:, :)
      integer :: variables, states, lag, row, status
      real(dp) :: infinity

      out%b = fit%b
      if (.not. valid_fit_layout(fit)) then
         out%info = 1
         return
      end if
      variables = size(fit%alpha, 1)
      states = variables*(fit%lag_order + 1)
      allocate(companion(states, states))
      companion = 0.0_dp
      companion(:variables, :variables) = identity_matrix(variables) + fit%pi
      if (fit%lag_order > 0) then
         companion(:variables, :variables) = &
            companion(:variables, :variables) + fit%gamma(:, :variables)
         do lag = 2, fit%lag_order
            companion(:variables, (lag - 1)*variables + 1:lag*variables) = &
               fit%gamma(:, (lag - 1)*variables + 1:lag*variables) - &
               fit%gamma(:, (lag - 2)*variables + 1:(lag - 1)*variables)
         end do
         companion(:variables, fit%lag_order*variables + 1:states) = &
            -fit%gamma(:, (fit%lag_order - 1)*variables + 1: &
            fit%lag_order*variables)
         do row = 1, variables*fit%lag_order
            companion(variables + row, row) = 1.0_dp
         end do
      end if
      call general_eigenvalues(companion, out%eigenvalues, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      allocate(out%roots(states))
      infinity = ieee_value(0.0_dp, ieee_positive_inf)
      do row = 1, states
         if (abs(out%eigenvalues(row)) <= 100.0_dp*epsilon(1.0_dp)) then
            out%roots(row) = cmplx(infinity, 0.0_dp, dp)
         else
            out%roots(row) = 1.0_dp/out%eigenvalues(row)
         end if
      end do
      call sort_complex_modulus(out%roots)
      out%minimum_root_modulus = minval(abs(out%roots))
      out%maximum_inverse_root_modulus = maxval(abs(out%eigenvalues))
      out%outside_unit_circle = out%minimum_root_modulus > 1.0_dp
   end function fcvar_characteristic_roots

   pure logical function valid_fit_layout(fit) result(valid)
      !! Check allocation and dimensions required by FCVAR postestimation.
      type(fcvar_fit_t), intent(in) :: fit !! Previously fitted model.
      integer :: variables

      valid = allocated(fit%alpha) .and. allocated(fit%beta) .and. &
         allocated(fit%pi) .and. allocated(fit%rho) .and. &
         allocated(fit%gamma) .and. allocated(fit%level) .and. &
         allocated(fit%unrestricted_constant)
      if (.not. valid) return
      variables = size(fit%alpha, 1)
      valid = variables > 0 .and. fit%lag_order >= 0 .and. &
         fit%rank >= 0 .and. fit%rank <= variables .and. &
         all(shape(fit%alpha) == [variables, fit%rank]) .and. &
         all(shape(fit%beta) == [variables, fit%rank]) .and. &
         all(shape(fit%pi) == [variables, variables]) .and. &
         size(fit%gamma, 1) == variables .and. &
         size(fit%gamma, 2) == variables*fit%lag_order .and. &
         (size(fit%level) == 0 .or. size(fit%level) == variables) .and. &
         (size(fit%rho) == 0 .or. size(fit%rho) == fit%rank) .and. &
         (size(fit%unrestricted_constant) == 0 .or. &
         size(fit%unrestricted_constant) == variables)
      if (.not. valid) return
      valid = ieee_is_finite(fit%d) .and. ieee_is_finite(fit%b) .and. &
         all(ieee_is_finite(fit%alpha)) .and. &
         all(ieee_is_finite(fit%beta)) .and. all(ieee_is_finite(fit%pi)) .and. &
         all(ieee_is_finite(fit%level)) .and. &
         all(ieee_is_finite(fit%rho)) .and. &
         all(ieee_is_finite(fit%gamma)) .and. &
         all(ieee_is_finite(fit%unrestricted_constant))
   end function valid_fit_layout

   pure integer function inference_parameter_count(fit, equal_orders) &
      result(count)
      !! Count identified mean parameters represented in the Hessian vector.
      type(fcvar_fit_t), intent(in) :: fit !! Previously fitted model.
      logical, intent(in) :: equal_orders !! Flag controlling equal orders.
      integer :: variables

      variables = size(fit%alpha, 1)
      count = merge(1, 2, equal_orders) + size(fit%level) + size(fit%alpha) + &
         max(0, variables - fit%rank)*fit%rank + size(fit%rho) + &
         size(fit%gamma) + size(fit%unrestricted_constant)
   end function inference_parameter_count

   pure subroutine map_standard_errors(out, fit)
      !! Map packed standard errors back to coefficient array shapes.
      type(fcvar_standard_errors_t), intent(inout) :: out !! Procedure result, updated in place.
      type(fcvar_fit_t), intent(in) :: fit !! Previously fitted model.
      integer :: variables, index, block

      variables = size(fit%alpha, 1)
      allocate(out%alpha, mold=fit%alpha)
      allocate(out%beta, mold=fit%beta)
      allocate(out%rho, mold=fit%rho)
      allocate(out%gamma, mold=fit%gamma)
      allocate(out%level, mold=fit%level)
      allocate(out%unrestricted_constant, mold=fit%unrestricted_constant)
      out%alpha = 0.0_dp
      out%beta = 0.0_dp
      out%rho = 0.0_dp
      out%gamma = 0.0_dp
      out%level = 0.0_dp
      out%unrestricted_constant = 0.0_dp
      index = 1
      out%d = out%standard_error(index)
      index = index + 1
      if (out%equal_orders) then
         out%b = out%d
      else
         out%b = out%standard_error(index)
         index = index + 1
      end if
      block = size(out%level)
      if (block > 0) then
         out%level = out%standard_error(index:index + block - 1)
         index = index + block
      end if
      block = size(out%alpha)
      if (block > 0) then
         out%alpha = reshape(out%standard_error(index:index + block - 1), &
            shape(out%alpha))
         index = index + block
      end if
      block = max(0, variables - fit%rank)*fit%rank
      if (block > 0) then
         out%beta(fit%rank + 1:, :) = reshape( &
            out%standard_error(index:index + block - 1), &
            [variables - fit%rank, fit%rank])
         index = index + block
      end if
      block = size(out%rho)
      if (block > 0) then
         out%rho = out%standard_error(index:index + block - 1)
         index = index + block
      end if
      block = size(out%gamma)
      if (block > 0) then
         out%gamma = reshape(out%standard_error(index:index + block - 1), &
            shape(out%gamma))
         index = index + block
      end if
      block = size(out%unrestricted_constant)
      if (block > 0) out%unrestricted_constant = &
         out%standard_error(index:index + block - 1)
   end subroutine map_standard_errors

   pure subroutine sort_complex_modulus(values)
      !! Sort complex values by decreasing modulus.
      complex(dp), intent(inout) :: values(:) !! Input values, updated in place.
      complex(dp) :: selected
      integer :: i, j

      do i = 2, size(values)
         selected = values(i)
         j = i - 1
         do while (j >= 1)
            if (abs(values(j)) >= abs(selected)) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = selected
      end do
   end subroutine sort_complex_modulus

   pure function fcvar_restricted_estimate(x, lag_order, rank, bounds, &
      restrictions, options, equal_orders, constrain_b, grid_points, &
      max_iterations) result(out)
      !! Estimate an FCVAR model subject to exact linear coefficient restrictions.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2, 2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_restrictions_t), intent(in) :: restrictions !! Restrictions.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      integer, intent(in), optional :: grid_points !! Grid points.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_restricted_estimation_t) :: out
      type(fcvar_options_t) :: selected
      type(fcvar_estimation_t) :: unrestricted
      type(optimization_result_t) :: optimized
      type(fcvar_fit_t) :: coefficients
      type(fcvar_switching_estimation_t) :: switching
      real(dp), allocatable :: parameters(:), particular(:), basis(:, :)
      real(dp), allocatable :: coordinates(:)
      logical :: equal, constrained, switching_applicable
      integer :: points, limit, status, parameter_count

      selected = fcvar_options_t()
      if (present(options)) selected = options
      equal = .true.
      if (present(equal_orders)) equal = equal_orders
      constrained = .true.
      if (present(constrain_b)) constrained = constrain_b
      points = 11
      if (present(grid_points)) points = grid_points
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      unrestricted = estimate_rank_model(x, lag_order, rank, bounds, selected, &
         equal, constrained, points, limit)
      if (unrestricted%info /= 0) then
         out%info = 10 + unrestricted%info
         return
      end if
      out%unrestricted_fit = unrestricted%fit
      parameters = fcvar_pack_parameters(out%unrestricted_fit, equal)
      parameter_count = size(parameters)
      call build_constraint_system(out%unrestricted_fit, equal, restrictions, &
         out%constraint_matrix, out%constraint_value, status)
      if (status /= 0 .or. size(out%constraint_matrix, 1) < 1) then
         out%info = 20 + max(1, status)
         return
      end if
      call constraint_geometry(out%constraint_matrix, out%constraint_value, &
         parameters, particular, basis, out%restriction_rank, status)
      if (status /= 0 .or. out%restriction_rank < 1) then
         out%info = 30 + max(1, status)
         return
      end if
      if (size(basis, 2) == 0) then
         coefficients = fcvar_unpack_parameters(particular, &
            out%unrestricted_fit, equal)
         out%converged = .true.
      else
         coordinates = matmul(transpose(basis), parameters - particular)
         optimized = bfgs_minimize_fd(objective, coordinates, limit, &
            1.0e-6_dp, 1.0e-5_dp)
         out%iterations = optimized%iterations
         out%optimizer_info = optimized%info
         out%converged = optimized%converged
         if (.not. allocated(optimized%parameters) .or. &
            .not. ieee_is_finite(optimized%objective) .or. &
            optimized%objective >= 0.2_dp*huge(1.0_dp)) then
            out%info = 40 + max(1, optimized%info)
            return
         end if
         coefficients = fcvar_unpack_parameters(particular + &
            matmul(basis, optimized%parameters), out%unrestricted_fit, equal)
      end if
      if (coefficients%info /= 0) then
         out%info = 50 + coefficients%info
         return
      end if
      switching_applicable = rank > 0 .and. .not. selected%level_parameter .and. &
         .not. allocated(restrictions%parameter_matrix) .and. &
         .not. allocated(restrictions%order_matrix) .and. &
         .not. allocated(restrictions%level_matrix) .and. &
         (allocated(restrictions%alpha_matrix) .or. &
         allocated(restrictions%beta_matrix))
      if (switching_applicable) then
         switching = fcvar_estimate_restricted_fixed(x, lag_order, rank, &
            coefficients%d, coefficients%b, restrictions, selected, limit, &
            1.0e-8_dp)
         if (switching%info /= 0) then
            out%info = 55 + switching%info
            return
         end if
         out%restricted_fit = switching%fit
         out%iterations = out%iterations + switching%iterations
         out%converged = out%converged .and. switching%converged
      else
         out%restricted_fit = refresh_fcvar_fit(x, coefficients, selected)
      end if
      if (out%restricted_fit%info /= 0) then
         out%info = 60 + out%restricted_fit%info
         return
      end if
      out%test = fcvar_likelihood_ratio( &
         out%unrestricted_fit%log_likelihood, &
         out%restricted_fit%log_likelihood, parameter_count, &
         parameter_count - out%restriction_rank)
      if (out%test%info /= 0) out%info = 70 + out%test%info

   contains

      pure function objective(coordinates_value) result(value)
         !! Return negative likelihood in exact restriction-null-space coordinates.
         real(dp), intent(in) :: coordinates_value(:) !! Coordinates value.
         real(dp) :: value, likelihood
         type(fcvar_fit_t) :: adjusted

         adjusted = fcvar_unpack_parameters(particular + &
            matmul(basis, coordinates_value), out%unrestricted_fit, equal)
         if (adjusted%info /= 0 .or. adjusted%d < bounds(1, 1) .or. &
            adjusted%d > bounds(2, 1) .or. adjusted%b < bounds(1, 2) .or. &
            adjusted%b > bounds(2, 2) .or. &
            (constrained .and. adjusted%b > adjusted%d)) then
            value = 0.25_dp*huge(1.0_dp)
            return
         end if
         likelihood = fcvar_full_log_likelihood(x, adjusted, selected)
         if (likelihood <= -0.5_dp*huge(1.0_dp)) then
            value = 0.25_dp*huge(1.0_dp)
         else
            value = -likelihood
         end if
      end function objective

   end function fcvar_restricted_estimate

   pure function fcvar_bootstrap_hypothesis_from_signs(x, lag_order, rank, &
      bounds, restrictions, signs, options, equal_orders, constrain_b, &
      grid_points, max_iterations) result(out)
      !! Bootstrap a linear-restriction LR test using supplied wild signs.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2, 2) !! Bounds.
      real(dp), intent(in) :: signs(:, :) !! Signs.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_restrictions_t), intent(in) :: restrictions !! Restrictions.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      integer, intent(in), optional :: grid_points !! Grid points.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_bootstrap_hypothesis_t) :: out
      type(fcvar_options_t) :: selected
      type(fcvar_restricted_estimation_t) :: estimation
      type(fcvar_path_t) :: path
      real(dp), allocatable :: bootstrap_sample(:, :)
      logical :: equal, constrained
      integer :: periods, replications, replication, points, limit, variables

      selected = fcvar_options_t()
      if (present(options)) selected = options
      equal = .true.
      if (present(equal_orders)) equal = equal_orders
      constrained = .true.
      if (present(constrain_b)) constrained = constrain_b
      points = 11
      if (present(grid_points)) points = grid_points
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      periods = size(x, 1) - selected%initial_values
      replications = size(signs, 2)
      variables = size(x, 2)
      out%replications = replications
      if (periods < 1 .or. replications < 1 .or. &
         size(signs, 1) /= periods .or. lag_order < 0 .or. &
         lag_order + 1 > size(x, 1) .or. variables < 1 .or. &
         .not. all(ieee_is_finite(signs)) .or. &
         any(abs(abs(signs) - 1.0_dp) > 10.0_dp*epsilon(1.0_dp))) then
         out%info = 1
         return
      end if
      estimation = fcvar_restricted_estimate(x, lag_order, rank, bounds, &
         restrictions, selected, equal, constrained, points, limit)
      if (estimation%info /= 0) then
         out%info = 10 + estimation%info
         return
      end if
      out%unrestricted_fit = estimation%unrestricted_fit
      out%restricted_fit = estimation%restricted_fit
      out%observed_statistic = estimation%test%statistic
      out%restriction_rank = estimation%restriction_rank
      allocate(out%statistic(replications))
      do replication = 1, replications
         path = fcvar_bootstrap_from_signs(x(:lag_order + 1, :), &
            out%restricted_fit, out%restricted_fit%residuals, &
            signs(:, replication), selected)
         if (path%info /= 0) then
            out%info = 100 + replication
            return
         end if
         allocate(bootstrap_sample(lag_order + 1 + periods, variables))
         bootstrap_sample(:lag_order + 1, :) = x(:lag_order + 1, :)
         bootstrap_sample(lag_order + 2:, :) = path%series
         estimation = fcvar_restricted_estimate(bootstrap_sample, lag_order, &
            rank, bounds, restrictions, selected, equal, constrained, points, &
            limit)
         if (estimation%info /= 0) then
            out%info = 200 + replication
            return
         end if
         out%statistic(replication) = estimation%test%statistic
         deallocate(bootstrap_sample)
      end do
      call sort(out%statistic)
      out%exceedances = count(out%statistic > out%observed_statistic)
      out%p_value = real(out%exceedances, dp)/real(replications, dp)
   end function fcvar_bootstrap_hypothesis_from_signs

   function fcvar_bootstrap_hypothesis(x, lag_order, rank, bounds, &
      restrictions, replications, options, equal_orders, constrain_b, &
      grid_points, max_iterations) result(out)
      !! Draw shared-RNG wild signs and bootstrap a linear-restriction LR test.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2, 2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      integer, intent(in) :: replications !! Replications.
      type(fcvar_restrictions_t), intent(in) :: restrictions !! Restrictions.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in), optional :: equal_orders !! Flag controlling equal orders.
      logical, intent(in), optional :: constrain_b !! Flag controlling constrain b.
      integer, intent(in), optional :: grid_points !! Grid points.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fcvar_bootstrap_hypothesis_t) :: out
      type(fcvar_options_t) :: selected
      real(dp), allocatable :: signs(:, :)
      integer :: periods, replication, time

      selected = fcvar_options_t()
      if (present(options)) selected = options
      periods = size(x, 1) - selected%initial_values
      if (periods < 1 .or. replications < 1) then
         out%info = 1
         return
      end if
      allocate(signs(periods, replications))
      do replication = 1, replications
         do time = 1, periods
            signs(time, replication) = merge(1.0_dp, -1.0_dp, &
               random_uniform() > 0.5_dp)
         end do
      end do
      out = fcvar_bootstrap_hypothesis_from_signs(x, lag_order, rank, bounds, &
         restrictions, signs, selected, equal_orders, constrain_b, grid_points, &
         max_iterations)
   end function fcvar_bootstrap_hypothesis

   pure function switching_covariance(s00, s01, s11, alpha, beta) result(omega)
      !! Evaluate the concentrated covariance for alpha and beta candidates.
      real(dp), intent(in) :: s00(:, :) !! S00.
      real(dp), intent(in) :: s01(:, :) !! S01.
      real(dp), intent(in) :: s11(:, :) !! S11.
      real(dp), intent(in) :: alpha(:, :) !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: beta(:, :) !! Regression or model coefficients.
      real(dp), allocatable :: omega(:, :)
      real(dp), allocatable :: cross(:, :)

      cross = matmul(alpha, transpose(beta))
      omega = s00 - matmul(s01, transpose(cross)) - &
         matmul(cross, transpose(s01)) + &
         matmul(cross, matmul(s11, transpose(cross)))
      omega = symmetrize(omega)
   end function switching_covariance

   pure subroutine constrained_quadratic_projection(unrestricted, weight, &
      matrix, value, adjusted, info)
      !! Project a quadratic optimum onto affine linear restrictions.
      real(dp), intent(in) :: unrestricted(:) !! Unrestricted.
      real(dp), intent(in) :: weight(:, :) !! Weight.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), intent(in) :: value(:) !! Input value.
      real(dp), allocatable, intent(out) :: adjusted(:) !! Adjusted.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: inverse(:, :), gram(:, :), gram_inverse(:, :)
      real(dp) :: tolerance
      integer :: parameters, restrictions_count

      parameters = size(unrestricted)
      restrictions_count = size(matrix, 1)
      info = 0
      if (any(shape(weight) /= [parameters, parameters]) .or. &
         size(matrix, 2) /= parameters .or. &
         size(value) /= restrictions_count) then
         allocate(adjusted(0))
         info = 1
         return
      end if
      if (restrictions_count == 0) then
         adjusted = unrestricted
         return
      end if
      call invert_matrix(symmetrize(weight), inverse, info)
      if (info /= 0) return
      gram = matmul(matrix, matmul(inverse, transpose(matrix)))
      allocate(gram_inverse(restrictions_count, restrictions_count))
      tolerance = 1.0e-11_dp*max(1.0_dp, maxval(abs(gram)))
      call symmetric_pseudoinverse(gram, gram_inverse, info, tolerance)
      if (info /= 0) return
      adjusted = unrestricted - matmul(inverse, matmul(transpose(matrix), &
         matmul(gram_inverse, matmul(matrix, unrestricted) - value)))
      tolerance = 1.0e-8_dp*max(1.0_dp, maxval(abs(value)), &
         maxval(abs(matmul(matrix, adjusted))))
      if (maxval(abs(matmul(matrix, adjusted) - value)) > tolerance) info = 2
   end subroutine constrained_quadratic_projection

   pure subroutine build_constraint_system(fit, equal_orders, restrictions, &
      matrix, value, info)
      !! Map convenient FCVAR restrictions into the packed parameter vector.
      type(fcvar_fit_t), intent(in) :: fit !! Previously fitted model.
      logical, intent(in) :: equal_orders !! Flag controlling equal orders.
      type(fcvar_restrictions_t), intent(in) :: restrictions !! Restrictions.
      real(dp), allocatable, intent(out) :: matrix(:, :) !! Input matrix.
      real(dp), allocatable, intent(out) :: value(:) !! Input value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      integer :: rows, parameter_count, order_count, alpha_count, beta_rows
      integer :: target, source, restriction, variable, relation, packed
      integer :: level_start, alpha_start, beta_start, rho_start, variables
      real(dp) :: coefficient, fixed_value

      info = 0
      rows = 0
      parameter_count = inference_parameter_count(fit, equal_orders)
      variables = size(fit%alpha, 1)
      order_count = merge(1, 2, equal_orders)
      alpha_count = size(fit%alpha)
      beta_rows = variables + merge(1, 0, size(fit%rho) > 0)
      if (allocated(restrictions%parameter_matrix) .neqv. &
         allocated(restrictions%parameter_value)) info = 1
      if (allocated(restrictions%order_matrix) .neqv. &
         allocated(restrictions%order_value)) info = 1
      if (allocated(restrictions%alpha_matrix) .neqv. &
         allocated(restrictions%alpha_value)) info = 1
      if (allocated(restrictions%beta_matrix) .neqv. &
         allocated(restrictions%beta_value)) info = 1
      if (allocated(restrictions%level_matrix) .neqv. &
         allocated(restrictions%level_value)) info = 1
      if (allocated(restrictions%order_inequality_matrix) .or. &
         allocated(restrictions%order_inequality_value)) info = 4
      if (info /= 0) then
         allocate(matrix(0, parameter_count), value(0))
         return
      end if
      if (allocated(restrictions%parameter_matrix)) then
         if (size(restrictions%parameter_matrix, 2) /= parameter_count .or. &
            size(restrictions%parameter_matrix, 1) /= &
            size(restrictions%parameter_value)) info = 2
         rows = rows + size(restrictions%parameter_matrix, 1)
      end if
      if (allocated(restrictions%order_matrix)) then
         if (size(restrictions%order_matrix, 2) /= 2 .or. &
            size(restrictions%order_matrix, 1) /= &
            size(restrictions%order_value)) info = 2
         rows = rows + size(restrictions%order_matrix, 1)
      end if
      if (allocated(restrictions%level_matrix)) then
         if (size(restrictions%level_matrix, 2) /= size(fit%level) .or. &
            size(restrictions%level_matrix, 1) /= &
            size(restrictions%level_value)) info = 2
         rows = rows + size(restrictions%level_matrix, 1)
      end if
      if (allocated(restrictions%alpha_matrix)) then
         if (size(restrictions%alpha_matrix, 2) /= alpha_count .or. &
            size(restrictions%alpha_matrix, 1) /= &
            size(restrictions%alpha_value)) info = 2
         rows = rows + size(restrictions%alpha_matrix, 1)
      end if
      if (allocated(restrictions%beta_matrix)) then
         if (size(restrictions%beta_matrix, 2) /= beta_rows*fit%rank .or. &
            size(restrictions%beta_matrix, 1) /= &
            size(restrictions%beta_value)) info = 2
         rows = rows + size(restrictions%beta_matrix, 1)
      end if
      allocate(matrix(rows, parameter_count), value(rows))
      matrix = 0.0_dp
      value = 0.0_dp
      if (info /= 0 .or. rows == 0) return
      if (.not. all_constraints_finite(restrictions)) then
         info = 3
         return
      end if
      target = 0
      if (allocated(restrictions%parameter_matrix)) then
         source = size(restrictions%parameter_matrix, 1)
         matrix(target + 1:target + source, :) = restrictions%parameter_matrix
         value(target + 1:target + source) = restrictions%parameter_value
         target = target + source
      end if
      if (allocated(restrictions%order_matrix)) then
         source = size(restrictions%order_matrix, 1)
         if (equal_orders) then
            matrix(target + 1:target + source, 1) = &
               restrictions%order_matrix(:, 1) + &
               restrictions%order_matrix(:, 2)
         else
            matrix(target + 1:target + source, 1:2) = &
               restrictions%order_matrix
         end if
         value(target + 1:target + source) = restrictions%order_value
         target = target + source
      end if
      level_start = order_count
      if (allocated(restrictions%level_matrix)) then
         source = size(restrictions%level_matrix, 1)
         matrix(target + 1:target + source, &
            level_start + 1:level_start + size(fit%level)) = &
            restrictions%level_matrix
         value(target + 1:target + source) = restrictions%level_value
         target = target + source
      end if
      alpha_start = level_start + size(fit%level)
      if (allocated(restrictions%alpha_matrix)) then
         source = size(restrictions%alpha_matrix, 1)
         matrix(target + 1:target + source, &
            alpha_start + 1:alpha_start + alpha_count) = &
            restrictions%alpha_matrix
         value(target + 1:target + source) = restrictions%alpha_value
         target = target + source
      end if
      beta_start = alpha_start + alpha_count
      rho_start = beta_start + (variables - fit%rank)*fit%rank
      if (allocated(restrictions%beta_matrix)) then
         source = size(restrictions%beta_matrix, 1)
         value(target + 1:target + source) = restrictions%beta_value
         do restriction = 1, source
            do relation = 1, fit%rank
               do variable = 1, beta_rows
                  coefficient = restrictions%beta_matrix(restriction, &
                     (relation - 1)*beta_rows + variable)
                  if (variable <= fit%rank) then
                     fixed_value = merge(1.0_dp, 0.0_dp, variable == relation)
                     value(target + restriction) = &
                        value(target + restriction) - coefficient*fixed_value
                  else if (variable <= variables) then
                     packed = beta_start + (relation - 1)* &
                        (variables - fit%rank) + variable - fit%rank
                     matrix(target + restriction, packed) = &
                        matrix(target + restriction, packed) + coefficient
                  else
                     packed = rho_start + relation
                     matrix(target + restriction, packed) = &
                        matrix(target + restriction, packed) + coefficient
                  end if
               end do
            end do
         end do
      end if
   end subroutine build_constraint_system

   pure logical function all_constraints_finite(restrictions) result(finite)
      !! Check every allocated restriction matrix and right-hand side.
      type(fcvar_restrictions_t), intent(in) :: restrictions !! Restrictions.

      finite = .true.
      if (allocated(restrictions%parameter_matrix)) finite = finite .and. &
         all(ieee_is_finite(restrictions%parameter_matrix)) .and. &
         all(ieee_is_finite(restrictions%parameter_value))
      if (allocated(restrictions%order_matrix)) finite = finite .and. &
         all(ieee_is_finite(restrictions%order_matrix)) .and. &
         all(ieee_is_finite(restrictions%order_value))
      if (allocated(restrictions%order_inequality_matrix) .and. &
         allocated(restrictions%order_inequality_value)) then
         finite = finite .and. all(ieee_is_finite( &
            restrictions%order_inequality_matrix)) .and. &
            all(ieee_is_finite(restrictions%order_inequality_value))
      else if (allocated(restrictions%order_inequality_matrix) .or. &
         allocated(restrictions%order_inequality_value)) then
         finite = .false.
      end if
      if (allocated(restrictions%alpha_matrix)) finite = finite .and. &
         all(ieee_is_finite(restrictions%alpha_matrix)) .and. &
         all(ieee_is_finite(restrictions%alpha_value))
      if (allocated(restrictions%beta_matrix)) finite = finite .and. &
         all(ieee_is_finite(restrictions%beta_matrix)) .and. &
         all(ieee_is_finite(restrictions%beta_value))
      if (allocated(restrictions%level_matrix)) finite = finite .and. &
         all(ieee_is_finite(restrictions%level_matrix)) .and. &
         all(ieee_is_finite(restrictions%level_value))
   end function all_constraints_finite

   pure subroutine constraint_geometry(matrix, value, initial, particular, &
      basis, rank, info)
      !! Project to the affine constraint set and construct its null space.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), intent(in) :: value(:) !! Input value.
      real(dp), intent(in) :: initial(:) !! Initial value.
      real(dp), allocatable, intent(out) :: particular(:) !! Particular.
      real(dp), allocatable, intent(out) :: basis(:, :) !! Basis.
      integer, intent(out) :: rank !! Matrix or cointegration rank.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: gram(:, :), inverse(:, :), cross(:, :)
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      real(dp) :: tolerance, residual_tolerance
      integer :: rows, columns

      rows = size(matrix, 1)
      columns = size(matrix, 2)
      info = 0
      rank = 0
      if (rows < 1 .or. size(value) /= rows .or. size(initial) /= columns) then
         allocate(particular(0), basis(0, 0))
         info = 1
         return
      end if
      gram = matmul(matrix, transpose(matrix))
      allocate(inverse(rows, rows))
      tolerance = 1.0e-11_dp*max(1.0_dp, maxval(abs(gram)))
      call symmetric_pseudoinverse(gram, inverse, info, tolerance)
      if (info /= 0) return
      particular = initial - matmul(transpose(matrix), &
         matmul(inverse, matmul(matrix, initial) - value))
      residual_tolerance = 1.0e-8_dp*max(1.0_dp, maxval(abs(value)), &
         maxval(abs(matmul(matrix, particular))))
      if (maxval(abs(matmul(matrix, particular) - value)) > &
         residual_tolerance) then
         info = 2
         return
      end if
      cross = matmul(transpose(matrix), matrix)
      call symmetric_eigen(symmetrize(cross), eigenvalues, eigenvectors, info)
      if (info /= 0) return
      tolerance = 1.0e-10_dp*max(1.0_dp, maxval(abs(eigenvalues)))
      rank = count(eigenvalues > tolerance)
      allocate(basis(columns, columns - rank))
      if (rank < columns) basis = eigenvectors(:, rank + 1:)
   end subroutine constraint_geometry

   pure subroutine maximize_equal_order(x, lag_order, rank, bounds, options, &
      best_order, best_likelihood)
      !! Maximize the likelihood on the equality line d equals b.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: bounds(2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      real(dp), intent(out) :: best_order !! Best order.
      real(dp), intent(out) :: best_likelihood !! Best likelihood.
      real(dp) :: lower, upper, left, right, left_value, right_value, value
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      integer :: iteration

      lower = bounds(1)
      upper = bounds(2)
      left = upper - ratio*(upper - lower)
      right = lower + ratio*(upper - lower)
      left_value = fcvar_log_likelihood_fixed(x, lag_order, rank, left, left, &
         options)
      right_value = fcvar_log_likelihood_fixed(x, lag_order, rank, right, right, &
         options)
      do iteration = 1, 60
         if (left_value < right_value) then
            lower = left
            left = right
            left_value = right_value
            right = lower + ratio*(upper - lower)
            right_value = fcvar_log_likelihood_fixed(x, lag_order, rank, right, &
               right, options)
         else
            upper = right
            right = left
            right_value = left_value
            left = upper - ratio*(upper - lower)
            left_value = fcvar_log_likelihood_fixed(x, lag_order, rank, left, &
               left, options)
         end if
      end do
      if (left_value >= right_value) then
         best_order = left
         best_likelihood = left_value
      else
         best_order = right
         best_likelihood = right_value
      end if
      value = fcvar_log_likelihood_fixed(x, lag_order, rank, bounds(1), &
         bounds(1), options)
      if (value > best_likelihood) then
         best_order = bounds(1)
         best_likelihood = value
      end if
      value = fcvar_log_likelihood_fixed(x, lag_order, rank, bounds(2), &
         bounds(2), options)
      if (value > best_likelihood) then
         best_order = bounds(2)
         best_likelihood = value
      end if
   end subroutine maximize_equal_order

   pure subroutine maximize_b_order(x, lag_order, rank, d, bounds, options, &
      constrained, best_b, best_likelihood)
      !! Maximize the likelihood over b for a fixed d, including both bounds.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: bounds(2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in) :: constrained !! Flag controlling constrained.
      real(dp), intent(out) :: best_b !! Best b.
      real(dp), intent(out) :: best_likelihood !! Best likelihood.
      real(dp) :: lower, upper, left, right, left_value, right_value, value
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      integer :: iteration

      lower = bounds(1)
      upper = bounds(2)
      if (constrained) upper = min(upper, d)
      if (upper <= lower) then
         best_b = upper
         best_likelihood = fcvar_log_likelihood_fixed(x, lag_order, rank, d, &
            best_b, options)
         return
      end if
      left = upper - ratio*(upper - lower)
      right = lower + ratio*(upper - lower)
      left_value = fcvar_log_likelihood_fixed(x, lag_order, rank, d, left, options)
      right_value = fcvar_log_likelihood_fixed(x, lag_order, rank, d, right, options)
      do iteration = 1, 60
         if (left_value < right_value) then
            lower = left
            left = right
            left_value = right_value
            right = lower + ratio*(upper - lower)
            right_value = fcvar_log_likelihood_fixed(x, lag_order, rank, d, &
               right, options)
         else
            upper = right
            right = left
            right_value = left_value
            left = upper - ratio*(upper - lower)
            left_value = fcvar_log_likelihood_fixed(x, lag_order, rank, d, &
               left, options)
         end if
      end do
      if (left_value >= right_value) then
         best_b = left
         best_likelihood = left_value
      else
         best_b = right
         best_likelihood = right_value
      end if
      value = fcvar_log_likelihood_fixed(x, lag_order, rank, d, bounds(1), &
         options)
      if (value > best_likelihood) then
         best_b = bounds(1)
         best_likelihood = value
      end if
      value = fcvar_log_likelihood_fixed(x, lag_order, rank, d, &
         min(bounds(2), merge(d, bounds(2), constrained)), options)
      if (value > best_likelihood) then
         best_b = min(bounds(2), merge(d, bounds(2), constrained))
         best_likelihood = value
      end if
   end subroutine maximize_b_order

   pure subroutine maximize_d_order(x, lag_order, rank, b, bounds, options, &
      constrained, best_d, best_likelihood)
      !! Maximize the likelihood over d for a fixed b, including both bounds.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: b !! B.
      real(dp), intent(in) :: bounds(2) !! Bounds.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in), optional :: options !! Algorithm options.
      logical, intent(in) :: constrained !! Flag controlling constrained.
      real(dp), intent(out) :: best_d !! Best d.
      real(dp), intent(out) :: best_likelihood !! Best likelihood.
      real(dp) :: lower, upper, left, right, left_value, right_value, value
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      integer :: iteration

      lower = bounds(1)
      if (constrained) lower = max(lower, b)
      upper = bounds(2)
      if (upper <= lower) then
         best_d = lower
         best_likelihood = fcvar_log_likelihood_fixed(x, lag_order, rank, &
            best_d, b, options)
         return
      end if
      left = upper - ratio*(upper - lower)
      right = lower + ratio*(upper - lower)
      left_value = fcvar_log_likelihood_fixed(x, lag_order, rank, left, b, options)
      right_value = fcvar_log_likelihood_fixed(x, lag_order, rank, right, b, options)
      do iteration = 1, 60
         if (left_value < right_value) then
            lower = left
            left = right
            left_value = right_value
            right = lower + ratio*(upper - lower)
            right_value = fcvar_log_likelihood_fixed(x, lag_order, rank, right, &
               b, options)
         else
            upper = right
            right = left
            right_value = left_value
            left = upper - ratio*(upper - lower)
            left_value = fcvar_log_likelihood_fixed(x, lag_order, rank, left, &
               b, options)
         end if
      end do
      if (left_value >= right_value) then
         best_d = left
         best_likelihood = left_value
      else
         best_d = right
         best_likelihood = right_value
      end if
      value = fcvar_log_likelihood_fixed(x, lag_order, rank, &
         max(bounds(1), merge(b, bounds(1), constrained)), b, options)
      if (value > best_likelihood) then
         best_d = max(bounds(1), merge(b, bounds(1), constrained))
         best_likelihood = value
      end if
      value = fcvar_log_likelihood_fixed(x, lag_order, rank, bounds(2), b, &
         options)
      if (value > best_likelihood) then
         best_d = bounds(2)
         best_likelihood = value
      end if
   end subroutine maximize_d_order

   pure logical function orders_feasible(orders, matrix, value, tolerance) &
      result(feasible)
      !! Test an order pair against affine greater-than-or-equal constraints.
      real(dp), intent(in) :: orders(:) !! Orders.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), intent(in) :: value(:) !! Input value.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.

      feasible = size(orders) == 2 .and. size(matrix, 2) == 2 .and. &
         size(value) == size(matrix, 1)
      if (.not. feasible) return
      feasible = all(matmul(matrix, orders) >= value - tolerance)
   end function orders_feasible

   pure subroutine feasible_line_interval(coefficient, value, lower, upper, info)
      !! Intersect scalar affine inequalities into one closed interval.
      real(dp), intent(in) :: coefficient(:) !! Coefficient.
      real(dp), intent(in) :: value(:) !! Input value.
      real(dp), intent(out) :: lower !! Lower.
      real(dp), intent(out) :: upper !! Upper.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: tolerance
      integer :: constraint

      lower = -huge(1.0_dp)
      upper = huge(1.0_dp)
      info = 0
      if (size(coefficient) /= size(value)) then
         info = 1
         return
      end if
      tolerance = 1.0e-12_dp
      do constraint = 1, size(value)
         if (coefficient(constraint) > tolerance) then
            lower = max(lower, value(constraint)/coefficient(constraint))
         else if (coefficient(constraint) < -tolerance) then
            upper = min(upper, value(constraint)/coefficient(constraint))
         else if (value(constraint) > tolerance) then
            info = 2
            return
         end if
      end do
      if (lower > upper + 1.0e-10_dp*max(1.0_dp, abs(lower), abs(upper))) then
         info = 3
      else if (lower > upper) then
         lower = 0.5_dp*(lower + upper)
         upper = lower
      end if
   end subroutine feasible_line_interval

   pure subroutine boundary_line_interval(point, direction, matrix, value, &
      lower, upper, info)
      !! Find the feasible coordinate interval on an affine boundary line.
      real(dp), intent(in) :: point(2) !! Point.
      real(dp), intent(in) :: direction(2) !! Direction.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), intent(in) :: value(:) !! Input value.
      real(dp), intent(out) :: lower !! Lower.
      real(dp), intent(out) :: upper !! Upper.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: coefficient(:), reduced_value(:)

      coefficient = matmul(matrix, direction)
      reduced_value = value - matmul(matrix, point)
      call feasible_line_interval(coefficient, reduced_value, lower, upper, info)
   end subroutine boundary_line_interval

   pure subroutine maximize_fcvar_order_line(x, lag_order, rank, options, &
      line_lower, line_upper, point, direction, orders, likelihood)
      !! Maximize the FCVAR likelihood over one closed affine line segment.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: line_lower !! Line lower.
      real(dp), intent(in) :: line_upper !! Line upper.
      real(dp), intent(in) :: point(2) !! Point.
      real(dp), intent(in) :: direction(2) !! Direction.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in) :: options !! Algorithm options.
      real(dp), allocatable, intent(out) :: orders(:) !! Orders.
      real(dp), intent(out) :: likelihood !! Likelihood.
      real(dp) :: left, right, left_value, right_value, lower, upper
      real(dp) :: endpoint_value
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      integer :: iteration

      lower = line_lower
      upper = line_upper
      left = upper - ratio*(upper - lower)
      right = lower + ratio*(upper - lower)
      left_value = fcvar_order_line_likelihood(x, lag_order, rank, options, &
         point, direction, left)
      right_value = fcvar_order_line_likelihood(x, lag_order, rank, options, &
         point, direction, right)
      do iteration = 1, 60
         if (left_value < right_value) then
            lower = left
            left = right
            left_value = right_value
            right = lower + ratio*(upper - lower)
            right_value = fcvar_order_line_likelihood(x, lag_order, rank, &
               options, point, direction, right)
         else
            upper = right
            right = left
            right_value = left_value
            left = upper - ratio*(upper - lower)
            left_value = fcvar_order_line_likelihood(x, lag_order, rank, &
               options, point, direction, left)
         end if
      end do
      if (left_value >= right_value) then
         orders = point + left*direction
         likelihood = left_value
      else
         orders = point + right*direction
         likelihood = right_value
      end if
      endpoint_value = fcvar_order_line_likelihood(x, lag_order, rank, options, &
         point, direction, line_lower)
      if (endpoint_value > likelihood) then
         orders = point + line_lower*direction
         likelihood = endpoint_value
      end if
      endpoint_value = fcvar_order_line_likelihood(x, lag_order, rank, options, &
         point, direction, line_upper)
      if (endpoint_value > likelihood) then
         orders = point + line_upper*direction
         likelihood = endpoint_value
      end if
   end subroutine maximize_fcvar_order_line

   pure function fcvar_order_line_likelihood(x, lag_order, rank, options, &
      point, direction, coordinate) result(value)
      !! Evaluate the fixed-order likelihood at one affine-line coordinate.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: point(2) !! Point.
      real(dp), intent(in) :: direction(2) !! Direction.
      real(dp), intent(in) :: coordinate !! Coordinate.
      integer, intent(in) :: lag_order !! Model lag order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      type(fcvar_options_t), intent(in) :: options !! Algorithm options.
      real(dp) :: value
      real(dp) :: orders(2)

      orders = point + coordinate*direction
      value = fcvar_log_likelihood_fixed(x, lag_order, rank, orders(1), &
         orders(2), options)
   end function fcvar_order_line_likelihood

   pure function regular_grid(bounds, points) result(grid)
      !! Construct an equally spaced grid including both bounds.
      real(dp), intent(in) :: bounds(2) !! Bounds.
      integer, intent(in) :: points !! Points.
      real(dp), allocatable :: grid(:)
      integer :: index

      allocate(grid(points))
      do index = 1, points
         grid(index) = bounds(1) + real(index - 1, dp)* &
            (bounds(2) - bounds(1))/real(points - 1, dp)
      end do
   end function regular_grid

   pure subroutine encode_orders(d, b, d_bounds, b_bounds, constrained, &
      coordinates)
      !! Map bounded fractional orders to unconstrained optimizer coordinates.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: b !! B.
      real(dp), intent(in) :: d_bounds(2) !! D bounds.
      real(dp), intent(in) :: b_bounds(2) !! B bounds.
      logical, intent(in) :: constrained !! Flag controlling constrained.
      real(dp), intent(out) :: coordinates(2) !! Coordinates.
      real(dp) :: d_fraction, b_fraction, upper_b, tolerance

      tolerance = 1.0e-3_dp
      d_fraction = (d - d_bounds(1))/(d_bounds(2) - d_bounds(1))
      d_fraction = min(1.0_dp - tolerance, max(tolerance, d_fraction))
      coordinates(1) = log(d_fraction/(1.0_dp - d_fraction))
      upper_b = b_bounds(2)
      if (constrained) upper_b = min(upper_b, d)
      if (upper_b <= b_bounds(1)) then
         b_fraction = tolerance
      else
         b_fraction = (b - b_bounds(1))/(upper_b - b_bounds(1))
         b_fraction = min(1.0_dp - tolerance, max(tolerance, b_fraction))
      end if
      coordinates(2) = log(b_fraction/(1.0_dp - b_fraction))
   end subroutine encode_orders

   pure subroutine decode_orders(coordinates, d_bounds, b_bounds, constrained, &
      d, b)
      !! Map unconstrained optimizer coordinates into the bounded parameter space.
      real(dp), intent(in) :: coordinates(:) !! Coordinates.
      real(dp), intent(in) :: d_bounds(2) !! D bounds.
      real(dp), intent(in) :: b_bounds(2) !! B bounds.
      logical, intent(in) :: constrained !! Flag controlling constrained.
      real(dp), intent(out) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(out) :: b !! B.
      real(dp) :: d_fraction, b_fraction, upper_b

      d_fraction = stable_logistic(coordinates(1))
      d = d_bounds(1) + (d_bounds(2) - d_bounds(1))*d_fraction
      upper_b = b_bounds(2)
      if (constrained) upper_b = min(upper_b, d)
      if (upper_b <= b_bounds(1)) then
         b = upper_b
      else
         b_fraction = stable_logistic(coordinates(2))
         b = b_bounds(1) + (upper_b - b_bounds(1))*b_fraction
      end if
   end subroutine decode_orders

   pure elemental function stable_logistic(value) result(transformed)
      !! Evaluate a logistic transform without exponential overflow.
      real(dp), intent(in) :: value !! Input value.
      real(dp) :: transformed

      if (value >= 0.0_dp) then
         transformed = 1.0_dp/(1.0_dp + exp(-value))
      else
         transformed = exp(value)/(1.0_dp + exp(value))
      end if
   end function stable_logistic

   pure subroutine covariance_log_determinant(covariance, log_determinant, info)
      !! Compute a positive-definite covariance log determinant by Cholesky.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix.
      real(dp), intent(out) :: log_determinant !! Log determinant.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: lower(:, :)
      integer :: index

      call cholesky_lower(symmetrize(covariance), lower, info)
      if (info /= 0) then
         log_determinant = huge(1.0_dp)
         return
      end if
      log_determinant = 0.0_dp
      do index = 1, size(lower, 1)
         log_determinant = log_determinant + 2.0_dp*log(lower(index, index))
      end do
   end subroutine covariance_log_determinant

end module fcvar_mod
