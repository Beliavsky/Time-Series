! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R mAr package.
! Multivariate autoregressive analysis translated from the CRAN mAr package.
module mar_mod
   use kind_mod, only: dp
   use time_series_linalg_mod, only: cholesky_lower, identity_matrix, &
      invert_matrix, symmetric_eigen, general_eigenvalues
   use time_series_random_mod, only: random_standard_normal_matrix
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value
   use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf
   implicit none
   private

   type, public :: mar_fit_t
      ! Stabilized VAR estimates, residuals, covariance, and Schwarz criterion.
      real(dp), allocatable :: intercept(:), coefficients(:, :)
      real(dp), allocatable :: ar(:, :, :), covariance(:, :), residuals(:, :)
      real(dp) :: sbc = 0.0_dp
      integer :: order = 0
      integer :: info = 0
   end type mar_fit_t

   type, public :: mar_modes_t
      ! Companion roots, periods, damping times, and phase-normalized modes.
      complex(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: periods(:), damping_times(:)
      logical :: stable = .false.
      integer :: info = 0
   end type mar_modes_t

   type, public :: mar_pca_fit_t
      ! PCA-space stabilized VAR fit and modes mapped to original variables.
      type(mar_fit_t) :: model
      type(mar_modes_t) :: modes
      real(dp), allocatable :: fraction_variance(:)
      real(dp), allocatable :: scores(:, :), loadings(:, :)
      complex(dp), allocatable :: eigenvectors(:, :)
      integer :: components = 0
      integer :: info = 0
   end type mar_pca_fit_t

   type, public :: mar_simulation_t
      ! Simulated VAR observations and the innovations used after burn-in.
      real(dp), allocatable :: series(:, :), innovations(:, :)
      integer :: info = 0
   end type mar_simulation_t

   public :: mar_estimate, mar_eigenmodes, mar_pca
   public :: mar_simulate_from_innovations, mar_simulate

contains

   pure function mar_estimate(x, order) result(out)
      ! Estimate a VAR by mAr's augmented, column-scaled QR algorithm.
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: order
      type(mar_fit_t) :: out
      real(dp), allocatable :: augmented(:, :), design_response(:, :)
      real(dp), allocatable :: q_factor(:, :), r_factor(:, :), scale(:)
      real(dp), allocatable :: r11(:, :), r12(:, :), r22(:, :), solution(:, :)
      real(dp) :: delta, ratio, log_determinant
      integer :: n, dimension, observations, predictors, columns
      integer :: row, lag, status

      n = size(x, 1)
      dimension = size(x, 2)
      observations = n - order
      predictors = dimension*order + 1
      columns = predictors + dimension
      if (dimension < 1 .or. order < 1 .or. observations <= predictors .or. &
         .not. all(ieee_is_finite(x))) then
         out%info = 1
         return
      end if
      allocate(design_response(observations, columns))
      design_response(:, 1) = 1.0_dp
      do row = 1, observations
         do lag = 1, order
            design_response(row, 2 + (lag - 1)*dimension: &
               1 + lag*dimension) = x(order + row - lag, :)
         end do
         design_response(row, predictors + 1:columns) = x(order + row, :)
      end do
      delta = real(columns*columns + columns + 1, dp)*epsilon(1.0_dp)
      allocate(scale(columns))
      scale = sqrt(delta)*sqrt(sum(design_response**2, dim=1))
      if (any(scale <= tiny(1.0_dp))) then
         out%info = 2
         return
      end if
      allocate(augmented(observations + columns, columns))
      augmented = 0.0_dp
      augmented(:observations, :) = design_response
      do row = 1, columns
         augmented(observations + row, row) = scale(row)
      end do
      call thin_qr(augmented, q_factor, r_factor, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      r11 = r_factor(:predictors, :predictors)
      r12 = r_factor(:predictors, predictors + 1:columns)
      r22 = r_factor(predictors + 1:columns, predictors + 1:columns)
      ratio = maxval(scale(2:columns))/scale(1)
      r11(:, 1) = ratio*r11(:, 1)
      call solve_upper_matrix(r11, r12, solution, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      solution = transpose(solution)
      allocate(out%intercept(dimension))
      allocate(out%coefficients(dimension, dimension*order))
      allocate(out%ar(dimension, dimension, order))
      allocate(out%covariance(dimension, dimension))
      allocate(out%residuals(observations, dimension))
      out%intercept = ratio*solution(:, 1)
      out%coefficients = solution(:, 2:predictors)
      do lag = 1, order
         out%ar(:, :, lag) = out%coefficients(:, &
            1 + (lag - 1)*dimension:lag*dimension)
      end do
      out%covariance = matmul(transpose(r22), r22)/ &
         real(observations - predictors, dp)
      do row = 1, observations
         out%residuals(row, :) = x(order + row, :) - out%intercept
         do lag = 1, order
            out%residuals(row, :) = out%residuals(row, :) - &
               matmul(out%ar(:, :, lag), x(order + row - lag, :))
         end do
      end do
      log_determinant = 2.0_dp*sum(log(abs([(r22(row, row), &
         row=1, dimension)])))
      out%sbc = log_determinant/real(dimension, dp) - &
         log(real(observations, dp))* &
         real(observations - predictors, dp)/real(observations, dp)
      out%order = order
   end function mar_estimate

   pure function mar_eigenmodes(coefficients) result(out)
      ! Compute phase-normalized companion eigenmodes and decay diagnostics.
      real(dp), intent(in) :: coefficients(:, :)
      type(mar_modes_t) :: out
      real(dp), allocatable :: companion(:, :), real_part(:), imaginary_part(:)
      complex(dp), allocatable :: vector(:)
      real(dp) :: magnitude, phase, norm_real, norm_imaginary, tolerance
      real(dp) :: pi_value
      integer :: dimension, order, states, mode, row, status

      dimension = size(coefficients, 1)
      if (dimension < 1 .or. mod(size(coefficients, 2), dimension) /= 0 .or. &
         .not. all(ieee_is_finite(coefficients))) then
         out%info = 1
         return
      end if
      order = size(coefficients, 2)/dimension
      if (order < 1) then
         out%info = 1
         return
      end if
      states = dimension*order
      allocate(companion(states, states))
      companion = 0.0_dp
      companion(:dimension, :) = coefficients
      if (order > 1) then
         do row = dimension + 1, states
            companion(row, row - dimension) = 1.0_dp
         end do
      end if
      call general_eigenvalues(companion, out%eigenvalues, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      allocate(out%eigenvectors(dimension, states))
      allocate(out%periods(states), out%damping_times(states))
      allocate(real_part(states), imaginary_part(states), vector(states))
      tolerance = 1000.0_dp*epsilon(1.0_dp)
      pi_value = acos(-1.0_dp)
      do mode = 1, states
         call inverse_eigenvector(companion, out%eigenvalues(mode), vector, &
            status)
         if (status /= 0) then
            out%info = 20 + mode
            return
         end if
         real_part = real(vector, dp)
         imaginary_part = aimag(vector)
         phase = 0.5_dp*atan2(2.0_dp*dot_product(real_part, imaginary_part), &
            dot_product(imaginary_part, imaginary_part) - &
            dot_product(real_part, real_part))
         norm_real = sqrt(sum((cos(phase)*real_part - &
            sin(phase)*imaginary_part)**2))
         norm_imaginary = sqrt(sum((sin(phase)*real_part + &
            cos(phase)*imaginary_part)**2))
         if (norm_imaginary > norm_real .and. phase < 0.0_dp) then
            phase = phase - 0.5_dp*pi_value
         else if (norm_imaginary > norm_real .and. phase > 0.0_dp) then
            phase = phase + 0.5_dp*pi_value
         end if
         vector = vector*exp(cmplx(0.0_dp, phase, dp))
         out%eigenvectors(:, mode) = vector(states - dimension + 1:states)
         magnitude = abs(out%eigenvalues(mode))
         if (magnitude <= tiny(1.0_dp)) then
            out%damping_times(mode) = 0.0_dp
         else if (abs(magnitude - 1.0_dp) <= tolerance) then
            out%damping_times(mode) = ieee_value(0.0_dp, ieee_positive_inf)
         else
            out%damping_times(mode) = -2.0_dp/log(magnitude*magnitude)
         end if
         if (abs(aimag(out%eigenvalues(mode))) <= &
            tolerance*max(1.0_dp, magnitude)) then
            if (real(out%eigenvalues(mode), dp) >= 0.0_dp) then
               out%periods(mode) = ieee_value(0.0_dp, ieee_positive_inf)
            else
               out%periods(mode) = 2.0_dp
            end if
         else
            out%periods(mode) = 2.0_dp*pi_value/abs(atan2( &
               aimag(out%eigenvalues(mode)), real(out%eigenvalues(mode), dp)))
         end if
      end do
      call sort_mode_diagnostics(out%periods, out%damping_times)
      out%stable = all(abs(out%eigenvalues) <= 1.0_dp + tolerance)
   end function mar_eigenmodes

   pure function mar_pca(x, order, components) result(out)
      ! Fit the stabilized mAr model in a standardized PCA subspace.
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: order, components
      type(mar_pca_fit_t) :: out
      real(dp), allocatable :: centered(:, :), standardized(:, :)
      real(dp), allocatable :: covariance(:, :), values(:), vectors(:, :)
      real(dp), allocatable :: mean(:), scale(:)
      real(dp) :: total
      integer :: n, variables, variable, status

      n = size(x, 1)
      variables = size(x, 2)
      if (n < 3 .or. variables < 2 .or. order < 1 .or. components <= 1 .or. &
         components > min(n, variables) .or. .not. all(ieee_is_finite(x))) then
         out%info = 1
         return
      end if
      allocate(mean(variables), scale(variables), centered(n, variables))
      allocate(standardized(n, variables))
      mean = sum(x, dim=1)/real(n, dp)
      centered = x - spread(mean, 1, n)
      do variable = 1, variables
         scale(variable) = sqrt(sum(centered(:, variable)**2)/real(n - 1, dp))
         if (scale(variable) <= sqrt(epsilon(1.0_dp))) then
            out%info = 2
            return
         end if
         standardized(:, variable) = centered(:, variable)/scale(variable)
      end do
      covariance = matmul(transpose(standardized), standardized)
      call symmetric_eigen(covariance, values, vectors, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      values = max(values, 0.0_dp)
      total = sum(values)
      if (total <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      out%components = components
      out%loadings = vectors(:, :components)
      out%scores = matmul(standardized, out%loadings)
      allocate(out%fraction_variance(components))
      do variable = 1, components
         out%fraction_variance(variable) = sum(values(:variable))/total
      end do
      out%model = mar_estimate(out%scores, order)
      if (out%model%info /= 0) then
         out%info = 100 + out%model%info
         return
      end if
      out%modes = mar_eigenmodes(out%model%coefficients)
      if (out%modes%info /= 0) then
         out%info = 200 + out%modes%info
         return
      end if
      out%eigenvectors = matmul(cmplx(out%loadings, 0.0_dp, dp), &
         out%modes%eigenvectors)
   end function mar_pca

   pure function mar_simulate_from_innovations(intercept, ar, innovations, &
      burnin) result(out)
      ! Simulate a VAR from supplied innovations after a fixed burn-in.
      real(dp), intent(in) :: intercept(:), ar(:, :, :), innovations(:, :)
      integer, intent(in) :: burnin
      type(mar_simulation_t) :: out
      real(dp), allocatable :: work(:, :), mean_state(:), inverse(:, :)
      real(dp), allocatable :: equilibrium_matrix(:, :)
      integer :: dimension, order, total, observations, time, lag, status

      dimension = size(intercept)
      order = size(ar, 3)
      total = size(innovations, 1)
      observations = total - burnin
      if (dimension < 1 .or. order < 1 .or. burnin < 0 .or. &
         observations < 1 .or. size(ar, 1) /= dimension .or. &
         size(ar, 2) /= dimension .or. &
         size(innovations, 2) /= dimension .or. &
         .not. all(ieee_is_finite(intercept)) .or. &
         .not. all(ieee_is_finite(ar)) .or. &
         .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      equilibrium_matrix = identity_matrix(dimension) - sum(ar, dim=3)
      call invert_matrix(equilibrium_matrix, inverse, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      mean_state = matmul(inverse, intercept)
      allocate(work(order + total, dimension))
      work(:order, :) = spread(mean_state, 1, order)
      do time = order + 1, order + total
         work(time, :) = intercept + innovations(time - order, :)
         do lag = 1, order
            work(time, :) = work(time, :) + &
               matmul(ar(:, :, lag), work(time - lag, :))
         end do
      end do
      allocate(out%series(observations, dimension))
      allocate(out%innovations(observations, dimension))
      out%series = work(order + burnin + 1:order + total, :)
      out%innovations = innovations(burnin + 1:total, :)
   end function mar_simulate_from_innovations

   function mar_simulate(intercept, ar, covariance, observations, burnin) &
      result(out)
      ! Simulate a Gaussian VAR using the shared random stream.
      real(dp), intent(in) :: intercept(:), ar(:, :, :), covariance(:, :)
      integer, intent(in) :: observations
      integer, intent(in), optional :: burnin
      type(mar_simulation_t) :: out
      real(dp), allocatable :: lower(:, :), standard(:, :), innovations(:, :)
      integer :: dimension, discard, status

      dimension = size(intercept)
      discard = 1000
      if (present(burnin)) discard = burnin
      if (dimension < 1 .or. observations < 1 .or. discard < 0 .or. &
         any(shape(covariance) /= [dimension, dimension])) then
         out%info = 1
         return
      end if
      call cholesky_lower(covariance, lower, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      allocate(standard(dimension, observations + discard))
      allocate(innovations(observations + discard, dimension))
      call random_standard_normal_matrix(standard)
      innovations = transpose(matmul(lower, standard))
      out = mar_simulate_from_innovations(intercept, ar, innovations, discard)
   end function mar_simulate

   pure subroutine thin_qr(a, q, r, info)
      ! Compute a reorthogonalized thin QR factorization.
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: q(:, :), r(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: vector(:)
      real(dp) :: correction, norm_value
      integer :: rows, columns, column, previous

      rows = size(a, 1)
      columns = size(a, 2)
      allocate(q(rows, columns), r(columns, columns), vector(rows))
      q = 0.0_dp
      r = 0.0_dp
      info = 0
      do column = 1, columns
         vector = a(:, column)
         do previous = 1, column - 1
            r(previous, column) = dot_product(q(:, previous), vector)
            vector = vector - r(previous, column)*q(:, previous)
         end do
         do previous = 1, column - 1
            correction = dot_product(q(:, previous), vector)
            r(previous, column) = r(previous, column) + correction
            vector = vector - correction*q(:, previous)
         end do
         norm_value = sqrt(dot_product(vector, vector))
         if (norm_value <= tiny(1.0_dp)) then
            info = column
            return
         end if
         r(column, column) = norm_value
         q(:, column) = vector/norm_value
      end do
   end subroutine thin_qr

   pure subroutine solve_upper_matrix(upper, right, solution, info)
      ! Solve an upper-triangular system with multiple right-hand sides.
      real(dp), intent(in) :: upper(:, :), right(:, :)
      real(dp), allocatable, intent(out) :: solution(:, :)
      integer, intent(out) :: info
      integer :: n, row

      n = size(upper, 1)
      allocate(solution(n, size(right, 2)))
      solution = right
      info = 0
      do row = n, 1, -1
         if (abs(upper(row, row)) <= tiny(1.0_dp)) then
            info = row
            return
         end if
         if (row < n) solution(row, :) = solution(row, :) - &
            matmul(upper(row, row + 1:n), solution(row + 1:n, :))
         solution(row, :) = solution(row, :)/upper(row, row)
      end do
   end subroutine solve_upper_matrix

   pure subroutine inverse_eigenvector(a, eigenvalue, vector, info)
      ! Recover one right eigenvector by regularized inverse iteration.
      real(dp), intent(in) :: a(:, :)
      complex(dp), intent(in) :: eigenvalue
      complex(dp), intent(out) :: vector(:)
      integer, intent(out) :: info
      complex(dp), allocatable :: matrix(:, :), next(:)
      complex(dp) :: shifted_value
      real(dp) :: norm_value, delta
      integer :: n, index, iteration, solve_info

      n = size(a, 1)
      delta = sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(eigenvalue), &
         maxval(abs(a)))
      shifted_value = eigenvalue + cmplx(delta, 0.5_dp*delta, dp)
      allocate(matrix(n, n), next(n))
      matrix = cmplx(a, 0.0_dp, dp)
      do index = 1, n
         matrix(index, index) = matrix(index, index) - shifted_value
      end do
      vector = [(cmplx(1.0_dp + real(index, dp)/real(n, dp), &
         real(index - 1, dp)/real(n, dp), dp), index=1, n)]
      norm_value = sqrt(real(dot_product(vector, vector), dp))
      vector = vector/norm_value
      do iteration = 1, 12
         call solve_complex_system(matrix, vector, next, solve_info)
         if (solve_info /= 0) then
            info = solve_info
            return
         end if
         norm_value = sqrt(real(dot_product(next, next), dp))
         if (norm_value <= tiny(1.0_dp)) then
            info = 1
            return
         end if
         vector = next/norm_value
      end do
      info = 0
   end subroutine inverse_eigenvector

   pure subroutine solve_complex_system(a, b, x, info)
      ! Solve a dense complex system by partial-pivot elimination.
      complex(dp), intent(in) :: a(:, :), b(:)
      complex(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      complex(dp), allocatable :: work(:, :), row(:)
      complex(dp) :: factor
      integer :: n, column, pivot, target

      n = size(b)
      allocate(work(n, n + 1), row(n + 1))
      work(:, :n) = a
      work(:, n + 1) = b
      info = 0
      do column = 1, n
         pivot = column - 1 + maxloc(abs(work(column:, column)), dim=1)
         if (abs(work(pivot, column)) <= tiny(1.0_dp)) then
            info = column
            return
         end if
         if (pivot /= column) then
            row = work(column, :)
            work(column, :) = work(pivot, :)
            work(pivot, :) = row
         end if
         do target = column + 1, n
            factor = work(target, column)/work(column, column)
            work(target, column:n + 1) = work(target, column:n + 1) - &
               factor*work(column, column:n + 1)
         end do
      end do
      x = cmplx(0.0_dp, 0.0_dp, dp)
      do column = n, 1, -1
         x(column) = work(column, n + 1)
         if (column < n) x(column) = x(column) - &
            sum(work(column, column + 1:n)*x(column + 1:n))
         x(column) = x(column)/work(column, column)
      end do
   end subroutine solve_complex_system

   pure subroutine sort_mode_diagnostics(periods, damping_times)
      ! Sort period/damping rows by decreasing damping time.
      real(dp), intent(inout) :: periods(:), damping_times(:)
      real(dp) :: period_value, damping_value
      integer :: index, position

      do index = 2, size(damping_times)
         period_value = periods(index)
         damping_value = damping_times(index)
         position = index - 1
         do while (position >= 1)
            if (damping_times(position) >= damping_value) exit
            damping_times(position + 1) = damping_times(position)
            periods(position + 1) = periods(position)
            position = position - 1
         end do
         damping_times(position + 1) = damping_value
         periods(position + 1) = period_value
      end do
   end subroutine sort_mode_diagnostics

end module mar_mod
