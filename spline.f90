! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original spline basis and penalized-regression infrastructure.
module spline_mod
   !! Reusable B-spline construction, penalization, fitting, and prediction.
   use kind_mod, only: dp
   use stats_mod, only: quantile, sorted
   use linalg_mod, only: invert_matrix
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   type, public :: spline_basis_t
      !! Open-knot B-spline basis definition.
      real(dp), allocatable :: knots(:)
      real(dp) :: lower_bound = 0.0_dp
      real(dp) :: upper_bound = 0.0_dp
      integer :: basis_count = 0
      integer :: degree = 3
      integer :: info = 0
   end type spline_basis_t

   type, public :: penalized_regression_t
      !! Penalized least-squares coefficients, diagnostics, and covariance.
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: effective_df = 0.0_dp
      real(dp) :: gcv = huge(1.0_dp)
      real(dp) :: residual_variance = 0.0_dp
      integer :: info = 0
   end type penalized_regression_t

   type, public :: spline_fit_t
      !! Penalized univariate spline fit and its reusable basis.
      type(spline_basis_t) :: basis
      type(penalized_regression_t) :: regression
      integer :: difference_order = 2
      integer :: info = 0
   end type spline_fit_t

   public :: spline_basis_create, spline_basis_values, spline_basis_matrix
   public :: spline_difference_penalty, penalized_regression_fit
   public :: spline_fit, spline_gcv_fit, spline_predict

contains

   pure function spline_basis_create(data, basis_count, degree) result(basis)
      !! Construct an open B-spline basis with interior sample-quantile knots.
      real(dp), intent(in) :: data(:) !! Data defining the knot range and quantiles.
      integer, intent(in) :: basis_count !! Number of requested basis functions.
      integer, intent(in), optional :: degree !! Nonnegative polynomial degree.
      type(spline_basis_t) :: basis
      real(dp), allocatable :: ordered(:)
      integer :: selected_degree, interior_count, knot, position

      selected_degree = 3
      if (present(degree)) selected_degree = degree
      if (size(data) < 2 .or. selected_degree < 0 .or. &
         basis_count < selected_degree + 1 .or. &
         .not. all(ieee_is_finite(data)) .or. maxval(data) <= minval(data)) then
         basis%info = 1
         return
      end if
      ordered = sorted(data)
      basis%lower_bound = ordered(1)
      basis%upper_bound = ordered(size(ordered))
      basis%basis_count = basis_count
      basis%degree = selected_degree
      interior_count = basis_count - selected_degree - 1
      allocate(basis%knots(basis_count + selected_degree + 1))
      basis%knots(:selected_degree + 1) = basis%lower_bound
      do knot = 1, interior_count
         position = selected_degree + 1 + knot
         basis%knots(position) = quantile(ordered, &
            real(knot, dp)/real(interior_count + 1, dp))
      end do
      basis%knots(selected_degree + interior_count + 2:) = basis%upper_bound
   end function spline_basis_create

   pure function spline_basis_values(basis, value) result(values)
      !! Evaluate all functions of a B-spline basis at one value.
      type(spline_basis_t), intent(in) :: basis !! B-spline basis definition.
      real(dp), intent(in) :: value !! Evaluation coordinate.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: previous(:), current(:)
      real(dp) :: coordinate, left_denominator, right_denominator
      integer :: work_count, function_index, recursion

      if (basis%info /= 0 .or. basis%basis_count < 1 .or. &
         .not. allocated(basis%knots)) then
         allocate(values(0))
         return
      end if
      allocate(values(basis%basis_count))
      values = 0.0_dp
      coordinate = max(basis%lower_bound, min(basis%upper_bound, value))
      if (coordinate >= basis%upper_bound) then
         values(basis%basis_count) = 1.0_dp
         return
      end if
      work_count = size(basis%knots) - 1
      allocate(previous(work_count), current(work_count))
      previous = 0.0_dp
      do function_index = 1, work_count
         if (coordinate >= basis%knots(function_index) .and. &
            coordinate < basis%knots(function_index + 1)) &
            previous(function_index) = 1.0_dp
      end do
      do recursion = 1, basis%degree
         current = 0.0_dp
         do function_index = 1, work_count - recursion
            left_denominator = basis%knots(function_index + recursion) - &
               basis%knots(function_index)
            right_denominator = basis%knots(function_index + recursion + 1) - &
               basis%knots(function_index + 1)
            if (left_denominator > 0.0_dp) current(function_index) = &
               current(function_index) + (coordinate - &
               basis%knots(function_index))*previous(function_index)/left_denominator
            if (right_denominator > 0.0_dp) current(function_index) = &
               current(function_index) + (basis%knots(function_index + &
               recursion + 1) - coordinate)*previous(function_index + 1)/ &
               right_denominator
         end do
         previous = current
      end do
      values = previous(:basis%basis_count)
   end function spline_basis_values

   pure function spline_basis_matrix(basis, values) result(matrix)
      !! Evaluate a B-spline basis for a vector of coordinates.
      type(spline_basis_t), intent(in) :: basis !! B-spline basis definition.
      real(dp), intent(in) :: values(:) !! Evaluation coordinates.
      real(dp), allocatable :: matrix(:, :)
      integer :: row

      if (basis%info /= 0 .or. basis%basis_count < 1) then
         allocate(matrix(0, 0))
         return
      end if
      allocate(matrix(size(values), basis%basis_count))
      do row = 1, size(values)
         matrix(row, :) = spline_basis_values(basis, values(row))
      end do
   end function spline_basis_matrix

   pure function spline_difference_penalty(basis_count, difference_order) &
      result(penalty)
      !! Construct a coefficient finite-difference roughness penalty.
      integer, intent(in) :: basis_count !! Number of spline coefficients.
      integer, intent(in), optional :: difference_order !! Positive difference order.
      real(dp), allocatable :: penalty(:, :)
      real(dp), allocatable :: difference(:, :), next_difference(:, :)
      integer :: selected_order, row, recursion

      selected_order = 2
      if (present(difference_order)) selected_order = difference_order
      if (basis_count < 1 .or. selected_order < 1 .or. &
         selected_order >= basis_count) then
         allocate(penalty(0, 0))
         return
      end if
      allocate(difference(basis_count, basis_count))
      difference = 0.0_dp
      do row = 1, basis_count
         difference(row, row) = 1.0_dp
      end do
      do recursion = 1, selected_order
         allocate(next_difference(size(difference, 1) - 1, basis_count))
         do row = 1, size(next_difference, 1)
            next_difference(row, :) = difference(row + 1, :) - difference(row, :)
         end do
         call move_alloc(next_difference, difference)
      end do
      penalty = matmul(transpose(difference), difference)
   end function spline_difference_penalty

   pure function penalized_regression_fit(design, response, penalty, lambda) &
      result(fit)
      !! Fit a quadratic-penalty least-squares regression.
      real(dp), intent(in) :: design(:, :) !! Regression design matrix.
      real(dp), intent(in) :: response(:) !! Response vector.
      real(dp), intent(in) :: penalty(:, :) !! Symmetric coefficient penalty matrix.
      real(dp), intent(in) :: lambda !! Nonnegative penalty multiplier.
      type(penalized_regression_t) :: fit
      real(dp), allocatable :: cross_product(:, :), inverse(:, :), bread(:, :)
      real(dp) :: denominator
      integer :: columns, row, status

      columns = size(design, 2)
      if (size(design, 1) /= size(response) .or. columns < 1 .or. &
         any(shape(penalty) /= [columns, columns]) .or. lambda < 0.0_dp .or. &
         .not. all(ieee_is_finite(design)) .or. &
         .not. all(ieee_is_finite(response)) .or. &
         .not. all(ieee_is_finite(penalty))) then
         fit%info = 1
         return
      end if
      cross_product = matmul(transpose(design), design)
      bread = cross_product + lambda*penalty
      call invert_matrix(bread, inverse, status)
      if (status /= 0) then
         fit%info = 2
         return
      end if
      fit%coefficients = matmul(inverse, matmul(transpose(design), response))
      fit%fitted = matmul(design, fit%coefficients)
      fit%residuals = response - fit%fitted
      fit%rss = sum(fit%residuals*fit%residuals)
      fit%effective_df = 0.0_dp
      bread = matmul(inverse, cross_product)
      do row = 1, columns
         fit%effective_df = fit%effective_df + bread(row, row)
      end do
      denominator = real(size(response), dp) - fit%effective_df
      if (denominator <= tiny(1.0_dp)) then
         fit%info = 3
         return
      end if
      fit%lambda = lambda
      fit%gcv = real(size(response), dp)*fit%rss/(denominator*denominator)
      fit%residual_variance = fit%rss/denominator
      fit%covariance = fit%residual_variance*matmul( &
         matmul(inverse, cross_product), transpose(inverse))
   end function penalized_regression_fit

   pure function spline_fit(x, y, basis_count, lambda, degree, &
      difference_order) result(fit)
      !! Fit one penalized B-spline for a fixed smoothing parameter.
      real(dp), intent(in) :: x(:) !! Predictor coordinates.
      real(dp), intent(in) :: y(:) !! Response values.
      integer, intent(in) :: basis_count !! Number of B-spline basis functions.
      real(dp), intent(in) :: lambda !! Nonnegative smoothing parameter.
      integer, intent(in), optional :: degree !! Nonnegative spline degree.
      integer, intent(in), optional :: difference_order !! Positive penalty difference order.
      type(spline_fit_t) :: fit
      real(dp), allocatable :: design(:, :), penalty(:, :)
      integer :: selected_degree, selected_difference

      selected_degree = 3
      selected_difference = 2
      if (present(degree)) selected_degree = degree
      if (present(difference_order)) selected_difference = difference_order
      if (size(x) /= size(y) .or. .not. all(ieee_is_finite(y))) then
         fit%info = 1
         return
      end if
      fit%basis = spline_basis_create(x, basis_count, selected_degree)
      if (fit%basis%info /= 0) then
         fit%info = fit%basis%info
         return
      end if
      design = spline_basis_matrix(fit%basis, x)
      penalty = spline_difference_penalty(basis_count, selected_difference)
      fit%regression = penalized_regression_fit(design, y, penalty, lambda)
      fit%difference_order = selected_difference
      fit%info = fit%regression%info
   end function spline_fit

   pure function spline_gcv_fit(x, y, basis_count, lambda_values, degree, &
      difference_order) result(fit)
      !! Select a penalized B-spline smoothing parameter by generalized cross-validation.
      real(dp), intent(in) :: x(:) !! Predictor coordinates.
      real(dp), intent(in) :: y(:) !! Response values.
      integer, intent(in) :: basis_count !! Number of B-spline basis functions.
      real(dp), intent(in), optional :: lambda_values(:) !! Candidate nonnegative multipliers.
      integer, intent(in), optional :: degree !! Nonnegative spline degree.
      integer, intent(in), optional :: difference_order !! Positive penalty difference order.
      type(spline_fit_t) :: fit
      type(spline_fit_t) :: candidate
      real(dp), allocatable :: candidates(:)
      real(dp) :: best_gcv
      integer :: selected_degree, selected_difference, index

      selected_degree = 3
      selected_difference = 2
      if (present(degree)) selected_degree = degree
      if (present(difference_order)) selected_difference = difference_order
      if (present(lambda_values)) then
         if (size(lambda_values) < 1 .or. any(lambda_values < 0.0_dp) .or. &
            .not. all(ieee_is_finite(lambda_values))) then
            fit%info = 1
            return
         end if
         candidates = lambda_values
      else
         allocate(candidates(61))
         do index = 1, size(candidates)
            candidates(index) = exp(log(1.0e-6_dp) + real(index - 1, dp)* &
               (log(1.0e6_dp) - log(1.0e-6_dp))/real(size(candidates) - 1, dp))
         end do
      end if
      best_gcv = huge(1.0_dp)
      do index = 1, size(candidates)
         candidate = spline_fit(x, y, basis_count, candidates(index), &
            selected_degree, selected_difference)
         if (candidate%info == 0 .and. candidate%regression%gcv < best_gcv) then
            best_gcv = candidate%regression%gcv
            fit = candidate
         end if
      end do
      if (.not. ieee_is_finite(best_gcv)) fit%info = 2
   end function spline_gcv_fit

   pure function spline_predict(fit, xnew) result(prediction)
      !! Predict a fitted penalized spline at new coordinates.
      type(spline_fit_t), intent(in) :: fit !! Fitted penalized spline.
      real(dp), intent(in) :: xnew(:) !! New predictor coordinates.
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: design(:, :)

      if (fit%info /= 0 .or. .not. allocated(fit%regression%coefficients)) then
         allocate(prediction(0))
         return
      end if
      design = spline_basis_matrix(fit%basis, xnew)
      prediction = matmul(design, fit%regression%coefficients)
   end function spline_predict

end module spline_mod
