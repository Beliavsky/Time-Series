! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original linear-algebra infrastructure for this Fortran library.
! Shared dense linear-algebra helpers for small and medium time-series models.
module linalg_mod
   use kind_mod, only: dp
   implicit none
   private
   public :: invert_matrix, inverse_logdet, cholesky_lower, symmetric_eigen
   public :: cholesky_lower_semidefinite
   public :: general_eigenvalues
   public :: symmetric_pseudoinverse
   public :: symmetrize, outer_product, diagonal_part, diagonal_matrix
   public :: identity_matrix, matrix_rank
   public :: solve_matrix, solve_upper_matrix, solve_complex_system
   public :: kronecker_product

contains

   pure function symmetrize(a) result(s)
      !! Remove numerical asymmetry from a square matrix.
      real(dp), intent(in) :: a(:, :) !! A.
      real(dp) :: s(size(a, 1), size(a, 2))
      s = 0.5_dp*(a + transpose(a))
   end function symmetrize

   pure function outer_product(a, b) result(product)
      !! Form the matrix outer product of two vectors.
      real(dp), intent(in) :: a(:) !! A.
      real(dp), intent(in) :: b(:) !! B.
      real(dp) :: product(size(a), size(b))
      product = spread(a, 2, size(b))*spread(b, 1, size(a))
   end function outer_product

   pure function diagonal_part(a) result(diagonal)
      !! Return a matrix containing only the input diagonal.
      real(dp), intent(in) :: a(:, :) !! A.
      real(dp) :: diagonal(size(a, 1), size(a, 2))
      integer :: i
      diagonal = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         diagonal(i, i) = a(i, i)
      end do
   end function diagonal_part

   pure function diagonal_matrix(diagonal) result(matrix)
      !! Construct a square matrix with the supplied diagonal.
      real(dp), intent(in) :: diagonal(:) !! Diagonal elements.
      real(dp) :: matrix(size(diagonal), size(diagonal))
      integer :: i

      matrix = 0.0_dp
      do i = 1, size(diagonal)
         matrix(i, i) = diagonal(i)
      end do
   end function diagonal_matrix

   pure function identity_matrix(n) result(identity)
      !! Construct a square identity matrix.
      integer, intent(in) :: n !! Number of observations or elements.
      real(dp) :: identity(n, n)
      integer :: i
      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
   end function identity_matrix

   pure function kronecker_product(left, right) result(product)
      !! Form the Kronecker product of two real matrices.
      real(dp), intent(in) :: left(:, :) !! Left.
      real(dp), intent(in) :: right(:, :) !! Right.
      real(dp) :: product(size(left, 1)*size(right, 1), size(left, 2)*size(right, 2))
      integer :: row, column, right_rows, right_columns

      right_rows = size(right, 1)
      right_columns = size(right, 2)
      do column = 1, size(left, 2)
         do row = 1, size(left, 1)
            product((row - 1)*right_rows + 1:row*right_rows, &
               (column - 1)*right_columns + 1:column*right_columns) = left(row, column)*right
         end do
      end do
   end function kronecker_product

   pure subroutine solve_matrix(matrix, right, solution, info)
      !! Solve a real dense system with multiple right-hand sides.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), intent(in) :: right(:, :) !! Right.
      real(dp), intent(out) :: solution(:, :) !! Solution.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: augmented(size(matrix, 1), size(matrix, 2) + size(right, 2))
      real(dp) :: held(size(augmented, 2)), pivot_value
      integer :: i, j, pivot, n

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. size(right, 1) /= n .or. &
         size(solution, 1) /= n .or. size(solution, 2) /= size(right, 2)) then
         info = 1
         return
      end if
      augmented(:, :n) = matrix
      augmented(:, n + 1:) = right
      info = 0
      do i = 1, n
         pivot = i - 1 + maxloc(abs(augmented(i:, i)), dim=1)
         if (abs(augmented(pivot, i)) <= tiny(1.0_dp)) then
            info = i
            solution = 0.0_dp
            return
         end if
         if (pivot /= i) then
            held = augmented(i, :)
            augmented(i, :) = augmented(pivot, :)
            augmented(pivot, :) = held
         end if
         pivot_value = augmented(i, i)
         augmented(i, :) = augmented(i, :)/pivot_value
         do j = 1, n
            if (j /= i) augmented(j, :) = augmented(j, :) - augmented(j, i)*augmented(i, :)
         end do
      end do
      solution = augmented(:, n + 1:)
   end subroutine solve_matrix

   pure subroutine solve_upper_matrix(upper, right, solution, info)
      !! Solve an upper-triangular system with multiple right-hand sides.
      real(dp), intent(in) :: upper(:, :) !! Upper.
      real(dp), intent(in) :: right(:, :) !! Right.
      real(dp), allocatable, intent(out) :: solution(:, :) !! Solution.
      integer, intent(out) :: info !! Status code; zero indicates success.
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

   pure subroutine solve_complex_system(a, b, x, info)
      !! Solve a dense complex system by partial-pivot elimination.
      complex(dp), intent(in) :: a(:, :) !! A.
      complex(dp), intent(in) :: b(:) !! B.
      complex(dp), intent(out) :: x(:) !! Input data or predictor values.
      integer, intent(out) :: info !! Status code; zero indicates success.
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

   pure integer function matrix_rank(a, tolerance) result(rank)
      !! Estimate matrix rank by tolerance-scaled Gaussian elimination.
      real(dp), intent(in) :: a(:, :) !! A.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      real(dp), allocatable :: work(:, :), row(:)
      integer :: column, pivot_row, active, nrow, ncol
      nrow = size(a, 1)
      ncol = size(a, 2)
      allocate(work(nrow, ncol), row(ncol))
      work = a
      rank = 0
      active = 1
      do column = 1, ncol
         if (active > nrow) exit
         pivot_row = active - 1 + maxloc(abs(work(active:, column)), dim=1)
         if (abs(work(pivot_row, column)) <= tolerance) cycle
         row = work(active, :)
         work(active, :) = work(pivot_row, :)
         work(pivot_row, :) = row
         work(active + 1:, :) = work(active + 1:, :) - &
            spread(work(active + 1:, column)/work(active, column), 2, ncol)* &
            spread(work(active, :), 1, nrow - active)
         rank = rank + 1
         active = active + 1
      end do
   end function matrix_rank

   pure subroutine invert_matrix(a, inverse, info)
      !! Invert a square matrix using partial-pivot Gauss-Jordan elimination.
      real(dp), intent(in) :: a(:, :) !! A.
      real(dp), allocatable, intent(out) :: inverse(:, :) !! Inverse.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: work(:, :), row(:)
      integer :: i, j, pivot_row, n
      n = size(a, 1)
      allocate(work(n, 2*n), row(2*n), inverse(n, n))
      work = 0.0_dp
      work(:, :n) = a
      do i = 1, n
         work(i, n + i) = 1.0_dp
      end do
      info = 0
      do i = 1, n
         pivot_row = i - 1 + maxloc(abs(work(i:, i)), dim=1)
         if (abs(work(pivot_row, i)) <= 100.0_dp*epsilon(1.0_dp)) then
            info = i
            return
         end if
         if (pivot_row /= i) then
            row = work(i, :)
            work(i, :) = work(pivot_row, :)
            work(pivot_row, :) = row
         end if
         work(i, :) = work(i, :)/work(i, i)
         do j = 1, n
            if (j == i) cycle
            work(j, :) = work(j, :) - work(j, i)*work(i, :)
         end do
      end do
      inverse = work(:, n + 1:)
   end subroutine invert_matrix

   pure subroutine inverse_logdet(a, inverse, logdet, info, tolerance)
      !! Invert a positive-definite matrix and compute its log determinant.
      real(dp), intent(in) :: a(:, :) !! A.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(out) :: inverse(:, :) !! Inverse.
      real(dp), intent(out) :: logdet !! Logdet.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: work(:, :), row(:)
      real(dp) :: pivot
      integer :: i, j, n, pivot_row
      n = size(a, 1)
      allocate(work(n, 2*n), row(2*n))
      work = 0.0_dp
      work(:, :n) = a
      do i = 1, n
         work(i, n + i) = 1.0_dp
      end do
      logdet = 0.0_dp
      info = 0
      do i = 1, n
         pivot_row = i - 1 + maxloc(abs(work(i:n, i)), dim=1)
         pivot = work(pivot_row, i)
         if (abs(pivot) <= tolerance) then
            info = i
            inverse = 0.0_dp
            return
         end if
         if (pivot_row /= i) then
            row = work(i, :)
            work(i, :) = work(pivot_row, :)
            work(pivot_row, :) = row
         end if
         logdet = logdet + log(abs(work(i, i)))
         work(i, :) = work(i, :)/work(i, i)
         do j = 1, n
            if (j == i) cycle
            work(j, :) = work(j, :) - work(j, i)*work(i, :)
         end do
      end do
      inverse = work(:, n + 1:)
   end subroutine inverse_logdet

   pure subroutine cholesky_lower(a, lower, info)
      !! Compute an unpivoted lower Cholesky factor.
      real(dp), intent(in) :: a(:, :) !! A.
      real(dp), allocatable, intent(out) :: lower(:, :) !! Lower.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: value
      integer :: i, j, n
      n = size(a, 1)
      allocate(lower(n, n))
      lower = 0.0_dp
      info = 0
      do i = 1, n
         do j = 1, i
            value = a(i, j)
            if (j > 1) value = value - dot_product(lower(i, :j - 1), lower(j, :j - 1))
            if (i == j) then
               if (value <= 100.0_dp*epsilon(1.0_dp)) then
                  info = i
                  return
               end if
               lower(i, j) = sqrt(value)
            else
               lower(i, j) = value/lower(j, j)
            end if
         end do
      end do
   end subroutine cholesky_lower

   pure subroutine cholesky_lower_semidefinite(matrix, lower, info)
      !! Compute a lower Cholesky factor allowing zero diagonal pivots.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), allocatable, intent(out) :: lower(:, :) !! Lower.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: value
      integer :: i, j

      allocate(lower(size(matrix, 1), size(matrix, 2)))
      lower = 0.0_dp
      info = 0
      do i = 1, size(matrix, 1)
         do j = 1, i
            value = matrix(i, j) - dot_product(lower(i, 1:j - 1), lower(j, 1:j - 1))
            if (i == j) then
               if (value < -sqrt(epsilon(1.0_dp))) then
                  info = i
                  return
               end if
               lower(i, j) = sqrt(max(value, 0.0_dp))
            else if (lower(j, j) > tiny(1.0_dp)) then
               lower(i, j) = value/lower(j, j)
            end if
         end do
      end do
   end subroutine cholesky_lower_semidefinite

   pure subroutine symmetric_eigen(a, values, vectors, info)
      !! Compute descending eigenpairs of a symmetric matrix by Jacobi rotations.
      real(dp), intent(in) :: a(:, :) !! A.
      real(dp), allocatable, intent(out) :: values(:) !! Input values.
      real(dp), allocatable, intent(out) :: vectors(:, :) !! Vectors.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: work(:, :), swap(:)
      real(dp) :: app, aqq, apq, tau, tangent, c, s, temp, threshold
      integer :: i, j, p, q, n, iteration, max_iterations, selected
      n = size(a, 1)
      allocate(work(n, n), vectors(n, n), values(n), swap(n))
      work = a
      vectors = 0.0_dp
      do i = 1, n
         vectors(i, i) = 1.0_dp
      end do
      threshold = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))
      max_iterations = 100*n*n
      info = 0
      do iteration = 1, max_iterations
         apq = 0.0_dp
         p = 1
         q = min(2, n)
         do i = 1, n - 1
            do j = i + 1, n
               if (abs(work(i, j)) > abs(apq)) then
                  apq = work(i, j)
                  p = i
                  q = j
               end if
            end do
         end do
         if (abs(apq) <= threshold) exit
         app = work(p, p)
         aqq = work(q, q)
         tau = (aqq - app)/(2.0_dp*apq)
         tangent = sign(1.0_dp, tau)/(abs(tau) + sqrt(1.0_dp + tau*tau))
         c = 1.0_dp/sqrt(1.0_dp + tangent*tangent)
         s = tangent*c
         do i = 1, n
            if (i == p .or. i == q) cycle
            temp = work(i, p)
            work(i, p) = c*temp - s*work(i, q)
            work(p, i) = work(i, p)
            work(i, q) = s*temp + c*work(i, q)
            work(q, i) = work(i, q)
         end do
         work(p, p) = app - tangent*apq
         work(q, q) = aqq + tangent*apq
         work(p, q) = 0.0_dp
         work(q, p) = 0.0_dp
         do i = 1, n
            temp = vectors(i, p)
            vectors(i, p) = c*temp - s*vectors(i, q)
            vectors(i, q) = s*temp + c*vectors(i, q)
         end do
      end do
      if (iteration > max_iterations) then
         info = 1
         return
      end if
      values = [(work(i, i), i=1, n)]
      do i = 1, n - 1
         selected = i - 1 + maxloc(values(i:), dim=1)
         if (selected /= i) then
            temp = values(i)
            values(i) = values(selected)
            values(selected) = temp
            swap = vectors(:, i)
            vectors(:, i) = vectors(:, selected)
            vectors(:, selected) = swap
         end if
      end do
   end subroutine symmetric_eigen

   pure subroutine general_eigenvalues(a, values, info)
      !! Compute general real-matrix eigenvalues by shifted complex QR deflation.
      real(dp), intent(in) :: a(:, :) !! A.
      complex(dp), allocatable, intent(out) :: values(:) !! Input values.
      integer, intent(out) :: info !! Status code; zero indicates success.
      complex(dp), allocatable :: work(:, :), q(:, :), r(:, :), shifted(:, :)
      real(dp), allocatable :: hessenberg(:, :)
      complex(dp) :: trace_value, determinant, discriminant, shift
      real(dp) :: tolerance, scale
      integer :: n, active, iteration, maximum_iterations, index, qr_info

      n = size(a, 1)
      allocate(values(n), work(n, n))
      call upper_hessenberg(a, hessenberg)
      work = cmplx(hessenberg, 0.0_dp, dp)
      tolerance = 1000.0_dp*epsilon(1.0_dp)
      maximum_iterations = max(2000, 500*n*n)
      active = n
      info = 0
      do while (active > 1)
         iteration = 0
         do
            scale = max(1.0_dp, abs(work(active - 1, active - 1)) + &
               abs(work(active, active)))
            if (abs(work(active, active - 1)) <= tolerance*scale) exit
            iteration = iteration + 1
            if (iteration > maximum_iterations) then
               info = active
               return
            end if
            trace_value = work(active - 1, active - 1) + work(active, active)
            determinant = work(active - 1, active - 1)*work(active, active) - &
               work(active - 1, active)*work(active, active - 1)
            discriminant = sqrt(trace_value*trace_value - 4.0_dp*determinant)
            shift = 0.5_dp*(trace_value + discriminant)
            if (abs(shift - work(active, active)) > &
               abs(0.5_dp*(trace_value - discriminant) - &
               work(active, active))) then
               shift = 0.5_dp*(trace_value - discriminant)
            end if
            allocate(shifted(active, active))
            shifted = work(:active, :active)
            do index = 1, active
               shifted(index, index) = shifted(index, index) - shift
            end do
            call complex_qr(shifted, q, r, qr_info)
            deallocate(shifted)
            if (qr_info /= 0) then
               info = 100 + qr_info
               return
            end if
            work(:active, :active) = matmul(r, q)
            do index = 1, active
               work(index, index) = work(index, index) + shift
            end do
            deallocate(q, r)
         end do
         values(active) = work(active, active)
         active = active - 1
      end do
      if (active == 1) values(1) = work(1, 1)
   end subroutine general_eigenvalues

   pure subroutine upper_hessenberg(a, hessenberg)
      !! Reduce a real square matrix by Householder similarity transformations.
      real(dp), intent(in) :: a(:, :) !! A.
      real(dp), allocatable, intent(out) :: hessenberg(:, :) !! Hessenberg.
      real(dp), allocatable :: vector(:), product(:)
      real(dp) :: norm_value, reflector_norm
      integer :: n, column, rows

      n = size(a, 1)
      allocate(hessenberg(n, n))
      hessenberg = a
      do column = 1, n - 2
         rows = n - column
         allocate(vector(rows))
         vector = hessenberg(column + 1:, column)
         norm_value = sqrt(dot_product(vector, vector))
         if (norm_value <= tiny(1.0_dp)) then
            deallocate(vector)
            cycle
         end if
         vector(1) = vector(1) + sign(norm_value, vector(1))
         reflector_norm = sqrt(dot_product(vector, vector))
         if (reflector_norm <= tiny(1.0_dp)) then
            deallocate(vector)
            cycle
         end if
         vector = vector/reflector_norm
         allocate(product(n - column + 1))
         product = matmul(vector, hessenberg(column + 1:, column:))
         hessenberg(column + 1:, column:) = &
            hessenberg(column + 1:, column:) - &
            2.0_dp*spread(vector, 2, n - column + 1)* &
            spread(product, 1, rows)
         deallocate(product)
         allocate(product(n))
         product = matmul(hessenberg(:, column + 1:), vector)
         hessenberg(:, column + 1:) = hessenberg(:, column + 1:) - &
            2.0_dp*spread(product, 2, rows)*spread(vector, 1, n)
         hessenberg(column + 2:, column) = 0.0_dp
         deallocate(product, vector)
      end do
   end subroutine upper_hessenberg

   pure subroutine complex_qr(a, q, r, info)
      !! Compute a square complex QR factorization by modified Gram-Schmidt.
      complex(dp), intent(in) :: a(:, :) !! A.
      complex(dp), allocatable, intent(out) :: q(:, :) !! Model order, dimension, or parameter.
      complex(dp), allocatable, intent(out) :: r(:, :) !! R.
      integer, intent(out) :: info !! Status code; zero indicates success.
      complex(dp), allocatable :: vector(:)
      real(dp) :: norm_value
      integer :: n, column, previous

      n = size(a, 1)
      allocate(q(n, n), r(n, n), vector(n))
      q = cmplx(0.0_dp, 0.0_dp, dp)
      r = cmplx(0.0_dp, 0.0_dp, dp)
      info = 0
      do column = 1, n
         vector = a(:, column)
         do previous = 1, column - 1
            r(previous, column) = dot_product(q(:, previous), vector)
            vector = vector - r(previous, column)*q(:, previous)
         end do
         norm_value = sqrt(real(dot_product(vector, vector), dp))
         if (norm_value <= tiny(1.0_dp)) then
            info = column
            return
         end if
         r(column, column) = cmplx(norm_value, 0.0_dp, dp)
         q(:, column) = vector/norm_value
      end do
   end subroutine complex_qr

   pure subroutine symmetric_pseudoinverse(matrix, inverse, info, tolerance)
      !! Form a spectral pseudoinverse of a positive-semidefinite matrix.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), intent(out) :: inverse(:, :) !! Inverse.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      real(dp), allocatable :: values(:), vectors(:, :), reciprocal(:)
      real(dp) :: threshold
      integer :: n

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. any(shape(inverse) /= [n, n])) then
         info = 1
         return
      end if
      call symmetric_eigen(symmetrize(matrix), values, vectors, info)
      if (info /= 0) return
      threshold = 1.0e-12_dp*max(1.0_dp, maxval(abs(values)))
      if (present(tolerance)) threshold = tolerance
      if (threshold < 0.0_dp .or. &
         minval(values) < -100.0_dp*threshold) then
         info = 2
         return
      end if
      allocate(reciprocal(n))
      reciprocal = 0.0_dp
      where (values > threshold) reciprocal = 1.0_dp/values
      inverse = matmul(vectors*spread(reciprocal, 1, n), transpose(vectors))
      inverse = symmetrize(inverse)
      info = 0
   end subroutine symmetric_pseudoinverse
end module linalg_mod
