! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R KFAS package.
! Modern Fortran interface for Gaussian state-space algorithms from KFAS.
module kfas_mod
   use kind_mod, only: dp
   use linalg_mod, only: symmetrize, outer_product, diagonal_part, identity_matrix, matrix_rank, &
      inverse_logdet
   use linalg_mod, only: symmetric_pseudoinverse
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   real(dp), parameter :: log_two_pi = log(2.0_dp*acos(-1.0_dp))

   type, public :: ssm_model_t
      ! Linear Gaussian model y=Z*alpha+epsilon, alpha(+1)=T*alpha+R*eta.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: z(:, :, :), h(:, :, :)
      real(dp), allocatable :: transition(:, :, :), r(:, :, :), q(:, :, :)
      real(dp), allocatable :: a1(:), p1(:, :)
      real(dp), allocatable :: p1inf(:, :)
      logical, allocatable :: missing(:, :)
   end type

   type, public :: kfs_filter_t
      ! Predicted and filtered moments and observation innovations.
      real(dp), allocatable :: a_pred(:, :), p_pred(:, :, :)
      real(dp), allocatable :: a_filt(:, :), p_filt(:, :, :)
      real(dp), allocatable :: innovation(:, :), innovation_cov(:, :, :)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: observations = 0
      integer :: info = 0
   end type

   type, public :: kfs_disturbance_t
      ! Conditional means of observation and state disturbances.
      real(dp), allocatable :: observation(:, :), state(:, :)
      integer :: info = 0
   end type

   type, public :: kfs_smoother_t
      ! Smoothed state moments conditional on the complete observation series.
      real(dp), allocatable :: state(:, :), covariance(:, :, :)
      real(dp), allocatable :: lag_one_covariance(:, :, :)
      real(dp), allocatable :: conditional_matrix(:, :, :)
      real(dp), allocatable :: conditional_covariance(:, :, :)
      integer :: info = 0
   end type

   type, public :: kfs_conditional_t
      ! Forward conditional state regression matrices and covariances.
      real(dp), allocatable :: matrix(:, :, :), covariance(:, :, :)
      integer :: info = 0
   end type

   type, public :: kfs_prediction_t
      ! Future observation and state moments.
      real(dp), allocatable :: mean(:, :), covariance(:, :, :)
      real(dp), allocatable :: state(:, :), state_covariance(:, :, :)
      integer :: info = 0
   end type

   type, public :: kfs_loglik_result_t
      ! Gaussian log likelihood and its filter status code.
      real(dp) :: value = 0.0_dp
      integer :: info = 0
   end type

   public :: kfs_filter, kfs_filter_diffuse, kfs_smooth, kfs_loglik, kfs_predict
   public :: kfs_fast_smooth, kfs_conditional_covariance
   public :: kfs_standardized_innovations, kfs_disturbance_smooth
   public :: make_local_level, make_local_linear_trend, validate_ssm

contains

   pure function make_local_level(y, observation_variance, level_variance, a1, p1) result(model)
      !! Construct a univariate local-level state-space model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: level_variance !! Level variance.
      real(dp), intent(in), optional :: a1 !! A1.
      real(dp), intent(in), optional :: p1 !! P1.
      type(ssm_model_t) :: model

      allocate(model%y(size(y), 1), model%z(1, 1, 1), model%h(1, 1, 1))
      allocate(model%transition(1, 1, 1), model%r(1, 1, 1), model%q(1, 1, 1))
      allocate(model%a1(1), model%p1(1, 1), model%missing(size(y), 1))
      allocate(model%p1inf(1, 1))
      model%y(:, 1) = y
      model%z = 1.0_dp
      model%h = observation_variance
      model%transition = 1.0_dp
      model%r = 1.0_dp
      model%q = level_variance
      model%a1 = y(1)
      model%p1 = observation_variance + level_variance
      model%p1inf = 0.0_dp
      if (present(a1)) model%a1 = a1
      if (present(p1)) model%p1 = p1
      model%missing(:, 1) = .not. ieee_is_finite(y)
   end function make_local_level

   pure function make_local_linear_trend(y, observation_variance, level_variance, slope_variance, a1, p1) result(model)
      !! Construct a univariate local linear-trend state-space model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: level_variance !! Level variance.
      real(dp), intent(in) :: slope_variance !! Slope variance.
      real(dp), intent(in), optional :: a1(2) !! A1.
      real(dp), intent(in), optional :: p1(2, 2) !! P1.
      type(ssm_model_t) :: model

      allocate(model%y(size(y), 1), model%z(1, 2, 1), model%h(1, 1, 1))
      allocate(model%transition(2, 2, 1), model%r(2, 2, 1), model%q(2, 2, 1))
      allocate(model%a1(2), model%p1(2, 2), model%missing(size(y), 1))
      allocate(model%p1inf(2, 2))
      model%y(:, 1) = y
      model%z = 0.0_dp
      model%z(1, 1, 1) = 1.0_dp
      model%h = observation_variance
      model%transition = 0.0_dp
      model%transition(:, :, 1) = reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [2, 2])
      model%r = 0.0_dp
      model%r(1, 1, 1) = 1.0_dp
      model%r(2, 2, 1) = 1.0_dp
      model%q = 0.0_dp
      model%q(1, 1, 1) = level_variance
      model%q(2, 2, 1) = slope_variance
      model%a1 = [y(1), 0.0_dp]
      model%p1 = 0.0_dp
      model%p1(1, 1) = observation_variance + level_variance
      model%p1(2, 2) = max(slope_variance, epsilon(1.0_dp))
      model%p1inf = 0.0_dp
      if (present(a1)) model%a1 = a1
      if (present(p1)) model%p1 = p1
      model%missing(:, 1) = .not. ieee_is_finite(y)
   end function make_local_linear_trend

   pure function kfs_filter_diffuse(model, tolerance) result(out)
      !! Run KFAS exact diffuse filtering for diagonal observation covariance.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(kfs_filter_t) :: out
      real(dp), allocatable :: pinf(:, :), k0(:), kinf(:), z(:), noise(:, :)
      real(dp) :: f0, finf, v, tol
      integer :: t, i, n, p, m, rank_inf

      call validate_ssm(model, out%info)
      if (out%info /= 0) return
      if (.not. allocated(model%p1inf)) then
         out = kfs_filter(model, tolerance)
         return
      end if
      n = size(model%y, 1)
      p = size(model%y, 2)
      m = size(model%a1)
      tol = 100.0_dp*epsilon(1.0_dp)
      if (present(tolerance)) tol = tolerance
      do t = 1, size(model%h, 3)
         if (maxval(abs(model%h(:, :, t) - diagonal_part(model%h(:, :, t)))) > tol) then
            out%info = 2
            return
         end if
      end do
      allocate(out%a_pred(m, n), out%p_pred(m, m, n))
      allocate(out%a_filt(m, n), out%p_filt(m, m, n))
      allocate(out%innovation(n, p), out%innovation_cov(p, p, n))
      allocate(pinf(m, m), k0(m), kinf(m), z(m))
      out%innovation = ieee_value(0.0_dp, ieee_quiet_nan)
      out%innovation_cov = ieee_value(0.0_dp, ieee_quiet_nan)
      out%a_pred(:, 1) = model%a1
      out%p_pred(:, :, 1) = symmetrize(model%p1)
      pinf = symmetrize(model%p1inf)
      rank_inf = matrix_rank(pinf, tol)
      do t = 1, n
         out%a_filt(:, t) = out%a_pred(:, t)
         out%p_filt(:, :, t) = out%p_pred(:, :, t)
         do i = 1, p
            if (.not. is_observed(model, t, i)) cycle
            z = model%z(i, :, time_index(model%z, t))
            k0 = matmul(out%p_filt(:, :, t), z)
            kinf = matmul(pinf, z)
            f0 = dot_product(z, k0) + model%h(i, i, time_index(model%h, t))
            finf = dot_product(z, kinf)
            v = model%y(t, i) - dot_product(z, out%a_filt(:, t))
            out%innovation(t, i) = v
            if (finf > tol) then
               out%innovation_cov(i, i, t) = finf
               out%a_filt(:, t) = out%a_filt(:, t) + kinf*v/finf
               out%p_filt(:, :, t) = out%p_filt(:, :, t) + outer_product(kinf, kinf)*f0/(finf*finf) &
                  - (outer_product(k0, kinf) + outer_product(kinf, k0))/finf
               pinf = pinf - outer_product(kinf, kinf)/finf
               pinf = symmetrize(pinf)
               rank_inf = max(0, rank_inf - 1)
               out%log_likelihood = out%log_likelihood - 0.5_dp*log(finf)
            else if (f0 > tol) then
               out%innovation_cov(i, i, t) = f0
               out%a_filt(:, t) = out%a_filt(:, t) + k0*v/f0
               out%p_filt(:, :, t) = out%p_filt(:, :, t) - outer_product(k0, k0)/f0
               out%log_likelihood = out%log_likelihood - 0.5_dp*(log_two_pi + log(f0) + v*v/f0)
            end if
            out%p_filt(:, :, t) = symmetrize(out%p_filt(:, :, t))
            out%observations = out%observations + 1
         end do
         if (t < n) then
            out%a_pred(:, t + 1) = matmul(model%transition(:, :, time_index(model%transition, t)), out%a_filt(:, t))
            noise = matmul(matmul(model%r(:, :, time_index(model%r, t)), model%q(:, :, time_index(model%q, t))), &
               transpose(model%r(:, :, time_index(model%r, t))))
            out%p_pred(:, :, t + 1) = matmul(matmul(model%transition(:, :, time_index(model%transition, t)), &
               out%p_filt(:, :, t)), transpose(model%transition(:, :, time_index(model%transition, t)))) + noise
            pinf = matmul(matmul(model%transition(:, :, time_index(model%transition, t)), pinf), &
               transpose(model%transition(:, :, time_index(model%transition, t))))
            pinf = symmetrize(pinf)
         end if
      end do
      if (rank_inf > 0) out%info = 3
   end function kfs_filter_diffuse

   pure subroutine validate_ssm(model, info)
      !! Check state-space dimensions and supported time variation.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      integer, intent(out) :: info !! Status code; zero indicates success.
      integer :: n, p, m, rr

      info = 1
      if (.not. allocated(model%y)) return
      if (.not. allocated(model%z) .or. .not. allocated(model%h)) return
      if (.not. allocated(model%transition) .or. .not. allocated(model%r)) return
      if (.not. allocated(model%q) .or. .not. allocated(model%a1)) return
      if (.not. allocated(model%p1)) return
      n = size(model%y, 1)
      p = size(model%y, 2)
      m = size(model%a1)
      rr = size(model%r, 2)
      if (n < 1 .or. p < 1 .or. m < 1 .or. rr < 1) return
      if (size(model%z, 1) /= p .or. size(model%z, 2) /= m) return
      if (size(model%h, 1) /= p .or. size(model%h, 2) /= p) return
      if (size(model%transition, 1) /= m .or. size(model%transition, 2) /= m) return
      if (size(model%r, 1) /= m) return
      if (size(model%q, 1) /= rr .or. size(model%q, 2) /= rr) return
      if (any(shape(model%p1) /= [m, m])) return
      if (allocated(model%p1inf)) then
         if (any(shape(model%p1inf) /= [m, m])) return
      end if
      if (.not. valid_time_extent(size(model%z, 3), n)) return
      if (.not. valid_time_extent(size(model%h, 3), n)) return
      if (.not. valid_time_extent(size(model%transition, 3), n)) return
      if (.not. valid_time_extent(size(model%r, 3), n)) return
      if (.not. valid_time_extent(size(model%q, 3), n)) return
      if (allocated(model%missing)) then
         if (any(shape(model%missing) /= shape(model%y))) return
      end if
      info = 0
   end subroutine validate_ssm

   pure function kfs_filter(model, tolerance) result(out)
      !! Run the multivariate Gaussian Kalman filter and evaluate its likelihood.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(kfs_filter_t) :: out
      real(dp), allocatable :: zobs(:, :), hobs(:, :), f(:, :), finv(:, :), gain(:, :)
      real(dp), allocatable :: v(:), yobs(:), pwork(:, :), noise(:, :)
      real(dp) :: logdet, tol
      integer, allocatable :: observed(:)
      integer :: t, i, j, n, p, m, no, info

      call validate_ssm(model, out%info)
      if (out%info /= 0) return
      n = size(model%y, 1)
      p = size(model%y, 2)
      m = size(model%a1)
      tol = 100.0_dp*epsilon(1.0_dp)
      if (present(tolerance)) tol = tolerance
      allocate(out%a_pred(m, n), out%p_pred(m, m, n))
      allocate(out%a_filt(m, n), out%p_filt(m, m, n))
      allocate(out%innovation(n, p), out%innovation_cov(p, p, n))
      out%innovation = ieee_value(0.0_dp, ieee_quiet_nan)
      out%innovation_cov = ieee_value(0.0_dp, ieee_quiet_nan)
      out%a_pred(:, 1) = model%a1
      out%p_pred(:, :, 1) = symmetrize(model%p1)
      out%log_likelihood = 0.0_dp

      do t = 1, n
         observed = pack([(i, i=1, p)], [(is_observed(model, t, i), i=1, p)])
         no = size(observed)
         if (no > 0) then
            allocate(zobs(no, m), hobs(no, no), yobs(no), v(no), f(no, no), finv(no, no), gain(m, no))
            do i = 1, no
               yobs(i) = model%y(t, observed(i))
               zobs(i, :) = model%z(observed(i), :, time_index(model%z, t))
               do j = 1, no
                  hobs(i, j) = model%h(observed(i), observed(j), time_index(model%h, t))
               end do
            end do
            v = yobs - matmul(zobs, out%a_pred(:, t))
            f = matmul(matmul(zobs, out%p_pred(:, :, t)), transpose(zobs)) + hobs
            call inverse_logdet(f, finv, logdet, info, tol)
            if (info /= 0) then
               out%info = 10 + t
               return
            end if
            gain = matmul(matmul(out%p_pred(:, :, t), transpose(zobs)), finv)
            out%a_filt(:, t) = out%a_pred(:, t) + matmul(gain, v)
            pwork = identity_matrix(m) - matmul(gain, zobs)
            out%p_filt(:, :, t) = &
               matmul(matmul(pwork, out%p_pred(:, :, t)), transpose(pwork)) + &
               matmul(matmul(gain, hobs), transpose(gain))
            out%p_filt(:, :, t) = symmetrize(out%p_filt(:, :, t))
            do i = 1, no
               out%innovation(t, observed(i)) = v(i)
               do j = 1, no
                  out%innovation_cov(observed(i), observed(j), t) = f(i, j)
               end do
            end do
            out%log_likelihood = out%log_likelihood - 0.5_dp*(no*log_two_pi + logdet + dot_product(v, matmul(finv, v)))
            out%observations = out%observations + no
            deallocate(zobs, hobs, yobs, v, f, finv, gain)
         else
            out%a_filt(:, t) = out%a_pred(:, t)
            out%p_filt(:, :, t) = out%p_pred(:, :, t)
         end if
         if (t < n) then
            out%a_pred(:, t + 1) = matmul(model%transition(:, :, time_index(model%transition, t)), out%a_filt(:, t))
            pwork = matmul(matmul(model%transition(:, :, time_index(model%transition, t)), out%p_filt(:, :, t)), &
               transpose(model%transition(:, :, time_index(model%transition, t))))
            noise = matmul(matmul(model%r(:, :, time_index(model%r, t)), model%q(:, :, time_index(model%q, t))), &
               transpose(model%r(:, :, time_index(model%r, t))))
            out%p_pred(:, :, t + 1) = symmetrize(pwork + noise)
         end if
      end do
   end function kfs_filter

   pure function kfs_smooth(model, filtered, tolerance) result(out)
      !! Smooth states, lag-one covariances, and forward conditionals.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      type(kfs_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(kfs_smoother_t) :: out
      real(dp), allocatable :: pinv(:, :), gain(:, :)
      type(kfs_conditional_t) :: conditional
      integer :: t, n, m, info

      if (filtered%info /= 0) then
         out%info = filtered%info
         return
      end if
      n = size(model%y, 1)
      m = size(model%a1)
      allocate(out%state(m, n), out%covariance(m, m, n), pinv(m, m), gain(m, m))
      allocate(out%lag_one_covariance(m, m, max(0, n - 1)))
      out%state(:, n) = filtered%a_filt(:, n)
      out%covariance(:, :, n) = filtered%p_filt(:, :, n)
      do t = n - 1, 1, -1
         call symmetric_pseudoinverse(filtered%p_pred(:, :, t + 1), pinv, &
            info, tolerance)
         if (info /= 0) then
            out%info = 20 + t
            return
         end if
         gain = matmul(matmul(filtered%p_filt(:, :, t), &
            transpose(model%transition(:, :, time_index(model%transition, t)))), pinv)
         out%state(:, t) = filtered%a_filt(:, t) + matmul(gain, out%state(:, t + 1) - filtered%a_pred(:, t + 1))
         out%covariance(:, :, t) = filtered%p_filt(:, :, t) + matmul(matmul(gain, &
            out%covariance(:, :, t + 1) - filtered%p_pred(:, :, t + 1)), transpose(gain))
         out%covariance(:, :, t) = symmetrize(out%covariance(:, :, t))
         out%lag_one_covariance(:, :, t) = &
            matmul(gain, out%covariance(:, :, t + 1))
      end do
      conditional = kfs_conditional_covariance(out%covariance, &
         out%lag_one_covariance, tolerance)
      if (conditional%info /= 0) then
         out%info = 100 + conditional%info
         return
      end if
      out%conditional_matrix = conditional%matrix
      out%conditional_covariance = conditional%covariance
   end function kfs_smooth

   pure function kfs_fast_smooth(model, filtered, tolerance) result(out)
      !! Smooth state means without storing smoothed covariance arrays.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      type(kfs_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(kfs_smoother_t) :: out
      real(dp), allocatable :: inverse(:, :), gain(:, :)
      integer :: state, times, time, info

      if (filtered%info /= 0) then
         out%info = filtered%info
         return
      end if
      state = size(model%a1)
      times = size(model%y, 1)
      allocate(out%state(state, times), inverse(state, state))
      allocate(gain(state, state))
      out%state(:, times) = filtered%a_filt(:, times)
      do time = times - 1, 1, -1
         call symmetric_pseudoinverse(filtered%p_pred(:, :, time + 1), &
            inverse, info, tolerance)
         if (info /= 0) then
            out%info = 20 + time
            return
         end if
         gain = matmul(matmul(filtered%p_filt(:, :, time), &
            transpose(model%transition(:, :, &
            time_index(model%transition, time)))), inverse)
         out%state(:, time) = filtered%a_filt(:, time) + &
            matmul(gain, out%state(:, time + 1) - &
            filtered%a_pred(:, time + 1))
      end do
   end function kfs_fast_smooth

   pure function kfs_conditional_covariance(covariance, lag_one_covariance, &
      tolerance) result(out)
      !! Convert smoothed joint covariances to forward conditionals.
      real(dp), intent(in) :: covariance(:, :, :) !! Covariance matrix.
      real(dp), intent(in) :: lag_one_covariance(:, :, :) !! Lag one covariance matrix.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(kfs_conditional_t) :: out
      real(dp), allocatable :: inverse(:, :)
      integer :: state, times, time, info

      state = size(covariance, 1)
      times = size(covariance, 3)
      if (state < 1 .or. times < 1 .or. size(covariance, 2) /= state .or. &
         any(shape(lag_one_covariance) /= &
         [state, state, max(0, times - 1)]) .or. &
         .not. all(ieee_is_finite(covariance)) .or. &
         .not. all(ieee_is_finite(lag_one_covariance))) then
         out%info = 1
         return
      end if
      allocate(out%matrix(state, state, times))
      allocate(out%covariance(state, state, times), inverse(state, state))
      out%matrix = 0.0_dp
      out%covariance(:, :, 1) = symmetrize(covariance(:, :, 1))
      do time = 2, times
         call symmetric_pseudoinverse(covariance(:, :, time - 1), inverse, &
            info, tolerance)
         if (info /= 0) then
            out%info = 10 + time
            return
         end if
         out%matrix(:, :, time) = &
            matmul(transpose(lag_one_covariance(:, :, time - 1)), inverse)
         out%covariance(:, :, time) = covariance(:, :, time) - &
            matmul(out%matrix(:, :, time), &
            lag_one_covariance(:, :, time - 1))
         out%covariance(:, :, time) = &
            symmetrize(out%covariance(:, :, time))
      end do
   end function kfs_conditional_covariance

   pure function kfs_loglik(model, tolerance) result(out)
      !! Return the Gaussian log likelihood and filter status.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(kfs_loglik_result_t) :: out
      type(kfs_filter_t) :: filtered

      filtered = kfs_filter(model, tolerance)
      out%value = filtered%log_likelihood
      out%info = filtered%info
   end function kfs_loglik

   pure function kfs_predict(model, filtered, h) result(out)
      !! Propagate filtered state moments and observation moments h steps ahead.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      type(kfs_filter_t), intent(in) :: filtered !! Filtered.
      integer, intent(in) :: h !! H.
      type(kfs_prediction_t) :: out
      real(dp), allocatable :: pwork(:, :), noise(:, :)
      integer :: i, n, p, m, tt, tz, th, tr, tq

      if (filtered%info /= 0 .or. h < 1) then
         out%info = 1
         return
      end if
      n = size(model%y, 1)
      p = size(model%y, 2)
      m = size(model%a1)
      allocate(out%mean(h, p), out%covariance(p, p, h))
      allocate(out%state(m, h), out%state_covariance(m, m, h))
      out%state(:, 1) = filtered%a_filt(:, n)
      out%state_covariance(:, :, 1) = filtered%p_filt(:, :, n)
      do i = 1, h
         tt = future_index(model%transition, n + i - 1)
         tr = future_index(model%r, n + i - 1)
         tq = future_index(model%q, n + i - 1)
         out%state(:, i) = matmul(model%transition(:, :, tt), out%state(:, i))
         pwork = matmul(matmul(model%transition(:, :, tt), out%state_covariance(:, :, i)), &
            transpose(model%transition(:, :, tt)))
         noise = matmul(matmul(model%r(:, :, tr), model%q(:, :, tq)), transpose(model%r(:, :, tr)))
         out%state_covariance(:, :, i) = symmetrize(pwork + noise)
         tz = future_index(model%z, n + i)
         th = future_index(model%h, n + i)
         out%mean(i, :) = matmul(model%z(:, :, tz), out%state(:, i))
         out%covariance(:, :, i) = matmul(matmul(model%z(:, :, tz), out%state_covariance(:, :, i)), &
            transpose(model%z(:, :, tz))) + model%h(:, :, th)
         if (i < h) then
            out%state(:, i + 1) = out%state(:, i)
            out%state_covariance(:, :, i + 1) = out%state_covariance(:, :, i)
         end if
      end do
   end function kfs_predict

   pure function kfs_standardized_innovations(filtered) result(residuals)
      !! Standardize one-step innovations by their marginal standard deviations.
      type(kfs_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), allocatable :: residuals(:, :)
      integer :: t, i

      allocate(residuals(size(filtered%innovation, 1), size(filtered%innovation, 2)))
      residuals = ieee_value(0.0_dp, ieee_quiet_nan)
      do t = 1, size(residuals, 1)
         do i = 1, size(residuals, 2)
            if (.not. ieee_is_finite(filtered%innovation(t, i))) cycle
            if (filtered%innovation_cov(i, i, t) <= 0.0_dp) cycle
            residuals(t, i) = filtered%innovation(t, i)/sqrt(filtered%innovation_cov(i, i, t))
         end do
      end do
   end function kfs_standardized_innovations

   pure function kfs_disturbance_smooth(model, filtered, tolerance) result(out)
      !! Compute conditional means of Gaussian observation and state disturbances.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      type(kfs_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(kfs_disturbance_t) :: out
      real(dp), allocatable :: rnext(:), rtrans(:), zobs(:, :), hobs(:, :)
      real(dp), allocatable :: f(:, :), finv(:, :), v(:), kraw(:, :), u(:), eps(:)
      real(dp) :: logdet, tol
      integer, allocatable :: observed(:)
      integer :: t, i, j, n, p, m, rr, no, info

      if (filtered%info /= 0) then
         out%info = filtered%info
         return
      end if
      n = size(model%y, 1)
      p = size(model%y, 2)
      m = size(model%a1)
      rr = size(model%q, 1)
      tol = 100.0_dp*epsilon(1.0_dp)
      if (present(tolerance)) tol = tolerance
      allocate(out%observation(n, p), out%state(n, rr), rnext(m), rtrans(m))
      out%observation = ieee_value(0.0_dp, ieee_quiet_nan)
      out%state = 0.0_dp
      rnext = 0.0_dp
      do t = n, 1, -1
         out%state(t, :) = matmul(model%q(:, :, time_index(model%q, t)), &
            matmul(transpose(model%r(:, :, time_index(model%r, t))), rnext))
         rtrans = matmul(transpose(model%transition(:, :, time_index(model%transition, t))), rnext)
         observed = pack([(i, i=1, p)], [(is_observed(model, t, i), i=1, p)])
         no = size(observed)
         if (no > 0) then
            allocate(zobs(no, m), hobs(no, no), f(no, no), finv(no, no), v(no), kraw(m, no), u(no), eps(no))
            do i = 1, no
               zobs(i, :) = model%z(observed(i), :, time_index(model%z, t))
               v(i) = filtered%innovation(t, observed(i))
               do j = 1, no
                  hobs(i, j) = model%h(observed(i), observed(j), time_index(model%h, t))
                  f(i, j) = filtered%innovation_cov(observed(i), observed(j), t)
               end do
            end do
            call inverse_logdet(f, finv, logdet, info, tol)
            if (info /= 0) then
               out%info = 30 + t
               return
            end if
            kraw = matmul(filtered%p_pred(:, :, t), transpose(zobs))
            u = matmul(finv, v - matmul(transpose(kraw), rtrans))
            eps = matmul(hobs, u)
            do i = 1, no
               out%observation(t, observed(i)) = eps(i)
            end do
            rnext = matmul(transpose(zobs), matmul(finv, v)) + &
               matmul(transpose(identity_matrix(m) - matmul(matmul(kraw, finv), zobs)), rtrans)
            deallocate(zobs, hobs, f, finv, v, kraw, u, eps)
         else
            rnext = rtrans
         end if
      end do
   end function kfs_disturbance_smooth

   pure elemental logical function valid_time_extent(extent, n)
      !! Report whether a system matrix is constant or available for every time.
      integer, intent(in) :: extent !! Extent.
      integer, intent(in) :: n !! Number of observations or elements.
      valid_time_extent = extent == 1 .or. extent >= n
   end function valid_time_extent

   pure integer function time_index(x, t) result(index)
      !! Select the constant or time-varying slice for an in-sample time.
      real(dp), intent(in) :: x(:, :, :) !! Input data or predictor values.
      integer, intent(in) :: t !! T.
      index = min(t, size(x, 3))
   end function time_index

   pure integer function future_index(x, t) result(index)
      !! Select an available matrix slice, holding the final slice fixed if needed.
      real(dp), intent(in) :: x(:, :, :) !! Input data or predictor values.
      integer, intent(in) :: t !! T.
      index = min(t, size(x, 3))
   end function future_index

   pure logical function is_observed(model, t, i)
      !! Test whether one observation should enter the measurement update.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: t !! T.
      integer, intent(in) :: i !! I.
      is_observed = ieee_is_finite(model%y(t, i))
      if (allocated(model%missing)) is_observed = is_observed .and. .not. model%missing(t, i)
   end function is_observed






end module kfas_mod
