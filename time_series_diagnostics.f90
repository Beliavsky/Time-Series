! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Shared diagnostics adapted from translated GPL time-series packages.
! Shared residual diagnostic tests for time-series models.
module time_series_diagnostics_mod
   use kind_mod, only: dp
   use forecast_mod, only: acf_values, pacf_values
   use itsmr_mod, only: regularized_gamma_q
   use time_series_linalg_mod, only: invert_matrix
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: box_test_box_pierce = 1
   integer, parameter, public :: box_test_ljung_box = 2
   integer, parameter, public :: box_test_monti = 3
   integer, parameter, public :: residual_raw = 0
   integer, parameter, public :: residual_squared = 1
   integer, parameter, public :: residual_log_squared = 2
   integer, parameter, public :: residual_absolute = 3

   type, public :: weighted_box_test_t
      ! Portmanteau statistic, reference distribution, and selected variant.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: gamma_shape = 0.0_dp
      real(dp) :: gamma_scale = 0.0_dp
      integer :: observations = 0
      integer :: lag = 0
      integer :: fitted_parameters = 0
      integer :: degrees_of_freedom = 0
      integer :: test_type = box_test_box_pierce
      integer :: transform = residual_raw
      integer :: info = 0
      logical :: weighted = .true.
   end type weighted_box_test_t

   type, public :: multivariate_white_noise_test_t
      ! Univariate and multivariate residual serial-correlation tests.
      real(dp), allocatable :: q_statistic(:), q_p_value(:)
      real(dp), allocatable :: lm_statistic(:), lm_p_value(:)
      real(dp) :: multivariate_q_statistic = 0.0_dp
      real(dp) :: multivariate_q_p_value = 1.0_dp
      integer :: observations = 0
      integer :: variables = 0
      integer :: lag = 0
      integer :: info = 0
   end type multivariate_white_noise_test_t

   public :: weighted_box_test
   public :: multivariate_white_noise_test

contains

   pure function weighted_box_test(residuals, lag, test_type, &
      fitted_parameter_count, weighted, transform) result(out)
      ! Compute weighted or classical Box-Pierce, Ljung-Box, and Monti tests.
      real(dp), intent(in) :: residuals(:)
      integer, intent(in) :: lag
      integer, intent(in), optional :: test_type, fitted_parameter_count, transform
      logical, intent(in), optional :: weighted
      type(weighted_box_test_t) :: out
      real(dp), allocatable :: working(:), correlation(:), observed(:)
      real(dp) :: denominator, sample_mean, weight
      integer :: i, requested_fit

      out%observations = size(residuals)
      out%lag = lag
      if (present(test_type)) out%test_type = test_type
      if (present(fitted_parameter_count)) &
         out%fitted_parameters = fitted_parameter_count
      if (present(weighted)) out%weighted = weighted
      if (present(transform)) out%transform = transform
      requested_fit = out%fitted_parameters
      if (lag < 1 .or. lag >= out%observations .or. requested_fit < 0 .or. &
         lag < 3*requested_fit - 1 .or. &
         out%test_type < box_test_box_pierce .or. &
         out%test_type > box_test_monti .or. &
         out%transform < residual_raw .or. &
         out%transform > residual_absolute) then
         out%info = 1
         return
      end if
      if (.not. all(ieee_is_finite(residuals))) then
         out%info = 2
         return
      end if
      working = residuals
      select case (out%transform)
      case (residual_squared, residual_log_squared)
         working = working**2
         if (out%transform == residual_log_squared) then
            if (any(working <= 0.0_dp)) then
               out%info = 2
               return
            end if
            working = log(working)
         end if
         out%fitted_parameters = 0
      case (residual_absolute)
         working = abs(working)
         out%fitted_parameters = 0
      end select
      if (.not. all(ieee_is_finite(working))) then
         out%info = 2
         return
      end if
      sample_mean = sum(working)/real(out%observations, dp)
      if (sum((working - sample_mean)**2) <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      if (out%test_type == box_test_monti) then
         observed = pacf_values(working, lag)
      else
         correlation = acf_values(working, lag)
         observed = correlation(2:lag + 1)
      end if
      if (.not. all(ieee_is_finite(observed))) then
         out%info = 2
         return
      end if
      do i = 1, lag
         weight = 1.0_dp
         if (out%weighted) weight = real(lag - i + 1, dp)/real(lag, dp)
         if (out%test_type == box_test_box_pierce) then
            out%statistic = out%statistic + weight*observed(i)**2
         else
            out%statistic = out%statistic + weight*observed(i)**2/ &
               real(out%observations - i, dp)
         end if
      end do
      if (out%test_type == box_test_box_pierce) then
         out%statistic = real(out%observations, dp)*out%statistic
      else
         out%statistic = real(out%observations, dp)* &
            real(out%observations + 2, dp)*out%statistic
      end if
      out%degrees_of_freedom = lag - out%fitted_parameters
      if (out%weighted) then
         denominator = real(2*lag*lag + 3*lag + 1 - &
            6*lag*out%fitted_parameters, dp)
         if (denominator <= 0.0_dp) then
            out%info = 1
            return
         end if
         out%gamma_shape = 0.75_dp*real((lag + 1)**2*lag, dp)/denominator
         out%gamma_scale = (2.0_dp/3.0_dp)*denominator/real(lag*(lag + 1), dp)
      else
         if (out%degrees_of_freedom <= 0) then
            out%info = 1
            return
         end if
         out%gamma_shape = 0.5_dp*real(out%degrees_of_freedom, dp)
         out%gamma_scale = 2.0_dp
      end if
      out%p_value = regularized_gamma_q(out%gamma_shape, &
         out%statistic/out%gamma_scale)
      out%p_value = max(0.0_dp, min(1.0_dp, out%p_value))
   end function weighted_box_test

   pure function multivariate_white_noise_test(residuals, lag) result(out)
      ! Compute FCVAR-style Q and heteroskedasticity-robust LM tests.
      real(dp), intent(in) :: residuals(:, :)
      integer, intent(in) :: lag
      type(multivariate_white_noise_test_t) :: out
      real(dp) :: statistic, probability
      integer :: variable, status

      out%observations = size(residuals, 1)
      out%variables = size(residuals, 2)
      out%lag = lag
      allocate(out%q_statistic(out%variables), out%q_p_value(out%variables))
      allocate(out%lm_statistic(out%variables), out%lm_p_value(out%variables))
      out%q_statistic = 0.0_dp
      out%q_p_value = 1.0_dp
      out%lm_statistic = 0.0_dp
      out%lm_p_value = 1.0_dp
      if (out%observations < 2 .or. out%variables < 1 .or. lag < 1 .or. &
         lag >= out%observations .or. .not. all(ieee_is_finite(residuals))) then
         out%info = 1
         return
      end if
      do variable = 1, out%variables
         call multivariate_q_statistic(residuals(:, variable:variable), lag, &
            statistic, probability, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%q_statistic(variable) = statistic
         out%q_p_value(variable) = probability
         call robust_lm_statistic(residuals(:, variable), lag, statistic, &
            probability, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         out%lm_statistic(variable) = statistic
         out%lm_p_value(variable) = probability
      end do
      call multivariate_q_statistic(residuals, lag, statistic, probability, status)
      if (status /= 0) then
         out%info = 30 + status
         return
      end if
      out%multivariate_q_statistic = statistic
      out%multivariate_q_p_value = probability
   end function multivariate_white_noise_test

   pure subroutine multivariate_q_statistic(residuals, lag, statistic, &
      p_value, info)
      ! Compute the multivariate Ljung-Box statistic of Luetkepohl.
      real(dp), intent(in) :: residuals(:, :)
      integer, intent(in) :: lag
      real(dp), intent(out) :: statistic, p_value
      integer, intent(out) :: info
      real(dp), allocatable :: covariance(:, :), inverse(:, :), cross(:, :)
      real(dp), allocatable :: product(:, :)
      integer :: observations, variables, current, i

      observations = size(residuals, 1)
      variables = size(residuals, 2)
      covariance = matmul(transpose(residuals), residuals)/real(observations, dp)
      call invert_matrix(covariance, inverse, info)
      statistic = 0.0_dp
      p_value = 1.0_dp
      if (info /= 0) return
      do current = 1, lag
         cross = matmul(transpose(residuals(current + 1:, :)), &
            residuals(:observations - current, :))/ &
            real(observations - current, dp)
         product = matmul(matmul(transpose(cross), inverse), &
            matmul(cross, inverse))
         do i = 1, variables
            statistic = statistic + product(i, i)/real(observations - current, dp)
         end do
      end do
      statistic = statistic*real(observations*(observations + 2), dp)
      p_value = regularized_gamma_q(0.5_dp*real(variables*variables*lag, dp), &
         0.5_dp*statistic)
      p_value = max(0.0_dp, min(1.0_dp, p_value))
   end subroutine multivariate_q_statistic

   pure subroutine robust_lm_statistic(residuals, lag, statistic, p_value, info)
      ! Compute FCVAR's heteroskedasticity-consistent serial-correlation LM test.
      real(dp), intent(in) :: residuals(:)
      integer, intent(in) :: lag
      real(dp), intent(out) :: statistic, p_value
      integer, intent(out) :: info
      real(dp), allocatable :: centered(:), scores(:, :), covariance(:, :)
      real(dp), allocatable :: inverse(:, :), score_mean(:)
      integer :: observations, rows, current

      observations = size(residuals)
      rows = observations - lag
      allocate(centered(observations), scores(rows, lag), score_mean(lag))
      centered = residuals - sum(residuals)/real(observations, dp)
      do current = 1, lag
         scores(:, current) = centered(lag + 1:)* &
            centered(lag + 1 - current:observations - current)
      end do
      score_mean = sum(scores, dim=1)/real(rows, dp)
      scores = scores - spread(score_mean, 1, rows)
      covariance = matmul(transpose(scores), scores)/real(observations, dp)
      call invert_matrix(covariance, inverse, info)
      statistic = 0.0_dp
      p_value = 1.0_dp
      if (info /= 0) return
      statistic = real(observations, dp)*dot_product(score_mean, &
         matmul(inverse, score_mean))
      p_value = regularized_gamma_q(0.5_dp*real(lag, dp), 0.5_dp*statistic)
      p_value = max(0.0_dp, min(1.0_dp, p_value))
   end subroutine robust_lm_statistic

end module time_series_diagnostics_mod
