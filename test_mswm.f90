! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Numerical tests for algorithms translated from R MSwM.
program test_mswm
   use kind_mod, only: dp
   use mswm_mod
   use random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   integer, parameter :: observations = 600
   real(dp) :: response(observations), design(observations, 2)
   real(dp) :: initial_coefficients(2, 2), initial_sd(2), initial_transition(2, 2)
   logical :: switching(2)
   type(mswm_filter_t) :: filtered
   type(mswm_fit_t) :: fit
   type(mswm_fit_t) :: ar_fit
   type(mswm_glm_fit_t) :: glm_fit, diagnostic_fit, glm_ar_fit
   type(mswm_ar_data_t) :: ar_data
   type(mswm_inference_t) :: inference
   type(mswm_intervals_t) :: intervals
   type(mswm_multistart_t) :: multistart, random_multistart, repeated_multistart
   type(mswm_glm_multistart_t) :: glm_multistart
   type(mswm_glm_multistart_t) :: glm_random_multistart
   type(mswm_glm_multistart_t) :: glm_repeated_multistart
   integer, allocatable :: state(:)
   real(dp), allocatable :: residual(:)
   real(dp) :: count_response(observations)
   real(dp) :: binary_response(observations)
   real(dp) :: gamma_response(observations)
   real(dp) :: glm_coefficients(2, 2), probability
   real(dp) :: start_coefficients(2, 2, 2)
   real(dp) :: start_deviation(2, 2), start_transition(2, 2, 2)
   real(dp) :: glm_start_coefficients(2, 2, 2)
   real(dp) :: glm_start_transition(2, 2, 2)
   integer :: time, active_state, unit

   design(:, 1) = 1.0_dp
   do time = 1, observations
      design(time, 2) = sin(0.037_dp*real(time, dp))
      active_state = merge(1, 2, mod((time - 1)/75, 2) == 0)
      if (active_state == 1) then
         response(time) = -1.2_dp + 0.5_dp*design(time, 2) + &
            0.18_dp*sin(1.71_dp*real(time, dp))
      else
         response(time) = 1.0_dp + 0.5_dp*design(time, 2) + &
            0.28_dp*cos(1.33_dp*real(time, dp))
      end if
   end do
   switching = [.true., .false.]
   initial_coefficients = reshape([-1.0_dp, 1.0_dp, 0.4_dp, 0.4_dp], [2, 2])
   initial_sd = [0.25_dp, 0.35_dp]
   initial_transition = reshape([0.97_dp, 0.03_dp, 0.03_dp, 0.97_dp], [2, 2])

   ar_data = mswm_ar_data(response, 2)
   call check(ar_data%info == 0, 'autoregressive data construction')
   call check(size(ar_data%response) == observations - 2 .and. &
      all(shape(ar_data%design) == [observations - 2, 3]), &
      'autoregressive data dimensions')
   call check(abs(ar_data%design(1, 2) - response(2)) < 1.0e-12_dp .and. &
      abs(ar_data%design(1, 3) - response(1)) < 1.0e-12_dp, &
      'autoregressive lag alignment')
   ar_fit = mswm_gaussian_ar_fit(response, 2, 1, &
      [.false., .false., .false.], .false., max_iterations=100, &
      tolerance=1.0e-8_dp)
   call check(ar_fit%info == 0, 'Gaussian autoregressive switching wrapper')
   call check(ar_fit%lag_order == 2 .and. &
      maxval(abs(ar_fit%terminal_response - response(observations - 1:))) < &
      1.0e-12_dp, 'Gaussian autoregressive terminal history')

   filtered = mswm_gaussian_filter(response, design, initial_coefficients, &
      initial_sd, initial_transition, [0.5_dp, 0.5_dp])
   call check(filtered%info == 0, 'Hamilton filter status')
   call check(maxval(abs(sum(filtered%filtered_probability, dim=2) - 1.0_dp)) < &
      1.0e-12_dp, 'filtered probability normalization')

   fit = mswm_gaussian_fit(response, design, 2, switching, .true., &
      initial_coefficients, initial_sd, initial_transition, 200, 1.0e-9_dp)
   call check(fit%info == 0, 'Gaussian switching fit status')
   call check(fit%converged, 'Gaussian switching fit convergence')
   call check(maxval(abs(sum(fit%transition, dim=2) - 1.0_dp)) < 1.0e-12_dp, &
      'transition probability normalization')
   call check(maxval(abs(sum(fit%smoothed_probability, dim=2) - 1.0_dp)) < &
      1.0e-12_dp, 'smoothed probability normalization')
   call check(abs(fit%coefficients(1, 1) + 1.2_dp) < 0.12_dp, &
      'first regime intercept')
   call check(abs(fit%coefficients(2, 1) - 1.0_dp) < 0.12_dp, &
      'second regime intercept')
   call check(maxval(abs(fit%coefficients(:, 2) - 0.5_dp)) < 0.08_dp, &
      'shared slope estimate')
   call check(fit%transition(1, 1) > 0.90_dp .and. &
      fit%transition(2, 2) > 0.90_dp, 'persistent transition estimates')

   state = mswm_states(fit)
   residual = mswm_residuals(fit)
   call check(size(state) == observations, 'decoded state size')
   call check(size(residual) == observations, 'weighted residual size')
   call check(count(state == 1) > 200 .and. count(state == 2) > 200, &
      'both decoded regimes represented')

   inference = mswm_gaussian_inference(fit, response, design, 2.0e-4_dp)
   call check(inference%info == 0, 'Gaussian switching Hessian inference')
   call check(all(ieee_is_finite(inference%coefficient_standard_error)), &
      'finite coefficient standard errors')
   call check(all(inference%coefficient_standard_error >= 0.0_dp), &
      'nonnegative coefficient standard errors')
   intervals = mswm_gaussian_intervals(fit, inference, 0.90_dp)
   call check(intervals%info == 0, 'Gaussian switching intervals')
   call check(all(intervals%coefficient_lower <= fit%coefficients .and. &
      fit%coefficients <= intervals%coefficient_upper), &
      'coefficient interval containment')

   start_coefficients(:, :, 1) = initial_coefficients
   start_coefficients(:, :, 2) = initial_coefficients
   start_coefficients(:, 1, 2) = [-0.7_dp, 0.7_dp]
   start_deviation(:, 1) = initial_sd
   start_deviation(:, 2) = [0.4_dp, 0.4_dp]
   start_transition(:, :, 1) = initial_transition
   start_transition(:, :, 2) = reshape([ &
      0.90_dp, 0.10_dp, 0.10_dp, 0.90_dp], [2, 2])
   multistart = mswm_gaussian_multistart(response, design, switching, .true., &
      start_coefficients, start_deviation, start_transition, 150, 1.0e-8_dp)
   call check(multistart%info == 0 .and. multistart%successful == 2, &
      'supplied Gaussian multistart')
   call check(multistart%best >= 1 .and. multistart%best <= 2, &
      'supplied Gaussian best start')

   call set_random_seed(7219)
   random_multistart = mswm_gaussian_random_starts(response, design, 2, &
      switching, .true., 3, 150, 1.0e-8_dp)
   call set_random_seed(7219)
   repeated_multistart = mswm_gaussian_random_starts(response, design, 2, &
      switching, .true., 3, 150, 1.0e-8_dp)
   call check(random_multistart%info == 0, 'random Gaussian multistart')
   call check(random_multistart%best == repeated_multistart%best, &
      'random Gaussian multistart best reproducibility')
   call check(abs(random_multistart%fits(random_multistart%best)%log_likelihood - &
      repeated_multistart%fits(repeated_multistart%best)%log_likelihood) < &
      1.0e-12_dp, 'random Gaussian multistart likelihood reproducibility')

   glm_coefficients = reshape([0.20_dp, 1.40_dp, 0.25_dp, 0.25_dp], [2, 2])
   do time = 1, observations
      active_state = merge(1, 2, mod((time - 1)/75, 2) == 0)
      probability = exp(glm_coefficients(active_state, 1) + &
         glm_coefficients(active_state, 2)*design(time, 2))
      count_response(time) = real(max(0, nint(probability + &
         0.65_dp*sqrt(probability)*sin(1.17_dp*real(time, dp)))), dp)
      probability = 1.0_dp/(1.0_dp + exp(-merge(-1.4_dp, 1.4_dp, &
         active_state == 1)))
      binary_response(time) = merge(1.0_dp, 0.0_dp, &
         modulo(37*time, 100) < nint(100.0_dp*probability))
      gamma_response(time) = exp(merge(-0.3_dp, 0.5_dp, active_state == 1) + &
         0.15_dp*design(time, 2))*(1.0_dp + &
         0.12_dp*sin(0.83_dp*real(time, dp)))
   end do

   filtered = mswm_glm_filter(count_response, design, glm_coefficients, &
      initial_transition, [0.5_dp, 0.5_dp], 'poisson')
   call check(filtered%info == 0, 'Poisson Hamilton filter status')
   call check(maxval(abs(sum(filtered%filtered_probability, dim=2) - 1.0_dp)) < &
      1.0e-12_dp, 'Poisson filtered probability normalization')

   glm_fit = mswm_glm_fit(count_response, design, 2, switching, 'poisson', &
      glm_coefficients, initial_transition, 200, 1.0e-8_dp, 30)
   call check(glm_fit%info == 0, 'Poisson switching fit status')
   call check(glm_fit%converged, 'Poisson switching fit convergence')
   call check(abs(glm_fit%coefficients(1, 1) - 0.20_dp) < 0.35_dp, &
      'first Poisson regime intercept')
   call check(abs(glm_fit%coefficients(2, 1) - 1.40_dp) < 0.35_dp, &
      'second Poisson regime intercept')
   call check(maxval(abs(glm_fit%coefficients(:, 2) - &
      glm_fit%coefficients(1, 2))) < 1.0e-12_dp, 'shared Poisson slope')
   state = mswm_states(glm_fit)
   call check(size(state) == observations .and. all(state >= 1 .and. state <= 2), &
      'decoded Poisson states')
   inference = mswm_glm_inference(glm_fit, count_response, design, 2.0e-4_dp)
   call check(inference%info == 0, 'Poisson switching Hessian inference')
   call check(size(inference%deviation_standard_error) == 0, &
      'GLM inference has no Gaussian scale errors')
   call check(all(ieee_is_finite(inference%coefficient_standard_error)), &
      'finite Poisson coefficient standard errors')
   call check(all(inference%coefficient_standard_error >= 0.0_dp), &
      'nonnegative Poisson coefficient standard errors')
   intervals = mswm_glm_intervals(glm_fit, inference, 0.90_dp)
   call check(intervals%info == 0, 'Poisson switching intervals')
   call check(size(intervals%deviation_lower) == 0, &
      'GLM intervals have no Gaussian scale limits')
   call check(all(intervals%coefficient_lower <= glm_fit%coefficients .and. &
      glm_fit%coefficients <= intervals%coefficient_upper), &
      'Poisson coefficient interval containment')
   residual = mswm_residuals(glm_fit, 1)
   call check(maxval(abs(residual - glm_fit%residuals(:, 1)/ &
      sqrt(glm_fit%conditional_mean(:, 1)))) < 1.0e-12_dp, &
      'regime-specific Poisson residuals')
   residual = mswm_residuals(glm_fit)
   call check(maxval(abs(residual - sum(glm_fit%smoothed_probability* &
      glm_fit%residuals/sqrt(glm_fit%conditional_mean), dim=2))) < 1.0e-12_dp, &
      'probability-weighted Poisson residuals')

   diagnostic_fit%conditional_mean = reshape([0.25_dp, 0.80_dp], [2, 1])
   diagnostic_fit%residuals = reshape([0.75_dp, -0.80_dp], [2, 1])
   diagnostic_fit%smoothed_probability = reshape([1.0_dp, 1.0_dp], [2, 1])
   diagnostic_fit%family = 'binomial'
   residual = mswm_residuals(diagnostic_fit)
   call check(maxval(abs(residual - [0.75_dp/sqrt(0.1875_dp), &
      -0.80_dp/sqrt(0.16_dp)])) < 1.0e-12_dp, &
      'standardized binomial residuals')
   diagnostic_fit%conditional_mean = reshape([2.0_dp, 4.0_dp], [2, 1])
   diagnostic_fit%residuals = reshape([1.0_dp, -2.0_dp], [2, 1])
   diagnostic_fit%family = 'gamma_log'
   residual = mswm_residuals(diagnostic_fit)
   call check(maxval(abs(residual - [0.5_dp, -0.5_dp])) < 1.0e-12_dp, &
      'standardized Gamma residuals')

   glm_ar_fit = mswm_glm_ar_fit(count_response, 1, 1, &
      [.false., .false.], 'poisson', max_iterations=100, &
      tolerance=1.0e-8_dp, irls_iterations=30)
   call check(glm_ar_fit%info == 0, 'Poisson autoregressive switching wrapper')
   call check(glm_ar_fit%lag_order == 1 .and. &
      abs(glm_ar_fit%terminal_response(1) - count_response(observations)) < &
      1.0e-12_dp, 'Poisson autoregressive terminal history')

   glm_start_coefficients(:, :, 1) = glm_coefficients
   glm_start_coefficients(:, :, 2) = glm_coefficients
   glm_start_coefficients(:, 1, 2) = [0.0_dp, 1.1_dp]
   glm_start_transition(:, :, 1) = initial_transition
   glm_start_transition(:, :, 2) = reshape([ &
      0.90_dp, 0.10_dp, 0.10_dp, 0.90_dp], [2, 2])
   glm_multistart = mswm_glm_multistart(count_response, design, switching, &
      'poisson', glm_start_coefficients, glm_start_transition, 150, 1.0e-8_dp, 30)
   call check(glm_multistart%info == 0 .and. glm_multistart%successful == 2, &
      'supplied Poisson multistart')
   call set_random_seed(8137)
   glm_random_multistart = mswm_glm_random_starts(count_response, design, 2, &
      switching, 'poisson', 3, 150, 1.0e-8_dp, 30)
   call set_random_seed(8137)
   glm_repeated_multistart = mswm_glm_random_starts(count_response, design, 2, &
      switching, 'poisson', 3, 150, 1.0e-8_dp, 30)
   call check(glm_random_multistart%info == 0, 'random Poisson multistart')
   call check(glm_random_multistart%best == glm_repeated_multistart%best, &
      'random Poisson multistart best reproducibility')
   call check(abs(glm_random_multistart%fits(glm_random_multistart%best)% &
      log_likelihood - glm_repeated_multistart%fits( &
      glm_repeated_multistart%best)%log_likelihood) < 1.0e-12_dp, &
      'random Poisson multistart likelihood reproducibility')

   filtered = mswm_glm_filter(binary_response, design, &
      reshape([-1.4_dp, 1.4_dp, 0.0_dp, 0.0_dp], [2, 2]), &
      initial_transition, [0.5_dp, 0.5_dp], 'binomial')
   call check(filtered%info == 0, 'binomial Hamilton filter status')
   call check(maxval(abs(sum(filtered%filtered_probability, dim=2) - 1.0_dp)) < &
      1.0e-12_dp, 'binomial filtered probability normalization')

   filtered = mswm_glm_filter(gamma_response, design, &
      reshape([-0.3_dp, 0.5_dp, 0.15_dp, 0.15_dp], [2, 2]), &
      initial_transition, [0.5_dp, 0.5_dp], 'gamma_log')
   call check(filtered%info == 0, 'Gamma Hamilton filter status')
   call check(maxval(abs(sum(filtered%filtered_probability, dim=2) - 1.0_dp)) < &
      1.0e-12_dp, 'Gamma filtered probability normalization')

   open(newunit=unit, status='scratch', action='write')
   call display(filtered, unit)
   call display(fit, unit)
   call display(glm_fit, unit)
   call display(ar_data, unit)
   call display(inference, unit)
   call display(intervals, unit)
   call display(multistart, unit)
   call display(glm_multistart, unit)
   close(unit)

   print '(a)', 'MSwM tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_mswm
