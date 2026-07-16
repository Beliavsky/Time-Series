! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from astsa with MIT-licensed SV adaptations.
! Numerical algorithms and compatibility adapters from the R package astsa.
! The SV filtering structure is adapted from GARCH-BFGS/sv.f90 under the MIT license.
module astsa_mod
   use kind_mod, only: dp
   use time_series_linalg_mod, only: symmetric_eigen, inverse_logdet, invert_matrix, symmetrize, identity_matrix
   use time_series_stats_mod, only: ols_fit, yule_walker_result_t, yule_walker_fit
   use time_series_utils_mod, only: inverse_standard_normal
   use time_series_random_mod, only: random_standard_normal_matrix, random_uniform, random_gamma, &
      multivariate_normal_from_standard, random_multivariate_normal
   use kfas_mod, only: ssm_model_t, kfs_filter_t, kfs_smoother_t, kfs_filter, kfs_filter_diffuse, kfs_smooth
   use forecast_mod, only: acf_values, ccf_values
   use time_series_optimization_mod, only: optimization_result_t, bfgs_minimize_fd, &
      finite_difference_hessian
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private

   type, public :: astsa_filter_t
      ! astsa-style predicted and filtered state moments and innovations.
      real(dp), allocatable :: xp(:, :), pp(:, :, :), xf(:, :), pf(:, :, :)
      real(dp), allocatable :: innovation(:, :), innovation_cov(:, :, :)
      real(dp) :: likelihood = 0.0_dp
      integer :: info = 0
   end type

   type, public :: astsa_smoother_t
      ! astsa-style smoothed states together with forward-filter output.
      type(astsa_filter_t) :: filter
      real(dp), allocatable :: xs(:, :), ps(:, :, :), x0n(:), p0n(:, :)
      real(dp), allocatable :: j(:, :, :), j0(:, :)
      integer :: info = 0
   end type

   type, public :: astsa_em_result_t
      ! EM estimates for an astsa linear Gaussian state-space model.
      real(dp), allocatable :: transition(:, :), q(:, :), r(:, :), mu0(:), sigma0(:, :)
      real(dp), allocatable :: likelihood(:)
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type

   type, public :: astsa_ffbs_result_t
      ! One conditional state trajectory drawn by forward-filter backward-sampling.
      real(dp), allocatable :: states(:, :), initial_state(:)
      integer :: info = 0
   end type

   type, public :: astsa_spectrum_t
      ! Graphics-independent univariate or multivariate spectral estimates.
      real(dp), allocatable :: frequency(:), spectrum(:, :)
      complex(dp), allocatable :: spectral_matrix(:, :, :)
      real(dp), allocatable :: coherence(:, :, :), phase(:, :, :)
      real(dp) :: bandwidth = 0.0_dp
      real(dp) :: degrees_freedom = 0.0_dp
      real(dp) :: sampling_frequency = 1.0_dp
      integer :: original_size = 0
      integer :: used_size = 0
      integer :: info = 0
   end type

   type, public :: astsa_matrix_power_result_t
      ! Symmetric matrix power and its eigendecomposition status.
      real(dp), allocatable :: value(:, :)
      integer :: info = 0
   end type

   type, public :: astsa_simulation_t
      ! Simulated SARIMA series, seasonal frequency, and status code.
      real(dp), allocatable :: series(:)
      integer :: frequency = 1
      integer :: info = 0
   end type

   type :: astsa_root_check_t
      ! Internal polynomial stability result.
      logical :: outside_unit_circle = .true.
      logical :: converged = .true.
   end type

   type :: astsa_polynomial_roots_t
      ! Internal complex polynomial roots and convergence status.
      complex(dp), allocatable :: roots(:)
      logical :: converged = .true.
   end type

   type, public :: astsa_arma_diagnostic_t
      ! ARMA roots, stability flags, and approximate common-factor diagnostics.
      complex(dp), allocatable :: ar_roots(:), ma_roots(:)
      complex(dp), allocatable :: seasonal_ar_roots(:), seasonal_ma_roots(:)
      real(dp) :: redundancy_tolerance = 0.1_dp
      integer :: seasonal_period = 1
      integer :: info = 0
      logical :: causal = .true.
      logical :: invertible = .true.
      logical :: redundant = .false.
      logical :: seasonal_redundant = .false.
      logical :: redundancy_checked = .false.
   end type

   type, public :: astsa_prewhite_result_t
      ! Selected AR whitening model, aligned series, and cross-correlations.
      real(dp), allocatable :: ar(:), first(:), second(:), correlation(:)
      integer, allocatable :: lags(:)
      integer :: order = 0
      integer :: differences = 0
      integer :: info = 0
   end type

   type, public :: astsa_sarima_likelihood_t
      ! Conditional SARIMA likelihood, innovations, fitted values, and criteria.
      real(dp), allocatable :: transformed(:), fitted(:), residuals(:)
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: aicc = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: observations = 0
      integer :: parameters = 0
      integer :: info = 0
   end type

   type, public :: astsa_sarima_fit_t
      ! Estimated SARIMA parameters, likelihood output, and optimizer status.
      real(dp), allocatable :: coefficients(:), covariance(:, :), standard_errors(:), statistics(:)
      type(astsa_sarima_likelihood_t) :: likelihood
      integer :: iterations = 0
      integer :: info = 0
      integer :: p = 0
      integer :: d = 0
      integer :: q = 0
      integer :: seasonal_p = 0
      integer :: seasonal_difference = 0
      integer :: seasonal_q = 0
      integer :: season = 1
      integer :: regression_count = 0
      logical :: includes_intercept = .false.
      logical :: includes_drift = .false.
      logical :: uses_exact_likelihood = .false.
      logical :: converged = .false.
   end type

   type, public :: astsa_sarima_forecast_t
      ! SARIMA forecast means, standard errors, and symmetric intervals.
      real(dp), allocatable :: mean(:), standard_error(:), lower(:), upper(:)
      real(dp) :: interval_multiplier = 2.0_dp
      integer :: info = 0
   end type

   type, public :: astsa_sarima_diagnostics_t
      ! Standardized residual, ACF, Ljung-Box, and normal Q-Q diagnostics.
      real(dp), allocatable :: standardized_residuals(:), residual_acf(:)
      real(dp), allocatable :: ljung_box(:), p_values(:), qq_theoretical(:), qq_sample(:)
      integer, allocatable :: lags(:), degrees_of_freedom(:)
      integer :: info = 0
   end type

   type, public :: astsa_lag_regression_t
      ! Frequency-domain lag coefficients and aligned regression output.
      integer, allocatable :: lags(:), selected_lags(:)
      real(dp), allocatable :: coefficients(:), selected_coefficients(:)
      real(dp), allocatable :: output(:), fitted(:), residuals(:)
      real(dp) :: intercept = 0.0_dp
      real(dp) :: mse = 0.0_dp
      integer :: info = 0
   end type

   type, public :: astsa_signal_extraction_t
      ! Two-sided ideal-band filter, filtered series, and frequency responses.
      integer, allocatable :: lags(:)
      real(dp), allocatable :: coefficients(:), filtered(:), frequency(:)
      real(dp), allocatable :: desired_response(:), attained_response(:)
      integer :: valid_start = 0
      integer :: valid_end = 0
      integer :: info = 0
   end type

   type, public :: astsa_stochastic_regression_t
      ! Spectral regression powers, impulse coefficients, F statistics, and coherence.
      real(dp), allocatable :: frequency(:), full_power(:), reduced_power(:)
      real(dp), allocatable :: coefficients(:, :), f_statistic(:), coherence(:)
      integer, allocatable :: lags(:)
      real(dp) :: critical_f = 0.0_dp
      real(dp) :: critical_coherence = 0.0_dp
      integer :: numerator_df = 0
      integer :: denominator_df = 0
      integer :: info = 0
   end type

   type, public :: astsa_sv_filter_t
      ! Mixture-filter predictions, variances, responsibilities, and objective.
      real(dp), allocatable :: predicted_log_variance(:), prediction_variance(:)
      real(dp), allocatable :: component_zero_probability(:), component_one_probability(:)
      real(dp) :: negative_log_likelihood = 0.0_dp
      integer :: info = 0
   end type

   type, public :: astsa_sv_fit_t
      ! Estimated astsa stochastic-volatility mixture parameters and inference.
      real(dp), allocatable :: coefficients(:), covariance(:, :), standard_errors(:)
      type(astsa_sv_filter_t) :: filter
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: feedback = .false.
      logical :: leverage = .false.
   end type

   type, public :: astsa_sv_particle_filter_t
      ! Conditional particles, weights, ancestors, selected trajectory, and likelihood.
      real(dp), allocatable :: particles(:, :), weights(:, :), trajectory(:)
      integer, allocatable :: ancestors(:, :)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: info = 0
   end type

   type, public :: astsa_sv_mcmc_t
      ! Retained stochastic-volatility parameter and latent-state draws.
      real(dp), allocatable :: phi(:), state_sd(:), observation_scale(:), latent(:, :), effective_size(:)
      real(dp) :: acceptance_rate = 0.0_dp
      integer :: info = 0
   end type

   type, public :: astsa_ssm_result_t
      ! Scalar state-space estimates, filtering, smoothing, and inference.
      real(dp), allocatable :: coefficients(:), covariance(:, :), standard_errors(:)
      real(dp), allocatable :: predicted_state(:), predicted_variance(:)
      real(dp), allocatable :: filtered_state(:), filtered_variance(:)
      real(dp), allocatable :: smoothed_state(:), smoothed_variance(:)
      real(dp) :: negative_log_likelihood = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: fixed_phi = .false.
   end type

   type, public :: astsa_ar_bootstrap_t
      ! Original and bootstrap AR coefficients, simulated series, and summaries.
      real(dp), allocatable :: original_coefficients(:), coefficient_draws(:, :)
      real(dp), allocatable :: simulated_series(:, :), coefficient_means(:), quantiles(:, :)
      real(dp), allocatable :: probabilities(:)
      integer :: info = 0
   end type

   type, public :: astsa_ar_mcmc_t
      ! Retained Bayesian AR coefficient and innovation-scale posterior draws.
      real(dp), allocatable :: coefficient_draws(:, :), innovation_sd(:)
      real(dp), allocatable :: means(:), standard_deviations(:), quantiles(:, :)
      real(dp), allocatable :: probabilities(:), effective_size(:)
      integer :: info = 0
   end type

   type, public :: astsa_spectrum_ic_t
      ! AR order criteria, selected model, and theoretical spectrum.
      integer, allocatable :: orders(:)
      real(dp), allocatable :: aic(:), bic(:), coefficients(:), frequency(:), spectrum(:)
      real(dp) :: innovation_variance = 0.0_dp
      integer :: selected_order = 0
      integer :: info = 0
   end type

   type, public :: astsa_linearity_test_t
      ! Block-bispectrum nonlinearity statistics and p-values.
      real(dp), allocatable :: frequency(:), statistic(:, :), p_values(:, :)
      complex(dp), allocatable :: bispectrum(:, :)
      real(dp) :: noncentrality = 0.0_dp
      integer :: block_length = 0
      integer :: block_count = 0
      integer :: info = 0
   end type

   public :: astsa_kfilter, astsa_ksmooth, astsa_kfilter_correlated, astsa_ksmooth_correlated
   public :: astsa_em
   public :: astsa_ffbs, astsa_ffbs_draws
   public :: arma_spectrum, mv_periodogram
   public :: sarima_sim, sarima_sim_from_innovations
   public :: arma_check
   public :: pre_white
   public :: sarima_likelihood, sarima_exact_likelihood
   public :: sarima_fit
   public :: sarima_forecast
   public :: sarima_diagnostics
   public :: lag_reg
   public :: signal_extract
   public :: stochastic_regression
   public :: sv_filter, sv_mle
   public :: sv_particle_filter_draws, sv_particle_filter
   public :: sv_mcmc
   public :: ssm_fit
   public :: ar_boot_draws, ar_boot
   public :: ar_mcmc_draws, ar_mcmc
   public :: spectrum_ic
   public :: test_linearity
   public :: poly_mul, arma_to_ar, symmetric_matrix_power, fdr_cutoff

contains

   pure function astsa_kfilter(y, observation, mu0, sigma0, transition, sq, sr) result(out)
      ! Run astsa Kfilter version 1 by adapting the model to kfas_mod.
      real(dp), intent(in) :: y(:, :), observation(:, :, :)
      real(dp), intent(in) :: mu0(:), sigma0(:, :), transition(:, :), sq(:, :), sr(:, :)
      type(astsa_filter_t) :: out
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      real(dp), allocatable :: process_cov(:, :), observation_cov(:, :)
      integer :: i, j, t, n, p, q, process_rank

      n = size(y, 1)
      q = size(y, 2)
      p = size(mu0)
      process_rank = size(sq, 2)
      if (size(observation, 1) /= q .or. size(observation, 2) /= p) then
         out%info = 1
         return
      end if
      if (size(observation, 3) /= 1 .and. size(observation, 3) < n) then
         out%info = 2
         return
      end if
      if (any(shape(sigma0) /= [p, p]) .or. any(shape(transition) /= [p, p])) then
         out%info = 3
         return
      end if
      if (size(sq, 1) /= p .or. size(sr, 1) /= q) then
         out%info = 4
         return
      end if
      process_cov = matmul(sq, transpose(sq))
      observation_cov = matmul(sr, transpose(sr))
      allocate(model%y(n, q), model%z(q, p, size(observation, 3)), model%h(q, q, 1))
      allocate(model%transition(p, p, 1), model%r(p, process_rank, 1))
      allocate(model%q(process_rank, process_rank, 1), model%a1(p), model%p1(p, p), model%p1inf(p, p))
      allocate(model%missing(n, q))
      model%y = y
      model%z = observation
      do t = 1, n
         do j = 1, q
            if (.not. ieee_is_finite(model%y(t, j))) model%y(t, j) = 0.0_dp
         end do
      end do
      do t = 1, size(model%z, 3)
         do j = 1, p
            do i = 1, q
               if (.not. ieee_is_finite(model%z(i, j, t))) model%z(i, j, t) = 0.0_dp
            end do
         end do
      end do
      model%h(:, :, 1) = observation_cov
      model%transition(:, :, 1) = transition
      model%r(:, :, 1) = sq
      model%q(:, :, 1) = 0.0_dp
      do i = 1, process_rank
         model%q(i, i, 1) = 1.0_dp
      end do
      model%a1 = matmul(transition, mu0)
      model%p1 = matmul(matmul(transition, sigma0), transpose(transition)) + process_cov
      model%p1inf = 0.0_dp
      model%missing = .false.
      filtered = kfs_filter(model)
      out%info = filtered%info
      if (out%info /= 0) return
      out%xp = filtered%a_pred
      out%pp = filtered%p_pred
      out%xf = filtered%a_filt
      out%pf = filtered%p_filt
      out%innovation = filtered%innovation
      out%innovation_cov = filtered%innovation_cov
      out%likelihood = -filtered%log_likelihood &
         - 0.5_dp*real(filtered%observations, dp)*log(2.0_dp*acos(-1.0_dp))
   end function astsa_kfilter

   pure function astsa_ksmooth(y, observation, mu0, sigma0, transition, sq, sr) result(out)
      ! Run astsa Ksmooth version 1 using the shared KFAS filter and smoother.
      real(dp), intent(in) :: y(:, :), observation(:, :, :)
      real(dp), intent(in) :: mu0(:), sigma0(:, :), transition(:, :), sq(:, :), sr(:, :)
      type(astsa_smoother_t) :: out
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      real(dp), allocatable :: process_cov(:, :), observation_cov(:, :), inverse(:, :), gain(:, :)
      real(dp) :: logdet
      integer :: i, n, p, q, rank, info

      out%filter = astsa_kfilter(y, observation, mu0, sigma0, transition, sq, sr)
      if (out%filter%info /= 0) then
         out%info = out%filter%info
         return
      end if
      n = size(y, 1)
      q = size(y, 2)
      p = size(mu0)
      rank = size(sq, 2)
      process_cov = matmul(sq, transpose(sq))
      observation_cov = matmul(sr, transpose(sr))
      allocate(model%y(n, q), model%z(q, p, size(observation, 3)), model%h(q, q, 1))
      allocate(model%transition(p, p, 1), model%r(p, rank, 1), model%q(rank, rank, 1))
      allocate(model%a1(p), model%p1(p, p), model%p1inf(p, p), model%missing(n, q))
      model%y = y
      where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
      model%z = observation
      where (.not. ieee_is_finite(model%z)) model%z = 0.0_dp
      model%h(:, :, 1) = observation_cov
      model%transition(:, :, 1) = transition
      model%r(:, :, 1) = sq
      model%q = 0.0_dp
      do i = 1, rank
         model%q(i, i, 1) = 1.0_dp
      end do
      model%a1 = matmul(transition, mu0)
      model%p1 = matmul(matmul(transition, sigma0), transpose(transition)) + process_cov
      model%p1inf = 0.0_dp
      model%missing = .false.
      filtered = kfs_filter(model)
      smoothed = kfs_smooth(model, filtered)
      out%info = smoothed%info
      if (out%info /= 0) return
      out%xs = smoothed%state
      out%ps = smoothed%covariance
      allocate(inverse(p, p), gain(p, p), out%x0n(p), out%p0n(p, p), out%j(p, p, n), out%j0(p, p))
      out%j = 0.0_dp
      do i = 1, n - 1
         call local_inverse_logdet(filtered%p_pred(:, :, i + 1), inverse, logdet, info)
         if (info /= 0) then
            out%info = 20 + info
            return
         end if
         out%j(:, :, i) = matmul(matmul(filtered%p_filt(:, :, i), transpose(transition)), inverse)
      end do
      call local_inverse_logdet(filtered%p_pred(:, :, 1), inverse, logdet, info)
      if (info /= 0) then
         out%info = 10 + info
         return
      end if
      gain = matmul(matmul(sigma0, transpose(transition)), inverse)
      out%j0 = gain
      out%x0n = mu0 + matmul(gain, out%xs(:, 1) - filtered%a_pred(:, 1))
      out%p0n = sigma0 + matmul(matmul(gain, out%ps(:, :, 1) - filtered%p_pred(:, :, 1)), transpose(gain))
   end function astsa_ksmooth

   pure function astsa_kfilter_correlated(y, observation, mu0, sigma0, transition, sq, sr, correlation) result(out)
      ! Run astsa Kfilter version 2 with correlated state and observation errors.
      real(dp), intent(in) :: y(:, :), observation(:, :, :)
      real(dp), intent(in) :: mu0(:), sigma0(:, :), transition(:, :), sq(:, :), sr(:, :)
      real(dp), intent(in) :: correlation(:, :)
      type(astsa_filter_t) :: out
      real(dp), allocatable :: process_cov(:, :), observation_cov(:, :), cross_cov(:, :)
      real(dp), allocatable :: inverse(:, :), gain(:, :), pz(:, :), innovation(:), z(:, :)
      real(dp) :: logdet
      integer :: t, n, p, q, info, slice

      n = size(y, 1)
      q = size(y, 2)
      p = size(mu0)
      if (size(observation, 1) /= q .or. size(observation, 2) /= p) then
         out%info = 1
         return
      end if
      if (size(observation, 3) /= 1 .and. size(observation, 3) < n) then
         out%info = 2
         return
      end if
      if (any(shape(sigma0) /= [p, p]) .or. any(shape(transition) /= [p, p])) then
         out%info = 3
         return
      end if
      if (size(correlation, 1) /= size(sq, 2) .or. size(correlation, 2) /= size(sr, 2)) then
         out%info = 4
         return
      end if
      process_cov = matmul(sq, transpose(sq))
      observation_cov = matmul(sr, transpose(sr))
      cross_cov = matmul(sq, matmul(correlation, transpose(sr)))
      allocate(out%xp(p, n), out%pp(p, p, n), out%xf(p, n), out%pf(p, p, n))
      allocate(out%innovation(n, q), out%innovation_cov(q, q, n))
      allocate(inverse(q, q), gain(p, q), pz(p, q), innovation(q), z(q, p))
      out%xp(:, 1) = matmul(transition, mu0)
      out%pp(:, :, 1) = symmetrize(matmul(matmul(transition, sigma0), transpose(transition)) + process_cov)
      out%likelihood = 0.0_dp
      do t = 1, n
         slice = min(t, size(observation, 3))
         z = observation(:, :, slice)
         where (.not. ieee_is_finite(z)) z = 0.0_dp
         innovation = y(t, :) - matmul(z, out%xp(:, t))
         where (.not. ieee_is_finite(innovation)) innovation = 0.0_dp
         pz = matmul(out%pp(:, :, t), transpose(z))
         out%innovation_cov(:, :, t) = symmetrize(matmul(z, pz) + observation_cov)
         call inverse_logdet(out%innovation_cov(:, :, t), inverse, logdet, info, &
            100.0_dp*epsilon(1.0_dp))
         if (info /= 0) then
            out%info = 10 + t
            return
         end if
         out%innovation(t, :) = innovation
         out%xf(:, t) = out%xp(:, t) + matmul(pz, matmul(inverse, innovation))
         out%pf(:, :, t) = symmetrize(out%pp(:, :, t) - matmul(pz, matmul(inverse, transpose(pz))))
         out%likelihood = out%likelihood + 0.5_dp*(logdet + dot_product(innovation, matmul(inverse, innovation)))
         if (t == n) cycle
         gain = matmul(matmul(transition, pz) + cross_cov, inverse)
         out%xp(:, t + 1) = matmul(transition, out%xp(:, t)) + matmul(gain, innovation)
         out%pp(:, :, t + 1) = symmetrize(matmul(matmul(transition, out%pp(:, :, t)), transpose(transition)) &
            + process_cov - matmul(gain, matmul(out%innovation_cov(:, :, t), transpose(gain))))
      end do
   end function astsa_kfilter_correlated

   pure function astsa_ksmooth_correlated(y, observation, mu0, sigma0, transition, sq, sr, correlation) result(out)
      ! Run astsa Ksmooth version 2 after correlated-error forward filtering.
      real(dp), intent(in) :: y(:, :), observation(:, :, :)
      real(dp), intent(in) :: mu0(:), sigma0(:, :), transition(:, :), sq(:, :), sr(:, :)
      real(dp), intent(in) :: correlation(:, :)
      type(astsa_smoother_t) :: out
      real(dp), allocatable :: inverse(:, :), gain(:, :)
      real(dp) :: logdet
      integer :: t, n, p, info

      out%filter = astsa_kfilter_correlated(y, observation, mu0, sigma0, transition, sq, sr, correlation)
      if (out%filter%info /= 0) then
         out%info = out%filter%info
         return
      end if
      n = size(y, 1)
      p = size(mu0)
      allocate(out%xs(p, n), out%ps(p, p, n), inverse(p, p), gain(p, p), out%j(p, p, n), out%j0(p, p))
      out%j = 0.0_dp
      out%xs(:, n) = out%filter%xf(:, n)
      out%ps(:, :, n) = out%filter%pf(:, :, n)
      do t = n - 1, 1, -1
         call inverse_logdet(out%filter%pp(:, :, t + 1), inverse, logdet, info, &
            100.0_dp*epsilon(1.0_dp))
         if (info /= 0) then
            out%info = 10 + t
            return
         end if
         gain = matmul(matmul(out%filter%pf(:, :, t), transpose(transition)), inverse)
         out%j(:, :, t) = gain
         out%xs(:, t) = out%filter%xf(:, t) + matmul(gain, out%xs(:, t + 1) - out%filter%xp(:, t + 1))
         out%ps(:, :, t) = symmetrize(out%filter%pf(:, :, t) + matmul(matmul(gain, &
            out%ps(:, :, t + 1) - out%filter%pp(:, :, t + 1)), transpose(gain)))
      end do
      allocate(out%x0n(p), out%p0n(p, p))
      call inverse_logdet(out%filter%pp(:, :, 1), inverse, logdet, info, &
         100.0_dp*epsilon(1.0_dp))
      gain = matmul(matmul(sigma0, transpose(transition)), inverse)
      out%j0 = gain
      out%x0n = mu0 + matmul(gain, out%xs(:, 1) - out%filter%xp(:, 1))
      out%p0n = symmetrize(sigma0 + matmul(matmul(gain, &
         out%ps(:, :, 1) - out%filter%pp(:, :, 1)), transpose(gain)))
   end function astsa_ksmooth_correlated

   pure function astsa_em(y, observation, mu0, sigma0, transition, q_covariance, r_covariance, &
      max_iterations, tolerance) result(out)
      ! Estimate a no-input Gaussian state-space model by astsa's EM recursions.
      real(dp), intent(in) :: y(:, :), observation(:, :, :)
      real(dp), intent(in) :: mu0(:), sigma0(:, :), transition(:, :)
      real(dp), intent(in) :: q_covariance(:, :), r_covariance(:, :)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(astsa_em_result_t) :: out
      type(astsa_smoother_t) :: smooth
      type(astsa_matrix_power_result_t) :: power_result
      real(dp), allocatable :: sq(:, :), sr(:, :), pcs(:, :, :), inverse(:, :)
      real(dp), allocatable :: s11(:, :), s10(:, :), s00(:, :), rsum(:, :), residual(:), z(:, :)
      real(dp) :: tol, convergence
      integer :: limit, iteration, t, i, n, p, obs_dim, info, slice

      limit = 100
      if (present(max_iterations)) limit = max_iterations
      tol = 0.0001_dp
      if (present(tolerance)) tol = tolerance
      n = size(y, 1)
      obs_dim = size(y, 2)
      p = size(mu0)
      if (limit < 1 .or. tol <= 0.0_dp .or. n < 2) then
         out%info = 1
         return
      end if
      out%transition = transition
      out%q = q_covariance
      out%r = r_covariance
      out%mu0 = mu0
      out%sigma0 = sigma0
      allocate(out%likelihood(limit), pcs(p, p, n), inverse(p, p))
      allocate(s11(p, p), s10(p, p), s00(p, p), rsum(obs_dim, obs_dim), residual(obs_dim), z(obs_dim, p))
      out%likelihood = 0.0_dp
      convergence = 1.0_dp + tol
      do iteration = 1, limit
         power_result = symmetric_matrix_power(out%q, 0.5_dp)
         if (power_result%info /= 0) then
            out%info = 10 + power_result%info
            return
         end if
         sq = power_result%value
         power_result = symmetric_matrix_power(out%r, 0.5_dp)
         if (power_result%info /= 0) then
            out%info = 20 + power_result%info
            return
         end if
         sr = power_result%value
         smooth = astsa_ksmooth(y, observation, out%mu0, out%sigma0, out%transition, sq, sr)
         if (smooth%info /= 0) then
            out%info = 30 + smooth%info
            return
         end if
         out%likelihood(iteration) = smooth%filter%likelihood
         if (iteration > 1) then
            convergence = (out%likelihood(iteration - 1) - out%likelihood(iteration)) &
               /abs(out%likelihood(iteration - 1))
            if (convergence < -100.0_dp*epsilon(1.0_dp)) then
               out%info = 2
               return
            end if
            if (abs(convergence) < tol) then
               out%converged = .true.
               exit
            end if
         end if
         call lag_one_covariances(smooth, observation, out%transition, pcs, info)
         if (info /= 0) then
            out%info = 40 + info
            return
         end if
         s11 = vector_outer(smooth%xs(:, 1), smooth%xs(:, 1)) + smooth%ps(:, :, 1)
         s10 = vector_outer(smooth%xs(:, 1), smooth%x0n) + pcs(:, :, 1)
         s00 = vector_outer(smooth%x0n, smooth%x0n) + smooth%p0n
         slice = min(1, size(observation, 3))
         z = observation(:, :, slice)
         residual = y(1, :) - matmul(z, smooth%xs(:, 1))
         rsum = vector_outer(residual, residual) + matmul(z, matmul(smooth%ps(:, :, 1), transpose(z)))
         do t = 2, n
            s11 = s11 + vector_outer(smooth%xs(:, t), smooth%xs(:, t)) + smooth%ps(:, :, t)
            s10 = s10 + vector_outer(smooth%xs(:, t), smooth%xs(:, t - 1)) + pcs(:, :, t)
            s00 = s00 + vector_outer(smooth%xs(:, t - 1), smooth%xs(:, t - 1)) + smooth%ps(:, :, t - 1)
            slice = min(t, size(observation, 3))
            z = observation(:, :, slice)
            residual = y(t, :) - matmul(z, smooth%xs(:, t))
            rsum = rsum + vector_outer(residual, residual) + &
               matmul(z, matmul(smooth%ps(:, :, t), transpose(z)))
         end do
         call invert_matrix(s00, inverse, info)
         if (info /= 0) then
            out%info = 50 + info
            return
         end if
         out%transition = matmul(s10, inverse)
         out%q = symmetrize((s11 - matmul(out%transition, transpose(s10)))/real(n, dp))
         out%r = 0.0_dp
         do i = 1, obs_dim
            out%r(i, i) = rsum(i, i)/real(n, dp)
         end do
         out%mu0 = smooth%x0n
         out%sigma0 = symmetrize(smooth%p0n)
      end do
      out%iterations = min(iteration, limit)
      out%likelihood = out%likelihood(:out%iterations)
   end function astsa_em

   pure function astsa_ffbs_draws(y, observation, mu0, sigma0, transition, sq, sr, normal_draws) result(out)
      ! Draw an FFBS trajectory using caller-supplied independent standard normals.
      real(dp), intent(in) :: y(:, :), observation(:, :, :)
      real(dp), intent(in) :: mu0(:), sigma0(:, :), transition(:, :), sq(:, :), sr(:, :)
      real(dp), intent(in) :: normal_draws(:, :)
      type(astsa_ffbs_result_t) :: out
      type(astsa_filter_t) :: filtered
      real(dp), allocatable :: inverse(:, :), gain(:, :), mean(:), covariance(:, :)
      real(dp) :: logdet
      integer :: t, n, p, info

      n = size(y, 1)
      p = size(mu0)
      if (any(shape(normal_draws) /= [p, n + 1])) then
         out%info = 1
         return
      end if
      filtered = astsa_kfilter(y, observation, mu0, sigma0, transition, sq, sr)
      if (filtered%info /= 0) then
         out%info = filtered%info
         return
      end if
      allocate(out%states(p, n), out%initial_state(p), inverse(p, p), gain(p, p), mean(p), covariance(p, p))
      call multivariate_normal_from_standard(filtered%xf(:, n), filtered%pf(:, :, n), normal_draws(:, n), &
         out%states(:, n), info)
      if (info /= 0) then
         out%info = 10 + info
         return
      end if
      do t = n - 1, 1, -1
         call inverse_logdet(filtered%pp(:, :, t + 1), inverse, logdet, info, &
            100.0_dp*epsilon(1.0_dp))
         if (info /= 0) then
            out%info = 20 + t
            return
         end if
         gain = matmul(matmul(filtered%pf(:, :, t), transpose(transition)), inverse)
         mean = filtered%xf(:, t) + matmul(gain, out%states(:, t + 1) - filtered%xp(:, t + 1))
         covariance = symmetrize(filtered%pf(:, :, t) - &
            matmul(gain, matmul(filtered%pp(:, :, t + 1), transpose(gain))))
         call multivariate_normal_from_standard(mean, covariance, normal_draws(:, t), out%states(:, t), info)
         if (info /= 0) then
            out%info = 30 + t
            return
         end if
      end do
      call inverse_logdet(filtered%pp(:, :, 1), inverse, logdet, info, &
         100.0_dp*epsilon(1.0_dp))
      if (info /= 0) then
         out%info = 40 + info
         return
      end if
      gain = matmul(matmul(sigma0, transpose(transition)), inverse)
      mean = mu0 + matmul(gain, out%states(:, 1) - filtered%xp(:, 1))
      covariance = symmetrize(sigma0 + matmul(gain, matmul(filtered%pp(:, :, 1), transpose(gain))))
      call multivariate_normal_from_standard(mean, covariance, normal_draws(:, n + 1), out%initial_state, info)
      if (info /= 0) out%info = 50 + info
   end function astsa_ffbs_draws

   function astsa_ffbs(y, observation, mu0, sigma0, transition, sq, sr) result(out)
      ! Draw an FFBS trajectory using the Fortran intrinsic random-number generator.
      real(dp), intent(in) :: y(:, :), observation(:, :, :)
      real(dp), intent(in) :: mu0(:), sigma0(:, :), transition(:, :), sq(:, :), sr(:, :)
      type(astsa_ffbs_result_t) :: out
      real(dp), allocatable :: draws(:, :)

      allocate(draws(size(mu0), size(y, 1) + 1))
      call random_standard_normal_matrix(draws)
      out = astsa_ffbs_draws(y, observation, mu0, sigma0, transition, sq, sr, draws)
   end function astsa_ffbs

   pure function arma_spectrum(ar, ma, noise_variance, n_frequency, sampling_frequency) result(out)
      ! Compute astsa arma.spec's theoretical ARMA spectrum without plotting.
      real(dp), intent(in) :: ar(:), ma(:)
      real(dp), intent(in), optional :: noise_variance, sampling_frequency
      integer, intent(in), optional :: n_frequency
      type(astsa_spectrum_t) :: out
      real(dp) :: variance, sample_rate, angle
      complex(dp) :: ar_polynomial, ma_polynomial
      integer :: i, j, count

      variance = 1.0_dp
      if (present(noise_variance)) variance = noise_variance
      sample_rate = 1.0_dp
      if (present(sampling_frequency)) sample_rate = sampling_frequency
      count = 500
      if (present(n_frequency)) count = n_frequency
      if (count < 2 .or. variance < 0.0_dp .or. sample_rate <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%frequency(count), out%spectrum(count, 1))
      do i = 1, count
         out%frequency(i) = 0.5_dp*sample_rate*real(i - 1, dp)/real(count - 1, dp)
         ar_polynomial = cmplx(1.0_dp, 0.0_dp, dp)
         ma_polynomial = cmplx(1.0_dp, 0.0_dp, dp)
         do j = 1, size(ar)
            angle = -2.0_dp*acos(-1.0_dp)*out%frequency(i)*real(j, dp)/sample_rate
            ar_polynomial = ar_polynomial - ar(j)*cmplx(cos(angle), sin(angle), dp)
         end do
         do j = 1, size(ma)
            angle = -2.0_dp*acos(-1.0_dp)*out%frequency(i)*real(j, dp)/sample_rate
            ma_polynomial = ma_polynomial + ma(j)*cmplx(cos(angle), sin(angle), dp)
         end do
         out%spectrum(i, 1) = variance*abs(ma_polynomial)**2/abs(ar_polynomial)**2
      end do
      out%sampling_frequency = sample_rate
      out%original_size = count
      out%used_size = count
   end function arma_spectrum

   pure function mv_periodogram(x, sampling_frequency, demean, detrend, taper, pad, span) result(out)
      ! Compute an mvspec-style spectral matrix, coherence, and phase estimate.
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in), optional :: sampling_frequency, taper
      logical, intent(in), optional :: demean, detrend
      integer, intent(in), optional :: pad, span
      type(astsa_spectrum_t) :: out
      real(dp), allocatable :: work(:, :), weights(:), beta(:), se(:), residuals(:), design(:, :)
      complex(dp), allocatable :: transform(:, :), full_spectrum(:, :, :), smoothed(:, :, :)
      real(dp) :: sample_rate, taper_fraction, u2, u4, rss, angle, lh
      integer :: n0, n, series, frequencies, padding, half_span, i, j, k, t, info
      logical :: remove_mean, remove_trend

      n0 = size(x, 1)
      series = size(x, 2)
      sample_rate = 1.0_dp
      if (present(sampling_frequency)) sample_rate = sampling_frequency
      remove_mean = .false.
      if (present(demean)) remove_mean = demean
      remove_trend = .true.
      if (present(detrend)) remove_trend = detrend
      if (remove_mean) remove_trend = .false.
      taper_fraction = 0.0_dp
      if (present(taper)) taper_fraction = taper
      padding = 0
      if (present(pad)) padding = pad
      half_span = 0
      if (present(span)) half_span = span/2
      if (n0 < 3 .or. series < 1 .or. sample_rate <= 0.0_dp .or. &
         taper_fraction < 0.0_dp .or. taper_fraction >= 0.5_dp .or. padding < 0) then
         out%info = 1
         return
      end if
      n = n0*(padding + 1)
      allocate(work(n, series))
      work = 0.0_dp
      work(:n0, :) = x
      if (remove_mean) then
         do j = 1, series
            work(:n0, j) = work(:n0, j) - sum(work(:n0, j))/real(n0, dp)
         end do
      else if (remove_trend) then
         allocate(design(n0, 2))
         design(:, 1) = 1.0_dp
         design(:, 2) = [(real(i, dp), i=1, n0)]
         do j = 1, series
            call ols_fit(design, work(:n0, j), beta, se, residuals, rss, info)
            if (info /= 0) then
               out%info = 10 + info
               return
            end if
            work(:n0, j) = residuals
         end do
      end if
      call apply_split_taper(work(:n0, :), taper_fraction)
      allocate(transform(n, series))
      do k = 0, n - 1
         do j = 1, series
            transform(k + 1, j) = cmplx(0.0_dp, 0.0_dp, dp)
            do t = 0, n - 1
               angle = -2.0_dp*acos(-1.0_dp)*real(k*t, dp)/real(n, dp)
               transform(k + 1, j) = transform(k + 1, j) + &
                  work(t + 1, j)*cmplx(cos(angle), sin(angle), dp)
            end do
         end do
      end do
      allocate(full_spectrum(n, series, series))
      do i = 1, series
         do j = 1, series
            full_spectrum(:, i, j) = transform(:, i)*conjg(transform(:, j))/real(n0, dp)/sample_rate
            full_spectrum(1, i, j) = 0.5_dp*(full_spectrum(2, i, j) + full_spectrum(n, i, j))
         end do
      end do
      lh = 1.0_dp
      if (half_span > 0) then
         call modified_daniell_weights(half_span, weights)
         allocate(smoothed(n, series, series))
         do i = 1, series
            do j = 1, series
               do k = 1, n
                  smoothed(k, i, j) = cmplx(0.0_dp, 0.0_dp, dp)
                  do t = -half_span, half_span
                     smoothed(k, i, j) = smoothed(k, i, j) + &
                        weights(t + half_span + 1)*full_spectrum(1 + modulo(k - 1 + t, n), i, j)
                  end do
               end do
            end do
         end do
         full_spectrum = smoothed
         lh = 1.0_dp/sum(weights*weights)
      end if
      frequencies = n/2
      allocate(out%frequency(frequencies), out%spectrum(frequencies, series))
      allocate(out%spectral_matrix(series, series, frequencies))
      allocate(out%coherence(series, series, frequencies), out%phase(series, series, frequencies))
      do k = 1, frequencies
         out%frequency(k) = sample_rate*real(k, dp)/real(n, dp)
         do i = 1, series
            do j = 1, series
               out%spectral_matrix(i, j, k) = full_spectrum(k + 1, i, j)
            end do
            out%spectrum(k, i) = real(out%spectral_matrix(i, i, k), dp)
         end do
         do i = 1, series
            do j = 1, series
               if (out%spectrum(k, i)*out%spectrum(k, j) > 0.0_dp) then
                  out%coherence(i, j, k) = abs(out%spectral_matrix(i, j, k))**2 &
                     /(out%spectrum(k, i)*out%spectrum(k, j))
               else
                  out%coherence(i, j, k) = 0.0_dp
               end if
               out%phase(i, j, k) = atan2(aimag(out%spectral_matrix(i, j, k)), &
                  real(out%spectral_matrix(i, j, k), dp))
            end do
         end do
      end do
      u2 = 1.0_dp - 1.25_dp*taper_fraction
      u4 = 1.0_dp - 93.0_dp*taper_fraction/64.0_dp
      out%spectrum = out%spectrum/u2
      out%degrees_freedom = 2.0_dp*lh/(u4/(u2*u2))*real(n0, dp)/real(n, dp)
      out%bandwidth = lh*sample_rate/real(n, dp)
      out%sampling_frequency = sample_rate
      out%original_size = n0
      out%used_size = n
   end function mv_periodogram

   pure subroutine apply_split_taper(x, fraction)
      ! Apply R's split cosine-bell taper to both ends of a series matrix.
      real(dp), intent(inout) :: x(:, :)
      real(dp), intent(in) :: fraction
      real(dp) :: weight
      integer :: i, m, n
      n = size(x, 1)
      m = int(floor(real(n, dp)*fraction))
      if (m < 1) return
      do i = 1, m
         weight = 0.5_dp*(1.0_dp - cos(acos(-1.0_dp)* &
            (real(i, dp) - 0.5_dp)/real(m, dp)))
         x(i, :) = weight*x(i, :)
         x(n - i + 1, :) = weight*x(n - i + 1, :)
      end do
   end subroutine apply_split_taper

   pure subroutine modified_daniell_weights(m, weights)
      ! Construct the normalized modified-Daniell kernel used by mvspec.
      integer, intent(in) :: m
      real(dp), allocatable, intent(out) :: weights(:)
      allocate(weights(2*m + 1))
      weights = 1.0_dp/real(2*m, dp)
      weights(1) = 0.5_dp/real(2*m, dp)
      weights(2*m + 1) = weights(1)
   end subroutine modified_daniell_weights

   pure subroutine lag_one_covariances(smooth, observation, transition, pcs, info)
      ! Compute astsa's lag-one smoothed state covariance recursion.
      type(astsa_smoother_t), intent(in) :: smooth
      real(dp), intent(in) :: observation(:, :, :), transition(:, :)
      real(dp), intent(out) :: pcs(:, :, :)
      integer, intent(out) :: info
      real(dp), allocatable :: inverse(:, :), measurement_gain(:, :), z(:, :)
      real(dp) :: logdet
      integer :: n, p, t, slice

      n = size(smooth%xs, 2)
      p = size(smooth%xs, 1)
      allocate(inverse(size(smooth%filter%innovation_cov, 1), size(smooth%filter%innovation_cov, 2)))
      slice = min(n, size(observation, 3))
      z = observation(:, :, slice)
      call inverse_logdet(smooth%filter%innovation_cov(:, :, n), inverse, logdet, info, &
         100.0_dp*epsilon(1.0_dp))
      if (info /= 0) return
      measurement_gain = matmul(matmul(smooth%filter%pp(:, :, n), transpose(z)), inverse)
      pcs(:, :, n) = matmul(matmul(identity_matrix(p) - matmul(measurement_gain, z), transition), &
         smooth%filter%pf(:, :, n - 1))
      do t = n, 3, -1
         pcs(:, :, t - 1) = matmul(smooth%filter%pf(:, :, t - 1), transpose(smooth%j(:, :, t - 2))) &
            + matmul(smooth%j(:, :, t - 1), matmul(pcs(:, :, t) &
            - matmul(transition, smooth%filter%pf(:, :, t - 1)), transpose(smooth%j(:, :, t - 2))))
      end do
      pcs(:, :, 1) = matmul(smooth%filter%pf(:, :, 1), transpose(smooth%j0)) &
         + matmul(smooth%j(:, :, 1), matmul(pcs(:, :, 2) &
         - matmul(transition, smooth%filter%pf(:, :, 1)), transpose(smooth%j0)))
      info = 0
   end subroutine lag_one_covariances

   pure function vector_outer(a, b) result(product)
      ! Form an outer product for EM sufficient statistics.
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: product(size(a), size(b))
      product = spread(a, 2, size(b))*spread(b, 1, size(a))
   end function vector_outer

   pure function poly_mul(p, q) result(product)
      ! Multiply two coefficient vectors by polynomial convolution.
      real(dp), intent(in) :: p(:), q(:)
      real(dp), allocatable :: product(:)
      integer :: i, j
      allocate(product(size(p) + size(q) - 1))
      product = 0.0_dp
      do i = 1, size(p)
         do j = 1, size(q)
            product(i + j - 1) = product(i + j - 1) + p(i)*q(j)
         end do
      end do
   end function poly_mul

   pure function arma_to_ar(ar, ma, lag_max) result(coefficients)
      ! Convert an ARMA model to the astsa infinite-AR truncation.
      real(dp), intent(in) :: ar(:), ma(:)
      integer, intent(in) :: lag_max
      real(dp), allocatable :: coefficients(:), psi(:)
      integer :: i, j
      allocate(coefficients(lag_max), psi(0:lag_max))
      psi = 0.0_dp
      psi(0) = 1.0_dp
      do i = 1, lag_max
         if (i <= size(ar)) psi(i) = -ar(i)
         do j = 1, min(size(ma), i)
            psi(i) = psi(i) - ma(j)*psi(i - j)
         end do
      end do
      coefficients = psi(1:)
   end function arma_to_ar

   pure function symmetric_matrix_power(a, power) result(out)
      ! Return a real symmetric matrix power and decomposition status.
      real(dp), intent(in) :: a(:, :), power
      type(astsa_matrix_power_result_t) :: out
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      integer :: status, i, n
      n = size(a, 1)
      allocate(out%value(n, n))
      if (size(a, 2) /= n) then
         out%value = 0.0_dp
         out%info = 1
         return
      end if
      call symmetric_eigen(a, eigenvalues, eigenvectors, status)
      if (status /= 0 .or. (power < 0.0_dp .and. minval(abs(eigenvalues)) <= epsilon(1.0_dp))) then
         out%value = 0.0_dp
         out%info = 2
         return
      end if
      out%value = 0.0_dp
      do i = 1, n
         out%value = out%value + eigenvalues(i)**power* &
            matmul(reshape(eigenvectors(:, i), [n, 1]), reshape(eigenvectors(:, i), [1, n]))
      end do
      out%info = 0
   end function symmetric_matrix_power

   pure integer function fdr_cutoff(p_values, q_level) result(index)
      ! Return astsa FDR's original index of the largest rejected p-value.
      real(dp), intent(in) :: p_values(:)
      real(dp), intent(in), optional :: q_level
      real(dp), allocatable :: sorted(:)
      integer, allocatable :: order(:)
      real(dp) :: level, temporary_value
      integer :: i, j, n, accepted, temporary
      level = 0.05_dp
      if (present(q_level)) level = q_level
      n = size(p_values)
      allocate(sorted(n), order(n))
      sorted = p_values
      order = [(i, i=1, n)]
      do i = 1, n - 1
         do j = i + 1, n
            if (sorted(j) < sorted(i)) then
               temporary_value = sorted(i)
               sorted(i) = sorted(j)
               sorted(j) = temporary_value
               temporary = order(i)
               order(i) = order(j)
               order(j) = temporary
            end if
         end do
      end do
      level = 0.05_dp
      if (present(q_level)) level = q_level
      accepted = 0
      do i = 1, n
         if (sorted(i) <= level*real(i, dp)/real(n, dp)) accepted = i
      end do
      index = 0
      if (accepted > 0) index = order(accepted)
   end function fdr_cutoff

   pure function arma_check(ar, ma, sar, sma, season, redundancy_tolerance) result(out)
      ! Check ARMA causality, invertibility, and approximate common factors.
      real(dp), intent(in), optional :: ar(:), ma(:), sar(:), sma(:)
      integer, intent(in), optional :: season
      real(dp), intent(in), optional :: redundancy_tolerance
      type(astsa_arma_diagnostic_t) :: out
      type(astsa_polynomial_roots_t) :: roots
      real(dp), allocatable :: polynomial(:)
      integer :: p, q, ps, qs, period, i

      p = 0
      if (present(ar)) p = size(ar)
      q = 0
      if (present(ma)) q = size(ma)
      ps = 0
      if (present(sar)) ps = size(sar)
      qs = 0
      if (present(sma)) qs = size(sma)
      period = 1
      if (present(season)) period = season
      if ((ps > 0 .or. qs > 0) .and. .not. present(season)) period = 12
      out%seasonal_period = period
      if (present(redundancy_tolerance)) then
         if (redundancy_tolerance >= 0.0_dp) out%redundancy_tolerance = redundancy_tolerance
      end if
      if (period < 1 .or. (present(season) .and. (period <= p .or. period <= q))) then
         out%info = 1
         return
      end if

      allocate(polynomial(p + 1))
      polynomial = 0.0_dp
      polynomial(1) = 1.0_dp
      if (p > 0) polynomial(2:) = -ar
      roots = polynomial_roots(polynomial)
      out%ar_roots = roots%roots
      if (.not. roots%converged) then
         out%info = 2
         return
      end if
      out%causal = all(abs(out%ar_roots) > 1.0_dp)

      deallocate(polynomial)
      allocate(polynomial(q + 1))
      polynomial = 0.0_dp
      polynomial(1) = 1.0_dp
      if (q > 0) polynomial(2:) = ma
      roots = polynomial_roots(polynomial)
      out%ma_roots = roots%roots
      if (.not. roots%converged) then
         out%info = 3
         return
      end if
      out%invertible = all(abs(out%ma_roots) > 1.0_dp)

      deallocate(polynomial)
      allocate(polynomial(ps*period + 1))
      polynomial = 0.0_dp
      polynomial(1) = 1.0_dp
      do i = 1, ps
         polynomial(i*period + 1) = -sar(i)
      end do
      roots = polynomial_roots(polynomial)
      out%seasonal_ar_roots = roots%roots
      if (.not. roots%converged) then
         out%info = 4
         return
      end if
      out%causal = out%causal .and. all(abs(out%seasonal_ar_roots) > 1.0_dp)

      deallocate(polynomial)
      allocate(polynomial(qs*period + 1))
      polynomial = 0.0_dp
      polynomial(1) = 1.0_dp
      do i = 1, qs
         polynomial(i*period + 1) = sma(i)
      end do
      roots = polynomial_roots(polynomial)
      out%seasonal_ma_roots = roots%roots
      if (.not. roots%converged) then
         out%info = 5
         return
      end if
      out%invertible = out%invertible .and. all(abs(out%seasonal_ma_roots) > 1.0_dp)

      if (.not. out%causal .or. .not. out%invertible) return
      out%redundancy_checked = .true.
      out%redundant = inverse_roots_overlap(out%ar_roots, out%ma_roots, out%redundancy_tolerance)
      out%seasonal_redundant = inverse_roots_overlap(out%seasonal_ar_roots, &
         out%seasonal_ma_roots, out%redundancy_tolerance)
   end function arma_check

   pure function pre_white(series1, series2, differences, max_lag, order_max) result(out)
      ! Select an AR prewhitener and return aligned filtered series and CCF values.
      real(dp), intent(in) :: series1(:), series2(:)
      integer, intent(in), optional :: differences, max_lag, order_max
      type(astsa_prewhite_result_t) :: out
      type(yule_walker_result_t) :: fit, candidate
      real(dp), allocatable :: first(:), second(:), ccf(:)
      real(dp) :: first_mean
      integer :: d, maximum_order, lag, n, p, i, j

      d = 0
      if (present(differences)) d = differences
      if (size(series1) /= size(series2) .or. d < 0 .or. size(series1) - d < 3) then
         out%info = 1
         return
      end if
      first = difference_series(series1, d)
      second = difference_series(series2, d)
      n = size(first)
      maximum_order = min(30, ceiling(0.15_dp*real(n, dp)))
      if (present(order_max)) maximum_order = order_max
      maximum_order = min(maximum_order, n - 2)
      if (maximum_order < 0) then
         out%info = 2
         return
      end if
      fit = yule_walker_fit(first, 0)
      do p = 1, maximum_order
         candidate = yule_walker_fit(first, p)
         if (candidate%info == 0 .and. candidate%criterion < fit%criterion) fit = candidate
      end do
      if (fit%info /= 0) then
         out%info = 3
         return
      end if

      out%order = size(fit%coefficients)
      out%differences = d
      out%ar = fit%coefficients
      allocate(out%first(n - out%order), out%second(n - out%order))
      first_mean = sum(first)/real(n, dp)
      do i = out%order + 1, n
         out%first(i - out%order) = first(i) - first_mean
         out%second(i - out%order) = second(i)
         do j = 1, out%order
            out%first(i - out%order) = out%first(i - out%order) - out%ar(j)*(first(i - j) - first_mean)
            out%second(i - out%order) = out%second(i - out%order) - out%ar(j)*second(i - j)
         end do
      end do
      lag = min(50, floor(0.2_dp*real(n, dp)))
      if (present(max_lag)) lag = max_lag
      lag = min(lag, size(out%first) - 1)
      if (lag < 0) then
         out%info = 4
         return
      end if
      ccf = ccf_values(out%first, out%second, lag)
      allocate(out%correlation(2*lag + 1))
      out%correlation = ccf
      out%lags = [(i, i=-lag, lag)]
   end function pre_white

   pure function ar_boot_draws(series, order, residual_indices, probabilities) result(out)
      ! Bootstrap a fixed-order Yule-Walker AR model from supplied residual indices.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: order
      integer, intent(in) :: residual_indices(:, :)
      real(dp), intent(in), optional :: probabilities(:)
      type(astsa_ar_bootstrap_t) :: out
      type(yule_walker_result_t) :: original, fitted
      real(dp), allocatable :: centered(:), residuals(:), simulated_centered(:), sorted(:)
      real(dp) :: series_mean
      integer :: n, residual_count, bootstrap_count, replicate, t, lag, coefficient, probability_index

      n = size(series)
      residual_count = n - order
      bootstrap_count = size(residual_indices, 2)
      if (order < 1 .or. residual_count < 2 .or. bootstrap_count < 1 .or. &
         size(residual_indices, 1) /= residual_count) then
         out%info = 1
         return
      end if
      if (any(residual_indices < 1) .or. any(residual_indices > residual_count)) then
         out%info = 2
         return
      end if
      original = yule_walker_fit(series, order)
      if (original%info /= 0) then
         out%info = 10 + original%info
         return
      end if
      series_mean = sum(series)/real(n, dp)
      centered = series - series_mean
      allocate(residuals(residual_count))
      do t = order + 1, n
         residuals(t - order) = centered(t)
         do lag = 1, order
            residuals(t - order) = residuals(t - order) - original%coefficients(lag)*centered(t - lag)
         end do
      end do
      out%original_coefficients = original%coefficients
      allocate(out%coefficient_draws(bootstrap_count, order), out%simulated_series(n, bootstrap_count))
      allocate(simulated_centered(n))
      do replicate = 1, bootstrap_count
         simulated_centered(:order) = centered(:order)
         do t = order + 1, n
            simulated_centered(t) = residuals(residual_indices(t - order, replicate))
            do lag = 1, order
               simulated_centered(t) = simulated_centered(t) + &
                  original%coefficients(lag)*simulated_centered(t - lag)
            end do
         end do
         out%simulated_series(:, replicate) = simulated_centered + series_mean
         fitted = yule_walker_fit(out%simulated_series(:, replicate), order)
         if (fitted%info /= 0) then
            out%info = 20 + replicate
            return
         end if
         out%coefficient_draws(replicate, :) = fitted%coefficients
      end do
      allocate(out%coefficient_means(order))
      out%coefficient_means = sum(out%coefficient_draws, dim=1)/real(bootstrap_count, dp)
      if (present(probabilities)) then
         out%probabilities = probabilities
      else
         out%probabilities = [0.01_dp, 0.025_dp, 0.05_dp, 0.1_dp, 0.25_dp, 0.5_dp, &
            0.75_dp, 0.9_dp, 0.95_dp, 0.975_dp, 0.99_dp]
      end if
      if (any(out%probabilities < 0.0_dp) .or. any(out%probabilities > 1.0_dp)) then
         out%info = 3
         return
      end if
      allocate(out%quantiles(size(out%probabilities), order))
      do coefficient = 1, order
         sorted = sorted_values(out%coefficient_draws(:, coefficient))
         do probability_index = 1, size(out%probabilities)
            out%quantiles(probability_index, coefficient) = &
               sample_quantile(sorted, out%probabilities(probability_index))
         end do
      end do
   end function ar_boot_draws

   function ar_boot(series, order, bootstrap_count, probabilities) result(out)
      ! Bootstrap a fixed-order AR model using shared random residual indices.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: order, bootstrap_count
      real(dp), intent(in), optional :: probabilities(:)
      type(astsa_ar_bootstrap_t) :: out
      integer, allocatable :: indices(:, :)
      integer :: i, j, residual_count

      residual_count = size(series) - order
      if (residual_count < 1 .or. bootstrap_count < 1) then
         out%info = 1
         return
      end if
      allocate(indices(residual_count, bootstrap_count))
      do j = 1, bootstrap_count
         do i = 1, residual_count
            indices(i, j) = min(residual_count, 1 + int(random_uniform()*real(residual_count, dp)))
         end do
      end do
      out = ar_boot_draws(series, order, indices, probabilities)
   end function ar_boot

   pure real(dp) function sample_quantile(sorted, probability) result(value)
      ! Return an R type-7 sample quantile from ascending values.
      real(dp), intent(in) :: sorted(:), probability
      real(dp) :: position, fraction
      integer :: lower

      if (size(sorted) == 1) then
         value = sorted(1)
         return
      end if
      position = 1.0_dp + real(size(sorted) - 1, dp)*probability
      lower = floor(position)
      if (lower >= size(sorted)) then
         value = sorted(size(sorted))
      else
         fraction = position - real(lower, dp)
         value = (1.0_dp - fraction)*sorted(lower) + fraction*sorted(lower + 1)
      end if
   end function sample_quantile

   pure function ar_mcmc_draws(series, order, retained_draws, warmup, normal_draws, gamma_draws, &
      prior_variance, prior_shape, prior_rate, probabilities) result(out)
      ! Sample a conjugate Bayesian AR posterior from supplied Gaussian and gamma draws.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: order, retained_draws, warmup
      real(dp), intent(in) :: normal_draws(:, :), gamma_draws(:)
      real(dp), intent(in), optional :: prior_variance, prior_shape, prior_rate, probabilities(:)
      type(astsa_ar_mcmc_t) :: out
      real(dp), allocatable :: design(:, :), response(:), xtx(:, :), inverse(:, :), covariance(:, :)
      real(dp), allocatable :: mean(:), draw(:), residuals(:), coefficients(:, :), variance_draws(:)
      real(dp), allocatable :: combined(:, :), sorted(:)
      real(dp) :: coefficient_prior, prior_a, rate_prior, posterior_shape, rate
      integer :: n, observations, parameter_count, total, iteration, lag, status, column, probability_index

      n = size(series)
      observations = n - order
      parameter_count = order + 1
      total = retained_draws + warmup
      coefficient_prior = 50.0_dp
      if (present(prior_variance)) coefficient_prior = prior_variance
      prior_a = 1.0_dp
      if (present(prior_shape)) prior_a = prior_shape
      rate_prior = 2.0_dp
      if (present(prior_rate)) rate_prior = prior_rate
      if (order < 0 .or. observations < 2 .or. retained_draws < 1 .or. warmup < 0 .or. &
         coefficient_prior <= 0.0_dp .or. prior_a <= 0.0_dp .or. rate_prior <= 0.0_dp .or. &
         any(shape(normal_draws) /= [parameter_count, total - 1]) .or. size(gamma_draws) /= total - 1 .or. &
         any(gamma_draws <= 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(design(observations, parameter_count), response(observations))
      design(:, 1) = 1.0_dp
      response = series(order + 1:)
      do lag = 1, order
         design(:, lag + 1) = series(order + 1 - lag:n - lag)
      end do
      xtx = matmul(transpose(design), design)
      do column = 1, parameter_count
         xtx(column, column) = xtx(column, column) + 1.0_dp/coefficient_prior
      end do
      call invert_matrix(xtx, inverse, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      mean = matmul(inverse, matmul(transpose(design), response))
      allocate(coefficients(parameter_count, total), variance_draws(total), residuals(observations))
      allocate(draw(parameter_count))
      coefficients = 0.0_dp
      variance_draws(1) = 1.0_dp
      posterior_shape = 0.5_dp*real(observations, dp) + prior_a
      do iteration = 2, total
         covariance = variance_draws(iteration - 1)*inverse
         call multivariate_normal_from_standard(mean, covariance, normal_draws(:, iteration - 1), draw, status)
         if (status /= 0) then
            out%info = 20 + iteration
            return
         end if
         coefficients(:, iteration) = draw
         residuals = response - matmul(design, draw)
         rate = 0.5_dp*sum(residuals**2) + rate_prior + &
            0.5_dp*sum(draw**2)/coefficient_prior
         variance_draws(iteration) = rate/gamma_draws(iteration - 1)
      end do
      out%coefficient_draws = transpose(coefficients(:, warmup + 1:))
      out%innovation_sd = sqrt(variance_draws(warmup + 1:))
      allocate(combined(retained_draws, parameter_count + 1))
      combined(:, :parameter_count) = out%coefficient_draws
      combined(:, parameter_count + 1) = out%innovation_sd
      allocate(out%means(parameter_count + 1), out%standard_deviations(parameter_count + 1))
      allocate(out%effective_size(parameter_count + 1))
      out%means = sum(combined, dim=1)/real(retained_draws, dp)
      do column = 1, parameter_count + 1
         if (retained_draws > 1) then
            out%standard_deviations(column) = sqrt(sum((combined(:, column) - out%means(column))**2)/ &
               real(retained_draws - 1, dp))
         else
            out%standard_deviations(column) = 0.0_dp
         end if
         out%effective_size(column) = effective_sample_size(combined(:, column))
      end do
      if (present(probabilities)) then
         out%probabilities = probabilities
      else
         out%probabilities = [0.01_dp, 0.025_dp, 0.05_dp, 0.1_dp, 0.25_dp, 0.5_dp, &
            0.75_dp, 0.9_dp, 0.95_dp, 0.975_dp, 0.99_dp]
      end if
      if (any(out%probabilities < 0.0_dp) .or. any(out%probabilities > 1.0_dp)) then
         out%info = 2
         return
      end if
      allocate(out%quantiles(size(out%probabilities), parameter_count + 1))
      do column = 1, parameter_count + 1
         sorted = sorted_values(combined(:, column))
         do probability_index = 1, size(out%probabilities)
            out%quantiles(probability_index, column) = sample_quantile(sorted, out%probabilities(probability_index))
         end do
      end do
   end function ar_mcmc_draws

   function ar_mcmc(series, order, retained_draws, warmup, prior_variance, prior_shape, &
      prior_rate, probabilities) result(out)
      ! Sample a Bayesian AR posterior using centralized Gaussian and gamma generation.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: order, retained_draws, warmup
      real(dp), intent(in), optional :: prior_variance, prior_shape, prior_rate, probabilities(:)
      type(astsa_ar_mcmc_t) :: out
      real(dp), allocatable :: normal_draws(:, :), gamma_draws(:)
      real(dp) :: shape
      integer :: total, observations, iteration

      total = retained_draws + warmup
      observations = size(series) - order
      shape = 1.0_dp
      if (present(prior_shape)) shape = prior_shape
      allocate(normal_draws(order + 1, total - 1), gamma_draws(total - 1))
      call random_standard_normal_matrix(normal_draws)
      do iteration = 1, total - 1
         gamma_draws(iteration) = random_gamma(0.5_dp*real(observations, dp) + shape)
      end do
      out = ar_mcmc_draws(series, order, retained_draws, warmup, normal_draws, gamma_draws, &
         prior_variance, prior_shape, prior_rate, probabilities)
   end function ar_mcmc

   pure function spectrum_ic(series, order_max, use_bic, detrend, frequency_count) result(out)
      ! Select a Yule-Walker AR spectrum using relative AIC or BIC.
      real(dp), intent(in) :: series(:)
      integer, intent(in), optional :: order_max, frequency_count
      logical, intent(in), optional :: use_bic, detrend
      type(astsa_spectrum_ic_t) :: out
      type(yule_walker_result_t) :: fit, selected_fit
      type(astsa_spectrum_t) :: selected_spectrum
      real(dp), allocatable :: work(:), design(:, :), beta(:), standard_errors(:), residuals(:)
      real(dp) :: rss, minimum
      integer :: maximum_order, count_frequency, order, n, status
      logical :: select_bic, remove_trend

      n = size(series)
      if (n < 3) then
         out%info = 1
         return
      end if
      maximum_order = min(100, ceiling(0.1_dp*real(n, dp)))
      if (present(order_max)) maximum_order = order_max
      maximum_order = min(maximum_order, n - 2)
      count_frequency = 500
      if (present(frequency_count)) count_frequency = frequency_count
      select_bic = .false.
      if (present(use_bic)) select_bic = use_bic
      remove_trend = .true.
      if (present(detrend)) remove_trend = detrend
      if (maximum_order < 0 .or. count_frequency < 2) then
         out%info = 1
         return
      end if
      work = series
      if (remove_trend) then
         allocate(design(n, 2))
         design(:, 1) = 1.0_dp
         design(:, 2) = [(real(order, dp), order=1, n)]
         call ols_fit(design, series, beta, standard_errors, residuals, rss, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         work = residuals
      end if
      allocate(out%orders(0:maximum_order), out%aic(0:maximum_order), out%bic(0:maximum_order))
      do order = 0, maximum_order
         out%orders(order) = order
         fit = yule_walker_fit(work, order)
         if (fit%info == 0) then
            out%aic(order) = real(n, dp)*log(fit%variance) + 2.0_dp*real(order, dp)
            out%bic(order) = real(n, dp)*log(fit%variance) + log(real(n, dp))*real(order, dp)
         else
            out%aic(order) = huge(1.0_dp)
            out%bic(order) = huge(1.0_dp)
         end if
      end do
      minimum = minval(out%aic)
      where (out%aic < huge(1.0_dp)) out%aic = out%aic - minimum
      minimum = minval(out%bic)
      where (out%bic < huge(1.0_dp)) out%bic = out%bic - minimum
      if (select_bic) then
         out%selected_order = minloc(out%bic, dim=1) - 1
      else
         out%selected_order = minloc(out%aic, dim=1) - 1
      end if
      selected_fit = yule_walker_fit(work, out%selected_order)
      if (selected_fit%info /= 0) then
         out%info = 20 + selected_fit%info
         return
      end if
      out%coefficients = selected_fit%coefficients
      out%innovation_variance = selected_fit%variance
      selected_spectrum = arma_spectrum(out%coefficients, [real(dp) ::], &
         noise_variance=out%innovation_variance, n_frequency=count_frequency)
      out%frequency = selected_spectrum%frequency
      out%spectrum = selected_spectrum%spectrum(:, 1)
   end function spectrum_ic

   pure function test_linearity(series, detrend) result(out)
      ! Compute astsa test.linear's normalized block bispectrum and p-values.
      real(dp), intent(in) :: series(:)
      logical, intent(in), optional :: detrend
      type(astsa_linearity_test_t) :: out
      real(dp), allocatable :: work(:), design(:, :), beta(:), standard_errors(:), residuals(:), power(:)
      complex(dp), allocatable :: transform(:, :)
      real(dp) :: rss, angle, scale
      integer :: n, block_length, block_count, half_length, block, k, t, k1, k2, status
      logical :: remove_trend

      n = size(series)
      if (n < 8) then
         out%info = 1
         return
      end if
      remove_trend = .false.
      if (present(detrend)) remove_trend = detrend
      work = series
      if (remove_trend) then
         allocate(design(n, 2))
         design(:, 1) = 1.0_dp
         design(:, 2) = [(real(t, dp), t=1, n)]
         call ols_fit(design, series, beta, standard_errors, residuals, rss, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         work = residuals
      else
         work = work - sum(work)/real(n, dp)
      end if
      block_length = floor(real(n, dp)**0.49_dp)
      block_count = n/block_length
      half_length = (block_length - 1)/2
      if (half_length < 1 .or. block_count < 2) then
         out%info = 2
         return
      end if
      out%block_length = block_length
      out%block_count = block_count
      allocate(transform(block_length - 1, block_count), power(block_length - 1))
      do block = 1, block_count
         do k = 1, block_length - 1
            transform(k, block) = cmplx(0.0_dp, 0.0_dp, dp)
            do t = 0, block_length - 1
               angle = -2.0_dp*acos(-1.0_dp)*real(k*t, dp)/real(block_length, dp)
               transform(k, block) = transform(k, block) + &
                  work((block - 1)*block_length + t + 1)*cmplx(cos(angle), sin(angle), dp)/ &
                  sqrt(real(block_length, dp))
            end do
         end do
      end do
      do k = 1, block_length - 1
         power(k) = sum(abs(transform(k, :))**2)/real(n, dp)
      end do
      allocate(out%bispectrum(half_length, half_length), out%statistic(half_length, half_length))
      allocate(out%p_values(half_length, half_length), out%frequency(half_length))
      do k = 1, half_length
         out%frequency(k) = 0.5_dp*real(k - 1, dp)/real(max(1, half_length - 1), dp)
      end do
      do k1 = 1, half_length
         do k2 = 1, half_length
            out%bispectrum(k1, k2) = sum(transform(k1, :)*transform(k2, :)* &
               transform(block_length - (k1 + k2), :))/ &
               (sqrt(max(tiny(1.0_dp), power(k1)*power(k2)*power(k1 + k2)))*real(n, dp))
         end do
      end do
      scale = 2.0_dp*real(n, dp)**(-0.02_dp)
      out%statistic = scale*abs(out%bispectrum)**2
      out%noncentrality = sum(out%statistic)/real(size(out%statistic), dp)
      do k1 = 1, half_length
         do k2 = 1, half_length
            out%p_values(k1, k2) = noncentral_chisq2_q(out%statistic(k1, k2), out%noncentrality)
         end do
      end do
   end function test_linearity

   pure real(dp) function noncentral_chisq2_q(value, noncentrality) result(probability)
      ! Return a two-degree noncentral chi-squared upper-tail probability.
      real(dp), intent(in) :: value, noncentrality
      real(dp) :: weight, term, half_lambda
      integer :: index

      half_lambda = 0.5_dp*max(0.0_dp, noncentrality)
      weight = exp(-half_lambda)
      probability = weight*regularized_gamma_q(1.0_dp, 0.5_dp*value)
      do index = 1, 10000
         weight = weight*half_lambda/real(index, dp)
         term = weight*regularized_gamma_q(real(index + 1, dp), 0.5_dp*value)
         probability = probability + term
         if (term <= epsilon(1.0_dp)*max(probability, tiny(1.0_dp))) exit
      end do
      probability = max(0.0_dp, min(1.0_dp, probability))
   end function noncentral_chisq2_q

   pure function difference_series(series, differences) result(value)
      ! Apply repeated first differences to a univariate series.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: differences
      real(dp), allocatable :: value(:), work(:)
      integer :: i

      value = series
      do i = 1, differences
         work = value(2:) - value(:size(value) - 1)
         value = work
      end do
   end function difference_series

   pure function sarima_likelihood(series, ar, d, ma, sar, seasonal_difference, sma, season, &
      intercept, drift, regressors, regression_coefficients) result(out)
      ! Evaluate a conditional Gaussian SARIMA likelihood for supplied parameters.
      real(dp), intent(in) :: series(:)
      real(dp), intent(in), optional :: ar(:), ma(:), sar(:), sma(:)
      integer, intent(in), optional :: d, seasonal_difference, season
      real(dp), intent(in), optional :: intercept, drift
      real(dp), intent(in), optional :: regressors(:, :), regression_coefficients(:)
      type(astsa_sarima_likelihood_t) :: out
      type(astsa_arma_diagnostic_t) :: diagnostic
      real(dp), allocatable :: adjusted(:), ar_base(:), ma_base(:), ar_seasonal(:), ma_seasonal(:)
      real(dp), allocatable :: ar_polynomial(:), ma_polynomial(:), ar_coefficient(:), ma_coefficient(:)
      real(dp) :: residual_sum
      integer :: ordinary_d, seasonal_d, period, p, q, ps, qs, n, start, i, j, parameter_count

      ordinary_d = 0
      if (present(d)) ordinary_d = d
      seasonal_d = 0
      if (present(seasonal_difference)) seasonal_d = seasonal_difference
      period = 1
      if (present(season)) period = season
      p = 0
      if (present(ar)) p = size(ar)
      q = 0
      if (present(ma)) q = size(ma)
      ps = 0
      if (present(sar)) ps = size(sar)
      qs = 0
      if (present(sma)) qs = size(sma)
      if (size(series) < 2 .or. ordinary_d < 0 .or. seasonal_d < 0 .or. period < 1) then
         out%info = 1
         return
      end if
      if ((ps > 0 .or. qs > 0 .or. seasonal_d > 0) .and. .not. present(season)) then
         out%info = 2
         return
      end if
      if (present(regressors) .neqv. present(regression_coefficients)) then
         out%info = 3
         return
      end if
      if (present(regressors)) then
         if (size(regressors, 1) /= size(series) .or. &
            size(regressors, 2) /= size(regression_coefficients)) then
            out%info = 3
            return
         end if
      end if
      diagnostic = arma_check(ar, ma, sar, sma, season)
      if (diagnostic%info /= 0 .or. .not. diagnostic%causal .or. .not. diagnostic%invertible) then
         out%info = 4
         return
      end if

      adjusted = series
      if (present(intercept)) adjusted = adjusted - intercept
      if (present(drift)) then
         do i = 1, size(adjusted)
            adjusted(i) = adjusted(i) - drift*real(i, dp)
         end do
      end if
      if (present(regressors)) adjusted = adjusted - matmul(regressors, regression_coefficients)
      adjusted = difference_series(adjusted, ordinary_d)
      do i = 1, seasonal_d
         if (size(adjusted) <= period) then
            out%info = 5
            return
         end if
         adjusted = adjusted(period + 1:) - adjusted(:size(adjusted) - period)
      end do
      n = size(adjusted)
      if (n < 1) then
         out%info = 5
         return
      end if

      allocate(ar_base(p + 1), ma_base(q + 1), ar_seasonal(ps*period + 1), ma_seasonal(qs*period + 1))
      ar_base = 0.0_dp
      ma_base = 0.0_dp
      ar_seasonal = 0.0_dp
      ma_seasonal = 0.0_dp
      ar_base(1) = 1.0_dp
      ma_base(1) = 1.0_dp
      ar_seasonal(1) = 1.0_dp
      ma_seasonal(1) = 1.0_dp
      if (p > 0) ar_base(2:) = -ar
      if (q > 0) ma_base(2:) = ma
      do i = 1, ps
         ar_seasonal(i*period + 1) = -sar(i)
      end do
      do i = 1, qs
         ma_seasonal(i*period + 1) = sma(i)
      end do
      ar_polynomial = poly_mul(ar_base, ar_seasonal)
      ma_polynomial = poly_mul(ma_base, ma_seasonal)
      ar_coefficient = -ar_polynomial(2:)
      ma_coefficient = ma_polynomial(2:)
      out%transformed = adjusted
      allocate(out%fitted(n), out%residuals(n))
      out%fitted = 0.0_dp
      out%residuals = 0.0_dp
      do i = 1, n
         do j = 1, min(size(ar_coefficient), i - 1)
            out%fitted(i) = out%fitted(i) + ar_coefficient(j)*adjusted(i - j)
         end do
         do j = 1, min(size(ma_coefficient), i - 1)
            out%fitted(i) = out%fitted(i) + ma_coefficient(j)*out%residuals(i - j)
         end do
         out%residuals(i) = adjusted(i) - out%fitted(i)
      end do
      start = max(size(ar_coefficient), size(ma_coefficient)) + 1
      if (start > n) then
         out%info = 6
         return
      end if
      out%observations = n - start + 1
      residual_sum = sum(out%residuals(start:)**2)
      out%sigma2 = residual_sum/real(out%observations, dp)
      if (out%sigma2 <= tiny(1.0_dp)) then
         out%info = 7
         return
      end if
      out%log_likelihood = -0.5_dp*real(out%observations, dp)* &
         (log(2.0_dp*acos(-1.0_dp)*out%sigma2) + 1.0_dp)
      parameter_count = p + q + ps + qs + 1
      if (present(intercept)) parameter_count = parameter_count + 1
      if (present(drift)) parameter_count = parameter_count + 1
      if (present(regression_coefficients)) parameter_count = parameter_count + size(regression_coefficients)
      out%parameters = parameter_count
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameter_count, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(out%observations, dp))*real(parameter_count, dp)
      out%aicc = out%aic
      if (out%observations > parameter_count + 1) then
         out%aicc = out%aic + 2.0_dp*real(parameter_count*(parameter_count + 1), dp)/ &
            real(out%observations - parameter_count - 1, dp)
      end if
   end function sarima_likelihood

   pure function sarima_exact_likelihood(series, ar, ma, sar, sma, season, intercept, drift, &
      regressors, regression_coefficients, d, seasonal_difference) result(out)
      ! Evaluate an exact stationary Gaussian SARIMA likelihood by Kalman filtering.
      real(dp), intent(in) :: series(:)
      real(dp), intent(in), optional :: ar(:), ma(:), sar(:), sma(:)
      integer, intent(in), optional :: season
      real(dp), intent(in), optional :: intercept, drift
      real(dp), intent(in), optional :: regressors(:, :), regression_coefficients(:)
      integer, intent(in), optional :: d, seasonal_difference
      type(astsa_sarima_likelihood_t) :: out
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(astsa_arma_diagnostic_t) :: diagnostic
      real(dp), allocatable :: ar_values(:), ma_values(:), sar_values(:), sma_values(:)
      real(dp), allocatable :: ar_polynomial(:), ma_polynomial(:), covariance(:, :), next_covariance(:, :)
      real(dp), allocatable :: adjusted(:), stationary_transition(:, :), stationary_loading(:)
      real(dp), allocatable :: stationary_covariance(:, :), difference_polynomial(:), factor(:)
      real(dp) :: scale_sum, log_variance_sum, convergence
      integer :: p, q, ps, qs, period, state_size, stationary_size, diffuse_count, ordinary_d, seasonal_d, i, n

      p = 0
      if (present(ar)) p = size(ar)
      q = 0
      if (present(ma)) q = size(ma)
      ps = 0
      if (present(sar)) ps = size(sar)
      qs = 0
      if (present(sma)) qs = size(sma)
      period = 1
      if (present(season)) period = season
      ordinary_d = 0
      if (present(d)) ordinary_d = d
      seasonal_d = 0
      if (present(seasonal_difference)) seasonal_d = seasonal_difference
      if (ordinary_d < 0 .or. seasonal_d < 0) then
         out%info = 1
         return
      end if
      allocate(ar_values(p), ma_values(q), sar_values(ps), sma_values(qs))
      if (p > 0) ar_values = ar
      if (q > 0) ma_values = ma
      if (ps > 0) sar_values = sar
      if (qs > 0) sma_values = sma
      if ((ps > 0 .or. qs > 0) .and. .not. present(season)) then
         out%info = 1
         return
      end if
      if (present(regressors) .neqv. present(regression_coefficients)) then
         out%info = 2
         return
      end if
      if (present(regressors)) then
         if (size(regressors, 1) /= size(series) .or. &
            size(regressors, 2) /= size(regression_coefficients)) then
            out%info = 2
            return
         end if
      end if
      if (ps > 0 .or. qs > 0) then
         diagnostic = arma_check(ar_values, ma_values, sar_values, sma_values, period)
      else
         diagnostic = arma_check(ar=ar_values, ma=ma_values)
      end if
      if (diagnostic%info /= 0 .or. .not. diagnostic%causal .or. .not. diagnostic%invertible) then
         out%info = 3
         return
      end if
      adjusted = series
      if (present(intercept)) adjusted = adjusted - intercept
      if (present(drift)) then
         do i = 1, size(adjusted)
            adjusted(i) = adjusted(i) - drift*real(i, dp)
         end do
      end if
      if (present(regressors)) adjusted = adjusted - matmul(regressors, regression_coefficients)
      call build_sarima_polynomials(ar_values, ma_values, sar_values, sma_values, period, &
         ar_polynomial, ma_polynomial)
      stationary_size = max(1, max(size(ar_polynomial) - 1, size(ma_polynomial)))
      difference_polynomial = [1.0_dp]
      do i = 1, ordinary_d
         difference_polynomial = poly_mul(difference_polynomial, [1.0_dp, -1.0_dp])
      end do
      allocate(factor(period + 1))
      factor = 0.0_dp
      factor(1) = 1.0_dp
      factor(period + 1) = -1.0_dp
      do i = 1, seasonal_d
         difference_polynomial = poly_mul(difference_polynomial, factor)
      end do
      diffuse_count = size(difference_polynomial) - 1
      state_size = stationary_size + diffuse_count
      n = size(series)
      allocate(model%y(n, 1), model%z(1, state_size, 1), model%h(1, 1, 1))
      allocate(model%transition(state_size, state_size, 1), model%r(state_size, 1, 1), model%q(1, 1, 1))
      allocate(model%a1(state_size), model%p1(state_size, state_size), model%p1inf(state_size, state_size))
      allocate(model%missing(n, 1), covariance(state_size, state_size), next_covariance(state_size, state_size))
      model%y(:, 1) = adjusted
      model%z = 0.0_dp
      model%z(1, 1, 1) = 1.0_dp
      model%h = 0.0_dp
      model%transition = 0.0_dp
      if (size(ar_polynomial) > 1) model%transition(:size(ar_polynomial) - 1, 1, 1) = -ar_polynomial(2:)
      do i = 1, stationary_size - 1
         model%transition(i, i + 1, 1) = 1.0_dp
      end do
      model%r = 0.0_dp
      model%r(1, 1, 1) = 1.0_dp
      if (size(ma_polynomial) > 1) model%r(2:size(ma_polynomial), 1, 1) = ma_polynomial(2:)
      model%q = 1.0_dp
      model%a1 = 0.0_dp
      model%p1inf = 0.0_dp
      model%missing = .false.
      covariance = 0.0_dp
      do i = 1, 10000
         next_covariance = matmul(matmul(model%transition(:, :, 1), covariance), &
            transpose(model%transition(:, :, 1))) + matmul(model%r(:, :, 1), transpose(model%r(:, :, 1)))
         convergence = maxval(abs(next_covariance - covariance))
         covariance = next_covariance
         if (convergence <= 100.0_dp*epsilon(1.0_dp)) exit
      end do
      if (convergence > 1.0e-8_dp) then
         out%info = 4
         return
      end if
      if (diffuse_count > 0) then
         stationary_transition = model%transition(:stationary_size, :stationary_size, 1)
         stationary_loading = model%r(:stationary_size, 1, 1)
         stationary_covariance = covariance(:stationary_size, :stationary_size)
         model%transition = 0.0_dp
         model%r = 0.0_dp
         model%z = 0.0_dp
         model%transition(:stationary_size, :stationary_size, 1) = stationary_transition
         model%r(:stationary_size, 1, 1) = stationary_loading
         model%transition(stationary_size + 1, :stationary_size, 1) = stationary_transition(1, :)
         model%transition(stationary_size + 1, stationary_size + 1:, 1) = -difference_polynomial(2:)
         model%r(stationary_size + 1, 1, 1) = stationary_loading(1)
         do i = 1, diffuse_count - 1
            model%transition(stationary_size + 1 + i, stationary_size + i, 1) = 1.0_dp
         end do
         model%z(1, stationary_size + 1, 1) = 1.0_dp
         model%p1 = 0.0_dp
         model%p1(:stationary_size, :stationary_size) = stationary_covariance
         model%p1inf = 0.0_dp
         do i = 1, diffuse_count
            model%p1inf(stationary_size + i, stationary_size + i) = 1.0_dp
         end do
         filtered = kfs_filter_diffuse(model)
      else
         model%p1 = covariance
         filtered = kfs_filter(model)
      end if
      if (filtered%info /= 0) then
         out%info = 10 + filtered%info
         return
      end if
      scale_sum = sum(filtered%innovation(diffuse_count + 1:, 1)**2/ &
         filtered%innovation_cov(1, 1, diffuse_count + 1:))
      out%observations = n - diffuse_count
      out%sigma2 = scale_sum/real(out%observations, dp)
      log_variance_sum = sum(log(filtered%innovation_cov(1, 1, diffuse_count + 1:)))
      out%log_likelihood = -0.5_dp*(real(out%observations, dp)* &
         (log(2.0_dp*acos(-1.0_dp)*out%sigma2) + 1.0_dp) + &
         log_variance_sum)
      out%transformed = adjusted
      out%residuals = filtered%innovation(:, 1)
      out%fitted = adjusted - out%residuals
      out%parameters = p + q + ps + qs + 1
      if (present(intercept)) out%parameters = out%parameters + 1
      if (present(drift)) out%parameters = out%parameters + 1
      if (present(regression_coefficients)) out%parameters = out%parameters + size(regression_coefficients)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(out%parameters, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n, dp))*real(out%parameters, dp)
      out%aicc = out%aic
      if (n > out%parameters + 1) out%aicc = out%aic + &
         2.0_dp*real(out%parameters*(out%parameters + 1), dp)/real(n - out%parameters - 1, dp)
   end function sarima_exact_likelihood

   pure function sarima_fit(series, p, d, q, seasonal_p, seasonal_difference, seasonal_q, &
      season, initial, include_intercept, max_iterations, tolerance, include_drift, &
      regressors, estimated, transform_parameters, exact_likelihood) result(out)
      ! Estimate a conditional Gaussian SARIMA model with finite-difference BFGS.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: p, d, q, seasonal_p, seasonal_difference, seasonal_q, season
      real(dp), intent(in), optional :: initial(:)
      logical, intent(in), optional :: include_intercept, include_drift, estimated(:), transform_parameters
      logical, intent(in), optional :: exact_likelihood
      real(dp), intent(in), optional :: regressors(:, :)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(astsa_sarima_fit_t) :: out
      type(optimization_result_t) :: optimization
      logical :: use_intercept, use_drift, use_transform, use_exact
      logical, allocatable :: estimate_mask(:)
      integer, allocatable :: free_index(:)
      real(dp), allocatable :: free_initial(:), hessian(:, :), inverse(:, :), optimizer_parameters(:), jacobian(:, :)
      integer :: count, limit, regression_count, free_count, dynamic_count, i, status
      real(dp) :: gradient_tolerance, objective_value

      use_intercept = .false.
      if (present(include_intercept)) use_intercept = include_intercept
      use_drift = .false.
      if (present(include_drift)) use_drift = include_drift
      use_exact = .false.
      if (present(exact_likelihood)) use_exact = exact_likelihood
      regression_count = 0
      if (present(regressors)) regression_count = size(regressors, 2)
      out%p = p
      out%d = d
      out%q = q
      out%seasonal_p = seasonal_p
      out%seasonal_difference = seasonal_difference
      out%seasonal_q = seasonal_q
      out%season = season
      out%regression_count = regression_count
      out%includes_intercept = use_intercept
      out%includes_drift = use_drift
      out%uses_exact_likelihood = use_exact
      count = p + q + seasonal_p + seasonal_q + regression_count
      if (use_intercept) count = count + 1
      if (use_drift) count = count + 1
      if (p < 0 .or. d < 0 .or. q < 0 .or. seasonal_p < 0 .or. &
         seasonal_difference < 0 .or. seasonal_q < 0) then
         out%info = 1
         return
      end if
      if ((seasonal_p > 0 .or. seasonal_q > 0 .or. seasonal_difference > 0) .and. season < 2) then
         out%info = 1
         return
      end if
      allocate(out%coefficients(count))
      out%coefficients = 0.0_dp
      if (present(initial)) then
         if (size(initial) /= count) then
            out%info = 2
            return
         end if
         out%coefficients = initial
      end if
      allocate(estimate_mask(count))
      estimate_mask = .true.
      if (present(estimated)) then
         if (size(estimated) /= count) then
            out%info = 2
            return
         end if
         estimate_mask = estimated
      end if
      free_count = sum(merge(1, 0, estimate_mask))
      dynamic_count = p + q + seasonal_p + seasonal_q
      use_transform = .true.
      if (present(transform_parameters)) use_transform = transform_parameters
      if (dynamic_count > 0) use_transform = use_transform .and. all(estimate_mask(:dynamic_count))
      allocate(optimizer_parameters(count))
      optimizer_parameters = out%coefficients
      if (use_transform) optimizer_parameters = inverse_sarima_transform(optimizer_parameters)
      allocate(free_index(free_count), free_initial(free_count))
      free_index = pack([(i, i=1, count)], estimate_mask)
      free_initial = pack(optimizer_parameters, estimate_mask)
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      gradient_tolerance = 1.0e-6_dp
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_initial, limit, gradient_tolerance)
         optimizer_parameters(free_index) = optimization%parameters
         out%coefficients = reported_parameters(optimizer_parameters)
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         out%info = optimization%info
      else
         out%converged = .true.
      end if
      out%likelihood = evaluate_parameters(out%coefficients)
      if (out%likelihood%info /= 0 .and. out%info == 0) out%info = 10 + out%likelihood%info
      out%likelihood%parameters = count + 1
      objective_value = -2.0_dp*out%likelihood%log_likelihood
      out%likelihood%aic = objective_value + 2.0_dp*real(count + 1, dp)
      out%likelihood%bic = objective_value + log(real(out%likelihood%observations, dp))*real(count + 1, dp)
      out%likelihood%aicc = out%likelihood%aic
      if (out%likelihood%observations > count + 2) then
         out%likelihood%aicc = out%likelihood%aic + 2.0_dp*real((count + 1)*(count + 2), dp)/ &
            real(out%likelihood%observations - count - 2, dp)
      end if
      allocate(out%covariance(count, count), out%standard_errors(count), out%statistics(count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      out%statistics = 0.0_dp
      if (free_count > 0 .and. out%likelihood%info == 0) then
         free_initial = pack(optimizer_parameters, estimate_mask)
         hessian = finite_difference_hessian(objective, free_initial)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            jacobian = coefficient_jacobian(free_initial)
            out%covariance = matmul(matmul(jacobian, inverse), transpose(jacobian))
            do i = 1, count
               if (out%covariance(i, i) > 0.0_dp) out%standard_errors(i) = sqrt(out%covariance(i, i))
            end do
            where (out%standard_errors > 0.0_dp)
               out%statistics = out%coefficients/out%standard_errors
            end where
         else if (out%info == 0) then
            out%info = 20 + status
         end if
      end if

   contains

      pure function objective(parameters) result(value)
         ! Return the penalized negative profile log likelihood.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value
         type(astsa_sarima_likelihood_t) :: evaluated

         evaluated = evaluate_parameters(reported_parameters(expand_parameters(parameters)))
         if (evaluated%info == 0) then
            value = -evaluated%log_likelihood
         else
            value = 1.0e30_dp + dot_product(parameters, parameters)
         end if
      end function objective

      pure function expand_parameters(free_parameters) result(parameters)
         ! Insert estimated parameters into the full parameter vector.
         real(dp), intent(in) :: free_parameters(:)
         real(dp) :: parameters(count)

         parameters = optimizer_parameters
         parameters(free_index) = free_parameters
      end function expand_parameters

      pure function reported_parameters(parameters) result(reported)
         ! Map optimizer parameters to reported causal and invertible coefficients.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: reported(size(parameters))

         reported = parameters
         if (.not. use_transform) return
         if (p > 0) reported(1:p) = pacf_to_ar(parameters(1:p))
         if (q > 0) reported(p + 1:p + q) = -pacf_to_ar(parameters(p + 1:p + q))
         if (seasonal_p > 0) reported(p + q + 1:p + q + seasonal_p) = &
            pacf_to_ar(parameters(p + q + 1:p + q + seasonal_p))
         if (seasonal_q > 0) reported(p + q + seasonal_p + 1:dynamic_count) = &
            -pacf_to_ar(parameters(p + q + seasonal_p + 1:dynamic_count))
      end function reported_parameters

      pure function inverse_sarima_transform(parameters) result(unconstrained)
         ! Map reported causal and invertible coefficients to optimizer coordinates.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: unconstrained(size(parameters))

         unconstrained = parameters
         if (p > 0) unconstrained(1:p) = ar_to_unconstrained(parameters(1:p))
         if (q > 0) unconstrained(p + 1:p + q) = ar_to_unconstrained(-parameters(p + 1:p + q))
         if (seasonal_p > 0) unconstrained(p + q + 1:p + q + seasonal_p) = &
            ar_to_unconstrained(parameters(p + q + 1:p + q + seasonal_p))
         if (seasonal_q > 0) unconstrained(p + q + seasonal_p + 1:dynamic_count) = &
            ar_to_unconstrained(-parameters(p + q + seasonal_p + 1:dynamic_count))
      end function inverse_sarima_transform

      pure function coefficient_jacobian(free_parameters) result(jacobian_value)
         ! Numerically map optimizer covariance to the reported coefficient scale.
         real(dp), intent(in) :: free_parameters(:)
         real(dp) :: jacobian_value(count, size(free_parameters))
         real(dp) :: plus(size(free_parameters)), minus(size(free_parameters)), step
         integer :: column

         do column = 1, size(free_parameters)
            step = epsilon(1.0_dp)**0.25_dp*max(1.0_dp, abs(free_parameters(column)))
            plus = free_parameters
            minus = free_parameters
            plus(column) = plus(column) + step
            minus(column) = minus(column) - step
            jacobian_value(:, column) = (reported_parameters(expand_parameters(plus)) - &
               reported_parameters(expand_parameters(minus)))/(2.0_dp*step)
         end do
      end function coefficient_jacobian

      pure function evaluate_parameters(parameters) result(evaluated)
         ! Unpack an optimizer vector and evaluate its SARIMA likelihood.
         real(dp), intent(in) :: parameters(:)
         type(astsa_sarima_likelihood_t) :: evaluated
         real(dp), allocatable :: ar_values(:), ma_values(:), sar_values(:), sma_values(:), beta(:), xreg(:, :)
         real(dp) :: intercept_value, drift_value
         integer :: offset

         offset = 0
         ar_values = parameters(offset + 1:offset + p)
         offset = offset + p
         ma_values = parameters(offset + 1:offset + q)
         offset = offset + q
         sar_values = parameters(offset + 1:offset + seasonal_p)
         offset = offset + seasonal_p
         sma_values = parameters(offset + 1:offset + seasonal_q)
         offset = offset + seasonal_q
         intercept_value = 0.0_dp
         if (use_intercept) then
            intercept_value = parameters(offset + 1)
            offset = offset + 1
         end if
         drift_value = 0.0_dp
         if (use_drift) then
            drift_value = parameters(offset + 1)
            offset = offset + 1
         end if
         allocate(xreg(size(series), regression_count), beta(regression_count))
         if (regression_count > 0) then
            xreg = regressors
            beta = parameters(offset + 1:offset + regression_count)
         end if
         if (seasonal_p > 0 .or. seasonal_q > 0 .or. seasonal_difference > 0) then
            if (use_exact) then
               evaluated = sarima_exact_likelihood(series, ar_values, ma_values, sar_values, &
                  sma_values, season, intercept_value, drift_value, xreg, beta, d, seasonal_difference)
            else
               evaluated = sarima_likelihood(series, ar_values, d, ma_values, sar_values, &
                  seasonal_difference, sma_values, season, intercept_value, drift_value, xreg, beta)
            end if
         else
            if (use_exact) then
               evaluated = sarima_exact_likelihood(series, ar=ar_values, ma=ma_values, &
                  intercept=intercept_value, drift=drift_value, regressors=xreg, &
                  regression_coefficients=beta, d=d, seasonal_difference=seasonal_difference)
            else
               evaluated = sarima_likelihood(series, ar=ar_values, d=d, ma=ma_values, &
                  intercept=intercept_value, drift=drift_value, regressors=xreg, &
                  regression_coefficients=beta)
            end if
         end if
      end function evaluate_parameters
   end function sarima_fit

   pure function pacf_to_ar(unconstrained) result(coefficients)
      ! Map unconstrained partial autocorrelations to stationary AR coefficients.
      real(dp), intent(in) :: unconstrained(:)
      real(dp) :: coefficients(size(unconstrained)), previous(size(unconstrained)), reflection
      integer :: order, j

      coefficients = 0.0_dp
      do order = 1, size(unconstrained)
         reflection = tanh(unconstrained(order))
         previous = coefficients
         do j = 1, order - 1
            coefficients(j) = previous(j) - reflection*previous(order - j)
         end do
         coefficients(order) = reflection
      end do
   end function pacf_to_ar

   pure function ar_to_unconstrained(coefficients) result(unconstrained)
      ! Map stationary AR coefficients to unconstrained partial correlations.
      real(dp), intent(in) :: coefficients(:)
      real(dp) :: unconstrained(size(coefficients)), work(size(coefficients)), previous(size(coefficients))
      real(dp) :: reflection, denominator
      integer :: order, j

      work = coefficients
      unconstrained = 0.0_dp
      do order = size(coefficients), 2, -1
         reflection = max(-0.999999_dp, min(0.999999_dp, work(order)))
         unconstrained(order) = atanh(reflection)
         denominator = 1.0_dp - reflection*reflection
         previous = work
         do j = 1, order - 1
            work(j) = (previous(j) + reflection*previous(order - j))/denominator
         end do
      end do
      if (size(coefficients) > 0) then
         reflection = max(-0.999999_dp, min(0.999999_dp, work(1)))
         unconstrained(1) = atanh(reflection)
      end if
   end function ar_to_unconstrained

   pure function sarima_forecast(fit, series, horizon, regressors, new_regressors, &
      interval_multiplier) result(out)
      ! Forecast a fitted SARIMA model and return Gaussian prediction intervals.
      type(astsa_sarima_fit_t), intent(in) :: fit
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: horizon
      real(dp), intent(in), optional :: regressors(:, :), new_regressors(:, :)
      real(dp), intent(in), optional :: interval_multiplier
      type(astsa_sarima_forecast_t) :: out
      real(dp), allocatable :: ar(:), ma(:), sar(:), sma(:), beta(:)
      real(dp), allocatable :: ar_base(:), ma_base(:), ar_seasonal(:), ma_seasonal(:)
      real(dp), allocatable :: ar_polynomial(:), ma_polynomial(:), ar_coefficient(:), ma_coefficient(:)
      real(dp), allocatable :: difference_polynomial(:), factor(:), denominator(:), impulse(:)
      real(dp), allocatable :: adjusted(:), adjusted_all(:), transformed_all(:), residual_all(:)
      real(dp) :: intercept_value, drift_value, variance
      integer :: offset, n, nt, i, j, k, absolute_time

      if (fit%info /= 0 .or. horizon < 1 .or. size(series) < 1) then
         out%info = 1
         return
      end if
      if (fit%regression_count > 0) then
         if (.not. present(regressors) .or. .not. present(new_regressors)) then
            out%info = 2
            return
         end if
         if (size(regressors, 1) /= size(series) .or. size(regressors, 2) /= fit%regression_count .or. &
            size(new_regressors, 1) /= horizon .or. size(new_regressors, 2) /= fit%regression_count) then
            out%info = 2
            return
         end if
      end if
      offset = 0
      ar = fit%coefficients(offset + 1:offset + fit%p)
      offset = offset + fit%p
      ma = fit%coefficients(offset + 1:offset + fit%q)
      offset = offset + fit%q
      sar = fit%coefficients(offset + 1:offset + fit%seasonal_p)
      offset = offset + fit%seasonal_p
      sma = fit%coefficients(offset + 1:offset + fit%seasonal_q)
      offset = offset + fit%seasonal_q
      intercept_value = 0.0_dp
      if (fit%includes_intercept) then
         intercept_value = fit%coefficients(offset + 1)
         offset = offset + 1
      end if
      drift_value = 0.0_dp
      if (fit%includes_drift) then
         drift_value = fit%coefficients(offset + 1)
         offset = offset + 1
      end if
      beta = fit%coefficients(offset + 1:offset + fit%regression_count)

      n = size(series)
      adjusted = series - intercept_value
      do i = 1, n
         adjusted(i) = adjusted(i) - drift_value*real(i, dp)
      end do
      if (fit%regression_count > 0) adjusted = adjusted - matmul(regressors, beta)
      call build_sarima_polynomials(ar, ma, sar, sma, fit%season, ar_polynomial, ma_polynomial)
      ar_coefficient = -ar_polynomial(2:)
      ma_coefficient = ma_polynomial(2:)
      nt = size(fit%likelihood%transformed)
      allocate(transformed_all(nt + horizon), residual_all(nt + horizon))
      transformed_all(:nt) = fit%likelihood%transformed
      residual_all(:nt) = fit%likelihood%residuals
      residual_all(nt + 1:) = 0.0_dp
      do i = nt + 1, nt + horizon
         transformed_all(i) = 0.0_dp
         do j = 1, min(size(ar_coefficient), i - 1)
            transformed_all(i) = transformed_all(i) + ar_coefficient(j)*transformed_all(i - j)
         end do
         do j = 1, min(size(ma_coefficient), i - 1)
            transformed_all(i) = transformed_all(i) + ma_coefficient(j)*residual_all(i - j)
         end do
      end do

      difference_polynomial = [1.0_dp]
      do i = 1, fit%d
         difference_polynomial = poly_mul(difference_polynomial, [1.0_dp, -1.0_dp])
      end do
      allocate(factor(fit%season + 1))
      factor = 0.0_dp
      factor(1) = 1.0_dp
      factor(fit%season + 1) = -1.0_dp
      do i = 1, fit%seasonal_difference
         difference_polynomial = poly_mul(difference_polynomial, factor)
      end do
      allocate(adjusted_all(n + horizon))
      adjusted_all(:n) = adjusted
      do k = 1, horizon
         absolute_time = n + k
         adjusted_all(absolute_time) = transformed_all(nt + k)
         do j = 1, size(difference_polynomial) - 1
            adjusted_all(absolute_time) = adjusted_all(absolute_time) - &
               difference_polynomial(j + 1)*adjusted_all(absolute_time - j)
         end do
      end do

      allocate(out%mean(horizon), out%standard_error(horizon), out%lower(horizon), out%upper(horizon))
      out%interval_multiplier = 2.0_dp
      if (present(interval_multiplier)) out%interval_multiplier = interval_multiplier
      do k = 1, horizon
         out%mean(k) = adjusted_all(n + k) + intercept_value + drift_value*real(n + k, dp)
      end do
      if (fit%regression_count > 0) out%mean = out%mean + matmul(new_regressors, beta)

      denominator = poly_mul(ar_polynomial, difference_polynomial)
      ar_coefficient = -denominator(2:)
      allocate(impulse(0:horizon - 1))
      impulse = 0.0_dp
      impulse(0) = 1.0_dp
      do k = 1, horizon - 1
         if (k <= size(ma_coefficient)) impulse(k) = ma_coefficient(k)
         do j = 1, min(size(ar_coefficient), k)
            impulse(k) = impulse(k) + ar_coefficient(j)*impulse(k - j)
         end do
      end do
      variance = 0.0_dp
      do k = 1, horizon
         variance = variance + fit%likelihood%sigma2*impulse(k - 1)**2
         out%standard_error(k) = sqrt(max(0.0_dp, variance))
      end do
      out%lower = out%mean - out%interval_multiplier*out%standard_error
      out%upper = out%mean + out%interval_multiplier*out%standard_error
   end function sarima_forecast

   pure function sarima_diagnostics(fit, max_lag) result(out)
      ! Compute numerical residual diagnostics for a fitted SARIMA model.
      type(astsa_sarima_fit_t), intent(in) :: fit
      integer, intent(in), optional :: max_lag
      type(astsa_sarima_diagnostics_t) :: out
      real(dp), allocatable :: residuals(:), correlations(:)
      real(dp) :: statistic, probability
      integer :: n, start, limit, parameter_count, first_lag, count_lags, h, i

      if (fit%info /= 0 .or. fit%likelihood%sigma2 <= 0.0_dp .or. &
         .not. allocated(fit%likelihood%residuals)) then
         out%info = 1
         return
      end if
      n = fit%likelihood%observations
      start = size(fit%likelihood%residuals) - n + 1
      if (n < 3 .or. start < 1) then
         out%info = 1
         return
      end if
      residuals = fit%likelihood%residuals(start:)
      out%standardized_residuals = residuals/sqrt(fit%likelihood%sigma2)
      parameter_count = fit%p + fit%q + fit%seasonal_p + fit%seasonal_q
      limit = min(merge(20, 3*fit%season, fit%season < 7), 52)
      limit = max(limit, parameter_count + 8)
      if (present(max_lag)) limit = max_lag
      limit = min(limit, n - 1)
      first_lag = parameter_count + 1
      if (limit < first_lag) then
         out%info = 2
         return
      end if
      correlations = acf_values(residuals, limit)
      out%residual_acf = correlations
      count_lags = limit - first_lag + 1
      allocate(out%lags(count_lags), out%degrees_of_freedom(count_lags))
      allocate(out%ljung_box(count_lags), out%p_values(count_lags))
      statistic = 0.0_dp
      i = 0
      do h = 1, limit
         statistic = statistic + real(n*(n + 2), dp)*correlations(h + 1)**2/real(n - h, dp)
         if (h >= first_lag) then
            i = i + 1
            out%lags(i) = h
            out%degrees_of_freedom(i) = h - parameter_count
            out%ljung_box(i) = statistic
            probability = regularized_gamma_q(0.5_dp*real(out%degrees_of_freedom(i), dp), 0.5_dp*statistic)
            out%p_values(i) = max(0.0_dp, min(1.0_dp, probability))
         end if
      end do
      out%qq_sample = sorted_values(out%standardized_residuals)
      allocate(out%qq_theoretical(n))
      do i = 1, n
         if (n <= 10) then
            probability = (real(i, dp) - 0.375_dp)/(real(n, dp) + 0.25_dp)
         else
            probability = (real(i, dp) - 0.5_dp)/real(n, dp)
         end if
         out%qq_theoretical(i) = inverse_standard_normal(probability)
      end do
   end function sarima_diagnostics

   pure function lag_reg(input, output, m, span, threshold, inverse) result(out)
      ! Estimate astsa LagReg's two-sided transfer coefficients and aligned fit.
      real(dp), intent(in) :: input(:), output(:)
      integer, intent(in), optional :: m, span
      real(dp), intent(in), optional :: threshold
      logical, intent(in), optional :: inverse
      type(astsa_lag_regression_t) :: out
      type(astsa_spectrum_t) :: spectrum
      real(dp), allocatable :: data(:, :)
      complex(dp), allocatable :: transfer(:)
      real(dp) :: input_mean, output_mean, cutoff, angle
      integer :: n, grid, half_grid, spectral_span, i, j, s, index, selected
      integer :: first_time, last_time
      logical :: inverse_direction

      n = size(input)
      if (size(output) /= n .or. n < 4) then
         out%info = 1
         return
      end if
      grid = min(40, 2*(n/2))
      if (present(m)) grid = m
      spectral_span = 3
      if (present(span)) spectral_span = span
      cutoff = 0.0_dp
      if (present(threshold)) cutoff = threshold
      inverse_direction = .false.
      if (present(inverse)) inverse_direction = inverse
      if (grid < 2 .or. mod(grid, 2) /= 0 .or. mod(n, grid) /= 0 .or. &
         spectral_span < 1 .or. cutoff < 0.0_dp) then
         out%info = 2
         return
      end if
      input_mean = sum(input)/real(n, dp)
      output_mean = sum(output)/real(n, dp)
      allocate(data(n, 2))
      data(:, 1) = output - output_mean
      data(:, 2) = input - input_mean
      spectrum = mv_periodogram(data, demean=.false., detrend=.false., span=spectral_span)
      if (spectrum%info /= 0) then
         out%info = 10 + spectrum%info
         return
      end if
      half_grid = grid/2
      allocate(transfer(half_grid), out%lags(grid - 1), out%coefficients(grid - 1))
      do j = 1, half_grid
         index = (n/grid)*j
         if (spectrum%spectrum(index, 2) <= tiny(1.0_dp)) then
            out%info = 3
            return
         end if
         transfer(j) = spectrum%spectral_matrix(1, 2, index)/spectrum%spectrum(index, 2)
      end do
      do i = 1, grid - 1
         s = -half_grid + i
         out%lags(i) = s
         out%coefficients(i) = 0.0_dp
         do j = 1, half_grid
            angle = 2.0_dp*acos(-1.0_dp)*real(j*s, dp)/real(grid, dp)
            out%coefficients(i) = out%coefficients(i) + &
               2.0_dp/real(grid, dp)*real(cmplx(cos(angle), sin(angle), dp)*transfer(j), dp)
         end do
      end do
      if (inverse_direction) then
         selected = count(out%lags < 0 .and. abs(out%coefficients) >= cutoff)
         allocate(out%selected_lags(selected), out%selected_coefficients(selected))
         out%selected_lags = pack(out%lags, out%lags < 0 .and. abs(out%coefficients) >= cutoff)
         out%selected_coefficients = pack(out%coefficients, &
            out%lags < 0 .and. abs(out%coefficients) >= cutoff)
      else
         selected = count(out%lags >= 0 .and. abs(out%coefficients) >= cutoff)
         allocate(out%selected_lags(selected), out%selected_coefficients(selected))
         out%selected_lags = pack(out%lags, out%lags >= 0 .and. abs(out%coefficients) >= cutoff)
         out%selected_coefficients = pack(out%coefficients, &
            out%lags >= 0 .and. abs(out%coefficients) >= cutoff)
      end if
      if (selected < 1) then
         out%info = 4
         return
      end if
      first_time = 1 + max(0, maxval(out%selected_lags))
      last_time = n + min(0, minval(out%selected_lags))
      if (last_time < first_time) then
         out%info = 5
         return
      end if
      out%intercept = output_mean - input_mean*sum(out%selected_coefficients)
      allocate(out%output(last_time - first_time + 1), out%fitted(last_time - first_time + 1))
      allocate(out%residuals(last_time - first_time + 1))
      do i = first_time, last_time
         out%output(i - first_time + 1) = output(i)
         out%fitted(i - first_time + 1) = out%intercept
         do j = 1, selected
            out%fitted(i - first_time + 1) = out%fitted(i - first_time + 1) + &
               out%selected_coefficients(j)*input(i - out%selected_lags(j))
         end do
      end do
      out%residuals = out%output - out%fitted
      out%mse = sum(out%residuals**2)/real(size(out%residuals), dp)
   end function lag_reg

   pure function signal_extract(series, m, max_frequency, min_frequency) result(out)
      ! Apply astsa SigExtract's tapered ideal low-frequency band-pass filter.
      real(dp), intent(in) :: series(:)
      integer, intent(in), optional :: m
      real(dp), intent(in), optional :: max_frequency, min_frequency
      type(astsa_signal_extraction_t) :: out
      real(dp) :: upper, lower, angle, desired, taper
      complex(dp) :: response
      integer :: grid, half_grid, n, i, j, s, t

      n = size(series)
      grid = 50
      if (present(m)) grid = m
      upper = 0.05_dp
      if (present(max_frequency)) upper = max_frequency
      lower = 0.01_dp
      if (present(min_frequency)) lower = min_frequency
      if (upper <= 0.0_dp .or. upper > 0.5_dp .or. lower < 0.0_dp .or. lower >= upper) then
         out%info = 1
         return
      end if
      if (upper < 1.0_dp/real(grid, dp)) grid = 2*floor(1.25_dp/upper)
      grid = 2*(grid/2)
      if (grid < 2 .or. n < grid - 1) then
         out%info = 2
         return
      end if
      half_grid = grid/2
      allocate(out%lags(grid - 1), out%coefficients(grid - 1), out%filtered(n))
      do i = 1, grid - 1
         s = -half_grid + i
         out%lags(i) = s
         out%coefficients(i) = 0.0_dp
         do j = 1, half_grid
            desired = merge(1.0_dp, 0.0_dp, real(j, dp)/real(grid, dp) > lower .and. &
               real(j, dp)/real(grid, dp) < upper)
            angle = 2.0_dp*acos(-1.0_dp)*real(j*s, dp)/real(grid, dp)
            out%coefficients(i) = out%coefficients(i) + &
               2.0_dp/real(grid, dp)*desired*cos(angle)
         end do
         taper = 0.5_dp*(1.0_dp + cos(2.0_dp*acos(-1.0_dp)*real(s, dp)/real(grid - 1, dp)))
         out%coefficients(i) = out%coefficients(i)*taper
      end do
      out%filtered = ieee_value(0.0_dp, ieee_quiet_nan)
      out%valid_start = half_grid
      out%valid_end = n - half_grid + 1
      do t = out%valid_start, out%valid_end
         out%filtered(t) = 0.0_dp
         do i = 1, grid - 1
            out%filtered(t) = out%filtered(t) + out%coefficients(i)*series(t - out%lags(i))
         end do
      end do
      allocate(out%frequency(n/2), out%desired_response(n/2), out%attained_response(n/2))
      do j = 1, n/2
         out%frequency(j) = real(j, dp)/real(n, dp)
         out%desired_response(j) = merge(1.0_dp, 0.0_dp, out%frequency(j) > lower .and. &
            out%frequency(j) < upper)
         response = cmplx(0.0_dp, 0.0_dp, dp)
         do i = 1, grid - 1
            angle = -2.0_dp*acos(-1.0_dp)*out%frequency(j)*real(out%lags(i), dp)
            response = response + out%coefficients(i)*cmplx(cos(angle), sin(angle), dp)
         end do
         out%attained_response(j) = real(response, dp)
      end do
   end function signal_extract

   pure function stochastic_regression(data, full_columns, reduced_columns, response_column, &
      span, m, alpha) result(out)
      ! Compute astsa stoch.reg's full and reduced spectral regressions.
      real(dp), intent(in) :: data(:, :)
      integer, intent(in) :: full_columns(:), reduced_columns(:), response_column
      integer, intent(in) :: span, m
      real(dp), intent(in), optional :: alpha
      type(astsa_stochastic_regression_t) :: out
      type(astsa_spectrum_t) :: spectrum
      complex(dp), allocatable :: xx(:, :), inverse(:, :), xy(:), yx(:), transfer(:, :)
      real(dp) :: test_size, numerator_mean, denominator_mean, angle
      integer :: n, frequencies, q, q1, q2, i, j, k, s, index, status

      n = size(data, 1)
      q = size(full_columns)
      q1 = size(reduced_columns)
      q2 = q - q1
      test_size = 0.05_dp
      if (present(alpha)) test_size = alpha
      if (n < 4 .or. q < 1 .or. q2 < 1 .or. response_column < 1 .or. &
         response_column > size(data, 2) .or. any(full_columns < 1) .or. &
         any(full_columns > size(data, 2)) .or. any(reduced_columns < 1) .or. &
         any(reduced_columns > size(data, 2)) .or. span <= q .or. mod(span, 2) == 0 .or. &
         m < 2 .or. mod(m, 2) /= 0 .or. mod(n, m) /= 0 .or. test_size <= 0.0_dp .or. test_size >= 1.0_dp) then
         out%info = 1
         return
      end if
      spectrum = mv_periodogram(data, demean=.true., detrend=.false., span=span)
      if (spectrum%info /= 0) then
         out%info = 10 + spectrum%info
         return
      end if
      frequencies = size(spectrum%frequency)
      allocate(out%frequency(frequencies), out%full_power(frequencies), out%reduced_power(frequencies))
      allocate(out%f_statistic(frequencies), out%coherence(frequencies), transfer(q, frequencies))
      out%frequency = spectrum%frequency
      do k = 1, frequencies
         allocate(xx(q, q), xy(q), yx(q))
         xx = spectrum%spectral_matrix(full_columns, full_columns, k)
         xy = spectrum%spectral_matrix(full_columns, response_column, k)
         yx = spectrum%spectral_matrix(response_column, full_columns, k)
         call invert_complex_matrix(xx, inverse, status)
         if (status /= 0) then
            out%info = 20 + k
            return
         end if
         out%full_power(k) = max(0.0_dp, real(spectrum%spectral_matrix(response_column, &
            response_column, k) - sum(yx*matmul(inverse, xy)), dp))
         transfer(:, k) = matmul(transpose(inverse), yx)
         deallocate(xx, xy, yx, inverse)
         if (q1 == 0) then
            out%reduced_power(k) = spectrum%spectrum(k, response_column)
         else
            allocate(xx(q1, q1), xy(q1), yx(q1))
            xx = spectrum%spectral_matrix(reduced_columns, reduced_columns, k)
            xy = spectrum%spectral_matrix(reduced_columns, response_column, k)
            yx = spectrum%spectral_matrix(response_column, reduced_columns, k)
            call invert_complex_matrix(xx, inverse, status)
            if (status /= 0) then
               out%info = 30 + k
               return
            end if
            out%reduced_power(k) = max(0.0_dp, real(spectrum%spectral_matrix(response_column, &
               response_column, k) - sum(yx*matmul(inverse, xy)), dp))
            deallocate(xx, xy, yx, inverse)
         end if
      end do
      out%numerator_df = 2*q2
      out%denominator_df = 2*(span - q)
      do k = 1, frequencies
         numerator_mean = real(span, dp)*(out%reduced_power(k) - out%full_power(k))/ &
            real(out%numerator_df, dp)
         denominator_mean = real(span, dp)*out%full_power(k)/real(out%denominator_df, dp)
         if (denominator_mean > tiny(1.0_dp)) then
            out%f_statistic(k) = max(0.0_dp, numerator_mean/denominator_mean)
            out%coherence(k) = out%f_statistic(k)/(out%f_statistic(k) + &
               real(out%denominator_df, dp)/real(out%numerator_df, dp))
         else
            out%f_statistic(k) = huge(1.0_dp)
            out%coherence(k) = 1.0_dp
         end if
      end do
      out%critical_f = f_distribution_quantile(1.0_dp - test_size, out%numerator_df, out%denominator_df)
      out%critical_coherence = out%critical_f/(out%critical_f + &
         real(out%denominator_df, dp)/real(out%numerator_df, dp))
      allocate(out%lags(m - 1), out%coefficients(m - 1, q))
      do i = 1, m - 1
         s = -m/2 + i
         out%lags(i) = s
         do j = 1, q
            out%coefficients(i, j) = 0.0_dp
            do k = 1, m/2
               index = (n/m)*k
               angle = 2.0_dp*acos(-1.0_dp)*real(k*s, dp)/real(m, dp)
               out%coefficients(i, j) = out%coefficients(i, j) + &
                  2.0_dp/real(m, dp)*real(cmplx(cos(angle), sin(angle), dp)*transfer(j, index), dp)
            end do
         end do
      end do
   end function stochastic_regression

   pure function sv_filter(returns, gamma, phi, state_sd, level, component_zero_sd, &
      component_one_mean, component_one_sd, rho) result(out)
      ! Run astsa SV.mle's two-component Gaussian-mixture volatility filter.
      real(dp), intent(in) :: returns(:), gamma, phi, state_sd, level
      real(dp), intent(in) :: component_zero_sd, component_one_mean, component_one_sd, rho
      type(astsa_sv_filter_t) :: out
      real(dp) :: state_variance, variance_zero, variance_one, phi_squared
      real(dp) :: covariance_zero, covariance_one, innovation_zero, innovation_one
      real(dp) :: gain_zero, gain_one, log_density_zero, log_density_one, maximum_log, denominator
      real(dp) :: updated_variance
      integer :: n, i

      n = size(returns)
      if (n < 2 .or. abs(phi) >= 1.0_dp .or. state_sd <= 0.0_dp .or. &
         component_zero_sd <= 0.0_dp .or. component_one_sd <= 0.0_dp .or. abs(rho) >= 1.0_dp) then
         out%info = 1
         return
      end if
      state_variance = state_sd**2
      variance_zero = component_zero_sd**2
      variance_one = component_one_sd**2
      phi_squared = phi**2
      covariance_zero = abs(state_sd*component_zero_sd)*rho
      covariance_one = abs(state_sd*component_one_sd)*rho
      allocate(out%predicted_log_variance(n), out%prediction_variance(n))
      allocate(out%component_zero_probability(n), out%component_one_probability(n))
      out%predicted_log_variance = 0.0_dp
      out%prediction_variance = 0.0_dp
      out%component_zero_probability = 0.5_dp
      out%component_one_probability = 0.5_dp
      out%prediction_variance(1) = phi_squared + state_variance
      do i = 2, n
         innovation_one = log(max(returns(i)**2, 1.0e-16_dp)) - &
            out%predicted_log_variance(i - 1) - component_one_mean - level
         innovation_zero = log(max(returns(i)**2, 1.0e-16_dp)) - &
            out%predicted_log_variance(i - 1) - level
         gain_one = (phi*out%prediction_variance(i - 1) + covariance_one)/ &
            (out%prediction_variance(i - 1) + variance_one)
         gain_zero = (phi*out%prediction_variance(i - 1) + covariance_zero)/ &
            (out%prediction_variance(i - 1) + variance_zero)
         log_density_one = -0.5_dp*(log(2.0_dp*acos(-1.0_dp)*(out%prediction_variance(i - 1) + &
            variance_one)) + innovation_one**2/(out%prediction_variance(i - 1) + variance_one))
         log_density_zero = -0.5_dp*(log(2.0_dp*acos(-1.0_dp)*(out%prediction_variance(i - 1) + &
            variance_zero)) + innovation_zero**2/(out%prediction_variance(i - 1) + variance_zero))
         maximum_log = max(log_density_zero, log_density_one)
         denominator = exp(log_density_zero - maximum_log) + exp(log_density_one - maximum_log)
         out%component_zero_probability(i) = exp(log_density_zero - maximum_log)/denominator
         out%component_one_probability(i) = exp(log_density_one - maximum_log)/denominator
         out%negative_log_likelihood = out%negative_log_likelihood - maximum_log - log(0.5_dp*denominator)
         out%predicted_log_variance(i) = gamma*returns(i - 1) + phi*out%predicted_log_variance(i - 1) + &
            out%component_zero_probability(i)*gain_zero*innovation_zero + &
            out%component_one_probability(i)*gain_one*innovation_one
         updated_variance = phi_squared*out%prediction_variance(i - 1) + state_variance - &
            out%component_zero_probability(i)*gain_zero**2*(out%prediction_variance(i - 1) + variance_zero) - &
            out%component_one_probability(i)*gain_one**2*(out%prediction_variance(i - 1) + variance_one)
         out%prediction_variance(i) = max(updated_variance, 1.0e-8_dp)
      end do
   end function sv_filter

   pure function sv_mle(returns, initial, feedback, leverage, max_iterations, tolerance) result(out)
      ! Estimate astsa's mixture stochastic-volatility model with BFGS.
      real(dp), intent(in) :: returns(:)
      real(dp), intent(in), optional :: initial(:)
      logical, intent(in), optional :: feedback, leverage
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(astsa_sv_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: optimizer_initial(:), hessian(:, :), inverse(:, :), jacobian(:, :)
      real(dp) :: gradient_tolerance
      integer :: count, limit, status, i

      out%feedback = .false.
      if (present(feedback)) out%feedback = feedback
      out%leverage = .false.
      if (present(leverage)) out%leverage = leverage
      if (out%leverage .and. .not. out%feedback) then
         out%info = 1
         return
      end if
      count = merge(7, 6, out%feedback)
      if (out%leverage) count = 8
      allocate(out%coefficients(count))
      if (out%feedback) then
         out%coefficients(:7) = [0.0_dp, 0.95_dp, 0.1_dp, &
            sum(log(max(returns**2, 1.0e-16_dp)))/real(size(returns), dp), &
            1.0_dp, -3.0_dp, 2.0_dp]
         if (out%leverage) out%coefficients(8) = 0.0_dp
      else
         out%coefficients = [0.95_dp, 0.1_dp, &
            sum(log(max(returns**2, 1.0e-16_dp)))/real(size(returns), dp), &
            1.0_dp, -3.0_dp, 2.0_dp]
      end if
      if (present(initial)) then
         if (size(initial) /= count) then
            out%info = 2
            return
         end if
         out%coefficients = initial
      end if
      optimizer_initial = sv_inverse_transform(out%coefficients)
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      gradient_tolerance = 1.0e-6_dp
      if (present(tolerance)) gradient_tolerance = tolerance
      optimization = bfgs_minimize_fd(objective, optimizer_initial, limit, gradient_tolerance)
      out%coefficients = sv_reported_parameters(optimization%parameters)
      out%iterations = optimization%iterations
      out%converged = optimization%converged
      out%info = optimization%info
      out%filter = evaluate(out%coefficients)
      if (out%filter%info /= 0 .and. out%info == 0) out%info = 10 + out%filter%info
      allocate(out%covariance(count, count), out%standard_errors(count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (out%filter%info == 0) then
         hessian = finite_difference_hessian(objective, optimization%parameters)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            jacobian = sv_parameter_jacobian(optimization%parameters)
            out%covariance = matmul(matmul(jacobian, inverse), transpose(jacobian))
            do i = 1, count
               if (out%covariance(i, i) > 0.0_dp) out%standard_errors(i) = sqrt(out%covariance(i, i))
            end do
         else if (out%info == 0) then
            out%info = 20 + status
         end if
      end if

   contains

      pure function objective(parameters) result(value)
         ! Return the mixture-filter negative log likelihood.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value
         type(astsa_sv_filter_t) :: filtered

         filtered = evaluate(sv_reported_parameters(parameters))
         if (filtered%info == 0) then
            value = filtered%negative_log_likelihood
         else
            value = 1.0e30_dp + dot_product(parameters, parameters)
         end if
      end function objective

      pure function evaluate(parameters) result(filtered)
         ! Unpack reported parameters into the mixture filter.
         real(dp), intent(in) :: parameters(:)
         type(astsa_sv_filter_t) :: filtered

         if (out%feedback) then
            if (out%leverage) then
               filtered = sv_filter(returns, parameters(1), parameters(2), parameters(3), parameters(4), &
                  parameters(5), parameters(6), parameters(7), parameters(8))
            else
               filtered = sv_filter(returns, parameters(1), parameters(2), parameters(3), parameters(4), &
                  parameters(5), parameters(6), parameters(7), 0.0_dp)
            end if
         else
            filtered = sv_filter(returns, 0.0_dp, parameters(1), parameters(2), parameters(3), &
               parameters(4), parameters(5), parameters(6), 0.0_dp)
         end if
      end function evaluate

      pure function sv_inverse_transform(parameters) result(unconstrained)
         ! Map reported volatility parameters to optimizer coordinates.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: unconstrained(size(parameters))
         integer :: phi_index, state_index, zero_index, one_index

         unconstrained = parameters
         phi_index = merge(2, 1, out%feedback)
         state_index = phi_index + 1
         zero_index = phi_index + 3
         one_index = phi_index + 5
         unconstrained(phi_index) = atanh(max(-0.999999_dp, min(0.999999_dp, parameters(phi_index))))
         unconstrained(state_index) = log(max(parameters(state_index), 1.0e-8_dp))
         unconstrained(zero_index) = log(max(parameters(zero_index), 1.0e-8_dp))
         unconstrained(one_index) = log(max(parameters(one_index), 1.0e-8_dp))
         if (out%leverage) unconstrained(8) = atanh(max(-0.999999_dp, min(0.999999_dp, parameters(8))))
      end function sv_inverse_transform

      pure function sv_reported_parameters(unconstrained) result(parameters)
         ! Map optimizer coordinates to constrained volatility parameters.
         real(dp), intent(in) :: unconstrained(:)
         real(dp) :: parameters(size(unconstrained))
         integer :: phi_index, state_index, zero_index, one_index

         parameters = unconstrained
         phi_index = merge(2, 1, out%feedback)
         state_index = phi_index + 1
         zero_index = phi_index + 3
         one_index = phi_index + 5
         parameters(phi_index) = tanh(unconstrained(phi_index))
         parameters(state_index) = exp(unconstrained(state_index))
         parameters(zero_index) = exp(unconstrained(zero_index))
         parameters(one_index) = exp(unconstrained(one_index))
         if (out%leverage) parameters(8) = tanh(unconstrained(8))
      end function sv_reported_parameters

      pure function sv_parameter_jacobian(parameters) result(jacobian_value)
         ! Compute the optimizer-to-reported parameter Jacobian numerically.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: jacobian_value(size(parameters), size(parameters))
         real(dp) :: plus(size(parameters)), minus(size(parameters)), step
         integer :: column

         do column = 1, size(parameters)
            step = epsilon(1.0_dp)**0.25_dp*max(1.0_dp, abs(parameters(column)))
            plus = parameters
            minus = parameters
            plus(column) = plus(column) + step
            minus(column) = minus(column) - step
            jacobian_value(:, column) = (sv_reported_parameters(plus) - &
               sv_reported_parameters(minus))/(2.0_dp*step)
         end do
      end function sv_parameter_jacobian
   end function sv_mle

   pure function sv_particle_filter_draws(observations, phi, state_variance, observation_scale, &
      conditioned_trajectory, normal_draws, resampling_uniforms, ancestor_uniforms, &
      terminal_uniform) result(out)
      ! Run a conditional particle filter with ancestor sampling from supplied draws.
      real(dp), intent(in) :: observations(:), phi, state_variance, observation_scale
      real(dp), intent(in) :: conditioned_trajectory(:), normal_draws(:, :)
      real(dp), intent(in) :: resampling_uniforms(:, :), ancestor_uniforms(:), terminal_uniform
      type(astsa_sv_particle_filter_t) :: out
      real(dp), allocatable :: log_weights(:), ancestor_weights(:), predicted(:)
      real(dp) :: maximum_log, total_weight, cumulative, uniform_value
      integer :: particles, periods, t, i, j, selected

      periods = size(observations)
      particles = size(normal_draws, 1)
      if (periods < 2 .or. particles < 2 .or. abs(phi) >= 1.0_dp .or. state_variance <= 0.0_dp .or. &
         observation_scale <= 0.0_dp .or. size(conditioned_trajectory) /= periods .or. &
         any(shape(normal_draws) /= [particles, periods - 1]) .or. &
         any(shape(resampling_uniforms) /= [particles, periods - 1]) .or. &
         size(ancestor_uniforms) /= periods - 1 .or. terminal_uniform < 0.0_dp .or. terminal_uniform > 1.0_dp) then
         out%info = 1
         return
      end if
      if (any(resampling_uniforms < 0.0_dp) .or. any(resampling_uniforms > 1.0_dp) .or. &
         any(ancestor_uniforms < 0.0_dp) .or. any(ancestor_uniforms > 1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%particles(particles, periods), out%weights(particles, periods))
      allocate(out%ancestors(particles, periods), out%trajectory(periods))
      allocate(log_weights(particles), ancestor_weights(particles), predicted(particles))
      out%particles = 0.0_dp
      out%ancestors = 0
      out%particles(particles, 1) = conditioned_trajectory(1)
      do t = 1, periods
         if (t > 1) then
            predicted = phi*out%particles(:, t - 1)
            do i = 1, particles
               out%ancestors(i, t) = categorical_index(out%weights(:, t - 1), &
                  resampling_uniforms(i, t - 1))
               out%particles(i, t) = predicted(out%ancestors(i, t)) + &
                  sqrt(state_variance)*normal_draws(i, t - 1)
            end do
            out%particles(particles, t) = conditioned_trajectory(t)
            do i = 1, particles
               log_weights(i) = log(max(out%weights(i, t - 1), tiny(1.0_dp))) - &
                  0.5_dp*(conditioned_trajectory(t) - predicted(i))**2/state_variance
            end do
            maximum_log = maxval(log_weights)
            ancestor_weights = exp(log_weights - maximum_log)
            ancestor_weights = ancestor_weights/sum(ancestor_weights)
            out%ancestors(particles, t) = categorical_index(ancestor_weights, ancestor_uniforms(t - 1))
         end if
         do i = 1, particles
            log_weights(i) = -log(observation_scale) - 0.5_dp*out%particles(i, t) - &
               0.5_dp*(observations(t)/observation_scale)**2/exp(out%particles(i, t))
         end do
         maximum_log = maxval(log_weights)
         out%weights(:, t) = exp(log_weights - maximum_log)
         total_weight = sum(out%weights(:, t))
         out%weights(:, t) = out%weights(:, t)/total_weight
         out%log_likelihood = out%log_likelihood + maximum_log + log(total_weight/real(particles, dp)) - &
            0.5_dp*log(2.0_dp*acos(-1.0_dp))
      end do
      selected = categorical_index(out%weights(:, periods), terminal_uniform)
      out%trajectory(periods) = out%particles(selected, periods)
      do t = periods, 2, -1
         selected = out%ancestors(selected, t)
         out%trajectory(t - 1) = out%particles(selected, t - 1)
      end do
   end function sv_particle_filter_draws

   function sv_particle_filter(observations, phi, state_variance, observation_scale, &
      conditioned_trajectory, particle_count) result(out)
      ! Run conditional particle filtering with shared-library random draws.
      real(dp), intent(in) :: observations(:), phi, state_variance, observation_scale
      real(dp), intent(in) :: conditioned_trajectory(:)
      integer, intent(in) :: particle_count
      type(astsa_sv_particle_filter_t) :: out
      real(dp), allocatable :: normal_draws(:, :), resampling_uniforms(:, :), ancestor_uniforms(:)
      real(dp) :: terminal_uniform
      integer :: i, j

      if (particle_count < 2 .or. size(observations) < 2) then
         out%info = 1
         return
      end if
      allocate(normal_draws(particle_count, size(observations) - 1))
      allocate(resampling_uniforms(particle_count, size(observations) - 1))
      allocate(ancestor_uniforms(size(observations) - 1))
      call random_standard_normal_matrix(normal_draws)
      do j = 1, size(observations) - 1
         ancestor_uniforms(j) = random_uniform()
         do i = 1, particle_count
            resampling_uniforms(i, j) = random_uniform()
         end do
      end do
      terminal_uniform = random_uniform()
      out = sv_particle_filter_draws(observations, phi, state_variance, observation_scale, &
         conditioned_trajectory, normal_draws, resampling_uniforms, ancestor_uniforms, terminal_uniform)
   end function sv_particle_filter

   pure integer function categorical_index(weights, uniform_value) result(index)
      ! Select one categorical outcome from normalized weights and a uniform draw.
      real(dp), intent(in) :: weights(:), uniform_value
      real(dp) :: cumulative
      integer :: i

      cumulative = 0.0_dp
      index = size(weights)
      do i = 1, size(weights)
         cumulative = cumulative + weights(i)
         if (uniform_value <= cumulative) then
            index = i
            return
         end if
      end do
   end function categorical_index

   function sv_mcmc(observations, retained_draws, burnin, particle_count, initial, tuning, &
      proposal_covariance) result(out)
      ! Run particle Gibbs with Metropolis updates for astsa's stochastic-volatility model.
      real(dp), intent(in) :: observations(:)
      integer, intent(in) :: retained_draws, burnin, particle_count
      real(dp), intent(in), optional :: initial(3), tuning, proposal_covariance(2, 2)
      type(astsa_sv_mcmc_t) :: out
      type(astsa_sv_particle_filter_t) :: filtered
      real(dp), allocatable :: phi_all(:), state_sd_all(:), scale_all(:), latent_all(:, :)
      real(dp) :: current(2), proposal(2), covariance(2, 2), base_covariance(2, 2)
      real(dp) :: step_scale, log_current, log_proposal, rate, shape, gamma_draw
      integer :: total, iteration, accepted, status, n

      n = size(observations)
      total = retained_draws + burnin
      if (n < 4 .or. retained_draws < 1 .or. burnin < 0 .or. particle_count < 2) then
         out%info = 1
         return
      end if
      allocate(phi_all(total), state_sd_all(total), scale_all(total), latent_all(total, n))
      phi_all(1) = 0.9_dp
      state_sd_all(1) = 0.5_dp
      scale_all(1) = 0.1_dp
      if (present(initial)) then
         phi_all(1) = initial(1)
         state_sd_all(1) = initial(2)
         scale_all(1) = initial(3)
      end if
      if (abs(phi_all(1)) >= 1.0_dp .or. state_sd_all(1) <= 0.0_dp .or. scale_all(1) <= 0.0_dp) then
         out%info = 2
         return
      end if
      step_scale = 0.03_dp
      if (present(tuning)) step_scale = tuning
      base_covariance = reshape([1.0_dp, -0.25_dp, -0.25_dp, 1.0_dp], [2, 2])
      if (present(proposal_covariance)) base_covariance = proposal_covariance
      covariance = step_scale*base_covariance
      latent_all(1, :) = 0.0_dp
      filtered = sv_particle_filter(observations, phi_all(1), state_sd_all(1)**2, &
         scale_all(1), latent_all(1, :), particle_count)
      if (filtered%info /= 0) then
         out%info = 10 + filtered%info
         return
      end if
      latent_all(1, :) = filtered%trajectory - sum(filtered%trajectory)/real(n, dp)
      accepted = 0
      do iteration = 2, total
         current = [phi_all(iteration - 1), state_sd_all(iteration - 1)]
         call random_multivariate_normal(current, covariance, proposal, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         log_current = sv_parameter_log_density(current, latent_all(iteration - 1, :))
         log_proposal = sv_parameter_log_density(proposal, latent_all(iteration - 1, :))
         if (abs(proposal(1)) < 1.0_dp .and. proposal(2) > 0.0_dp .and. &
            log(max(random_uniform(), tiny(1.0_dp))) < log_proposal - log_current) then
            phi_all(iteration) = proposal(1)
            state_sd_all(iteration) = proposal(2)
            accepted = accepted + 1
         else
            phi_all(iteration) = current(1)
            state_sd_all(iteration) = current(2)
         end if
         shape = 0.5_dp*real(n, dp) - 1.0_dp
         rate = 0.5_dp*sum(observations**2/exp(latent_all(iteration - 1, :)))
         gamma_draw = random_gamma(shape, 1.0_dp/max(rate, tiny(1.0_dp)))
         if (gamma_draw <= 0.0_dp) then
            out%info = 3
            return
         end if
         scale_all(iteration) = sqrt(1.0_dp/gamma_draw)
         filtered = sv_particle_filter(observations, phi_all(iteration), state_sd_all(iteration)**2, &
            scale_all(iteration), latent_all(iteration - 1, :), particle_count)
         if (filtered%info /= 0) then
            out%info = 30 + filtered%info
            return
         end if
         latent_all(iteration, :) = filtered%trajectory - sum(filtered%trajectory)/real(n, dp)
      end do
      out%phi = phi_all(burnin + 1:)
      out%state_sd = state_sd_all(burnin + 1:)
      out%observation_scale = scale_all(burnin + 1:)
      out%latent = latent_all(burnin + 1:, :)
      out%acceptance_rate = real(accepted, dp)/real(max(1, total - 1), dp)
      allocate(out%effective_size(3))
      out%effective_size = [effective_sample_size(out%phi), effective_sample_size(out%state_sd), &
         effective_sample_size(out%observation_scale)]
   end function sv_mcmc

   pure real(dp) function sv_parameter_log_density(parameters, latent) result(value)
      ! Return the latent-state likelihood plus astsa-style correlated Gaussian priors.
      real(dp), intent(in) :: parameters(2), latent(:)
      real(dp), parameter :: mean_phi = 0.9_dp, mean_sd = 0.5_dp
      real(dp), parameter :: sd_phi = 0.075_dp, sd_sd = 0.3_dp, correlation = -0.25_dp
      real(dp) :: phi, state_sd, prior, transition_sum

      phi = parameters(1)
      state_sd = parameters(2)
      if (abs(phi) >= 1.0_dp .or. state_sd <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      prior = ((phi - mean_phi)/sd_phi)**2 + ((state_sd - mean_sd)/sd_sd)**2 - &
         2.0_dp*correlation*(phi - mean_phi)*(state_sd - mean_sd)/(sd_phi*sd_sd)
      prior = -0.5_dp*prior/(1.0_dp - correlation**2)
      transition_sum = sum((latent(2:) - phi*latent(:size(latent) - 1))**2)
      value = 0.5_dp*log(1.0_dp - phi**2) - real(size(latent), dp)*log(state_sd) + prior - &
         0.5_dp*((1.0_dp - phi**2)*latent(1)**2 + transition_sum)/state_sd**2
   end function sv_parameter_log_density

   pure real(dp) function effective_sample_size(draws) result(value)
      ! Estimate effective sample size by the initial positive autocorrelation sequence.
      real(dp), intent(in) :: draws(:)
      real(dp) :: centered(size(draws)), variance, correlation_sum, correlation
      integer :: lag, n

      n = size(draws)
      if (n < 2) then
         value = real(n, dp)
         return
      end if
      centered = draws - sum(draws)/real(n, dp)
      variance = sum(centered**2)
      if (variance <= tiny(1.0_dp)) then
         value = 1.0_dp
         return
      end if
      correlation_sum = 0.0_dp
      do lag = 1, n - 1
         correlation = dot_product(centered(:n - lag), centered(lag + 1:))/variance
         if (correlation <= 0.0_dp) exit
         correlation_sum = correlation_sum + correlation
      end do
      value = min(real(n, dp), real(n, dp)/(1.0_dp + 2.0_dp*correlation_sum))
   end function effective_sample_size

   pure function ssm_fit(observations, measurement, phi, alpha, state_sd, observation_sd, &
      fix_phi, max_iterations, tolerance) result(out)
      ! Estimate astsa's scalar linear Gaussian state-space model.
      real(dp), intent(in) :: observations(:), measurement, phi, alpha, state_sd, observation_sd
      logical, intent(in), optional :: fix_phi
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(astsa_ssm_result_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: initial(:), hessian(:, :), inverse(:, :), jacobian(:, :)
      integer :: count, limit, status, i
      real(dp) :: gradient_tolerance

      out%fixed_phi = .false.
      if (present(fix_phi)) out%fixed_phi = fix_phi
      if (size(observations) < 20 .or. state_sd <= 0.0_dp .or. observation_sd <= 0.0_dp) then
         out%info = 1
         return
      end if
      count = merge(3, 4, out%fixed_phi)
      allocate(initial(count))
      if (out%fixed_phi) then
         initial = [alpha, log(state_sd), log(observation_sd)]
      else
         initial = [atanh(max(-0.999999_dp, min(0.999999_dp, phi))), alpha, &
            log(state_sd), log(observation_sd)]
      end if
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      gradient_tolerance = 1.0e-6_dp
      if (present(tolerance)) gradient_tolerance = tolerance
      optimization = bfgs_minimize_fd(objective, initial, limit, gradient_tolerance)
      out = evaluate(reported_parameters(optimization%parameters))
      out%iterations = optimization%iterations
      out%converged = optimization%converged
      out%fixed_phi = present_and_true(fix_phi)
      if (optimization%info /= 0) out%info = optimization%info
      hessian = finite_difference_hessian(objective, optimization%parameters)
      call invert_matrix(hessian, inverse, status)
      allocate(out%covariance(count, count), out%standard_errors(count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (status == 0) then
         jacobian = parameter_jacobian(optimization%parameters)
         out%covariance = matmul(matmul(jacobian, inverse), transpose(jacobian))
         do i = 1, count
            if (out%covariance(i, i) > 0.0_dp) out%standard_errors(i) = sqrt(out%covariance(i, i))
         end do
      else if (out%info == 0) then
         out%info = 20 + status
      end if

   contains

      pure function objective(parameters) result(value)
         ! Return the scalar state-space negative log likelihood.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value
         type(astsa_ssm_result_t) :: evaluated

         evaluated = evaluate(reported_parameters(parameters))
         if (evaluated%info == 0) then
            value = evaluated%negative_log_likelihood
         else
            value = 1.0e30_dp + dot_product(parameters, parameters)
         end if
      end function objective

      pure function reported_parameters(parameters) result(reported)
         ! Map optimizer coordinates to state-space parameters.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: reported(size(parameters))

         reported = parameters
         if (out%fixed_phi) then
            reported(2:3) = exp(parameters(2:3))
         else
            reported(1) = tanh(parameters(1))
            reported(3:4) = exp(parameters(3:4))
         end if
      end function reported_parameters

      pure function evaluate(parameters) result(evaluated)
         ! Filter and smooth the scalar model for reported parameters.
         real(dp), intent(in) :: parameters(:)
         type(astsa_ssm_result_t) :: evaluated
         real(dp) :: transition, state_intercept, process_sd, noise_sd
         real(dp) :: initial_mean, initial_variance, innovation, innovation_variance, gain, smoother_gain
         integer :: n, t

         if (out%fixed_phi) then
            transition = phi
            state_intercept = parameters(1)
            process_sd = parameters(2)
            noise_sd = parameters(3)
         else
            transition = parameters(1)
            state_intercept = parameters(2)
            process_sd = parameters(3)
            noise_sd = parameters(4)
         end if
         evaluated%coefficients = parameters
         n = size(observations)
         allocate(evaluated%predicted_state(n), evaluated%predicted_variance(n))
         allocate(evaluated%filtered_state(n), evaluated%filtered_variance(n))
         allocate(evaluated%smoothed_state(n), evaluated%smoothed_variance(n))
         initial_mean = sum(observations(:5))/5.0_dp
         initial_variance = sum((observations(:5) - initial_mean)**2)/4.0_dp
         do t = 1, n
            if (t == 1) then
               evaluated%predicted_state(t) = state_intercept + transition*initial_mean
               evaluated%predicted_variance(t) = transition**2*initial_variance + process_sd**2
            else
               evaluated%predicted_state(t) = state_intercept + transition*evaluated%filtered_state(t - 1)
               evaluated%predicted_variance(t) = transition**2*evaluated%filtered_variance(t - 1) + process_sd**2
            end if
            innovation_variance = measurement**2*evaluated%predicted_variance(t) + noise_sd**2
            if (innovation_variance <= tiny(1.0_dp)) then
               evaluated%info = 1
               return
            end if
            innovation = observations(t) - measurement*evaluated%predicted_state(t)
            gain = measurement*evaluated%predicted_variance(t)/innovation_variance
            evaluated%filtered_state(t) = evaluated%predicted_state(t) + gain*innovation
            evaluated%filtered_variance(t) = evaluated%predicted_variance(t)*(1.0_dp - gain*measurement)
            evaluated%negative_log_likelihood = evaluated%negative_log_likelihood + &
               0.5_dp*(log(innovation_variance) + innovation**2/innovation_variance)
         end do
         evaluated%smoothed_state(n) = evaluated%filtered_state(n)
         evaluated%smoothed_variance(n) = evaluated%filtered_variance(n)
         do t = n - 1, 1, -1
            smoother_gain = evaluated%filtered_variance(t)*transition/evaluated%predicted_variance(t + 1)
            evaluated%smoothed_state(t) = evaluated%filtered_state(t) + smoother_gain* &
               (evaluated%smoothed_state(t + 1) - evaluated%predicted_state(t + 1))
            evaluated%smoothed_variance(t) = evaluated%filtered_variance(t) + smoother_gain**2* &
               (evaluated%smoothed_variance(t + 1) - evaluated%predicted_variance(t + 1))
         end do
      end function evaluate

      pure function parameter_jacobian(parameters) result(jacobian_value)
         ! Compute the optimizer-to-reported state-space Jacobian.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: jacobian_value(size(parameters), size(parameters))
         real(dp) :: plus(size(parameters)), minus(size(parameters)), step
         integer :: column

         do column = 1, size(parameters)
            step = epsilon(1.0_dp)**0.25_dp*max(1.0_dp, abs(parameters(column)))
            plus = parameters
            minus = parameters
            plus(column) = plus(column) + step
            minus(column) = minus(column) - step
            jacobian_value(:, column) = (reported_parameters(plus) - reported_parameters(minus))/(2.0_dp*step)
         end do
      end function parameter_jacobian
   end function ssm_fit

   pure logical function present_and_true(value) result(answer)
      ! Return true only when an optional logical is present and true.
      logical, intent(in), optional :: value

      answer = .false.
      if (present(value)) answer = value
   end function present_and_true

   pure subroutine invert_complex_matrix(matrix, inverse, info)
      ! Invert a small complex matrix by pivoted Gauss-Jordan elimination.
      complex(dp), intent(in) :: matrix(:, :)
      complex(dp), allocatable, intent(out) :: inverse(:, :)
      integer, intent(out) :: info
      complex(dp), allocatable :: work(:, :), row(:)
      integer :: n, i, j, pivot

      n = size(matrix, 1)
      allocate(work(n, 2*n), row(2*n), inverse(n, n))
      work = cmplx(0.0_dp, 0.0_dp, dp)
      work(:, :n) = matrix
      do i = 1, n
         work(i, n + i) = cmplx(1.0_dp, 0.0_dp, dp)
      end do
      info = 0
      do i = 1, n
         pivot = i - 1 + maxloc(abs(work(i:, i)), dim=1)
         if (abs(work(pivot, i)) <= 100.0_dp*epsilon(1.0_dp)) then
            info = i
            return
         end if
         if (pivot /= i) then
            row = work(i, :)
            work(i, :) = work(pivot, :)
            work(pivot, :) = row
         end if
         work(i, :) = work(i, :)/work(i, i)
         do j = 1, n
            if (j /= i) work(j, :) = work(j, :) - work(j, i)*work(i, :)
         end do
      end do
      inverse = work(:, n + 1:)
   end subroutine invert_complex_matrix

   pure function sorted_values(values) result(sorted)
      ! Return values sorted in ascending order.
      real(dp), intent(in) :: values(:)
      real(dp), allocatable :: sorted(:)
      real(dp) :: temporary
      integer :: i, j

      sorted = values
      do i = 2, size(sorted)
         temporary = sorted(i)
         j = i - 1
         do while (j >= 1)
            if (sorted(j) <= temporary) exit
            sorted(j + 1) = sorted(j)
            j = j - 1
         end do
         sorted(j + 1) = temporary
      end do
   end function sorted_values

   pure real(dp) function regularized_gamma_q(shape, value) result(probability)
      ! Return the upper regularized incomplete gamma function.
      real(dp), intent(in) :: shape, value
      real(dp) :: term, sum_value, b, c, d, delta, log_factor
      integer :: i

      if (shape <= 0.0_dp .or. value < 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      if (value == 0.0_dp) then
         probability = 1.0_dp
         return
      end if
      log_factor = shape*log(value) - value - log_gamma(shape)
      if (value < shape + 1.0_dp) then
         term = 1.0_dp/shape
         sum_value = term
         do i = 1, 10000
            term = term*value/(shape + real(i, dp))
            sum_value = sum_value + term
            if (abs(term) <= epsilon(1.0_dp)*abs(sum_value)) exit
         end do
         probability = 1.0_dp - exp(log_factor)*sum_value
      else
         b = value + 1.0_dp - shape
         c = 1.0_dp/tiny(1.0_dp)
         d = 1.0_dp/b
         probability = d
         do i = 1, 10000
            term = -real(i, dp)*(real(i, dp) - shape)
            b = b + 2.0_dp
            d = term*d + b
            if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
            c = b + term/c
            if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
            d = 1.0_dp/d
            delta = d*c
            probability = probability*delta
            if (abs(delta - 1.0_dp) <= epsilon(1.0_dp)) exit
         end do
         probability = exp(log_factor)*probability
      end if
   end function regularized_gamma_q

   pure real(dp) function f_distribution_quantile(probability, numerator_df, denominator_df) result(value)
      ! Return an F-distribution quantile by bisection of its beta-form CDF.
      real(dp), intent(in) :: probability
      integer, intent(in) :: numerator_df, denominator_df
      real(dp) :: lower, upper, middle, cdf
      integer :: iteration

      lower = 0.0_dp
      upper = 1.0_dp
      do
         cdf = regularized_beta(real(numerator_df, dp)*upper/ &
            (real(numerator_df, dp)*upper + real(denominator_df, dp)), &
            0.5_dp*real(numerator_df, dp), 0.5_dp*real(denominator_df, dp))
         if (cdf >= probability .or. upper >= 1.0e12_dp) exit
         upper = 2.0_dp*upper
      end do
      do iteration = 1, 200
         middle = 0.5_dp*(lower + upper)
         cdf = regularized_beta(real(numerator_df, dp)*middle/ &
            (real(numerator_df, dp)*middle + real(denominator_df, dp)), &
            0.5_dp*real(numerator_df, dp), 0.5_dp*real(denominator_df, dp))
         if (cdf < probability) then
            lower = middle
         else
            upper = middle
         end if
         if (upper - lower <= 1.0e-10_dp*max(1.0_dp, middle)) exit
      end do
      value = 0.5_dp*(lower + upper)
   end function f_distribution_quantile

   pure real(dp) function regularized_beta(value, first_shape, second_shape) result(probability)
      ! Return the regularized incomplete beta function.
      real(dp), intent(in) :: value, first_shape, second_shape
      real(dp) :: factor

      if (value <= 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      if (value >= 1.0_dp) then
         probability = 1.0_dp
         return
      end if
      factor = exp(log_gamma(first_shape + second_shape) - log_gamma(first_shape) - &
         log_gamma(second_shape) + first_shape*log(value) + second_shape*log(1.0_dp - value))
      if (value < (first_shape + 1.0_dp)/(first_shape + second_shape + 2.0_dp)) then
         probability = factor*beta_continued_fraction(value, first_shape, second_shape)/first_shape
      else
         probability = 1.0_dp - factor*beta_continued_fraction(1.0_dp - value, &
            second_shape, first_shape)/second_shape
      end if
   end function regularized_beta

   pure real(dp) function beta_continued_fraction(value, first_shape, second_shape) result(fraction)
      ! Evaluate the incomplete-beta continued fraction.
      real(dp), intent(in) :: value, first_shape, second_shape
      real(dp) :: qab, qap, qam, c, d, h, aa, delta
      integer :: iteration, twice

      qab = first_shape + second_shape
      qap = first_shape + 1.0_dp
      qam = first_shape - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*value/qap
      if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
      d = 1.0_dp/d
      h = d
      do iteration = 1, 10000
         twice = 2*iteration
         aa = real(iteration, dp)*(second_shape - real(iteration, dp))*value/ &
            ((qam + real(twice, dp))*(first_shape + real(twice, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
         c = 1.0_dp + aa/c
         if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
         d = 1.0_dp/d
         h = h*d*c
         aa = -(first_shape + real(iteration, dp))*(qab + real(iteration, dp))*value/ &
            ((first_shape + real(twice, dp))*(qap + real(twice, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
         c = 1.0_dp + aa/c
         if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
         d = 1.0_dp/d
         delta = d*c
         h = h*delta
         if (abs(delta - 1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) exit
      end do
      fraction = h
   end function beta_continued_fraction

   pure subroutine build_sarima_polynomials(ar, ma, sar, sma, season, ar_polynomial, ma_polynomial)
      ! Combine ordinary and seasonal AR and MA polynomials.
      real(dp), intent(in) :: ar(:), ma(:), sar(:), sma(:)
      integer, intent(in) :: season
      real(dp), allocatable, intent(out) :: ar_polynomial(:), ma_polynomial(:)
      real(dp), allocatable :: base(:), seasonal(:)
      integer :: i

      allocate(base(size(ar) + 1), seasonal(size(sar)*season + 1))
      base = 0.0_dp
      seasonal = 0.0_dp
      base(1) = 1.0_dp
      seasonal(1) = 1.0_dp
      if (size(ar) > 0) base(2:) = -ar
      do i = 1, size(sar)
         seasonal(i*season + 1) = -sar(i)
      end do
      ar_polynomial = poly_mul(base, seasonal)
      deallocate(base, seasonal)
      allocate(base(size(ma) + 1), seasonal(size(sma)*season + 1))
      base = 0.0_dp
      seasonal = 0.0_dp
      base(1) = 1.0_dp
      seasonal(1) = 1.0_dp
      if (size(ma) > 0) base(2:) = ma
      do i = 1, size(sma)
         seasonal(i*season + 1) = sma(i)
      end do
      ma_polynomial = poly_mul(base, seasonal)
   end subroutine build_sarima_polynomials

   function sarima_sim(ar, d, ma, sar, seasonal_difference, sma, season, n, burnin) result(out)
      ! Simulate a SARIMA series using standard-normal innovations.
      real(dp), intent(in), optional :: ar(:), ma(:), sar(:), sma(:)
      integer, intent(in), optional :: d, seasonal_difference, season, burnin
      integer, intent(in) :: n
      type(astsa_simulation_t) :: out
      real(dp), allocatable :: draws(:, :)
      integer :: total, warmup

      warmup = sarima_default_burnin(ar, d, ma, sar, seasonal_difference, sma, season)
      if (present(burnin)) warmup = burnin
      if (n < 1 .or. warmup < 0) then
         out%info = 1
         return
      end if
      total = n + warmup
      allocate(draws(total, 1))
      call random_standard_normal_matrix(draws)
      out = sarima_sim_from_innovations(ar, d, ma, sar, seasonal_difference, sma, season, n, warmup, draws(:, 1))
   end function sarima_sim

   pure function sarima_sim_from_innovations(ar, d, ma, sar, seasonal_difference, sma, season, &
      n, burnin, innovations) result(out)
      ! Simulate a SARIMA series from caller-supplied innovations.
      real(dp), intent(in), optional :: ar(:), ma(:), sar(:), sma(:)
      integer, intent(in), optional :: d, seasonal_difference, season, burnin
      integer, intent(in) :: n
      real(dp), intent(in) :: innovations(:)
      type(astsa_simulation_t) :: out
      real(dp), allocatable :: ar_base(:), ma_base(:), ar_seasonal(:), ma_seasonal(:)
      real(dp), allocatable :: ar_polynomial(:), ma_polynomial(:), ar_coefficient(:), ma_coefficient(:), work(:)
      type(astsa_root_check_t) :: root_check
      integer :: ordinary_d, seasonal_d, period, warmup, total, p, q, ps, qs, i, j, t

      ordinary_d = 0
      if (present(d)) ordinary_d = d
      seasonal_d = 0
      if (present(seasonal_difference)) seasonal_d = seasonal_difference
      period = 1
      if (present(season)) period = season
      warmup = sarima_default_burnin(ar, d, ma, sar, seasonal_difference, sma, season)
      if (present(burnin)) warmup = burnin
      p = 0
      if (present(ar)) p = size(ar)
      q = 0
      if (present(ma)) q = size(ma)
      ps = 0
      if (present(sar)) ps = size(sar)
      qs = 0
      if (present(sma)) qs = size(sma)
      if (n < 1 .or. warmup < 0 .or. ordinary_d < 0 .or. seasonal_d < 0 .or. period < 1) then
         out%info = 1
         return
      end if
      if ((ps > 0 .or. qs > 0 .or. seasonal_d > 0) .and. .not. present(season)) then
         out%info = 2
         return
      end if
      if (present(season)) then
         if (period <= p .or. period <= q) then
            out%info = 4
            return
         end if
      end if
      total = n + warmup
      if (size(innovations) < total) then
         out%info = 3
         return
      end if

      allocate(ar_base(p + 1), ma_base(q + 1), ar_seasonal(ps*period + 1), ma_seasonal(qs*period + 1))
      ar_base = 0.0_dp
      ma_base = 0.0_dp
      ar_seasonal = 0.0_dp
      ma_seasonal = 0.0_dp
      ar_base(1) = 1.0_dp
      ma_base(1) = 1.0_dp
      ar_seasonal(1) = 1.0_dp
      ma_seasonal(1) = 1.0_dp
      if (p > 0) ar_base(2:) = -ar
      if (q > 0) ma_base(2:) = ma
      if (ps > 0) then
         do i = 1, ps
            ar_seasonal(i*period + 1) = -sar(i)
         end do
      end if
      if (qs > 0) then
         do i = 1, qs
            ma_seasonal(i*period + 1) = sma(i)
         end do
      end if
      ar_polynomial = poly_mul(ar_base, ar_seasonal)
      ma_polynomial = poly_mul(ma_base, ma_seasonal)
      root_check = polynomial_root_check(ar_polynomial)
      if (.not. root_check%converged) then
         out%info = 5
         return
      end if
      if (.not. root_check%outside_unit_circle) then
         out%info = 6
         return
      end if
      root_check = polynomial_root_check(ma_polynomial)
      if (.not. root_check%converged) then
         out%info = 7
         return
      end if
      if (.not. root_check%outside_unit_circle) then
         out%info = 8
         return
      end if
      ar_coefficient = -ar_polynomial(2:)
      ma_coefficient = ma_polynomial(2:)
      allocate(work(total))
      work = innovations(:total)
      do t = 1, total
         do j = 1, min(size(ma_coefficient), t - 1)
            work(t) = work(t) + ma_coefficient(j)*innovations(t - j)
         end do
         do j = 1, min(size(ar_coefficient), t - 1)
            work(t) = work(t) + ar_coefficient(j)*work(t - j)
         end do
      end do
      do i = 1, ordinary_d
         do t = 2, total
            work(t) = work(t) + work(t - 1)
         end do
      end do
      do i = 1, seasonal_d
         do t = period + 1, total
            work(t) = work(t) + work(t - period)
         end do
      end do
      out%series = work(warmup + 1:total)
      out%frequency = period
   end function sarima_sim_from_innovations

   pure integer function sarima_default_burnin(ar, d, ma, sar, seasonal_difference, sma, season) result(value)
      ! Return astsa's default SARIMA simulation burn-in length.
      real(dp), intent(in), optional :: ar(:), ma(:), sar(:), sma(:)
      integer, intent(in), optional :: d, seasonal_difference, season
      integer :: ordinary_d, seasonal_d, period, p, q, ps, qs

      ordinary_d = 0
      if (present(d)) ordinary_d = d
      seasonal_d = 0
      if (present(seasonal_difference)) seasonal_d = seasonal_difference
      period = 1
      if (present(season)) period = season
      p = 0
      if (present(ar)) p = size(ar)
      q = 0
      if (present(ma)) q = size(ma)
      ps = 0
      if (present(sar)) ps = size(sar)
      qs = 0
      if (present(sma)) qs = size(sma)
      value = 50 + ordinary_d + p + q
      if (present(season)) value = 50 + (seasonal_d + ps + qs)*period + ordinary_d + p + q
   end function sarima_default_burnin

   pure function polynomial_root_check(coefficients) result(out)
      ! Test whether every polynomial root lies outside the unit circle.
      real(dp), intent(in) :: coefficients(:)
      type(astsa_root_check_t) :: out
      type(astsa_polynomial_roots_t) :: roots

      roots = polynomial_roots(coefficients)
      out%converged = roots%converged
      if (roots%converged .and. size(roots%roots) > 0) then
         out%outside_unit_circle = minval(abs(roots%roots)) > 1.0_dp
      end if
   end function polynomial_root_check

   pure function polynomial_roots(coefficients) result(out)
      ! Compute complex roots of a real polynomial in ascending coefficient order.
      real(dp), intent(in) :: coefficients(:)
      type(astsa_polynomial_roots_t) :: out
      complex(dp), allocatable :: previous(:)
      complex(dp) :: numerator, denominator
      real(dp) :: radius, angle
      integer :: degree, iteration, i, j

      degree = size(coefficients) - 1
      do while (degree > 0 .and. abs(coefficients(degree + 1)) <= tiny(1.0_dp))
         degree = degree - 1
      end do
      allocate(out%roots(degree))
      if (degree < 1) return
      allocate(previous(degree))
      radius = 1.0_dp + maxval(abs(coefficients(1:degree)/coefficients(degree + 1)))
      do i = 1, degree
         angle = 2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/real(degree, dp)
         out%roots(i) = cmplx(radius*cos(angle), radius*sin(angle), dp)
      end do
      out%converged = .false.
      do iteration = 1, 2000
         previous = out%roots
         do i = 1, degree
            numerator = cmplx(coefficients(degree + 1), 0.0_dp, dp)
            do j = degree, 1, -1
               numerator = numerator*previous(i) + coefficients(j)
            end do
            denominator = cmplx(coefficients(degree + 1), 0.0_dp, dp)
            do j = 1, degree
               if (j /= i) denominator = denominator*(previous(i) - previous(j))
            end do
            if (abs(denominator) <= tiny(1.0_dp)) return
            out%roots(i) = previous(i) - numerator/denominator
         end do
         if (maxval(abs(out%roots - previous)) <= 100.0_dp*epsilon(1.0_dp)) then
            out%converged = .true.
            exit
         end if
      end do
   end function polynomial_roots

   pure logical function inverse_roots_overlap(first, second, tolerance) result(overlap)
      ! Test whether two root sets have inverse roots within a tolerance.
      complex(dp), intent(in) :: first(:), second(:)
      real(dp), intent(in) :: tolerance
      integer :: i, j

      overlap = .false.
      do i = 1, size(first)
         do j = 1, size(second)
            if (abs(1.0_dp/first(i) - 1.0_dp/second(j)) <= tolerance) then
               overlap = .true.
               return
            end if
         end do
      end do
   end function inverse_roots_overlap

   pure subroutine local_inverse_logdet(a, inverse, logdet, info)
      ! Call the shared inverse/log-determinant routine with library tolerance.
      use time_series_linalg_mod, only: inverse_logdet
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: inverse(:, :), logdet
      integer, intent(out) :: info
      call inverse_logdet(a, inverse, logdet, info, 100.0_dp*epsilon(1.0_dp))
   end subroutine local_inverse_logdet
end module astsa_mod
