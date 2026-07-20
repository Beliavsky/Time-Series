! SPDX-License-Identifier: GPL-2.0-only
! SPDX-FileComment: Algorithms translated from the R bayesforecast package.
! Native posterior kernels translated from the GPL-2 bayesforecast package.
module bayesforecast_mod
   use kind_mod, only: dp
   use stats_mod, only: sort
   use special_functions_mod, only: normal_log_density
   use random_mod, only: random_uniform, random_standard_normal, random_gamma
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: bf_prior_normal = 1
   integer, parameter, public :: bf_prior_beta = 2
   integer, parameter, public :: bf_prior_uniform = 3
   integer, parameter, public :: bf_prior_student = 4
   integer, parameter, public :: bf_prior_cauchy = 5
   integer, parameter, public :: bf_prior_inverse_gamma = 6
   integer, parameter, public :: bf_prior_inverse_chisq = 7
   integer, parameter, public :: bf_prior_jeffreys = 8
   integer, parameter, public :: bf_prior_gamma = 9
   integer, parameter, public :: bf_prior_exponential = 10
   integer, parameter, public :: bf_prior_chisq = 11
   integer, parameter, public :: bf_prior_laplace = 12

   type, public :: bf_filter_t
      ! Conditional means, innovations, scales, and pointwise likelihoods.
      real(dp), allocatable :: mean(:), residual(:), scale(:), log_likelihood(:)
      real(dp) :: total_log_likelihood = 0.0_dp
      integer :: info = 0
   end type

   type, public :: bf_ets_filter_t
      ! ETS states and its observation likelihood.
      type(bf_filter_t) :: observation
      real(dp), allocatable :: level(:), trend(:), seasonal(:)
      integer :: info = 0
   end type

   type, public :: bf_sv_filter_t
      ! Stochastic-volatility observation and latent-state likelihoods.
      type(bf_filter_t) :: observation
      real(dp), allocatable :: state_log_likelihood(:)
      real(dp) :: total_log_likelihood = 0.0_dp
      integer :: info = 0
   end type

   type, public :: bf_mcmc_t
      ! Random-walk Metropolis posterior draws and acceptance diagnostics.
      real(dp), allocatable :: draws(:, :), log_density(:)
      real(dp) :: acceptance_rate = 0.0_dp
      integer :: burnin = 0
      integer :: thin = 1
      integer :: info = 0
   end type

   type, public :: bf_information_criteria_t
      ! WAIC components computed from pointwise posterior log likelihoods.
      real(dp) :: lppd = 0.0_dp
      real(dp) :: effective_parameters = 0.0_dp
      real(dp) :: waic = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: aicc = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: info = 0
   end type

   type, public :: bf_predictive_t
      ! Posterior predictive draws, means, and equal-tail intervals.
      real(dp), allocatable :: draws(:, :), mean(:), interval(:, :)
      integer :: horizon = 0
      integer :: info = 0
   end type

   type, public :: bf_difference_t
      ! Differenced series and boundary values required for inversion.
      real(dp), allocatable :: values(:), ordinary_initial(:), seasonal_initial(:, :)
      integer :: ordinary_order = 0
      integer :: seasonal_order = 0
      integer :: period = 1
      integer :: info = 0
   end type

   abstract interface
      pure function bf_log_posterior_interface(parameters) result(value)
         !! Evaluate a log posterior at a parameter vector.
         import dp
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: value
      end function bf_log_posterior_interface
   end interface

   public :: bf_log_prior, bf_sarima_filter, bf_ets_filter
   public :: bf_garch_filter, bf_stochastic_volatility_filter
   public :: bf_metropolis_sample, bf_posterior_interval, bf_predictive_error, bf_waic
   public :: bf_log_prior_sum, bf_sarima_log_posterior
   public :: bf_log_posterior, bf_information_criteria
   public :: bf_bayes_factor
   public :: bf_difference, bf_inverse_difference, bf_sarima_predict
   public :: bf_ets_predict, bf_garch_predict, bf_stochastic_volatility_predict

contains

   pure function bf_log_prior(value, specification) result(log_density)
      !! Evaluate one bayesforecast prior encoded as location, scale, df, and code.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: specification(4) !! Specification.
      real(dp) :: log_density
      real(dp) :: location, scale, degrees, pi
      integer :: code

      location = specification(1)
      scale = specification(2)
      degrees = specification(3)
      code = nint(specification(4))
      pi = acos(-1.0_dp)
      log_density = -huge(1.0_dp)
      select case (code)
      case (bf_prior_normal)
         if (scale > 0.0_dp) log_density = -0.5_dp*log(2.0_dp*pi*scale**2) - &
            0.5_dp*((value - location)/scale)**2
      case (bf_prior_beta)
         if (value > 0.0_dp .and. value < 1.0_dp .and. location > 0.0_dp .and. scale > 0.0_dp) &
            log_density = (location - 1.0_dp)*log(value) + (scale - 1.0_dp)*log(1.0_dp - value) + &
            log_gamma(location + scale) - log_gamma(location) - log_gamma(scale)
      case (bf_prior_uniform)
         if (scale > location .and. value >= location .and. value <= scale) &
            log_density = -log(scale - location)
      case (bf_prior_student)
         if (scale > 0.0_dp .and. degrees > 0.0_dp) log_density = &
            log_gamma(0.5_dp*(degrees + 1.0_dp)) - log_gamma(0.5_dp*degrees) - &
            0.5_dp*log(degrees*pi) - log(scale) - 0.5_dp*(degrees + 1.0_dp)* &
            log(1.0_dp + ((value - location)/scale)**2/degrees)
      case (bf_prior_cauchy)
         if (scale > 0.0_dp) log_density = -log(pi*scale) - &
            log(1.0_dp + ((value - location)/scale)**2)
      case (bf_prior_inverse_gamma)
         if (value > 0.0_dp .and. location > 0.0_dp .and. scale > 0.0_dp) log_density = &
            location*log(scale) - log_gamma(location) - (location + 1.0_dp)*log(value) - scale/value
      case (bf_prior_inverse_chisq)
         if (value > 0.0_dp .and. degrees > 0.0_dp) log_density = &
            0.5_dp*degrees*log(0.5_dp*degrees) - log_gamma(0.5_dp*degrees) - &
            (0.5_dp*degrees + 1.0_dp)*log(value) - 0.5_dp*degrees/value
      case (bf_prior_jeffreys)
         if (value > 0.0_dp) log_density = -log(value)
      case (bf_prior_gamma)
         if (value > 0.0_dp .and. location > 0.0_dp .and. scale > 0.0_dp) log_density = &
            location*log(scale) - log_gamma(location) + (location - 1.0_dp)*log(value) - scale*value
      case (bf_prior_exponential)
         if (value >= 0.0_dp .and. scale > 0.0_dp) log_density = log(scale) - scale*value
      case (bf_prior_chisq)
         if (value > 0.0_dp .and. degrees > 0.0_dp) log_density = &
            (0.5_dp*degrees - 1.0_dp)*log(value) - 0.5_dp*value - &
            0.5_dp*degrees*log(2.0_dp) - log_gamma(0.5_dp*degrees)
      case (bf_prior_laplace)
         if (scale > 0.0_dp) log_density = -log(2.0_dp*scale) - abs(value - location)/scale
      end select
   end function bf_log_prior

   pure function bf_sarima_filter(series, sigma, intercept, ar, ma, seasonal_ar, &
      seasonal_ma, period, regressors, regression) result(out)
      !! Evaluate the package's additive seasonal ARIMA/ARIMAX recursion.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: sigma !! Scale parameter or standard deviation.
      real(dp), intent(in) :: intercept !! Model intercept.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: seasonal_ar(:) !! Seasonal autoregressive.
      real(dp), intent(in) :: seasonal_ma(:) !! Seasonal moving-average.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: regression(:) !! Regression.
      type(bf_filter_t) :: out
      integer :: n, t, lag

      n = size(series)
      if (n < 1 .or. sigma <= 0.0_dp .or. period < 1 .or. &
         (present(regressors) .neqv. present(regression))) then
         out%info = 1
         return
      end if
      if (present(regressors)) then
         if (size(regressors, 1) /= n .or. size(regressors, 2) /= size(regression)) then
            out%info = 1
            return
         end if
      end if
      allocate(out%mean(n), out%residual(n), out%scale(n), out%log_likelihood(n))
      out%mean = intercept
      if (present(regressors)) out%mean = out%mean + matmul(regressors, regression)
      out%residual = 0.0_dp
      out%scale = sigma
      do t = 1, n
         do lag = 1, min(size(ar), t - 1)
            out%mean(t) = out%mean(t) + ar(lag)*series(t - lag)
         end do
         do lag = 1, min(size(ma), t - 1)
            out%mean(t) = out%mean(t) + ma(lag)*out%residual(t - lag)
         end do
         do lag = 1, size(seasonal_ar)
            if (t > lag*period) out%mean(t) = out%mean(t) + &
               seasonal_ar(lag)*series(t - lag*period)
         end do
         do lag = 1, size(seasonal_ma)
            if (t > lag*period) out%mean(t) = out%mean(t) + &
               seasonal_ma(lag)*out%residual(t - lag*period)
         end do
         out%residual(t) = series(t) - out%mean(t)
         out%log_likelihood(t) = normal_log_density(out%residual(t), sigma)
      end do
      out%total_log_likelihood = sum(out%log_likelihood)
   end function bf_sarima_filter

   pure function bf_ets_filter(series, sigma, level_smoothing, initial_level, &
      trend_smoothing, initial_trend, damping, seasonal_smoothing, initial_seasonal, &
      regressors, regression, degrees_of_freedom) result(out)
      !! Evaluate local-level, Holt, damped, and seasonal ETS state recursions.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: sigma !! Scale parameter or standard deviation.
      real(dp), intent(in) :: level_smoothing !! Level smoothing.
      real(dp), intent(in) :: initial_level !! Initial level.
      real(dp), intent(in), optional :: trend_smoothing !! Trend smoothing.
      real(dp), intent(in), optional :: initial_trend !! Initial trend.
      real(dp), intent(in), optional :: damping !! Damping.
      real(dp), intent(in), optional :: seasonal_smoothing !! Seasonal smoothing.
      real(dp), intent(in), optional :: initial_seasonal(:) !! Initial seasonal.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: regression(:) !! Regression.
      real(dp), intent(in), optional :: degrees_of_freedom !! Degrees of freedom.
      type(bf_ets_filter_t) :: out
      logical :: has_trend, has_seasonal
      real(dp) :: trend_alpha, damp, seasonal_alpha
      integer :: n, period, t

      n = size(series)
      has_trend = present(trend_smoothing) .and. present(initial_trend)
      has_seasonal = present(seasonal_smoothing) .and. present(initial_seasonal)
      period = 1
      if (has_seasonal) period = size(initial_seasonal)
      if (n < 1 .or. sigma <= 0.0_dp .or. level_smoothing < 0.0_dp .or. &
         level_smoothing > 1.0_dp .or. (present(regressors) .neqv. present(regression))) then
         out%info = 1
         return
      end if
      trend_alpha = 0.0_dp
      damp = 1.0_dp
      seasonal_alpha = 0.0_dp
      if (has_trend) trend_alpha = trend_smoothing
      if (present(damping)) damp = damping
      if (has_seasonal) seasonal_alpha = seasonal_smoothing
      allocate(out%observation%mean(n), out%observation%residual(n), out%observation%scale(n))
      allocate(out%observation%log_likelihood(n), out%level(n), out%trend(n), out%seasonal(n))
      out%observation%mean = 0.0_dp
      if (present(regressors)) out%observation%mean = matmul(regressors, regression)
      out%level = initial_level
      out%trend = 0.0_dp
      out%seasonal = 0.0_dp
      if (has_trend) out%trend(1) = initial_trend
      if (has_seasonal) out%seasonal(:min(n, period)) = initial_seasonal(:min(n, period))
      out%observation%mean(1) = out%observation%mean(1) + initial_level
      out%observation%residual(1) = series(1) - out%observation%mean(1)
      do t = 2, n
         out%level(t) = level_smoothing*(series(t) - out%observation%mean(t)) + &
            (1.0_dp - level_smoothing)*out%level(t - 1)
         if (has_seasonal .and. t > period) out%level(t) = out%level(t) - &
            level_smoothing*out%seasonal(t - period)
         out%observation%mean(t) = out%observation%mean(t) + out%level(t - 1)
         if (has_trend) then
            out%trend(t - 1) = damp*out%trend(t - 1)
            out%trend(t) = trend_alpha*(out%level(t) - out%level(t - 1)) + &
               (1.0_dp - trend_alpha)*out%trend(t - 1)
            out%observation%mean(t) = out%observation%mean(t) + out%trend(t - 1)
         end if
         if (has_seasonal .and. t > period) then
            out%seasonal(t) = seasonal_alpha*(series(t) - out%level(t)) + &
               (1.0_dp - seasonal_alpha)*out%seasonal(t - period)
            if (has_trend) out%seasonal(t) = out%seasonal(t) - seasonal_alpha*out%trend(t)
            out%observation%mean(t) = out%observation%mean(t) + out%seasonal(t - period)
         end if
         out%observation%residual(t) = series(t) - out%observation%mean(t)
      end do
      out%observation%scale = sigma
      do t = 1, n
         if (present(degrees_of_freedom)) then
            out%observation%log_likelihood(t) = student_log_density(&
               out%observation%residual(t), sigma, degrees_of_freedom)
         else
            out%observation%log_likelihood(t) = normal_log_density(&
               out%observation%residual(t), sigma)
         end if
      end do
      out%observation%total_log_likelihood = sum(out%observation%log_likelihood)
   end function bf_ets_filter

   pure function bf_garch_filter(series, variance_constant, intercept, ar, ma, arch, garch, &
      mean_garch, asymmetry_type, asymmetry_scale, asymmetry_shape, regressors, regression, &
      degrees_of_freedom) result(out)
      !! Evaluate Gaussian or Student asymmetric ARMA-GARCH recursions.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: variance_constant !! Variance constant.
      real(dp), intent(in) :: intercept !! Model intercept.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: mean_garch(:) !! Mean garch.
      integer, intent(in) :: asymmetry_type !! Asymmetry type.
      real(dp), intent(in) :: asymmetry_scale !! Asymmetry scale.
      real(dp), intent(in) :: asymmetry_shape !! Asymmetry shape.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: regression(:) !! Regression.
      real(dp), intent(in), optional :: degrees_of_freedom !! Degrees of freedom.
      type(bf_filter_t) :: out
      real(dp), allocatable :: variance(:)
      real(dp) :: asymmetric
      integer :: n, t, lag

      n = size(series)
      if (n < 1 .or. variance_constant <= 0.0_dp .or. asymmetry_type < 0 .or. &
         asymmetry_type > 2 .or. (present(regressors) .neqv. present(regression))) then
         out%info = 1
         return
      end if
      allocate(out%mean(n), out%residual(n), out%scale(n), out%log_likelihood(n), variance(n))
      out%mean = intercept
      if (present(regressors)) out%mean = out%mean + matmul(regressors, regression)
      out%residual = 0.0_dp
      variance = variance_constant
      do t = 1, n
         do lag = 1, min(size(ar), t - 1)
            out%mean(t) = out%mean(t) + ar(lag)*series(t - lag)
         end do
         do lag = 1, min(size(ma), t - 1)
            out%mean(t) = out%mean(t) + ma(lag)*out%residual(t - lag)
         end do
         do lag = 1, min(size(arch), t - 1)
            variance(t) = variance(t) + arch(lag)*out%residual(t - lag)**2
         end do
         do lag = 1, min(size(garch), t - 1)
            variance(t) = variance(t) + garch(lag)*variance(t - lag)
         end do
         out%scale(t) = sqrt(max(variance(t), tiny(1.0_dp)))
         do lag = 1, min(size(mean_garch), t - 1)
            out%mean(t) = out%mean(t) + mean_garch(lag)*out%scale(t - lag)
         end do
         out%residual(t) = series(t) - out%mean(t)
         if (asymmetry_type > 0 .and. t > 1) then
            asymmetric = asymmetry_function(out%residual(t - 1), asymmetry_shape, asymmetry_type)
            out%scale(t) = out%scale(t) + asymmetry_scale*asymmetric
         end if
         if (present(degrees_of_freedom)) then
            out%log_likelihood(t) = student_log_density(out%residual(t), out%scale(t), &
               degrees_of_freedom)
         else
            out%log_likelihood(t) = normal_log_density(out%residual(t), out%scale(t))
         end if
      end do
      out%total_log_likelihood = sum(out%log_likelihood)
   end function bf_garch_filter

   pure function bf_stochastic_volatility_filter(series, log_variance, state_mean, persistence, &
      state_scale, intercept, ar, ma, regressors, regression) result(out)
      !! Evaluate the stochastic-volatility latent-state and observation densities.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: log_variance(:) !! Log variance.
      real(dp), intent(in) :: state_mean !! State mean.
      real(dp), intent(in) :: persistence !! Persistence.
      real(dp), intent(in) :: state_scale !! State scale.
      real(dp), intent(in) :: intercept !! Model intercept.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: regression(:) !! Regression.
      type(bf_sv_filter_t) :: out
      integer :: n, t, lag

      n = size(series)
      if (size(log_variance) /= n .or. n < 1 .or. state_scale <= 0.0_dp .or. &
         abs(persistence) >= 1.0_dp .or. (present(regressors) .neqv. present(regression))) then
         out%info = 1
         return
      end if
      allocate(out%observation%mean(n), out%observation%residual(n), &
         out%observation%scale(n), out%observation%log_likelihood(n))
      allocate(out%state_log_likelihood(n))
      out%observation%mean = intercept
      if (present(regressors)) out%observation%mean = out%observation%mean + &
         matmul(regressors, regression)
      out%observation%residual = 0.0_dp
      out%observation%scale = exp(0.5_dp*log_variance)
      do t = 1, n
         do lag = 1, min(size(ar), t - 1)
            out%observation%mean(t) = out%observation%mean(t) + ar(lag)*series(t - lag)
         end do
         do lag = 1, min(size(ma), t - 1)
            out%observation%mean(t) = out%observation%mean(t) + &
               ma(lag)*out%observation%residual(t - lag)
         end do
         out%observation%residual(t) = series(t) - out%observation%mean(t)
         out%observation%log_likelihood(t) = normal_log_density(&
            out%observation%residual(t), out%observation%scale(t))
      end do
      out%state_log_likelihood(1) = normal_log_density(log_variance(1) - state_mean, &
         state_scale/sqrt(1.0_dp - persistence**2))
      do t = 2, n
         out%state_log_likelihood(t) = normal_log_density(log_variance(t) - state_mean - &
            persistence*(log_variance(t - 1) - state_mean), state_scale)
      end do
      out%observation%total_log_likelihood = sum(out%observation%log_likelihood)
      out%total_log_likelihood = out%observation%total_log_likelihood + &
         sum(out%state_log_likelihood)
   end function bf_stochastic_volatility_filter

   pure elemental function student_log_density(value, scale, degrees) result(log_density)
      !! Evaluate a centered Student-t log density.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: scale !! Scale.
      real(dp), intent(in) :: degrees !! Degrees.
      real(dp) :: log_density

      log_density = -huge(1.0_dp)
      if (scale > 0.0_dp .and. degrees > 2.0_dp) log_density = &
         log_gamma(0.5_dp*(degrees + 1.0_dp)) - log_gamma(0.5_dp*degrees) - &
         0.5_dp*log(degrees*acos(-1.0_dp)) - log(scale) - &
         0.5_dp*(degrees + 1.0_dp)*log(1.0_dp + (value/scale)**2/degrees)
   end function student_log_density

   pure elemental function asymmetry_function(value, shape, asymmetry_type) result(result_value)
      !! Evaluate the package's logistic or exponential GARCH news function.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: shape !! Shape.
      integer, intent(in) :: asymmetry_type !! Asymmetry type.
      real(dp) :: result_value, weight

      weight = 1.0_dp
      if (asymmetry_type == 1) weight = 1.0_dp/(1.0_dp + exp(shape*value))
      if (asymmetry_type == 2) weight = 1.0_dp - exp(-shape*value**2)
      result_value = value**2*weight
   end function asymmetry_function

   function bf_metropolis_sample(log_posterior, initial, proposal_scale, draw_count, &
      burnin, thin) result(out)
      !! Draw posterior samples with component-scaled random-walk Metropolis.
      procedure(bf_log_posterior_interface) :: log_posterior !! Log posterior callback procedure.
      real(dp), intent(in) :: initial(:) !! Initial value.
      real(dp), intent(in) :: proposal_scale(:) !! Proposal scale.
      integer, intent(in) :: draw_count !! Number of draw.
      integer, intent(in) :: burnin !! Number of initial simulation draws to discard.
      integer, intent(in) :: thin !! Thin.
      type(bf_mcmc_t) :: out
      real(dp), allocatable :: current(:), proposal(:)
      real(dp) :: current_density, proposal_density, uniform_draw
      integer :: total, iteration, stored, accepted, i

      if (size(initial) < 1 .or. size(proposal_scale) /= size(initial) .or. &
         any(proposal_scale <= 0.0_dp) .or. draw_count < 1 .or. burnin < 0 .or. thin < 1) then
         out%info = 1
         return
      end if
      total = burnin + draw_count*thin
      allocate(current(size(initial)), proposal(size(initial)))
      allocate(out%draws(draw_count, size(initial)), out%log_density(draw_count))
      current = initial
      current_density = log_posterior(current)
      if (.not. ieee_is_finite(current_density)) then
         out%info = 2
         return
      end if
      stored = 0
      accepted = 0
      do iteration = 1, total
         do i = 1, size(initial)
            proposal(i) = current(i) + proposal_scale(i)*random_standard_normal()
         end do
         proposal_density = log_posterior(proposal)
         uniform_draw = max(random_uniform(), tiny(1.0_dp))
         if (ieee_is_finite(proposal_density)) then
            if (log(uniform_draw) < proposal_density - current_density) then
               current = proposal
               current_density = proposal_density
               accepted = accepted + 1
            end if
         end if
         if (iteration > burnin .and. mod(iteration - burnin, thin) == 0) then
            stored = stored + 1
            out%draws(stored, :) = current
            out%log_density(stored) = current_density
         end if
      end do
      out%acceptance_rate = real(accepted, dp)/real(total, dp)
      out%burnin = burnin
      out%thin = thin
   end function bf_metropolis_sample

   pure function bf_posterior_interval(draws, probability) result(interval)
      !! Return equal-tail posterior intervals for each sampled parameter.
      real(dp), intent(in) :: draws(:, :) !! Draws.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp), allocatable :: interval(:, :)
      real(dp), allocatable :: ordered(:)
      real(dp) :: lower_position, upper_position
      integer :: draw_count, parameter_count, parameter

      draw_count = size(draws, 1)
      parameter_count = size(draws, 2)
      if (draw_count < 1 .or. probability <= 0.0_dp .or. probability >= 1.0_dp) then
         allocate(interval(0, 0))
         return
      end if
      allocate(interval(parameter_count, 2), ordered(draw_count))
      lower_position = 0.5_dp*(1.0_dp - probability)*real(draw_count - 1, dp) + 1.0_dp
      upper_position = 0.5_dp*(1.0_dp + probability)*real(draw_count - 1, dp) + 1.0_dp
      do parameter = 1, parameter_count
         ordered = draws(:, parameter)
         call sort(ordered)
         interval(parameter, 1) = interpolated_order_statistic(ordered, lower_position)
         interval(parameter, 2) = interpolated_order_statistic(ordered, upper_position)
      end do
   end function bf_posterior_interval

   pure function bf_predictive_error(predictive_draws, observations) result(error)
      !! Subtract posterior predictive draws from observed values.
      real(dp), intent(in) :: predictive_draws(:, :) !! Predictive simulation draws.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), allocatable :: error(:, :)

      if (size(predictive_draws, 2) /= size(observations)) then
         allocate(error(0, 0))
         return
      end if
      error = spread(observations, 1, size(predictive_draws, 1)) - predictive_draws
   end function bf_predictive_error

   pure function bf_waic(pointwise_log_likelihood) result(out)
      !! Compute WAIC from draws by observation pointwise log likelihoods.
      real(dp), intent(in) :: pointwise_log_likelihood(:, :) !! Pointwise log likelihood.
      type(bf_information_criteria_t) :: out
      real(dp), allocatable :: shifted(:)
      real(dp) :: maximum, mean_value
      integer :: draws, observations, observation

      draws = size(pointwise_log_likelihood, 1)
      observations = size(pointwise_log_likelihood, 2)
      if (draws < 2 .or. observations < 1) then
         out%info = 1
         return
      end if
      do observation = 1, observations
         maximum = maxval(pointwise_log_likelihood(:, observation))
         shifted = exp(pointwise_log_likelihood(:, observation) - maximum)
         out%lppd = out%lppd + maximum + log(sum(shifted)/real(draws, dp))
         mean_value = sum(pointwise_log_likelihood(:, observation))/real(draws, dp)
         out%effective_parameters = out%effective_parameters + &
            sum((pointwise_log_likelihood(:, observation) - mean_value)**2)/real(draws - 1, dp)
      end do
      out%waic = -2.0_dp*(out%lppd - out%effective_parameters)
   end function bf_waic

   pure function bf_log_prior_sum(values, specifications) result(log_density)
      !! Sum independently encoded bayesforecast prior log densities.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp), intent(in) :: specifications(:, :) !! Specifications.
      real(dp) :: log_density
      integer :: i

      log_density = -huge(1.0_dp)
      if (size(specifications, 1) /= size(values) .or. size(specifications, 2) /= 4) return
      log_density = 0.0_dp
      do i = 1, size(values)
         log_density = log_density + bf_log_prior(values(i), specifications(i, :))
      end do
   end function bf_log_prior_sum

   pure function bf_log_posterior(log_likelihood, parameter_values, prior_specifications) &
      result(value)
      !! Combine any translated model likelihood with its independent priors.
      real(dp), intent(in) :: log_likelihood !! Log-likelihood value.
      real(dp), intent(in) :: parameter_values(:) !! Parameter values.
      real(dp), intent(in) :: prior_specifications(:, :) !! Prior specifications.
      real(dp) :: value

      value = log_likelihood + bf_log_prior_sum(parameter_values, prior_specifications)
   end function bf_log_posterior

   pure function bf_information_criteria(log_likelihood, parameter_count, observation_count) result(out)
      !! Compute log likelihood, AIC, corrected AIC, and BIC.
      real(dp), intent(in) :: log_likelihood !! Log-likelihood value.
      integer, intent(in) :: parameter_count !! Number of parameter.
      integer, intent(in) :: observation_count !! Observation count.
      type(bf_information_criteria_t) :: out

      if (parameter_count < 0 .or. observation_count < 1) then
         out%info = 1
         return
      end if
      out%log_likelihood = log_likelihood
      out%aic = -2.0_dp*log_likelihood + 2.0_dp*real(parameter_count, dp)
      out%aicc = out%aic
      if (observation_count > parameter_count + 1) out%aicc = out%aic + &
         2.0_dp*real(parameter_count*(parameter_count + 1), dp)/ &
         real(observation_count - parameter_count - 1, dp)
      out%bic = -2.0_dp*log_likelihood + log(real(observation_count, dp))* &
         real(parameter_count, dp)
   end function bf_information_criteria

   pure elemental function bf_bayes_factor(first_log_marginal, second_log_marginal, &
      log_scale) result(value)
      !! Compare two supplied log marginal likelihoods on log or natural scale.
      real(dp), intent(in) :: first_log_marginal !! First log marginal.
      real(dp), intent(in) :: second_log_marginal !! Second log marginal.
      logical, intent(in), optional :: log_scale !! Flag controlling log scale.
      real(dp) :: value
      logical :: return_log

      return_log = .false.
      if (present(log_scale)) return_log = log_scale
      value = first_log_marginal - second_log_marginal
      if (.not. return_log) value = exp(value)
   end function bf_bayes_factor

   pure function bf_sarima_log_posterior(series, sigma, intercept, ar, ma, seasonal_ar, &
      seasonal_ma, period, parameter_values, prior_specifications, regressors, regression) &
      result(log_posterior)
      !! Combine the SARIMA likelihood with independently specified parameter priors.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: sigma !! Scale parameter or standard deviation.
      real(dp), intent(in) :: intercept !! Model intercept.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: seasonal_ar(:) !! Seasonal autoregressive.
      real(dp), intent(in) :: seasonal_ma(:) !! Seasonal moving-average.
      real(dp), intent(in) :: parameter_values(:) !! Parameter values.
      real(dp), intent(in) :: prior_specifications(:, :) !! Prior specifications.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: regression(:) !! Regression.
      real(dp) :: log_posterior
      type(bf_filter_t) :: filtered

      if (present(regressors)) then
         filtered = bf_sarima_filter(series, sigma, intercept, ar, ma, seasonal_ar, &
            seasonal_ma, period, regressors, regression)
      else
         filtered = bf_sarima_filter(series, sigma, intercept, ar, ma, seasonal_ar, &
            seasonal_ma, period)
      end if
      if (filtered%info /= 0) then
         log_posterior = -huge(1.0_dp)
         return
      end if
      log_posterior = filtered%total_log_likelihood + &
         bf_log_prior_sum(parameter_values, prior_specifications)
   end function bf_sarima_log_posterior

   pure function bf_difference(series, ordinary_order, seasonal_order, period) result(out)
      !! Difference a series while retaining boundaries for exact inversion.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: ordinary_order !! Ordinary order.
      integer, intent(in) :: seasonal_order !! Seasonal order.
      integer, intent(in) :: period !! Seasonal period.
      type(bf_difference_t) :: out
      real(dp), allocatable :: working(:)
      integer :: i, n

      if (size(series) < 1 .or. ordinary_order < 0 .or. seasonal_order < 0 .or. &
         period < 1 .or. size(series) <= ordinary_order + seasonal_order*period) then
         out%info = 1
         return
      end if
      working = series
      allocate(out%ordinary_initial(ordinary_order))
      do i = 1, ordinary_order
         out%ordinary_initial(i) = working(size(working))
         working = working(2:) - working(:size(working) - 1)
      end do
      allocate(out%seasonal_initial(seasonal_order, period))
      do i = 1, seasonal_order
         n = size(working)
         out%seasonal_initial(i, :) = working(n - period + 1:n)
         working = working(period + 1:) - working(:n - period)
      end do
      out%values = working
      out%ordinary_order = ordinary_order
      out%seasonal_order = seasonal_order
      out%period = period
   end function bf_difference

   pure function bf_inverse_difference(values, difference) result(series)
      !! Undo forecast differencing using stored final observed boundaries.
      real(dp), intent(in) :: values(:) !! Input values.
      type(bf_difference_t), intent(in) :: difference !! Difference.
      real(dp), allocatable :: series(:), restored(:)
      integer :: i, t, period

      if (difference%info /= 0) then
         allocate(series(0))
         return
      end if
      restored = values
      period = difference%period
      do i = difference%seasonal_order, 1, -1
         do t = 1, size(restored)
            if (t <= period) then
               restored(t) = restored(t) + difference%seasonal_initial(i, t)
            else
               restored(t) = restored(t) + restored(t - period)
            end if
         end do
      end do
      do i = difference%ordinary_order, 1, -1
         restored(1) = restored(1) + difference%ordinary_initial(i)
         do t = 2, size(restored)
            restored(t) = restored(t) + restored(t - 1)
         end do
      end do
      series = restored
   end function bf_inverse_difference

   function bf_sarima_predict(series, sigma, intercept, ar, ma, seasonal_ar, seasonal_ma, &
      period, horizon, draw_count, probability, regressors, future_regressors, regression) result(out)
      !! Simulate SARIMA posterior predictive paths for fixed parameter values.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: sigma !! Scale parameter or standard deviation.
      real(dp), intent(in) :: intercept !! Model intercept.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: seasonal_ar(:) !! Seasonal autoregressive.
      real(dp), intent(in) :: seasonal_ma(:) !! Seasonal moving-average.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: draw_count !! Number of draw.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future regressors.
      real(dp), intent(in), optional :: regression(:) !! Regression.
      type(bf_predictive_t) :: out
      type(bf_filter_t) :: history_filter
      real(dp), allocatable :: extended(:), residual(:)
      real(dp) :: location
      integer :: n, draw, step, lag

      n = size(series)
      if (sigma <= 0.0_dp .or. horizon < 1 .or. draw_count < 1 .or. &
         probability <= 0.0_dp .or. probability >= 1.0_dp .or. &
         n < max(size(ar), size(ma), size(seasonal_ar)*period, size(seasonal_ma)*period) .or. &
         (present(regressors) .neqv. present(regression)) .or. &
         (present(future_regressors) .neqv. present(regression))) then
         out%info = 1
         return
      end if
      if (present(regressors)) then
         history_filter = bf_sarima_filter(series, sigma, intercept, ar, ma, seasonal_ar, &
            seasonal_ma, period, regressors, regression)
      else
         history_filter = bf_sarima_filter(series, sigma, intercept, ar, ma, seasonal_ar, &
            seasonal_ma, period)
      end if
      if (history_filter%info /= 0) then
         out%info = 10 + history_filter%info
         return
      end if
      if (present(future_regressors)) then
         if (size(future_regressors, 1) /= horizon .or. &
            size(future_regressors, 2) /= size(regression)) then
            out%info = 1
            return
         end if
      end if
      allocate(out%draws(draw_count, horizon), extended(n + horizon), residual(n + horizon))
      do draw = 1, draw_count
         extended(:n) = series
         residual(:n) = history_filter%residual
         do step = 1, horizon
            location = intercept
            if (present(future_regressors)) location = location + &
               dot_product(future_regressors(step, :), regression)
            do lag = 1, size(ar)
               location = location + ar(lag)*extended(n + step - lag)
            end do
            do lag = 1, size(ma)
               if (n + step - lag <= n) location = location + ma(lag)*residual(n + step - lag)
            end do
            do lag = 1, size(seasonal_ar)
               location = location + seasonal_ar(lag)*extended(n + step - lag*period)
            end do
            do lag = 1, size(seasonal_ma)
               if (n + step - lag*period <= n) location = location + &
                  seasonal_ma(lag)*residual(n + step - lag*period)
            end do
            residual(n + step) = sigma*random_standard_normal()
            extended(n + step) = location + residual(n + step)
            out%draws(draw, step) = extended(n + step)
         end do
      end do
      allocate(out%mean(horizon))
      out%mean = sum(out%draws, 1)/real(draw_count, dp)
      out%interval = bf_posterior_interval(out%draws, probability)
      out%horizon = horizon
   end function bf_sarima_predict

   function bf_ets_predict(series, sigma, level_smoothing, initial_level, horizon, &
      draw_count, probability, trend_smoothing, initial_trend, damping, &
      seasonal_smoothing, initial_seasonal, degrees_of_freedom, regressors, &
      future_regressors, regression) result(out)
      !! Simulate posterior predictive paths from a fixed ETS parameter draw.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: sigma !! Scale parameter or standard deviation.
      real(dp), intent(in) :: level_smoothing !! Level smoothing.
      real(dp), intent(in) :: initial_level !! Initial level.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: draw_count !! Number of draw.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp), intent(in), optional :: trend_smoothing !! Trend smoothing.
      real(dp), intent(in), optional :: initial_trend !! Initial trend.
      real(dp), intent(in), optional :: damping !! Damping.
      real(dp), intent(in), optional :: seasonal_smoothing !! Seasonal smoothing.
      real(dp), intent(in), optional :: initial_seasonal(:) !! Initial seasonal.
      real(dp), intent(in), optional :: degrees_of_freedom !! Degrees of freedom.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future regressors.
      real(dp), intent(in), optional :: regression(:) !! Regression.
      type(bf_predictive_t) :: out
      type(bf_ets_filter_t) :: history
      real(dp), allocatable :: seasonal_path(:)
      real(dp) :: level_value, trend_value, damped_trend, old_seasonal
      real(dp) :: location, innovation, generated, new_level, new_trend, new_seasonal
      real(dp) :: regression_value
      logical :: has_trend, has_seasonal, has_regression
      integer :: n, period, draw, step, time_index

      n = size(series)
      has_trend = present(trend_smoothing) .and. present(initial_trend)
      has_seasonal = present(seasonal_smoothing) .and. present(initial_seasonal)
      has_regression = present(regressors) .and. present(future_regressors) .and. present(regression)
      period = 1
      if (has_seasonal) period = size(initial_seasonal)
      if (horizon < 1 .or. draw_count < 1 .or. probability <= 0.0_dp .or. n < period .or. &
         probability >= 1.0_dp .or. &
         (present(trend_smoothing) .neqv. present(initial_trend)) .or. &
         (present(seasonal_smoothing) .neqv. present(initial_seasonal)) .or. &
         ((present(regressors) .or. present(future_regressors) .or. present(regression)) .and. &
         .not. has_regression)) then
         out%info = 1
         return
      end if
      if (has_regression) then
         if (size(regressors, 1) /= n .or. size(future_regressors, 1) /= horizon .or. &
            size(regressors, 2) /= size(regression) .or. &
            size(future_regressors, 2) /= size(regression)) then
            out%info = 1
            return
         end if
      end if
      if (has_regression) then
         history = bf_ets_filter(series, sigma, level_smoothing, initial_level, &
            trend_smoothing, initial_trend, damping, seasonal_smoothing, initial_seasonal, &
            regressors, regression, degrees_of_freedom)
      else
         history = bf_ets_filter(series, sigma, level_smoothing, initial_level, &
            trend_smoothing, initial_trend, damping, seasonal_smoothing, initial_seasonal, &
            degrees_of_freedom=degrees_of_freedom)
      end if
      if (history%info /= 0) then
         out%info = 10 + history%info
         return
      end if
      allocate(out%draws(draw_count, horizon), seasonal_path(n + horizon))
      do draw = 1, draw_count
         level_value = history%level(n)
         trend_value = history%trend(n)
         seasonal_path = 0.0_dp
         seasonal_path(:n) = history%seasonal
         do step = 1, horizon
            time_index = n + step
            location = level_value
            regression_value = 0.0_dp
            if (has_regression) regression_value = dot_product(future_regressors(step, :), regression)
            location = location + regression_value
            damped_trend = trend_value
            if (present(damping)) damped_trend = damping*trend_value
            if (has_trend) location = location + damped_trend
            old_seasonal = 0.0_dp
            if (has_seasonal) then
               old_seasonal = seasonal_path(time_index - period)
               location = location + old_seasonal
            end if
            innovation = sigma*random_standard_normal()
            if (present(degrees_of_freedom)) innovation = sigma*random_student(degrees_of_freedom)
            generated = location + innovation
            new_level = level_smoothing*(generated - regression_value) + &
               (1.0_dp - level_smoothing)*level_value
            if (has_seasonal) new_level = new_level - level_smoothing*old_seasonal
            new_trend = 0.0_dp
            if (has_trend) new_trend = trend_smoothing*(new_level - level_value) + &
               (1.0_dp - trend_smoothing)*damped_trend
            new_seasonal = 0.0_dp
            if (has_seasonal) then
               new_seasonal = seasonal_smoothing*(generated - new_level) + &
                  (1.0_dp - seasonal_smoothing)*old_seasonal
               if (has_trend) new_seasonal = new_seasonal - seasonal_smoothing*new_trend
               seasonal_path(time_index) = new_seasonal
            end if
            out%draws(draw, step) = generated
            level_value = new_level
            if (has_trend) trend_value = new_trend
         end do
      end do
      call finish_predictive(out, probability)
   end function bf_ets_predict

   function bf_garch_predict(series, variance_constant, intercept, ar, ma, arch, garch, &
      mean_garch, asymmetry_type, asymmetry_scale, asymmetry_shape, horizon, draw_count, &
      probability, degrees_of_freedom, regressors, future_regressors, regression) result(out)
      !! Simulate predictive paths from a fixed asymmetric ARMA-GARCH draw.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: variance_constant !! Variance constant.
      real(dp), intent(in) :: intercept !! Model intercept.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: mean_garch(:) !! Mean garch.
      integer, intent(in) :: asymmetry_type !! Asymmetry type.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: draw_count !! Number of draw.
      real(dp), intent(in) :: asymmetry_scale !! Asymmetry scale.
      real(dp), intent(in) :: asymmetry_shape !! Asymmetry shape.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp), intent(in), optional :: degrees_of_freedom !! Degrees of freedom.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future regressors.
      real(dp), intent(in), optional :: regression(:) !! Regression.
      type(bf_predictive_t) :: out
      type(bf_filter_t) :: history
      real(dp), allocatable :: extended(:), residual(:), scale(:)
      real(dp) :: location, variance_value, innovation
      logical :: has_regression
      integer :: n, draw, step, lag, index

      n = size(series)
      has_regression = present(regressors) .and. present(future_regressors) .and. present(regression)
      if (horizon < 1 .or. draw_count < 1 .or. probability <= 0.0_dp .or. &
         probability >= 1.0_dp .or. asymmetry_scale < 0.0_dp .or. asymmetry_shape < 0.0_dp .or. &
         n < max(size(ar), size(ma), size(arch), size(garch), size(mean_garch)) .or. &
         ((present(regressors) .or. present(future_regressors) .or. &
         present(regression)) .and. .not. has_regression)) then
         out%info = 1
         return
      end if
      if (has_regression) then
         history = bf_garch_filter(series, variance_constant, intercept, ar, ma, arch, garch, &
            mean_garch, asymmetry_type, asymmetry_scale, asymmetry_shape, regressors, regression, &
            degrees_of_freedom)
      else
         history = bf_garch_filter(series, variance_constant, intercept, ar, ma, arch, garch, &
            mean_garch, asymmetry_type, asymmetry_scale, asymmetry_shape, &
            degrees_of_freedom=degrees_of_freedom)
      end if
      if (history%info /= 0) then
         out%info = 10 + history%info
         return
      end if
      allocate(out%draws(draw_count, horizon), extended(n + horizon))
      allocate(residual(n + horizon), scale(n + horizon))
      do draw = 1, draw_count
         extended(:n) = series
         residual(:n) = history%residual
         scale(:n) = history%scale
         do step = 1, horizon
            index = n + step
            location = intercept
            if (has_regression) location = location + &
               dot_product(future_regressors(step, :), regression)
            do lag = 1, size(ar)
               location = location + ar(lag)*extended(index - lag)
            end do
            do lag = 1, size(ma)
               location = location + ma(lag)*residual(index - lag)
            end do
            variance_value = variance_constant
            do lag = 1, size(arch)
               variance_value = variance_value + arch(lag)*residual(index - lag)**2
            end do
            do lag = 1, size(garch)
               variance_value = variance_value + garch(lag)*scale(index - lag)**2
            end do
            scale(index) = sqrt(max(variance_value, tiny(1.0_dp)))
            do lag = 1, size(mean_garch)
               location = location + mean_garch(lag)*scale(index - lag)
            end do
            if (asymmetry_type > 0) scale(index) = scale(index) + asymmetry_scale* &
               asymmetry_function(residual(index - 1), asymmetry_shape, asymmetry_type)
            innovation = scale(index)*random_standard_normal()
            if (present(degrees_of_freedom)) innovation = &
               scale(index)*random_student(degrees_of_freedom)
            residual(index) = innovation
            extended(index) = location + innovation
            out%draws(draw, step) = extended(index)
         end do
      end do
      call finish_predictive(out, probability)
   end function bf_garch_predict

   function bf_stochastic_volatility_predict(series, log_variance, state_mean, persistence, &
      state_scale, intercept, ar, ma, horizon, draw_count, probability, regressors, &
      future_regressors, regression) result(out)
      !! Simulate stochastic-volatility latent states and observations forward.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: log_variance(:) !! Log variance.
      real(dp), intent(in) :: state_mean !! State mean.
      real(dp), intent(in) :: persistence !! Persistence.
      real(dp), intent(in) :: state_scale !! State scale.
      real(dp), intent(in) :: intercept !! Model intercept.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: probability !! Probability value.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: draw_count !! Number of draw.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future regressors.
      real(dp), intent(in), optional :: regression(:) !! Regression.
      type(bf_predictive_t) :: out
      type(bf_sv_filter_t) :: history
      real(dp), allocatable :: extended(:), residual(:)
      real(dp) :: latent, location, innovation
      logical :: has_regression
      integer :: n, draw, step, lag, index

      n = size(series)
      has_regression = present(regressors) .and. present(future_regressors) .and. present(regression)
      if (horizon < 1 .or. draw_count < 1 .or. probability <= 0.0_dp .or. &
         probability >= 1.0_dp .or. n < max(size(ar), size(ma)) .or. &
         ((present(regressors) .or. present(future_regressors) .or. &
         present(regression)) .and. .not. has_regression)) then
         out%info = 1
         return
      end if
      if (has_regression) then
         history = bf_stochastic_volatility_filter(series, log_variance, state_mean, persistence, &
            state_scale, intercept, ar, ma, regressors, regression)
      else
         history = bf_stochastic_volatility_filter(series, log_variance, state_mean, persistence, &
            state_scale, intercept, ar, ma)
      end if
      if (history%info /= 0) then
         out%info = 10 + history%info
         return
      end if
      allocate(out%draws(draw_count, horizon), extended(n + horizon), residual(n + horizon))
      do draw = 1, draw_count
         extended(:n) = series
         residual(:n) = history%observation%residual
         latent = log_variance(n)
         do step = 1, horizon
            index = n + step
            latent = state_mean + persistence*(latent - state_mean) + &
               state_scale*random_standard_normal()
            location = intercept
            if (has_regression) location = location + &
               dot_product(future_regressors(step, :), regression)
            do lag = 1, size(ar)
               location = location + ar(lag)*extended(index - lag)
            end do
            do lag = 1, size(ma)
               location = location + ma(lag)*residual(index - lag)
            end do
            innovation = exp(0.5_dp*latent)*random_standard_normal()
            residual(index) = innovation
            extended(index) = location + innovation
            out%draws(draw, step) = extended(index)
         end do
      end do
      call finish_predictive(out, probability)
   end function bf_stochastic_volatility_predict

   subroutine finish_predictive(out, probability)
      !! Complete predictive means, intervals, and horizon metadata.
      type(bf_predictive_t), intent(inout) :: out !! Procedure result, updated in place.
      real(dp), intent(in) :: probability !! Probability value.

      allocate(out%mean(size(out%draws, 2)))
      out%mean = sum(out%draws, 1)/real(size(out%draws, 1), dp)
      out%interval = bf_posterior_interval(out%draws, probability)
      out%horizon = size(out%draws, 2)
   end subroutine finish_predictive

   function random_student(degrees) result(value)
      !! Draw a standardized Student-t random variate.
      real(dp), intent(in) :: degrees !! Degrees.
      real(dp) :: value, chi_square

      chi_square = 2.0_dp*random_gamma(0.5_dp*degrees)
      value = random_standard_normal()/sqrt(chi_square/degrees)
   end function random_student

   pure function interpolated_order_statistic(values, position) result(value)
      !! Interpolate a sorted vector at a one-based fractional position.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp), intent(in) :: position !! Position.
      real(dp) :: value, fraction
      integer :: lower

      lower = max(1, min(size(values), floor(position)))
      if (lower == size(values)) then
         value = values(lower)
      else
         fraction = position - real(lower, dp)
         value = (1.0_dp - fraction)*values(lower) + fraction*values(lower + 1)
      end if
   end function interpolated_order_statistic
end module bayesforecast_mod
