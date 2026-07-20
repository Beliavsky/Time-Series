! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R TSSS package.
! Distinct numerical algorithms translated from the GPL-2 TSSS package.
module tsss_mod
   use kind_mod, only: dp
   use stats_mod, only: sort
   use linalg_mod, only: outer_product, solve_matrix, &
      cholesky_lower_semidefinite
   use random_mod, only: random_uniform, random_standard_normal, random_gamma
   use kfas_mod, only: ssm_model_t, kfs_filter_t, kfs_smoother_t, kfs_filter, kfs_smooth
   use forecast_mod, only: box_cox
   implicit none
   private

   integer, parameter, public :: tsss_system_gaussian = 0
   integer, parameter, public :: tsss_system_cauchy = 1
   integer, parameter, public :: tsss_system_mixture = 2
   integer, parameter, public :: tsss_initial_normal = 0
   integer, parameter, public :: tsss_initial_uniform = 1
   integer, parameter, public :: tsss_initial_cauchy = 2
   integer, parameter, public :: tsss_initial_zero = 3
   integer, parameter, public :: tsss_noise_gaussian = 1
   integer, parameter, public :: tsss_noise_pearson = 2
   integer, parameter, public :: tsss_noise_laplace = 3
   integer, parameter, public :: tsss_density_two_sided_exponential = 0
   integer, parameter, public :: tsss_density_gaussian = 1
   integer, parameter, public :: tsss_density_cauchy = 2
   integer, parameter, public :: tsss_density_pearson = 3
   integer, parameter, public :: tsss_density_exponential = 4
   integer, parameter, public :: tsss_density_chi_square = 5
   integer, parameter, public :: tsss_density_log_exponential = 6
   integer, parameter, public :: tsss_density_uniform = 7

   type, public :: tsss_particle_filter_t
      ! Fixed-lag particle smoothing quantiles, means, and log likelihood.
      real(dp), allocatable :: quantile(:, :), mean(:)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: particles = 0
      integer :: lag = 0
      integer :: info = 0
   end type

   type, public :: tsss_tvvar_t
      ! Smoothed log variance, variance, normalized observations, and diagnostics.
      real(dp), allocatable :: log_variance(:), variance(:), normalized(:)
      real(dp) :: system_variance = 0.0_dp
      real(dp) :: observation_variance = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      integer :: info = 0
   end type

   type, public :: tsss_tvar_t
      ! Smoothed time-varying AR and partial-autocorrelation coefficients.
      real(dp), allocatable :: ar(:, :), parcor(:, :)
      real(dp) :: system_variance = 0.0_dp
      real(dp) :: observation_variance = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      integer :: span = 0
      integer :: info = 0
   end type

   type, public :: tsss_tvspc_t
      ! Evolutionary log10 power spectrum on frequencies from zero to one half.
      real(dp), allocatable :: frequency(:), log_spectrum(:, :)
      integer :: info = 0
   end type

   type, public :: tsss_ngsmooth_t
      ! Grid densities, seven posterior quantiles, and marginal log likelihood.
      real(dp), allocatable :: grid(:), density(:, :), quantile(:, :)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: info = 0
   end type

   type, public :: tsss_simulation_t
      ! Simulated states and scalar observations from a state-space model.
      real(dp), allocatable :: state(:, :), observation(:)
      integer :: info = 0
   end type

   type, public :: tsss_structural_model_t
      ! Matrices and initial state for a TSSS structural time-series model.
      real(dp), allocatable :: transition(:, :), system_loading(:, :)
      real(dp), allocatable :: observation_loading(:), system_covariance(:, :), initial_state(:)
      real(dp) :: observation_variance = 0.0_dp
      integer :: trend_order = 0
      integer :: seasonal_order = 0
      integer :: period = 0
      integer :: ar_order = 0
      integer :: info = 0
   end type

   type, public :: tsss_tsmooth_t
      ! Predicted, filtered, and smoothed moments from the TSSS general model.
      real(dp), allocatable :: predicted_state(:, :), filtered_state(:, :), smoothed_state(:, :)
      real(dp), allocatable :: predicted_covariance(:, :, :), filtered_covariance(:, :, :)
      real(dp), allocatable :: smoothed_covariance(:, :, :), smoothed_variance(:, :)
      real(dp), allocatable :: prediction_error(:)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      integer :: observations = 0
      integer :: info = 0
   end type

   type, public :: tsss_density_t
      ! Density values evaluated on an equally spaced grid.
      real(dp), allocatable :: grid(:), density(:)
      integer :: model = 0
      integer :: info = 0
   end type

   type, public :: tsss_kl_t
      ! Four successive trapezoidal approximations to KL information and mass.
      integer :: intervals(4) = 0
      real(dp) :: spacing(4) = 0.0_dp
      real(dp) :: information(4) = 0.0_dp
      real(dp) :: reference_mass(4) = 0.0_dp
      integer :: info = 0
   end type

   type, public :: tsss_boxcox_t
      ! TSSS Box-Cox likelihood table and transform selected by adjusted AIC.
      real(dp) :: lambda(21) = 0.0_dp
      real(dp) :: adjusted_aic(21) = 0.0_dp, adjusted_log_likelihood(21) = 0.0_dp
      real(dp) :: aic(21) = 0.0_dp, log_likelihood(21) = 0.0_dp
      real(dp) :: mean(21) = 0.0_dp, variance(21) = 0.0_dp
      real(dp), allocatable :: transformed(:)
      real(dp) :: best_lambda = 0.0_dp
      integer :: info = 0
   end type

   public :: tsss_particle_filter, tsss_nonlinear_particle_filter
   public :: tsss_tvvar, tsss_tvar, tsss_tvspc
   public :: tsss_ngsmth, tsss_simssm, tsss_ngsim
   public :: tsss_structural_model
   public :: tsss_tsmooth
   public :: tsss_pdfunc, tsss_density_value, tsss_klinfo, tsss_boxcox

contains

   function tsss_particle_filter(observations, particle_count, system_model, smoothing_lag, &
      initial_distribution, observation_variance, system_variance, mixture_weight, &
      large_system_variance, initial_variance, lower_bound, upper_bound) result(out)
      !! Run the TSSS random-walk particle filter and fixed-lag smoother.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      integer, intent(in) :: particle_count !! Number of particle.
      integer, intent(in) :: system_model !! System model.
      integer, intent(in) :: smoothing_lag !! Smoothing lag.
      integer, intent(in) :: initial_distribution !! Initial distribution.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: system_variance !! System variance.
      real(dp), intent(in) :: mixture_weight !! Mixture weight.
      real(dp), intent(in) :: large_system_variance !! Large system variance.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: lower_bound !! Lower bound for bound.
      real(dp), intent(in) :: upper_bound !! Upper bound for bound.
      type(tsss_particle_filter_t) :: out
      real(dp), allocatable :: particles(:), predicted(:), weights(:), history(:, :), next_history(:, :)
      real(dp) :: innovation
      integer :: i

      if (size(observations) < 1 .or. particle_count < 2 .or. smoothing_lag < 0 .or. &
         system_model < 0 .or. system_model > 2 .or. initial_distribution < 0 .or. &
         initial_distribution > 3 .or. observation_variance <= 0.0_dp .or. &
         system_variance < 0.0_dp .or. large_system_variance < 0.0_dp .or. &
         initial_variance < 0.0_dp .or. mixture_weight < 0.0_dp .or. mixture_weight > 1.0_dp .or. &
         lower_bound >= upper_bound) then
         out%info = 1
         return
      end if
      allocate(particles(particle_count), predicted(particle_count), weights(particle_count))
      allocate(history(particle_count, 0:smoothing_lag), next_history(particle_count, 0:smoothing_lag))
      do i = 1, particle_count
         select case (initial_distribution)
         case (tsss_initial_normal)
            particles(i) = sqrt(initial_variance)*random_standard_normal()
         case (tsss_initial_uniform)
            particles(i) = lower_bound + (upper_bound - lower_bound)*random_uniform()
         case (tsss_initial_cauchy)
            particles(i) = sqrt(initial_variance)*tan(acos(-1.0_dp)*(random_uniform() - 0.5_dp))
         case default
            particles(i) = 0.0_dp
         end select
      end do
      history = spread(particles, 2, smoothing_lag + 1)
      call run_particle_filter(observations, particles, predicted, weights, history, next_history, &
         observation_variance, lower_bound, upper_bound, out, transition)
      out%particles = particle_count
      out%lag = smoothing_lag

   contains

      function transition(previous, time_index) result(value)
         !! Draw one random-walk system transition.
         real(dp), intent(in) :: previous !! Previous.
         integer, intent(in) :: time_index !! Index of time.
         real(dp) :: value

         select case (system_model)
         case (tsss_system_gaussian)
            innovation = sqrt(system_variance)*random_standard_normal()
         case (tsss_system_cauchy)
            innovation = sqrt(system_variance)*tan(acos(-1.0_dp)*(random_uniform() - 0.5_dp))
         case default
            if (random_uniform() <= mixture_weight) then
               innovation = sqrt(system_variance)*random_standard_normal()
            else
               innovation = sqrt(large_system_variance)*random_standard_normal()
            end if
         end select
         value = previous + innovation + 0.0_dp*real(time_index, dp)
      end function transition
   end function tsss_particle_filter

   function tsss_nonlinear_particle_filter(observations, particle_count, smoothing_lag, &
      observation_variance, system_variance, lower_bound, upper_bound) result(out)
      !! Run the nonlinear TSSS benchmark particle filter and smoother.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      integer, intent(in) :: particle_count !! Number of particle.
      integer, intent(in) :: smoothing_lag !! Smoothing lag.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: system_variance !! System variance.
      real(dp), intent(in) :: lower_bound !! Lower bound for bound.
      real(dp), intent(in) :: upper_bound !! Upper bound for bound.
      type(tsss_particle_filter_t) :: out
      real(dp), allocatable :: particles(:), predicted(:), weights(:), history(:, :), next_history(:, :)
      integer :: i

      if (size(observations) < 1 .or. particle_count < 2 .or. smoothing_lag < 0 .or. &
         observation_variance <= 0.0_dp .or. system_variance < 0.0_dp .or. &
         lower_bound >= upper_bound) then
         out%info = 1
         return
      end if
      allocate(particles(particle_count), predicted(particle_count), weights(particle_count))
      allocate(history(particle_count, 0:smoothing_lag), next_history(particle_count, 0:smoothing_lag))
      do i = 1, particle_count
         particles(i) = sqrt(5.0_dp)*random_standard_normal()
      end do
      history = spread(particles, 2, smoothing_lag + 1)
      call run_particle_filter(observations, particles, predicted, weights, history, next_history, &
         observation_variance, lower_bound, upper_bound, out, transition, nonlinear_observation=.true.)
      out%particles = particle_count
      out%lag = smoothing_lag

   contains

      function transition(previous, time_index) result(value)
         !! Draw the nonlinear benchmark system transition.
         real(dp), intent(in) :: previous !! Previous.
         integer, intent(in) :: time_index !! Index of time.
         real(dp) :: value

         value = 0.5_dp*previous + 25.0_dp*previous/(1.0_dp + previous**2) + &
            8.0_dp*cos(1.2_dp*real(time_index, dp)) + &
            sqrt(system_variance)*random_standard_normal()
      end function transition
   end function tsss_nonlinear_particle_filter

   subroutine run_particle_filter(observations, particles, predicted, weights, history, &
      next_history, observation_variance, lower_bound, upper_bound, out, transition, &
      nonlinear_observation)
      !! Apply weighting, systematic resampling, and fixed-lag summaries.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: lower_bound !! Lower bound for bound.
      real(dp), intent(in) :: upper_bound !! Upper bound for bound.
      real(dp), intent(inout) :: particles(:) !! Number of particles, updated in place.
      real(dp), intent(inout) :: predicted(:) !! Predicted values, updated in place.
      real(dp), intent(inout) :: weights(:) !! Observation or objective weights, updated in place.
      real(dp), intent(inout) :: history(:, 0:) !! History, updated in place.
      real(dp), intent(inout) :: next_history(:, 0:) !! Next history, updated in place.
      type(tsss_particle_filter_t), intent(out) :: out !! Procedure result.
      interface
         function transition(previous, time_index) result(value) !! State-transition callback procedure.
            !! Propagate one particle to the next time index.
            import dp
            real(dp), intent(in) :: previous !! Previous.
            integer, intent(in) :: time_index !! Index of time.
            real(dp) :: value
         end function transition
      end interface
      logical, intent(in), optional :: nonlinear_observation !! Flag controlling nonlinear observation.
      real(dp), allocatable :: log_weight(:), cumulative(:), summary_values(:)
      real(dp) :: maximum, total, target, observed_mean
      logical :: nonlinear
      integer :: n, particle_count, lag_count, time, i, j, ancestor, output_time

      n = size(observations)
      particle_count = size(particles)
      lag_count = ubound(history, 2)
      nonlinear = .false.
      if (present(nonlinear_observation)) nonlinear = nonlinear_observation
      allocate(out%quantile(n, 7), out%mean(n), log_weight(particle_count))
      allocate(cumulative(particle_count), summary_values(particle_count))
      out%quantile = 0.0_dp
      out%mean = 0.0_dp
      do time = 1, n
         do i = 1, particle_count
            predicted(i) = transition(particles(i), time)
            observed_mean = predicted(i)
            if (nonlinear) observed_mean = predicted(i)**2/20.0_dp
            log_weight(i) = -0.5_dp*(log(2.0_dp*acos(-1.0_dp)*observation_variance) + &
               (observations(time) - observed_mean)**2/observation_variance)
         end do
         maximum = maxval(log_weight)
         weights = exp(log_weight - maximum)
         total = sum(weights)
         if (total <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
         out%log_likelihood = out%log_likelihood + maximum + &
            log(total/real(particle_count, dp))
         cumulative(1) = weights(1)/total
         do i = 2, particle_count
            cumulative(i) = cumulative(i - 1) + weights(i)/total
         end do
         target = random_uniform()/real(particle_count, dp)
         ancestor = 1
         do i = 1, particle_count
            do while (ancestor < particle_count .and. cumulative(ancestor) < target)
               ancestor = ancestor + 1
            end do
            particles(i) = predicted(ancestor)
            next_history(i, 0) = particles(i)
            do j = 1, lag_count
               next_history(i, j) = history(ancestor, j - 1)
            end do
            target = target + 1.0_dp/real(particle_count, dp)
         end do
         history = next_history
         if (time > lag_count) then
            output_time = time - lag_count
            summary_values = history(:, lag_count)
            call summarize_particles(summary_values, lower_bound, upper_bound, &
               out%quantile(output_time, :), out%mean(output_time))
         end if
      end do
      do j = 0, min(lag_count - 1, n - 1)
         output_time = n - j
         summary_values = history(:, j)
         call summarize_particles(summary_values, lower_bound, upper_bound, &
            out%quantile(output_time, :), out%mean(output_time))
      end do
   end subroutine run_particle_filter

   subroutine summarize_particles(values, lower_bound, upper_bound, quantiles, mean_value)
      !! Return the seven TSSS smoothing quantiles and the particle mean.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp), intent(in) :: lower_bound !! Lower bound for bound.
      real(dp), intent(in) :: upper_bound !! Upper bound for bound.
      real(dp), intent(out) :: quantiles(7) !! Quantiles.
      real(dp), intent(out) :: mean_value !! Mean value.
      real(dp), allocatable :: ordered(:)
      real(dp), parameter :: probability(7) = [0.0013_dp, 0.0227_dp, 0.1587_dp, &
         0.5_dp, 0.8413_dp, 0.9773_dp, 0.9987_dp]
      real(dp) :: position, fraction
      integer :: i, lower

      ordered = min(upper_bound, max(lower_bound, values))
      call sort(ordered)
      mean_value = sum(values)/real(size(values), dp)
      do i = 1, 7
         position = 1.0_dp + probability(i)*real(size(values) - 1, dp)
         lower = max(1, min(size(values), floor(position)))
         if (lower == size(values)) then
            quantiles(i) = ordered(lower)
         else
            fraction = position - real(lower, dp)
            quantiles(i) = (1.0_dp - fraction)*ordered(lower) + fraction*ordered(lower + 1)
         end if
      end do
   end subroutine summarize_particles

   function tsss_tvvar(observations, trend_order, initial_system_variance, delta) result(out)
      !! Estimate a smoothly changing variance from paired squared observations.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      integer, intent(in) :: trend_order !! Trend order.
      real(dp), intent(in), optional :: initial_system_variance !! Initial system variance.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      type(tsss_tvvar_t) :: out
      real(dp), allocatable :: transformed(:), state(:, :)
      real(dp) :: candidate, step, best, likelihood, obs_variance
      integer :: i, j, n, search_count

      n = size(observations)/2
      if (n < max(3, trend_order + 1) .or. trend_order < 1 .or. trend_order > 3) then
         out%info = 1
         return
      end if
      allocate(transformed(n))
      do i = 1, n
         transformed(i) = log(max(0.5_dp*(observations(2*i - 1)**2 + observations(2*i)**2), &
            tiny(1.0_dp)))
      end do
      best = -huge(1.0_dp)
      search_count = 19
      do j = 1, search_count
         if (present(initial_system_variance) .and. present(delta)) then
            step = delta
            candidate = initial_system_variance + real(j - 9, dp)*step
         else
            candidate = 2.0_dp**(-j - merge(0, 5, trend_order == 1))
         end if
         if (candidate < 0.0_dp) cycle
         call smooth_local_level(transformed, trend_order, candidate, acos(-1.0_dp)**2/6.0_dp, &
            state, likelihood, obs_variance)
         if (likelihood > best) then
            best = likelihood
            out%system_variance = candidate
            out%observation_variance = obs_variance
         end if
      end do
      call smooth_local_level(transformed, trend_order, out%system_variance, &
         acos(-1.0_dp)**2/6.0_dp, state, out%log_likelihood, out%observation_variance)
      allocate(out%log_variance(n), out%variance(n), out%normalized(size(observations)))
      out%log_variance = state(1, :)
      out%variance = exp(out%log_variance + 0.57721_dp)
      do i = 1, size(observations)
         out%normalized(i) = observations(i)/sqrt(out%variance((i + 1)/2))
      end do
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(trend_order + 2, dp)
   end function tsss_tvvar

   function tsss_tvar(observations, ar_order, trend_order, span, &
      initial_system_variance, delta) result(out)
      !! Fit the TSSS locally observed, smoothly time-varying AR model.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: trend_order !! Trend order.
      integer, intent(in) :: span !! Span.
      real(dp), intent(in), optional :: initial_system_variance !! Initial system variance.
      real(dp), intent(in), optional :: delta !! Model increment or differencing parameter.
      type(tsss_tvar_t) :: out
      real(dp), allocatable :: response(:), regressors(:, :), state(:, :)
      real(dp) :: candidate, best, likelihood, obs_variance
      integer :: blocks, i, j, k

      blocks = size(observations)/span
      if (ar_order < 1 .or. trend_order < 1 .or. trend_order > 2 .or. span < 1 .or. &
         blocks < 2 .or. size(observations) <= ar_order) then
         out%info = 1
         return
      end if
      allocate(response(size(observations) - ar_order))
      allocate(regressors(ar_order, size(observations) - ar_order))
      do i = 1, size(response)
         k = i + ar_order
         response(i) = observations(k)
         do j = 1, ar_order
            regressors(j, i) = observations(k - j)
         end do
      end do
      best = -huge(1.0_dp)
      do j = 1, merge(19, 9, present(initial_system_variance) .and. present(delta))
         if (present(initial_system_variance) .and. present(delta)) then
            candidate = initial_system_variance + real(j - 9, dp)*delta
         else
            candidate = 10.0_dp**(-j - merge(0, 1, trend_order == 1))
         end if
         if (candidate < 0.0_dp) cycle
         call smooth_dynamic_ar(response, regressors, trend_order, span, candidate, state, &
            likelihood, obs_variance)
         if (likelihood > best) then
            best = likelihood
            out%system_variance = candidate
            out%observation_variance = obs_variance
         end if
      end do
      call smooth_dynamic_ar(response, regressors, trend_order, span, out%system_variance, state, &
         out%log_likelihood, out%observation_variance)
      allocate(out%ar(ar_order, blocks), out%parcor(ar_order, blocks))
      do i = 1, blocks
         k = min(size(state, 2), max(1, i*span - ar_order))
         do j = 1, ar_order
            out%ar(j, i) = state((j - 1)*trend_order + 1, k)
         end do
         out%parcor(:, i) = ar_to_parcor(out%ar(:, i))
      end do
      out%span = span
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(ar_order + 2, dp)
   end function tsss_tvar

   pure function tsss_tvspc(ar, observation_variance, span, frequency_count, variance) result(out)
      !! Compute the evolutionary AR spectrum, with optional variance correction.
      real(dp), intent(in) :: ar(:, :) !! Autoregressive coefficients.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      integer, intent(in) :: span !! Span.
      integer, intent(in) :: frequency_count !! Number of frequency.
      real(dp), intent(in), optional :: variance(:) !! Variance value or matrix.
      type(tsss_tvspc_t) :: out
      real(dp) :: angle, denominator, scale
      complex(dp) :: polynomial
      integer :: i, j, k, location

      if (size(ar, 1) < 1 .or. size(ar, 2) < 1 .or. observation_variance <= 0.0_dp .or. &
         span < 1 .or. frequency_count < 1) then
         out%info = 1
         return
      end if
      if (present(variance)) then
         if (size(variance) < size(ar, 2)*span) then
            out%info = 1
            return
         end if
      end if
      allocate(out%frequency(0:frequency_count), &
         out%log_spectrum(0:frequency_count, size(ar, 2)))
      do i = 0, frequency_count
         out%frequency(i) = 0.5_dp*real(i, dp)/real(frequency_count, dp)
         angle = 2.0_dp*acos(-1.0_dp)*out%frequency(i)
         do j = 1, size(ar, 2)
            polynomial = cmplx(1.0_dp, 0.0_dp, dp)
            do k = 1, size(ar, 1)
               polynomial = polynomial - ar(k, j)*exp(cmplx(0.0_dp, -angle*real(k, dp), dp))
            end do
            denominator = max(abs(polynomial)**2, tiny(1.0_dp))
            scale = 1.0_dp
            if (present(variance)) then
               location = max(1, min(size(variance), j*span - span/2))
               scale = max(variance(location), tiny(1.0_dp))
            end if
            out%log_spectrum(i, j) = log10(observation_variance/denominator) + log10(scale)
         end do
      end do
   end function tsss_tvspc

   pure function tsss_pdfunc(model, parameters, lower, upper, count) result(out)
      !! Evaluate one of the TSSS probability densities on a regular grid.
      integer, intent(in) :: model !! Model specification.
      integer, intent(in) :: count !! Count.
      real(dp), intent(in) :: parameters(3) !! Model parameter values.
      real(dp), intent(in) :: lower !! Lower.
      real(dp), intent(in) :: upper !! Upper.
      type(tsss_density_t) :: out
      integer :: i

      if (.not. valid_density_parameters(model, parameters) .or. count < 2 .or. lower >= upper) then
         out%info = 1
         return
      end if
      allocate(out%grid(count), out%density(count))
      do i = 1, count
         out%grid(i) = lower + (upper - lower)*real(i - 1, dp)/real(count - 1, dp)
         out%density(i) = tsss_density_value(out%grid(i), model, parameters)
      end do
      out%model = model
   end function tsss_pdfunc

   pure real(dp) function tsss_density_value(value, model, parameters) result(density)
      !! Return a scalar TSSS density value for the supplied model code.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: parameters(3) !! Model parameter values.
      integer, intent(in) :: model !! Model specification.
      real(dp) :: centered, dispersion, shape

      density = 0.0_dp
      if (.not. valid_density_parameters(model, parameters)) return
      select case (model)
      case (tsss_density_two_sided_exponential)
         density = 0.5_dp*parameters(1)*exp(-parameters(1)*abs(value))
      case (tsss_density_gaussian)
         density = exp(-0.5_dp*(value - parameters(1))**2/parameters(2))/ &
            sqrt(2.0_dp*acos(-1.0_dp)*parameters(2))
      case (tsss_density_cauchy)
         centered = value - parameters(1)
         density = sqrt(parameters(2))/(acos(-1.0_dp)*(parameters(2) + centered**2))
      case (tsss_density_pearson)
         centered = value - parameters(1)
         dispersion = parameters(2)
         shape = parameters(3)
         density = exp(log_gamma(shape) - log_gamma(shape - 0.5_dp))* &
            dispersion**(shape - 0.5_dp)/(sqrt(acos(-1.0_dp))* &
            (centered**2 + dispersion)**shape)
      case (tsss_density_exponential)
         if (value >= 0.0_dp) density = parameters(1)*exp(-parameters(1)*value)
      case (tsss_density_chi_square)
         if (value > 0.0_dp) then
            shape = 0.5_dp*parameters(1)
            density = exp(-0.5_dp*value + (shape - 1.0_dp)*log(0.5_dp*value) - &
               log(2.0_dp) - log_gamma(shape))
         end if
      case (tsss_density_log_exponential)
         centered = value - parameters(1)
         density = exp(centered - exp(centered))
      case (tsss_density_uniform)
         if (value > parameters(1) .and. value <= parameters(2)) then
            density = 1.0_dp/(parameters(2) - parameters(1))
         end if
      end select
   end function tsss_density_value

   pure function tsss_klinfo(reference_model, reference_parameters, model, parameters, &
      lower, upper) result(out)
      !! Integrate truncated KL information using the four TSSS grid refinements.
      integer, intent(in) :: reference_model !! Reference model.
      integer, intent(in) :: model !! Model specification.
      real(dp), intent(in) :: reference_parameters(3) !! Reference parameters.
      real(dp), intent(in) :: parameters(3) !! Model parameter values.
      real(dp), intent(in) :: lower !! Lower.
      real(dp), intent(in) :: upper !! Upper.
      type(tsss_kl_t) :: out
      real(dp) :: reference_density, model_density, weight, value
      integer :: level, i

      if (.not. valid_density_parameters(reference_model, reference_parameters) .or. &
         .not. valid_density_parameters(model, parameters) .or. lower >= upper) then
         out%info = 1
         return
      end if
      do level = 1, 4
         out%intervals(level) = max(1, int((upper - lower + 1.0e-5_dp)*2.0_dp**(level - 1)))
         out%spacing(level) = (upper - lower)/real(out%intervals(level), dp)
         do i = 0, out%intervals(level)
            value = lower + real(i, dp)*out%spacing(level)
            reference_density = tsss_density_value(value, reference_model, reference_parameters)
            model_density = tsss_density_value(value, model, parameters)
            weight = 1.0_dp
            if (i == 0 .or. i == out%intervals(level)) weight = 0.5_dp
            out%reference_mass(level) = out%reference_mass(level) + weight*reference_density
            if (reference_density > 0.0_dp) then
               if (model_density <= 0.0_dp) then
                  out%info = 2
                  out%information(level:4) = huge(1.0_dp)
                  return
               end if
               out%information(level) = out%information(level) + weight*reference_density* &
                  log(reference_density/model_density)
            end if
         end do
         out%reference_mass(level) = out%reference_mass(level)*out%spacing(level)
         out%information(level) = out%information(level)*out%spacing(level)
      end do
   end function tsss_klinfo

   pure function tsss_boxcox(observations) result(out)
      !! Select the TSSS Box-Cox transform on the fixed lambda grid from one to minus one.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      type(tsss_boxcox_t) :: out
      real(dp), allocatable :: transformed(:)
      real(dp) :: jacobian, best_aic
      integer :: i, best

      if (size(observations) < 2 .or. any(observations <= 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(transformed(size(observations)), out%transformed(size(observations)))
      best_aic = huge(1.0_dp)
      best = 1
      do i = 1, 21
         out%lambda(i) = 1.1_dp - 0.1_dp*real(i, dp)
         transformed = box_cox(observations, out%lambda(i))
         out%mean(i) = sum(transformed)/real(size(transformed), dp)
         out%variance(i) = sum((transformed - out%mean(i))**2)/real(size(transformed), dp)
         if (out%variance(i) <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
         out%log_likelihood(i) = -0.5_dp*real(size(transformed), dp)* &
            (log(2.0_dp*acos(-1.0_dp)*out%variance(i)) + 1.0_dp)
         out%aic(i) = -2.0_dp*out%log_likelihood(i) + 4.0_dp
         jacobian = (out%lambda(i) - 1.0_dp)*sum(log(observations))
         out%adjusted_log_likelihood(i) = out%log_likelihood(i) + jacobian
         out%adjusted_aic(i) = out%aic(i) - 2.0_dp*jacobian
         if (out%adjusted_aic(i) <= best_aic) then
            best_aic = out%adjusted_aic(i)
            best = i
            out%transformed = transformed
         end if
      end do
      out%best_lambda = out%lambda(best)
   end function tsss_boxcox

   function tsss_ngsmth(observations, system_noise, system_dispersion, system_shape, &
      observation_noise, observation_dispersion, observation_shape, intervals) result(out)
      !! Smooth a non-Gaussian random walk by TSSS grid-density recursions.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: system_dispersion !! System dispersion.
      real(dp), intent(in) :: system_shape !! System shape.
      real(dp), intent(in) :: observation_dispersion !! Observation dispersion.
      real(dp), intent(in) :: observation_shape !! Observation shape.
      integer, intent(in) :: system_noise !! System noise.
      integer, intent(in) :: observation_noise !! Observation noise.
      integer, intent(in) :: intervals !! Intervals.
      type(tsss_ngsmooth_t) :: out
      real(dp), allocatable :: filtered(:, :), predicted(:, :), kernel(:), work(:), ratio(:)
      real(dp) :: dx, center, spread, total
      integer :: i, j, k, n

      n = size(observations)
      if (n < 1 .or. intervals < 10 .or. system_dispersion <= 0.0_dp .or. &
         observation_dispersion <= 0.0_dp .or. system_noise < 1 .or. system_noise > 3 .or. &
         observation_noise < 1 .or. observation_noise > 3 .or. &
         (system_noise == tsss_noise_pearson .and. system_shape <= 0.5_dp) .or. &
         (observation_noise == tsss_noise_pearson .and. observation_shape <= 0.5_dp)) then
         out%info = 1
         return
      end if
      center = sum(observations)/real(n, dp)
      spread = sqrt(max(sum((observations - center)**2)/real(n, dp), tiny(1.0_dp)))
      allocate(out%grid(0:intervals), out%density(0:intervals, n), out%quantile(n, 7))
      allocate(filtered(0:intervals, n), predicted(0:intervals, n))
      allocate(kernel(-intervals:intervals), work(0:intervals), ratio(0:intervals))
      dx = (maxval(observations) - minval(observations) + 8.0_dp*spread)/real(intervals, dp)
      do i = 0, intervals
         out%grid(i) = minval(observations) - 4.0_dp*spread + real(i, dp)*dx
         work(i) = noise_density(out%grid(i) - center, tsss_noise_gaussian, spread**2, 1.0_dp)
      end do
      call normalize_density(work, dx, total)
      do i = -intervals, intervals
         kernel(i) = noise_density(real(i, dp)*dx, system_noise, system_dispersion, system_shape)
      end do
      do k = 1, n
         do i = 0, intervals
            predicted(i, k) = 0.0_dp
            do j = 0, intervals
               predicted(i, k) = predicted(i, k) + work(j)*kernel(i - j)*dx
            end do
            filtered(i, k) = predicted(i, k)*noise_density(observations(k) - out%grid(i), &
               observation_noise, observation_dispersion, observation_shape)
         end do
         call normalize_density(filtered(:, k), dx, total)
         if (total <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
         out%log_likelihood = out%log_likelihood + log(total)
         work = filtered(:, k)
      end do
      out%density(:, n) = filtered(:, n)
      do k = n - 1, 1, -1
         ratio = out%density(:, k + 1)/max(predicted(:, k + 1), tiny(1.0_dp))
         do i = 0, intervals
            work(i) = 0.0_dp
            do j = 0, intervals
               work(i) = work(i) + kernel(j - i)*ratio(j)*dx
            end do
            work(i) = filtered(i, k)*work(i)
         end do
         call normalize_density(work, dx, total)
         out%density(:, k) = work
      end do
      do k = 1, n
         call density_quantiles(out%grid, out%density(:, k), dx, out%quantile(k, :))
      end do
   end function tsss_ngsmth

   pure function tsss_structural_model(trend, seasonal_order, seasonal, ar_coefficients, &
      ar_initial, trend_variance, seasonal_variance, ar_variance, &
      observation_variance) result(out)
      !! Build TSSS trend, zero-sum seasonal, and AR component matrices.
      real(dp), intent(in), optional :: trend(:) !! Trend.
      real(dp), intent(in), optional :: seasonal(:) !! Seasonal.
      real(dp), intent(in), optional :: ar_coefficients(:) !! Autoregressive coefficients.
      real(dp), intent(in), optional :: ar_initial(:) !! Autoregressive initial.
      integer, intent(in) :: seasonal_order !! Seasonal order.
      real(dp), intent(in) :: trend_variance !! Trend variance.
      real(dp), intent(in) :: seasonal_variance !! Seasonal variance.
      real(dp), intent(in) :: ar_variance !! Autoregressive variance.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      type(tsss_structural_model_t) :: out
      integer :: ar_order, channel, dimension, i, offset, period, trend_order

      trend_order = 0
      if (present(trend)) trend_order = size(trend)
      ar_order = 0
      if (present(ar_coefficients)) ar_order = size(ar_coefficients)
      period = 0
      if (present(seasonal)) period = size(seasonal) + 1
      if (trend_order > 2 .or. seasonal_order < 0 .or. seasonal_order > 2 .or. &
         trend_variance < 0.0_dp .or. seasonal_variance < 0.0_dp .or. &
         ar_variance < 0.0_dp .or. observation_variance < 0.0_dp) then
         out%info = 1
         return
      end if
      if (seasonal_order > 0 .and. .not. present(seasonal)) then
         out%info = 1
         return
      end if
      if (seasonal_order > 0 .and. period < 2) then
         out%info = 1
         return
      end if
      if (seasonal_order == 0 .and. present(seasonal)) then
         out%info = 1
         return
      end if
      if (ar_order > 0 .and. .not. present(ar_initial)) then
         out%info = 1
         return
      end if
      if (present(ar_initial)) then
         if (size(ar_initial) /= ar_order) then
            out%info = 1
            return
         end if
      end if
      dimension = trend_order + seasonal_order*max(0, period - 1) + ar_order
      channel = merge(1, 0, trend_order > 0) + merge(1, 0, seasonal_order > 0) + &
         merge(1, 0, ar_order > 0)
      if (dimension < 1 .or. channel < 1) then
         out%info = 1
         return
      end if
      allocate(out%transition(dimension, dimension), out%system_loading(dimension, channel))
      allocate(out%observation_loading(dimension), out%system_covariance(channel, channel))
      allocate(out%initial_state(dimension))
      out%transition = 0.0_dp
      out%system_loading = 0.0_dp
      out%observation_loading = 0.0_dp
      out%system_covariance = 0.0_dp
      out%initial_state = 0.0_dp
      channel = 0
      offset = 0
      if (trend_order > 0) then
         channel = channel + 1
         if (trend_order == 1) then
            out%transition(1, 1) = 1.0_dp
         else
            out%transition(1, 1:2) = [2.0_dp, -1.0_dp]
            out%transition(2, 1) = 1.0_dp
         end if
         out%system_loading(1, channel) = 1.0_dp
         out%observation_loading(1) = 1.0_dp
         out%system_covariance(channel, channel) = trend_variance
         out%initial_state(1:trend_order) = trend
         offset = trend_order
      end if
      if (seasonal_order > 0) then
         channel = channel + 1
         out%transition(offset + 1, offset + 1:offset + period - 1) = -1.0_dp
         do i = 2, period - 1
            out%transition(offset + i, offset + i - 1) = 1.0_dp
         end do
         out%system_loading(offset + 1, channel) = 1.0_dp
         out%observation_loading(offset + 1) = 1.0_dp
         out%system_covariance(channel, channel) = seasonal_variance
         out%initial_state(offset + 1:offset + period - 1) = seasonal
         offset = offset + seasonal_order*(period - 1)
      end if
      if (ar_order > 0) then
         channel = channel + 1
         out%transition(offset + 1, offset + 1:offset + ar_order) = ar_coefficients
         do i = 2, ar_order
            out%transition(offset + i, offset + i - 1) = 1.0_dp
         end do
         out%system_loading(offset + 1, channel) = 1.0_dp
         out%observation_loading(offset + 1) = 1.0_dp
         out%system_covariance(channel, channel) = ar_variance
         out%initial_state(offset + 1:offset + ar_order) = ar_initial
      end if
      out%observation_variance = observation_variance
      out%trend_order = trend_order
      out%seasonal_order = seasonal_order
      out%period = period
      out%ar_order = ar_order
   end function tsss_structural_model

   pure function tsss_tsmooth(observations, transition, system_loading, observation_loading, &
      system_covariance, observation_variance, initial_state, initial_covariance, filter_end, &
      predict_end, observation_min, observation_max, missing_start, missing_count) result(out)
      !! Run the TSSS general Gaussian filter, smoother, interpolation, and prediction interface.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: system_loading(:, :) !! System loading.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: system_covariance(:, :) !! System covariance matrix.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in), optional :: filter_end !! Filter end.
      integer, intent(in), optional :: predict_end !! Predict end.
      integer, intent(in), optional :: missing_start(:) !! Missing start.
      integer, intent(in), optional :: missing_count(:) !! Number of missing.
      real(dp), intent(in), optional :: observation_min !! Observation min.
      real(dp), intent(in), optional :: observation_max !! Observation max.
      type(tsss_tsmooth_t) :: out
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      real(dp) :: lower, upper, observation_mean
      integer :: dimension, end_filter, end_prediction, i, j, last_missing

      dimension = size(initial_state)
      end_filter = size(observations)
      if (present(filter_end)) end_filter = filter_end
      end_prediction = size(observations)
      if (present(predict_end)) end_prediction = predict_end
      lower = -huge(1.0_dp)
      upper = huge(1.0_dp)
      if (present(observation_min)) lower = observation_min
      if (present(observation_max)) upper = observation_max
      if (size(observations) < 1 .or. dimension < 1 .or. end_filter < 1 .or. &
         end_filter > size(observations) .or. end_prediction < end_filter .or. &
         size(transition, 1) /= dimension .or. size(transition, 2) /= dimension .or. &
         size(system_loading, 1) /= dimension .or. &
         size(system_loading, 2) /= size(system_covariance, 1) .or. &
         size(system_covariance, 1) /= size(system_covariance, 2) .or. &
         size(observation_loading) /= dimension .or. &
         any(shape(initial_covariance) /= [dimension, dimension]) .or. &
         observation_variance < 0.0_dp .or. lower >= upper) then
         out%info = 1
         return
      end if
      if (present(missing_start) .neqv. present(missing_count)) then
         out%info = 1
         return
      end if
      if (present(missing_start)) then
         if (size(missing_start) /= size(missing_count) .or. any(missing_start < 1) .or. &
            any(missing_count < 0)) then
            out%info = 1
            return
         end if
         do i = 1, size(missing_start)
            last_missing = missing_start(i) + missing_count(i) - 1
            if (last_missing > size(observations)) then
               out%info = 1
               return
            end if
         end do
      end if
      observation_mean = sum(observations)/real(size(observations), dp)
      allocate(model%y(end_prediction, 1), model%z(1, dimension, 1), model%h(1, 1, 1))
      allocate(model%transition(dimension, dimension, 1))
      allocate(model%r(dimension, size(system_covariance, 1), 1))
      allocate(model%q(size(system_covariance, 1), size(system_covariance, 2), 1))
      allocate(model%a1(dimension), model%p1(dimension, dimension), model%p1inf(dimension, dimension))
      allocate(model%missing(end_prediction, 1))
      model%y = 0.0_dp
      model%y(1:size(observations), 1) = observations - observation_mean
      model%z(1, :, 1) = observation_loading
      model%h(1, 1, 1) = observation_variance
      model%transition(:, :, 1) = transition
      model%r(:, :, 1) = system_loading
      model%q(:, :, 1) = system_covariance
      model%a1 = initial_state
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing = .true.
      do i = 1, end_filter
         model%missing(i, 1) = observations(i) <= lower .or. observations(i) >= upper
      end do
      if (present(missing_start)) then
         do i = 1, size(missing_start)
            do j = missing_start(i), missing_start(i) + missing_count(i) - 1
               model%missing(j, 1) = .true.
            end do
         end do
      end if
      filtered = kfs_filter(model)
      if (filtered%info /= 0) then
         out%info = filtered%info
         return
      end if
      smoothed = kfs_smooth(model, filtered)
      if (smoothed%info /= 0) then
         out%info = smoothed%info
         return
      end if
      out%predicted_state = filtered%a_pred + observation_mean
      out%filtered_state = filtered%a_filt + observation_mean
      out%smoothed_state = smoothed%state + observation_mean
      out%predicted_covariance = filtered%p_pred
      out%filtered_covariance = filtered%p_filt
      out%smoothed_covariance = smoothed%covariance
      allocate(out%smoothed_variance(dimension, end_prediction), out%prediction_error(end_prediction))
      do i = 1, end_prediction
         do j = 1, dimension
            out%smoothed_variance(j, i) = smoothed%covariance(j, j, i)
         end do
         out%prediction_error(i) = 0.0_dp
         if (i <= size(observations)) then
            if (model%missing(i, 1)) then
               out%prediction_error(i) = observations(i) - &
                  dot_product(observation_loading, out%smoothed_state(:, i))
            end if
         end if
      end do
      out%log_likelihood = filtered%log_likelihood
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(dimension + 1, dp)
      out%observations = filtered%observations
   end function tsss_tsmooth

   function tsss_simssm(transition, system_loading, observation_loading, system_covariance, &
      observation_variance, initial_state, count) result(out)
      !! Simulate a Gaussian linear state-space model.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: system_loading(:, :) !! System loading.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: system_covariance(:, :) !! System covariance matrix.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      integer, intent(in) :: count !! Count.
      type(tsss_simulation_t) :: out

      out = simulate_state_space(transition, system_loading, observation_loading, &
         system_covariance, observation_variance, initial_state, count, &
         tsss_noise_gaussian, 1.0_dp, tsss_noise_gaussian, 1.0_dp)
   end function tsss_simssm

   function tsss_ngsim(transition, system_loading, observation_loading, system_covariance, &
      observation_variance, initial_state, count, system_noise, system_shape, &
      observation_noise, observation_shape) result(out)
      !! Simulate a linear state-space model with non-Gaussian innovations.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: system_loading(:, :) !! System loading.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: system_covariance(:, :) !! System covariance matrix.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      integer, intent(in) :: count !! Count.
      integer, intent(in) :: system_noise !! System noise.
      integer, intent(in) :: observation_noise !! Observation noise.
      real(dp), intent(in) :: system_shape !! System shape.
      real(dp), intent(in) :: observation_shape !! Observation shape.
      type(tsss_simulation_t) :: out

      out = simulate_state_space(transition, system_loading, observation_loading, &
         system_covariance, observation_variance, initial_state, count, system_noise, &
         system_shape, observation_noise, observation_shape)
   end function tsss_ngsim

   function simulate_state_space(transition, system_loading, observation_loading, &
      system_covariance, observation_variance, initial_state, count, system_noise, &
      system_shape, observation_noise, observation_shape) result(out)
      !! Implement shared Gaussian and non-Gaussian state simulation.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: system_loading(:, :) !! System loading.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: system_covariance(:, :) !! System covariance matrix.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      integer, intent(in) :: count !! Count.
      integer, intent(in) :: system_noise !! System noise.
      integer, intent(in) :: observation_noise !! Observation noise.
      real(dp), intent(in) :: system_shape !! System shape.
      real(dp), intent(in) :: observation_shape !! Observation shape.
      type(tsss_simulation_t) :: out
      real(dp), allocatable :: root(:, :), innovation(:), state(:)
      integer :: i, info, j

      if (count < 1 .or. size(transition, 1) /= size(transition, 2) .or. &
         size(initial_state) /= size(transition, 1) .or. &
         size(system_loading, 1) /= size(transition, 1) .or. &
         size(system_loading, 2) /= size(system_covariance, 1) .or. &
         size(system_covariance, 1) /= size(system_covariance, 2) .or. &
         size(observation_loading) /= size(initial_state) .or. observation_variance < 0.0_dp .or. &
         system_noise < 1 .or. system_noise > 3 .or. observation_noise < 1 .or. &
         observation_noise > 3) then
         out%info = 1
         return
      end if
      call cholesky_lower_semidefinite(system_covariance, root, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      allocate(out%state(size(initial_state), count), out%observation(count))
      allocate(innovation(size(system_covariance, 1)), state(size(initial_state)))
      state = initial_state
      do i = 1, count
         do j = 1, size(innovation)
            innovation(j) = noise_draw(system_noise, system_shape)
         end do
         state = matmul(transition, state) + matmul(system_loading, matmul(root, innovation))
         out%state(:, i) = state
         out%observation(i) = dot_product(observation_loading, state) + &
            sqrt(observation_variance)*noise_draw(observation_noise, observation_shape)
      end do
   end function simulate_state_space

   subroutine smooth_local_level(values, order, system_variance, fixed_variance, &
      state, log_likelihood, observation_variance)
      !! Smooth an integrated trend model with a scalar observation.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp), intent(in) :: system_variance !! System variance.
      real(dp), intent(in) :: fixed_variance !! Fixed variance.
      integer, intent(in) :: order !! Model or polynomial order.
      real(dp), allocatable, intent(out) :: state(:, :) !! State vector or state sequence.
      real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
      real(dp), intent(out) :: observation_variance !! Observation-error variance.
      real(dp), allocatable :: design(:, :)

      allocate(design(order, size(values)))
      design = 0.0_dp
      design(1, :) = 1.0_dp
      call smooth_dynamic_regression(values, design, order, system_variance, fixed_variance, &
         .false., state, log_likelihood, observation_variance)
   end subroutine smooth_local_level

   subroutine smooth_dynamic_ar(values, regressors, trend_order, span, system_variance, &
      state, log_likelihood, observation_variance)
      !! Smooth independent integrated trends for dynamic AR coefficients.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: system_variance !! System variance.
      integer, intent(in) :: trend_order !! Trend order.
      integer, intent(in) :: span !! Span.
      real(dp), allocatable, intent(out) :: state(:, :) !! State vector or state sequence.
      real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
      real(dp), intent(out) :: observation_variance !! Observation-error variance.
      real(dp), allocatable :: design(:, :)
      integer :: i, j

      allocate(design(size(regressors, 1)*trend_order, size(values)))
      design = 0.0_dp
      do i = 1, size(values)
         do j = 1, size(regressors, 1)
            design((j - 1)*trend_order + 1, i) = regressors(j, i)
         end do
      end do
      call smooth_dynamic_regression(values, design, trend_order, system_variance, 1.0_dp, &
         .true., state, log_likelihood, observation_variance, span)
   end subroutine smooth_dynamic_ar

   subroutine smooth_dynamic_regression(values, design, component_order, system_variance, &
      fixed_variance, estimate_variance, state, log_likelihood, observation_variance, transition_span)
      !! Apply a scalar-observation Kalman filter and fixed-interval smoother.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp), intent(in) :: design(:, :) !! Design.
      real(dp), intent(in) :: system_variance !! System variance.
      real(dp), intent(in) :: fixed_variance !! Fixed variance.
      integer, intent(in) :: component_order !! Component order.
      logical, intent(in) :: estimate_variance !! Whether to estimate the variance.
      integer, intent(in), optional :: transition_span !! Transition span.
      real(dp), allocatable, intent(out) :: state(:, :) !! State vector or state sequence.
      real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
      real(dp), intent(out) :: observation_variance !! Observation-error variance.
      real(dp), allocatable :: transition(:, :), filtered(:, :), predicted(:, :)
      real(dp), allocatable :: vf(:, :, :), vp(:, :, :), covariance(:, :), prediction_cov(:, :)
      real(dp), allocatable :: gain(:), projected(:), rhs(:, :), smoother_gain(:, :)
      real(dp) :: innovation, innovation_variance, sum_squares, log_determinant, pi
      integer :: dimension, components, evolve_every, i, j, time, info

      dimension = size(design, 1)
      components = dimension/component_order
      evolve_every = 1
      if (present(transition_span)) evolve_every = transition_span
      allocate(transition(dimension, dimension), filtered(dimension, size(values)))
      allocate(predicted(dimension, size(values)), vf(dimension, dimension, size(values)))
      allocate(vp(dimension, dimension, size(values)), covariance(dimension, dimension))
      allocate(prediction_cov(dimension, dimension), gain(dimension), projected(dimension))
      transition = 0.0_dp
      do i = 1, components
         do j = 1, component_order
            transition((i - 1)*component_order + 1, (i - 1)*component_order + j) = &
               (-1.0_dp)**real(j + 1, dp)*real(binomial(component_order, j), dp)
         end do
         do j = 2, component_order
            transition((i - 1)*component_order + j, (i - 1)*component_order + j - 1) = 1.0_dp
         end do
      end do
      filtered(:, 1) = 0.0_dp
      covariance = 0.0_dp
      do i = 1, dimension
         covariance(i, i) = 1.0e4_dp
      end do
      log_likelihood = 0.0_dp
      sum_squares = 0.0_dp
      log_determinant = 0.0_dp
      pi = acos(-1.0_dp)
      do time = 1, size(values)
         if (time == 1) then
            predicted(:, time) = 0.0_dp
            prediction_cov = covariance
         else if (mod(time - 1, evolve_every) == 0) then
            predicted(:, time) = matmul(transition, filtered(:, time - 1))
            prediction_cov = matmul(transition, matmul(covariance, transpose(transition)))
         else
            predicted(:, time) = filtered(:, time - 1)
            prediction_cov = covariance
         end if
         if (time == 1 .or. mod(time - 1, evolve_every) == 0) then
            do i = 1, components
               j = (i - 1)*component_order + 1
               prediction_cov(j, j) = prediction_cov(j, j) + system_variance
            end do
         end if
         projected = matmul(prediction_cov, design(:, time))
         innovation_variance = dot_product(design(:, time), projected) + fixed_variance
         innovation = values(time) - dot_product(design(:, time), predicted(:, time))
         gain = projected/innovation_variance
         filtered(:, time) = predicted(:, time) + gain*innovation
         covariance = prediction_cov - outer_product(gain, projected)
         covariance = 0.5_dp*(covariance + transpose(covariance))
         vp(:, :, time) = prediction_cov
         vf(:, :, time) = covariance
         sum_squares = sum_squares + innovation**2/innovation_variance
         log_determinant = log_determinant + log(innovation_variance)
         log_likelihood = log_likelihood - 0.5_dp*(log(2.0_dp*pi*innovation_variance) + &
            innovation**2/innovation_variance)
      end do
      observation_variance = fixed_variance
      if (estimate_variance) then
         observation_variance = max(sum_squares/real(size(values), dp), tiny(1.0_dp))
         log_likelihood = -0.5_dp*(real(size(values), dp)*(log(2.0_dp*pi*observation_variance) + &
            1.0_dp) + log_determinant)
      end if
      allocate(state(dimension, size(values)))
      state(:, size(values)) = filtered(:, size(values))
      allocate(rhs(dimension, dimension), smoother_gain(dimension, dimension))
      do time = size(values) - 1, 1, -1
         if (mod(time, evolve_every) == 0) then
            rhs = matmul(transition, vf(:, :, time))
         else
            rhs = vf(:, :, time)
         end if
         call solve_matrix(transpose(vp(:, :, time + 1)), rhs, smoother_gain, info)
         if (info /= 0) then
            state(:, time) = filtered(:, time)
         else
            smoother_gain = transpose(smoother_gain)
            state(:, time) = filtered(:, time) + matmul(smoother_gain, &
               state(:, time + 1) - predicted(:, time + 1))
         end if
      end do
   end subroutine smooth_dynamic_regression

   pure function ar_to_parcor(ar) result(parcor)
      !! Convert AR coefficients to partial autocorrelations by step-down recursion.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp) :: parcor(size(ar)), current(size(ar)), previous(size(ar)), denominator
      integer :: i, j

      current = ar
      parcor = 0.0_dp
      do i = size(ar), 1, -1
         parcor(i) = current(i)
         denominator = 1.0_dp - current(i)**2
         if (i > 1 .and. abs(denominator) > sqrt(epsilon(1.0_dp))) then
            previous = current
            do j = 1, i - 1
               current(j) = (previous(j) + current(i)*previous(i - j))/denominator
            end do
         end if
      end do
   end function ar_to_parcor

   pure integer function binomial(n, k) result(value)
      !! Return a small integer binomial coefficient.
      integer, intent(in) :: n !! Number of observations or elements.
      integer, intent(in) :: k !! K.
      integer :: i

      value = 1
      do i = 1, k
         value = value*(n - i + 1)/i
      end do
   end function binomial

   pure logical function valid_density_parameters(model, parameters) result(valid)
      !! Check the parameter domain for a supported TSSS density.
      integer, intent(in) :: model !! Model specification.
      real(dp), intent(in) :: parameters(3) !! Model parameter values.

      valid = model >= tsss_density_two_sided_exponential .and. &
         model <= tsss_density_uniform
      if (.not. valid) return
      select case (model)
      case (tsss_density_two_sided_exponential, tsss_density_exponential)
         valid = parameters(1) > 0.0_dp
      case (tsss_density_gaussian, tsss_density_cauchy)
         valid = parameters(2) > 0.0_dp
      case (tsss_density_pearson)
         valid = parameters(2) > 0.0_dp .and. parameters(3) > 0.5_dp
      case (tsss_density_chi_square)
         valid = parameters(1) > 0.0_dp
      case (tsss_density_uniform)
         valid = parameters(1) < parameters(2)
      end select
   end function valid_density_parameters

   pure real(dp) function noise_density(value, noise, dispersion, shape) result(density)
      !! Evaluate a centered Gaussian, Pearson, or Laplace density.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: dispersion !! Dispersion.
      real(dp), intent(in) :: shape !! Shape.
      integer, intent(in) :: noise !! Noise.
      real(dp) :: scale

      select case (noise)
      case (tsss_noise_gaussian)
         density = exp(-0.5_dp*value**2/dispersion)/sqrt(2.0_dp*acos(-1.0_dp)*dispersion)
      case (tsss_noise_pearson)
         scale = sqrt(dispersion)
         density = exp(log_gamma(shape) - log_gamma(shape - 0.5_dp))/ &
            (sqrt(acos(-1.0_dp))*scale)*(1.0_dp + (value/scale)**2)**(-shape)
      case default
         scale = sqrt(0.5_dp*dispersion)
         density = exp(-abs(value)/scale)/(2.0_dp*scale)
      end select
   end function noise_density

   pure subroutine normalize_density(density, spacing, integral)
      !! Normalize a grid density by the trapezoidal integral.
      real(dp), intent(inout) :: density(0:) !! Density, updated in place.
      real(dp), intent(in) :: spacing !! Spacing.
      real(dp), intent(out) :: integral !! Integral.

      integral = spacing*(sum(density) - 0.5_dp*(density(0) + density(ubound(density, 1))))
      if (integral > tiny(1.0_dp)) density = density/integral
   end subroutine normalize_density

   pure subroutine density_quantiles(grid, density, spacing, quantiles)
      !! Interpolate the seven TSSS posterior probability points.
      real(dp), intent(in) :: grid(0:) !! Grid.
      real(dp), intent(in) :: density(0:) !! Density.
      real(dp), intent(in) :: spacing !! Spacing.
      real(dp), intent(out) :: quantiles(7) !! Quantiles.
      real(dp), parameter :: probability(7) = [0.0013_dp, 0.0227_dp, 0.1587_dp, &
         0.5_dp, 0.8413_dp, 0.9773_dp, 0.9987_dp]
      real(dp) :: cumulative, previous
      integer :: i, j

      quantiles = grid(ubound(grid, 1))
      cumulative = 0.0_dp
      j = 1
      do i = 1, ubound(grid, 1)
         previous = cumulative
         cumulative = cumulative + 0.5_dp*spacing*(density(i - 1) + density(i))
         do while (j <= 7)
            if (cumulative < probability(j)) exit
            quantiles(j) = grid(i - 1) + spacing*(probability(j) - previous)/ &
               max(cumulative - previous, tiny(1.0_dp))
            j = j + 1
         end do
      end do
   end subroutine density_quantiles

   real(dp) function noise_draw(noise, shape) result(value)
      !! Draw a standardized Gaussian, Pearson, or Laplace innovation.
      integer, intent(in) :: noise !! Noise.
      real(dp), intent(in) :: shape !! Shape.
      real(dp) :: uniform, gamma_draw

      select case (noise)
      case (tsss_noise_gaussian)
         value = random_standard_normal()
      case (tsss_noise_pearson)
         if (shape <= 0.5_dp) then
            value = 0.0_dp
         else
            gamma_draw = random_gamma(shape - 0.5_dp)
            value = random_standard_normal()/sqrt(max(2.0_dp*gamma_draw, tiny(1.0_dp)))
         end if
      case default
         uniform = max(random_uniform(), tiny(1.0_dp))
         value = -sign(1.0_dp, uniform - 0.5_dp)*log(1.0_dp - 2.0_dp*abs(uniform - 0.5_dp))/sqrt(2.0_dp)
      end select
   end function noise_draw

end module tsss_mod
