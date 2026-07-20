! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R BAYSTAR package.
module baystar_mod
   !! Bayesian two-regime threshold autoregression translated from BAYSTAR.
   use kind_mod, only: dp
   use stats_mod, only: median, quantile, sorted, variance
   use linalg_mod, only: invert_matrix
   use random_mod, only: random_uniform, random_standard_normal, random_gamma, &
      random_multivariate_normal
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   type, public :: baystar_prior_t
      !! Gaussian coefficient precisions and inverse-gamma variance prior.
      real(dp), allocatable :: mean1(:)
      real(dp), allocatable :: precision1(:, :)
      real(dp), allocatable :: mean2(:)
      real(dp), allocatable :: precision2(:, :)
      real(dp) :: variance_degrees = 3.0_dp
      real(dp) :: variance_scale = 0.0_dp
      integer :: info = 0
   end type baystar_prior_t

   type, public :: baystar_coefficient_posterior_t
      !! Conditional Gaussian posterior for one regime's coefficients.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: covariance(:, :)
      integer :: observations = 0
      integer :: info = 0
   end type baystar_coefficient_posterior_t

   type, public :: baystar_variance_posterior_t
      !! Conditional inverse-gamma posterior represented by gamma shape and rate.
      real(dp) :: shape = 0.0_dp
      real(dp) :: rate = 0.0_dp
      integer :: observations = 0
      integer :: info = 0
   end type baystar_variance_posterior_t

   type, public :: baystar_summary_t
      !! Posterior mean, median, standard deviation, and central interval.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: median(:)
      real(dp), allocatable :: standard_deviation(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp) :: level = 0.95_dp
      integer :: info = 0
   end type baystar_summary_t

   type, public :: baystar_fit_t
      !! BAYSTAR MCMC draws, posterior summaries, residuals, and DIC.
      real(dp), allocatable :: coefficient1(:, :)
      real(dp), allocatable :: coefficient2(:, :)
      real(dp), allocatable :: innovation_variance(:, :)
      real(dp), allocatable :: threshold(:)
      real(dp), allocatable :: log_likelihood(:)
      real(dp), allocatable :: regime_mean(:, :)
      integer, allocatable :: delay(:)
      logical, allocatable :: threshold_accepted(:)
      real(dp), allocatable :: residuals(:)
      type(baystar_summary_t) :: summary
      real(dp) :: dic = huge(1.0_dp)
      real(dp) :: threshold_acceptance = 0.0_dp
      integer :: selected_delay = 0
      integer :: burnin = 0
      integer :: info = 0
   end type baystar_fit_t

   type, public :: baystar_simulation_t
      !! Simulated two-regime TAR series and retained innovations.
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: innovations(:)
      integer :: burnin = 0
      integer :: info = 0
   end type baystar_simulation_t

   public :: baystar_default_prior, baystar_log_likelihood
   public :: baystar_coefficient_posterior, baystar_variance_posterior
   public :: baystar_delay_probabilities, baystar_summary
   public :: baystar_simulate_from_innovations, baystar_simulate
   public :: baystar_fit

contains

   pure function baystar_default_prior(series, lags1, lags2, &
      include_constant) result(prior)
      !! Construct BAYSTAR's zero-mean coefficient and variance priors.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: lags1(:) !! Regime-one lag indices.
      integer, intent(in) :: lags2(:) !! Regime-two lag indices.
      logical, intent(in), optional :: include_constant !! Include regime intercepts.
      type(baystar_prior_t) :: prior
      integer :: constant, i

      constant = 1
      if (present(include_constant)) constant = merge(1, 0, include_constant)
      if (size(series) < 3 .or. size(lags1) < 1 .or. size(lags2) < 1 .or. &
         any(lags1 < 1) .or. any(lags2 < 1) .or. &
         .not. all(ieee_is_finite(series))) then
         prior%info = 1
         return
      end if
      allocate(prior%mean1(size(lags1) + constant), &
         prior%precision1(size(lags1) + constant, size(lags1) + constant))
      allocate(prior%mean2(size(lags2) + constant), &
         prior%precision2(size(lags2) + constant, size(lags2) + constant))
      prior%mean1 = 0.0_dp
      prior%mean2 = 0.0_dp
      prior%precision1 = 0.0_dp
      prior%precision2 = 0.0_dp
      do i = 1, size(prior%mean1)
         prior%precision1(i, i) = 0.1_dp
      end do
      do i = 1, size(prior%mean2)
         prior%precision2(i, i) = 0.1_dp
      end do
      prior%variance_scale = variance(series)/3.0_dp
   end function baystar_default_prior

   pure real(dp) function baystar_log_likelihood(series, coefficient1, &
      coefficient2, variance1, variance2, delay, threshold, lags1, lags2, &
      include_constant, threshold_variable) result(log_likelihood)
      !! Evaluate the BAYSTAR conditional Gaussian log likelihood.
      real(dp), intent(in) :: series(:) !! Observed time series.
      real(dp), intent(in) :: coefficient1(:) !! Regime-one coefficients.
      real(dp), intent(in) :: coefficient2(:) !! Regime-two coefficients.
      real(dp), intent(in) :: variance1 !! Regime-one innovation variance.
      real(dp), intent(in) :: variance2 !! Regime-two innovation variance.
      integer, intent(in) :: delay !! Threshold delay lag.
      real(dp), intent(in) :: threshold !! Regime threshold.
      integer, intent(in) :: lags1(:) !! Regime-one lag indices.
      integer, intent(in) :: lags2(:) !! Regime-two lag indices.
      logical, intent(in), optional :: include_constant !! Include regime intercepts.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.
      real(dp) :: prediction, residual
      integer :: start, time, count1, count2
      logical :: constant, low

      log_likelihood = -huge(1.0_dp)
      constant = .true.
      if (present(include_constant)) constant = include_constant
      if (.not. valid_specification(series, coefficient1, coefficient2, &
         variance1, variance2, delay, lags1, lags2, constant, &
         threshold_variable)) return
      start = baystar_start(lags1, lags2, constant)
      log_likelihood = 0.0_dp
      count1 = 0
      count2 = 0
      do time = start, size(series)
         low = threshold_value(series, time, delay, threshold_variable) <= threshold
         if (low) then
            prediction = regime_prediction(series, time, coefficient1, lags1, constant)
            residual = series(time) - prediction
            log_likelihood = log_likelihood - 0.5_dp*residual**2/variance1
            count1 = count1 + 1
         else
            prediction = regime_prediction(series, time, coefficient2, lags2, constant)
            residual = series(time) - prediction
            log_likelihood = log_likelihood - 0.5_dp*residual**2/variance2
            count2 = count2 + 1
         end if
      end do
      log_likelihood = log_likelihood - 0.5_dp*real(count1, dp)*log(variance1) - &
         0.5_dp*real(count2, dp)*log(variance2)
   end function baystar_log_likelihood

   pure function baystar_coefficient_posterior(regime, series, variance_value, &
      delay, threshold, lags1, lags2, prior_mean, prior_precision, &
      include_constant, threshold_variable) result(out)
      !! Compute one regime's conditional Gaussian coefficient posterior.
      integer, intent(in) :: regime !! Regime number, one or two.
      real(dp), intent(in) :: series(:) !! Observed time series.
      real(dp), intent(in) :: variance_value !! Current regime variance.
      integer, intent(in) :: delay !! Threshold delay lag.
      real(dp), intent(in) :: threshold !! Current threshold.
      integer, intent(in) :: lags1(:) !! Regime-one lag indices.
      integer, intent(in) :: lags2(:) !! Regime-two lag indices.
      real(dp), intent(in) :: prior_mean(:) !! Prior coefficient mean.
      real(dp), intent(in) :: prior_precision(:, :) !! Prior coefficient precision.
      logical, intent(in), optional :: include_constant !! Include regime intercepts.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.
      type(baystar_coefficient_posterior_t) :: out
      real(dp), allocatable :: design(:, :), response(:), precision(:, :), inverse(:, :)
      real(dp), allocatable :: right(:)
      integer, allocatable :: selected_lags(:)
      integer :: constant, start, time, row, columns, status
      logical :: low, select_row

      constant = 1
      if (present(include_constant)) constant = merge(1, 0, include_constant)
      if (regime == 1) then
         selected_lags = lags1
      else if (regime == 2) then
         selected_lags = lags2
      else
         out%info = 1
         return
      end if
      columns = size(selected_lags) + constant
      if (variance_value <= 0.0_dp .or. size(prior_mean) /= columns .or. &
         any(shape(prior_precision) /= [columns, columns])) then
         out%info = 1
         return
      end if
      start = baystar_start(lags1, lags2, constant == 1)
      out%observations = 0
      do time = start, size(series)
         low = threshold_value(series, time, delay, threshold_variable) <= threshold
         if ((regime == 1 .and. low) .or. (regime == 2 .and. .not. low)) &
            out%observations = out%observations + 1
      end do
      if (out%observations <= columns) then
         out%info = 2
         return
      end if
      allocate(design(out%observations, columns), response(out%observations))
      row = 0
      do time = start, size(series)
         low = threshold_value(series, time, delay, threshold_variable) <= threshold
         select_row = (regime == 1 .and. low) .or. (regime == 2 .and. .not. low)
         if (.not. select_row) cycle
         row = row + 1
         call fill_regressor(series, time, selected_lags, constant == 1, &
            design(row, :))
         response(row) = series(time)
      end do
      precision = matmul(transpose(design), design)/variance_value + prior_precision
      call invert_matrix(precision, inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      right = matmul(transpose(design), response)/variance_value + &
         matmul(prior_precision, prior_mean)
      out%mean = matmul(inverse, right)
      out%covariance = inverse
   end function baystar_coefficient_posterior

   pure function baystar_variance_posterior(regime, series, coefficient, &
      delay, threshold, lags1, lags2, prior_degrees, prior_scale, &
      include_constant, threshold_variable) result(out)
      !! Compute one regime's conditional inverse-gamma variance posterior.
      integer, intent(in) :: regime !! Regime number, one or two.
      real(dp), intent(in) :: series(:) !! Observed time series.
      real(dp), intent(in) :: coefficient(:) !! Current regime coefficients.
      integer, intent(in) :: delay !! Threshold delay lag.
      real(dp), intent(in) :: threshold !! Current threshold.
      integer, intent(in) :: lags1(:) !! Regime-one lag indices.
      integer, intent(in) :: lags2(:) !! Regime-two lag indices.
      real(dp), intent(in) :: prior_degrees !! Variance-prior degrees parameter.
      real(dp), intent(in) :: prior_scale !! Variance-prior scale parameter.
      logical, intent(in), optional :: include_constant !! Include regime intercepts.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.
      type(baystar_variance_posterior_t) :: out
      integer, allocatable :: selected_lags(:)
      real(dp) :: residual, sum_squares
      integer :: start, time
      logical :: constant, low, select_row

      constant = .true.
      if (present(include_constant)) constant = include_constant
      if (regime == 1) then
         selected_lags = lags1
      else if (regime == 2) then
         selected_lags = lags2
      else
         out%info = 1
         return
      end if
      if (size(coefficient) /= size(selected_lags) + merge(1, 0, constant) .or. &
         prior_degrees <= 0.0_dp .or. prior_scale <= 0.0_dp) then
         out%info = 1
         return
      end if
      start = baystar_start(lags1, lags2, constant)
      sum_squares = 0.0_dp
      do time = start, size(series)
         low = threshold_value(series, time, delay, threshold_variable) <= threshold
         select_row = (regime == 1 .and. low) .or. (regime == 2 .and. .not. low)
         if (.not. select_row) cycle
         residual = series(time) - &
            regime_prediction(series, time, coefficient, selected_lags, constant)
         sum_squares = sum_squares + residual**2
         out%observations = out%observations + 1
      end do
      if (out%observations < 1) then
         out%info = 2
         return
      end if
      out%shape = 0.5_dp*(prior_degrees + real(out%observations, dp))
      out%rate = 0.5_dp*(prior_degrees*prior_scale + sum_squares)
   end function baystar_variance_posterior

   pure function baystar_delay_probabilities(series, coefficient1, coefficient2, &
      variance1, variance2, threshold, lags1, lags2, maximum_delay, &
      include_constant, threshold_variable) result(probability)
      !! Compute BAYSTAR delay-lag probabilities including its decreasing prior.
      real(dp), intent(in) :: series(:) !! Observed time series.
      real(dp), intent(in) :: coefficient1(:) !! Regime-one coefficients.
      real(dp), intent(in) :: coefficient2(:) !! Regime-two coefficients.
      real(dp), intent(in) :: variance1 !! Regime-one innovation variance.
      real(dp), intent(in) :: variance2 !! Regime-two innovation variance.
      real(dp), intent(in) :: threshold !! Current threshold.
      integer, intent(in) :: lags1(:) !! Regime-one lag indices.
      integer, intent(in) :: lags2(:) !! Regime-two lag indices.
      integer, intent(in) :: maximum_delay !! Largest candidate delay.
      logical, intent(in), optional :: include_constant !! Include regime intercepts.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.
      real(dp), allocatable :: probability(:)
      real(dp), allocatable :: log_weight(:)
      integer :: delay

      if (maximum_delay < 1) then
         allocate(probability(0))
         return
      end if
      allocate(probability(maximum_delay), log_weight(maximum_delay))
      do delay = 1, maximum_delay
         log_weight(delay) = baystar_log_likelihood(series, coefficient1, &
            coefficient2, variance1, variance2, delay, threshold, lags1, lags2, &
            include_constant, threshold_variable) + &
            log(real(maximum_delay - delay + 1, dp))
      end do
      probability = exp(log_weight - maxval(log_weight))
      probability = probability/sum(probability)
   end function baystar_delay_probabilities

   pure function baystar_summary(draws, level) result(out)
      !! Summarize posterior draws by moments and central empirical intervals.
      real(dp), intent(in) :: draws(:, :) !! Draw-by-parameter matrix.
      real(dp), intent(in), optional :: level !! Central posterior interval level.
      type(baystar_summary_t) :: out
      real(dp), allocatable :: ordered(:)
      real(dp) :: selected_level
      integer :: parameter, draw_count

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      draw_count = size(draws, 1)
      if (draw_count < 2 .or. size(draws, 2) < 1 .or. selected_level <= 0.0_dp .or. &
         selected_level >= 1.0_dp .or. .not. all(ieee_is_finite(draws))) then
         out%info = 1
         return
      end if
      allocate(out%mean(size(draws, 2)), out%median(size(draws, 2)), &
         out%standard_deviation(size(draws, 2)), out%lower(size(draws, 2)), &
         out%upper(size(draws, 2)))
      do parameter = 1, size(draws, 2)
         ordered = sorted(draws(:, parameter))
         out%mean(parameter) = sum(draws(:, parameter))/real(draw_count, dp)
         out%median(parameter) = quantile(ordered, 0.5_dp)
         out%standard_deviation(parameter) = sqrt(sum((draws(:, parameter) - &
            out%mean(parameter))**2)/real(draw_count - 1, dp))
         out%lower(parameter) = quantile(ordered, 0.5_dp*(1.0_dp - selected_level))
         out%upper(parameter) = quantile(ordered, 0.5_dp*(1.0_dp + selected_level))
      end do
      out%level = selected_level
   end function baystar_summary

   pure function baystar_simulate_from_innovations(observations, coefficient1, &
      coefficient2, variance1, variance2, delay, threshold, lags1, lags2, &
      innovations, burnin) result(out)
      !! Simulate a two-regime TAR from supplied standard-normal innovations.
      integer, intent(in) :: observations !! Number of retained observations.
      real(dp), intent(in) :: coefficient1(:) !! Regime-one intercept and lag coefficients.
      real(dp), intent(in) :: coefficient2(:) !! Regime-two intercept and lag coefficients.
      real(dp), intent(in) :: variance1 !! Regime-one innovation variance.
      real(dp), intent(in) :: variance2 !! Regime-two innovation variance.
      integer, intent(in) :: delay !! Threshold delay lag.
      real(dp), intent(in) :: threshold !! Regime threshold.
      integer, intent(in) :: lags1(:) !! Regime-one lag indices.
      integer, intent(in) :: lags2(:) !! Regime-two lag indices.
      real(dp), intent(in) :: innovations(:) !! Standard-normal innovations including burn-in.
      integer, intent(in) :: burnin !! Number of initial observations to discard.
      type(baystar_simulation_t) :: out
      real(dp), allocatable :: work(:), errors(:)
      integer :: total, start, time
      logical :: low

      total = observations + burnin
      start = max(maxval(lags1), max(maxval(lags2), delay)) + 1
      if (observations < 1 .or. burnin < 0 .or. size(innovations) /= total .or. &
         size(coefficient1) /= size(lags1) + 1 .or. &
         size(coefficient2) /= size(lags2) + 1 .or. variance1 <= 0.0_dp .or. &
         variance2 <= 0.0_dp .or. total < start) then
         out%info = 1
         return
      end if
      allocate(work(total), errors(total))
      work = 0.0_dp
      errors = 0.0_dp
      do time = start, total
         low = work(time - delay) <= threshold
         if (low) then
            errors(time) = sqrt(variance1)*innovations(time)
            work(time) = regime_prediction(work, time, coefficient1, lags1, .true.) + &
               errors(time)
         else
            errors(time) = sqrt(variance2)*innovations(time)
            work(time) = regime_prediction(work, time, coefficient2, lags2, .true.) + &
               errors(time)
         end if
      end do
      out%series = work(burnin + 1:)
      out%innovations = errors(burnin + 1:)
      out%burnin = burnin
   end function baystar_simulate_from_innovations

   function baystar_simulate(observations, coefficient1, coefficient2, &
      variance1, variance2, delay, threshold, lags1, lags2, burnin) result(out)
      !! Simulate a two-regime TAR using the shared random stream.
      integer, intent(in) :: observations !! Number of retained observations.
      real(dp), intent(in) :: coefficient1(:) !! Regime-one intercept and lag coefficients.
      real(dp), intent(in) :: coefficient2(:) !! Regime-two intercept and lag coefficients.
      real(dp), intent(in) :: variance1 !! Regime-one innovation variance.
      real(dp), intent(in) :: variance2 !! Regime-two innovation variance.
      integer, intent(in) :: delay !! Threshold delay lag.
      real(dp), intent(in) :: threshold !! Regime threshold.
      integer, intent(in) :: lags1(:) !! Regime-one lag indices.
      integer, intent(in) :: lags2(:) !! Regime-two lag indices.
      integer, intent(in), optional :: burnin !! Number of initial observations to discard.
      type(baystar_simulation_t) :: out
      real(dp), allocatable :: innovations(:)
      integer :: discard, i

      discard = 1000
      if (present(burnin)) discard = burnin
      if (observations < 1 .or. discard < 0) then
         out%info = 1
         return
      end if
      allocate(innovations(observations + discard))
      do i = 1, size(innovations)
         innovations(i) = random_standard_normal()
      end do
      out = baystar_simulate_from_innovations(observations, coefficient1, &
         coefficient2, variance1, variance2, delay, threshold, lags1, lags2, &
         innovations, discard)
   end function baystar_simulate

   function baystar_fit(series, lags1, lags2, iterations, burnin, threshold_step, &
      maximum_delay, include_constant, threshold_variable, prior) result(out)
      !! Fit BAYSTAR's two-regime TAR by Gibbs and Metropolis-Hastings sampling.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: lags1(:) !! Regime-one lag indices.
      integer, intent(in) :: lags2(:) !! Regime-two lag indices.
      integer, intent(in) :: iterations !! Total MCMC iterations.
      integer, intent(in) :: burnin !! Discarded MCMC iterations.
      real(dp), intent(in) :: threshold_step !! Threshold random-walk standard deviation.
      integer, intent(in), optional :: maximum_delay !! Largest threshold delay.
      logical, intent(in), optional :: include_constant !! Include regime intercepts.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.
      type(baystar_prior_t), intent(in), optional :: prior !! Prior specification.
      type(baystar_fit_t) :: out
      type(baystar_prior_t) :: selected_prior
      type(baystar_coefficient_posterior_t) :: coefficient_posterior
      type(baystar_variance_posterior_t) :: variance_posterior
      real(dp), allocatable :: coefficient1(:), coefficient2(:), probability(:)
      real(dp), allocatable :: draw(:), retained(:, :), threshold_data(:)
      real(dp) :: variance1, variance2, threshold, proposal, old_likelihood
      real(dp) :: new_likelihood, lower_bound, upper_bound, uniform, cumulative
      integer :: delay, delay_limit, iteration, retained_count, selected, status
      integer :: parameter_count, mode_count, candidate_delay, start, time, row
      logical :: constant

      constant = .true.
      if (present(include_constant)) constant = include_constant
      delay_limit = 3
      if (present(maximum_delay)) delay_limit = maximum_delay
      if (iterations < 2 .or. burnin < 0 .or. burnin >= iterations .or. &
         threshold_step <= 0.0_dp .or. delay_limit < 1 .or. size(lags1) < 1 .or. &
         size(lags2) < 1 .or. any(lags1 < 1) .or. any(lags2 < 1) .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) < size(series) .or. &
            .not. all(ieee_is_finite(threshold_variable(:size(series))))) then
            out%info = 1
            return
         end if
         threshold_data = threshold_variable(:size(series))
      else
         threshold_data = series
      end if
      selected_prior = baystar_default_prior(series, lags1, lags2, constant)
      if (present(prior)) selected_prior = prior
      if (selected_prior%info /= 0 .or. selected_prior%variance_degrees <= 0.0_dp .or. &
         selected_prior%variance_scale <= 0.0_dp) then
         out%info = 2
         return
      end if
      allocate(coefficient1(size(lags1) + merge(1, 0, constant)), &
         coefficient2(size(lags2) + merge(1, 0, constant)))
      coefficient1 = 0.05_dp
      coefficient2 = 0.05_dp
      variance1 = 0.2_dp
      variance2 = 0.2_dp
      delay = 1
      threshold = median(threshold_data)
      lower_bound = quantile(sorted(threshold_data), 0.25_dp)
      upper_bound = quantile(sorted(threshold_data), 0.75_dp)
      allocate(out%coefficient1(size(coefficient1), iterations), &
         out%coefficient2(size(coefficient2), iterations), &
         out%innovation_variance(2, iterations), out%threshold(iterations), &
         out%delay(iterations), out%threshold_accepted(iterations), &
         out%log_likelihood(iterations), out%regime_mean(2, iterations))
      out%threshold_accepted = .false.
      do iteration = 1, iterations
         coefficient_posterior = baystar_coefficient_posterior(1, series, &
            variance1, delay, threshold, lags1, lags2, selected_prior%mean1, &
            selected_prior%precision1, constant, threshold_data)
         if (coefficient_posterior%info /= 0) then
            out%info = 10 + coefficient_posterior%info
            return
         end if
         if (allocated(draw)) deallocate(draw)
         allocate(draw(size(coefficient_posterior%mean)))
         call random_multivariate_normal(coefficient_posterior%mean, &
            coefficient_posterior%covariance, draw, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         coefficient1 = draw
         coefficient_posterior = baystar_coefficient_posterior(2, series, &
            variance2, delay, threshold, lags1, lags2, selected_prior%mean2, &
            selected_prior%precision2, constant, threshold_data)
         if (coefficient_posterior%info /= 0) then
            out%info = 30 + coefficient_posterior%info
            return
         end if
         if (allocated(draw)) deallocate(draw)
         allocate(draw(size(coefficient_posterior%mean)))
         call random_multivariate_normal(coefficient_posterior%mean, &
            coefficient_posterior%covariance, draw, status)
         if (status /= 0) then
            out%info = 40 + status
            return
         end if
         coefficient2 = draw
         variance_posterior = baystar_variance_posterior(1, series, coefficient1, &
            delay, threshold, lags1, lags2, selected_prior%variance_degrees, &
            selected_prior%variance_scale, constant, threshold_data)
         if (variance_posterior%info /= 0) then
            out%info = 50 + variance_posterior%info
            return
         end if
         variance1 = variance_posterior%rate/ &
            random_gamma(variance_posterior%shape)
         variance_posterior = baystar_variance_posterior(2, series, coefficient2, &
            delay, threshold, lags1, lags2, selected_prior%variance_degrees, &
            selected_prior%variance_scale, constant, threshold_data)
         if (variance_posterior%info /= 0) then
            out%info = 60 + variance_posterior%info
            return
         end if
         variance2 = variance_posterior%rate/ &
            random_gamma(variance_posterior%shape)
         probability = baystar_delay_probabilities(series, coefficient1, &
            coefficient2, variance1, variance2, threshold, lags1, lags2, &
            delay_limit, constant, threshold_data)
         uniform = random_uniform()
         cumulative = 0.0_dp
         delay = delay_limit
         do candidate_delay = 1, delay_limit
            cumulative = cumulative + probability(candidate_delay)
            if (uniform <= cumulative) then
               delay = candidate_delay
               exit
            end if
         end do
         do
            proposal = threshold + threshold_step*random_standard_normal()
            if (proposal >= lower_bound .and. proposal <= upper_bound) exit
         end do
         old_likelihood = baystar_log_likelihood(series, coefficient1, coefficient2, &
            variance1, variance2, delay, threshold, lags1, lags2, constant, &
            threshold_data)
         new_likelihood = baystar_log_likelihood(series, coefficient1, coefficient2, &
            variance1, variance2, delay, proposal, lags1, lags2, constant, &
            threshold_data)
         if (new_likelihood - old_likelihood > &
            log(max(random_uniform(), tiny(1.0_dp)))) then
            threshold = proposal
            out%threshold_accepted(iteration) = .true.
         end if
         out%coefficient1(:, iteration) = coefficient1
         out%coefficient2(:, iteration) = coefficient2
         out%innovation_variance(:, iteration) = [variance1, variance2]
         out%threshold(iteration) = threshold
         out%delay(iteration) = delay
         out%log_likelihood(iteration) = baystar_log_likelihood(series, &
            coefficient1, coefficient2, variance1, variance2, delay, threshold, &
            lags1, lags2, constant, threshold_data)
         if (constant) then
            out%regime_mean(1, iteration) = coefficient1(1)/ &
               (1.0_dp - sum(coefficient1(2:)))
            out%regime_mean(2, iteration) = coefficient2(1)/ &
               (1.0_dp - sum(coefficient2(2:)))
         end if
      end do
      retained_count = iterations - burnin
      parameter_count = size(coefficient1) + size(coefficient2) + 3 + &
         merge(2, 0, constant)
      allocate(retained(retained_count, parameter_count))
      retained(:, :size(coefficient1)) = transpose( &
         out%coefficient1(:, burnin + 1:))
      selected = size(coefficient1)
      retained(:, selected + 1:selected + size(coefficient2)) = transpose( &
         out%coefficient2(:, burnin + 1:))
      selected = selected + size(coefficient2)
      retained(:, selected + 1:selected + 2) = transpose( &
         out%innovation_variance(:, burnin + 1:))
      selected = selected + 2
      retained(:, selected + 1) = out%threshold(burnin + 1:)
      if (constant) retained(:, selected + 2:selected + 3) = transpose( &
         out%regime_mean(:, burnin + 1:))
      out%summary = baystar_summary(retained)
      mode_count = 0
      do candidate_delay = 1, delay_limit
         row = count(out%delay(burnin + 1:) == candidate_delay)
         if (row > mode_count) then
            mode_count = row
            out%selected_delay = candidate_delay
         end if
      end do
      new_likelihood = baystar_log_likelihood(series, &
         out%summary%mean(:size(coefficient1)), &
         out%summary%mean(size(coefficient1) + 1: &
         size(coefficient1) + size(coefficient2)), &
         out%summary%mean(size(coefficient1) + size(coefficient2) + 1), &
         out%summary%mean(size(coefficient1) + size(coefficient2) + 2), &
         out%selected_delay, &
         out%summary%mean(size(coefficient1) + size(coefficient2) + 3), &
         lags1, lags2, constant, threshold_data)
      out%dic = -4.0_dp*sum(out%log_likelihood(burnin + 1:))/ &
         real(retained_count, dp) + 2.0_dp*new_likelihood
      out%threshold_acceptance = real(count( &
         out%threshold_accepted), dp)/real(iterations, dp)
      start = baystar_start(lags1, lags2, constant)
      allocate(out%residuals(size(series) - start + 1))
      do time = start, size(series)
         if (threshold_value(series, time, out%selected_delay, threshold_data) <= &
            out%summary%mean(size(coefficient1) + size(coefficient2) + 3)) then
            out%residuals(time - start + 1) = series(time) - &
               regime_prediction(series, time, &
               out%summary%mean(:size(coefficient1)), lags1, constant)
         else
            out%residuals(time - start + 1) = series(time) - &
               regime_prediction(series, time, out%summary%mean( &
               size(coefficient1) + 1:size(coefficient1) + size(coefficient2)), &
               lags2, constant)
         end if
      end do
      out%burnin = burnin
   end function baystar_fit

   pure logical function valid_specification(series, coefficient1, coefficient2, &
      variance1, variance2, delay, lags1, lags2, constant, threshold_variable) &
      result(valid)
      !! Validate a two-regime TAR likelihood specification.
      real(dp), intent(in) :: series(:) !! Observed series.
      real(dp), intent(in) :: coefficient1(:) !! Regime-one coefficients.
      real(dp), intent(in) :: coefficient2(:) !! Regime-two coefficients.
      real(dp), intent(in) :: variance1 !! Regime-one variance.
      real(dp), intent(in) :: variance2 !! Regime-two variance.
      integer, intent(in) :: delay !! Threshold delay.
      integer, intent(in) :: lags1(:) !! Regime-one lags.
      integer, intent(in) :: lags2(:) !! Regime-two lags.
      logical, intent(in) :: constant !! Include intercepts.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.

      if (size(lags1) < 1 .or. size(lags2) < 1) then
         valid = .false.
         return
      end if
      valid = size(series) >= baystar_start(lags1, lags2, constant) .and. &
         size(coefficient1) == size(lags1) + merge(1, 0, constant) .and. &
         size(coefficient2) == size(lags2) + merge(1, 0, constant) .and. &
         variance1 > 0.0_dp .and. variance2 > 0.0_dp .and. delay >= 1 .and. &
         delay < baystar_start(lags1, lags2, constant)
      valid = valid .and. all(lags1 >= 1) .and. all(lags2 >= 1) .and. &
         all(ieee_is_finite(series))
      if (present(threshold_variable)) valid = valid .and. &
         size(threshold_variable) >= size(series) .and. &
         all(ieee_is_finite(threshold_variable(:size(series))))
   end function valid_specification

   pure integer function baystar_start(lags1, lags2, constant) result(start)
      !! Return BAYSTAR's first response index.
      integer, intent(in) :: lags1(:) !! Regime-one lags.
      integer, intent(in) :: lags2(:) !! Regime-two lags.
      logical, intent(in) :: constant !! Include intercepts.

      start = max(maxval(lags1), maxval(lags2)) + merge(2, 1, constant)
   end function baystar_start

   pure real(dp) function threshold_value(series, time, delay, &
      threshold_variable) result(value)
      !! Return the delayed internal or external threshold variable.
      real(dp), intent(in) :: series(:) !! Observed series.
      integer, intent(in) :: time !! Response time index.
      integer, intent(in) :: delay !! Threshold delay.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.

      if (present(threshold_variable)) then
         value = threshold_variable(time - delay)
      else
         value = series(time - delay)
      end if
   end function threshold_value

   pure real(dp) function regime_prediction(series, time, coefficient, lags, &
      constant) result(prediction)
      !! Evaluate one regime's sparse autoregressive mean.
      real(dp), intent(in) :: series(:) !! Observed or simulated series.
      integer, intent(in) :: time !! Response time index.
      real(dp), intent(in) :: coefficient(:) !! Regime coefficients.
      integer, intent(in) :: lags(:) !! Selected lag indices.
      logical, intent(in) :: constant !! Include an intercept.
      integer :: lag, offset

      prediction = 0.0_dp
      offset = 0
      if (constant) then
         prediction = coefficient(1)
         offset = 1
      end if
      do lag = 1, size(lags)
         prediction = prediction + coefficient(offset + lag)*series(time - lags(lag))
      end do
   end function regime_prediction

   pure subroutine fill_regressor(series, time, lags, constant, regressor)
      !! Fill one sparse autoregressive design row.
      real(dp), intent(in) :: series(:) !! Observed series.
      integer, intent(in) :: time !! Response time index.
      integer, intent(in) :: lags(:) !! Selected lag indices.
      logical, intent(in) :: constant !! Include an intercept.
      real(dp), intent(out) :: regressor(:) !! Filled design row.
      integer :: lag, offset

      offset = 0
      if (constant) then
         regressor(1) = 1.0_dp
         offset = 1
      end if
      do lag = 1, size(lags)
         regressor(offset + lag) = series(time - lags(lag))
      end do
   end subroutine fill_regressor

end module baystar_mod
