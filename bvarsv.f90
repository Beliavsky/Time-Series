! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Core algorithms translated from the R bvarsv package.
module bvarsv_mod
   !! Primiceri time-varying structural VAR state-simulation algorithms.
   use kind_mod, only: dp
   use linalg_mod, only: invert_matrix, cholesky_lower, &
      cholesky_lower_semidefinite, identity_matrix
   use random_mod, only: random_standard_normal, random_uniform, random_gamma
   implicit none
   private

   type, public :: bvarsv_state_draw_t
      !! Carter-Kohn state draw and Gaussian filtering log likelihood.
      real(dp), allocatable :: state(:, :)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: info = 0
   end type bvarsv_state_draw_t

   type, public :: bvarsv_data_t
      !! Prepared response and equation-stacked TVP-VAR design arrays.
      real(dp), allocatable :: response(:, :)
      real(dp), allocatable :: design(:, :, :)
      real(dp), allocatable :: series(:, :)
      integer :: variables = 0
      integer :: observations = 0
      integer :: states = 0
      integer :: lags = 0
      integer :: training_observations = 0
      integer :: info = 0
   end type bvarsv_data_t

   type, public :: bvarsv_volatility_draw_t
      !! KSC mixture indicators and random-walk log-volatility draw.
      real(dp), allocatable :: log_variance(:, :)
      real(dp), allocatable :: standard_deviation(:, :)
      integer, allocatable :: component(:, :)
      integer :: info = 0
   end type bvarsv_volatility_draw_t

   type, public :: bvarsv_contemporaneous_draw_t
      !! Time-varying packed contemporaneous coefficients.
      real(dp), allocatable :: alpha(:, :)
      integer :: info = 0
   end type bvarsv_contemporaneous_draw_t

   type, public :: bvarsv_covariance_draw_t
      !! Innovation covariance draw with a status code.
      real(dp), allocatable :: covariance(:, :)
      integer :: info = 0
   end type bvarsv_covariance_draw_t

   type, public :: bvarsv_prior_t
      !! Training-sample prior moments used by the Primiceri sampler.
      real(dp), allocatable :: coefficient_mean(:)
      real(dp), allocatable :: coefficient_covariance(:, :)
      real(dp), allocatable :: contemporaneous_mean(:)
      real(dp), allocatable :: contemporaneous_covariance(:, :)
      real(dp), allocatable :: log_variance_mean(:)
      integer :: training_observations = 0
      integer :: info = 0
   end type bvarsv_prior_t

   type, public :: bvarsv_draws_t
      !! Retained draws from the Primiceri TVP-SVAR Gibbs sampler.
      real(dp), allocatable :: coefficient(:, :, :)
      real(dp), allocatable :: contemporaneous(:, :, :)
      real(dp), allocatable :: log_variance(:, :, :)
      real(dp), allocatable :: covariance(:, :, :, :)
      real(dp), allocatable :: coefficient_state_covariance(:, :, :)
      real(dp), allocatable :: contemporaneous_state_covariance(:, :, :)
      real(dp), allocatable :: volatility_state_covariance(:, :, :)
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvarsv_draws_t

   type, public :: bvarsv_predictive_t
      !! Draw-wise recursive forecast means, covariances, and realizations.
      real(dp), allocatable :: mean(:, :, :)
      real(dp), allocatable :: covariance(:, :, :, :)
      real(dp), allocatable :: path(:, :, :)
      integer :: info = 0
   end type bvarsv_predictive_t

   type, public :: bvarsv_irf_t
      !! Draw-wise time-selected impulse-response matrices.
      real(dp), allocatable :: response(:, :, :, :)
      integer :: selected_time = 0
      integer :: scenario = 0
      integer :: info = 0
   end type bvarsv_irf_t

   type, public :: bvarsv_simulation_t
      !! Simulated TVP-VAR(1) observations and latent parameter paths.
      real(dp), allocatable :: series(:, :)
      real(dp), allocatable :: coefficient(:, :, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: contemporaneous(:, :)
      real(dp), allocatable :: log_variance(:, :)
      integer :: info = 0
   end type bvarsv_simulation_t

   type, public :: bvarsv_predictive_draws_t
      !! Selected predictive observations, means, and variances by draw.
      real(dp), allocatable :: observation(:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: variance(:)
      integer :: info = 0
   end type bvarsv_predictive_draws_t

   type, public :: bvarsv_parameter_path_t
      !! Selected parameter values with posterior draws in rows and time in columns.
      real(dp), allocatable :: value(:, :)
      integer :: info = 0
   end type bvarsv_parameter_path_t

   public :: bvarsv_carter_kohn
   public :: bvarsv_prepare, bvarsv_coefficient_update
   public :: bvarsv_ols_prior
   public :: bvarsv_gibbs
   public :: bvarsv_forecast
   public :: bvarsv_irf
   public :: bvarsv_predictive_density
   public :: bvarsv_predictive_draws, bvarsv_parameter_draws
   public :: bvarsv_simulate_var1
   public :: bvarsv_contemporaneous_matrix, bvarsv_covariance_from_state
   public :: bvarsv_log_volatility_update
   public :: bvarsv_contemporaneous_update
   public :: bvarsv_random_walk_covariance
   public :: bvarsv_contemporaneous_covariance_update

contains

   pure function bvarsv_predictive_draws(prediction, variable, horizon) result(out)
      !! Extract observation, mean, and variance draws for one forecast target.
      type(bvarsv_predictive_t), intent(in) :: prediction !! Draw-wise predictive moments and paths.
      integer, intent(in), optional :: variable !! One-based response-variable index.
      integer, intent(in), optional :: horizon !! One-based forecast-horizon index.
      type(bvarsv_predictive_draws_t) :: out
      integer :: selected_variable, selected_horizon, draws

      selected_variable = 1
      if (present(variable)) selected_variable = variable
      selected_horizon = 1
      if (present(horizon)) selected_horizon = horizon
      if (prediction%info /= 0 .or. .not. allocated(prediction%path) .or. &
         .not. allocated(prediction%mean) .or. &
         .not. allocated(prediction%covariance)) then
         out%info = 1
         return
      end if
      draws = size(prediction%path, 3)
      if (selected_variable < 1 .or. &
         selected_variable > size(prediction%path, 1) .or. &
         selected_variable > size(prediction%mean, 1) .or. &
         selected_variable > size(prediction%covariance, 1) .or. &
         selected_variable > size(prediction%covariance, 2) .or. &
         selected_horizon < 1 .or. &
         selected_horizon > size(prediction%path, 2) .or. &
         selected_horizon > size(prediction%mean, 2) .or. &
         selected_horizon > size(prediction%covariance, 3) .or. draws < 1 .or. &
         size(prediction%mean, 3) /= draws .or. &
         size(prediction%covariance, 4) /= draws) then
         out%info = 2
         return
      end if
      allocate(out%observation(draws), out%mean(draws), out%variance(draws))
      out%observation = prediction%path(selected_variable, selected_horizon, :)
      out%mean = prediction%mean(selected_variable, selected_horizon, :)
      out%variance = prediction%covariance(selected_variable, selected_variable, &
         selected_horizon, :)
   end function bvarsv_predictive_draws

   pure function bvarsv_parameter_draws(draws, data, parameter, row, column, &
      lag) result(out)
      !! Extract one intercept, lag, or covariance path from retained draws.
      type(bvarsv_draws_t), intent(in) :: draws !! Retained Primiceri posterior draws.
      type(bvarsv_data_t), intent(in) :: data !! Prepared model dimensions.
      character(len=*), intent(in) :: parameter !! Parameter category: intercept, lag, covariance, or vcv.
      integer, intent(in), optional :: row !! One-based response-equation or covariance row.
      integer, intent(in), optional :: column !! One-based predictor or covariance column.
      integer, intent(in), optional :: lag !! One-based lag index for lag coefficients.
      type(bvarsv_parameter_path_t) :: out
      integer :: selected_row, selected_column, selected_lag
      integer :: draw, time, position

      selected_row = 1
      if (present(row)) selected_row = row
      selected_column = 1
      if (present(column)) selected_column = column
      selected_lag = 1
      if (present(lag)) selected_lag = lag
      if (draws%info /= 0 .or. data%info /= 0 .or. &
         draws%retained_draws < 1 .or. data%observations < 1 .or. &
         selected_row < 1 .or. selected_row > data%variables .or. &
         selected_column < 1 .or. selected_column > data%variables) then
         out%info = 1
         return
      end if
      allocate(out%value(draws%retained_draws, data%observations))
      out%value = 0.0_dp
      select case (trim(parameter))
      case ('intercept')
         if (.not. allocated(draws%coefficient)) then
            out%info = 2
            return
         end if
         if (size(draws%coefficient, 1) < selected_row .or. &
            size(draws%coefficient, 2) /= data%observations .or. &
            size(draws%coefficient, 3) /= draws%retained_draws) then
            out%info = 2
            return
         end if
         do draw = 1, draws%retained_draws
            out%value(draw, :) = draws%coefficient(selected_row, :, draw)
         end do
      case ('lag')
         if (selected_lag < 1 .or. selected_lag > data%lags .or. &
            .not. allocated(draws%coefficient)) then
            out%info = 3
            return
         end if
         position = data%variables + (selected_lag - 1)*data%variables**2 + &
            (selected_row - 1)*data%variables + selected_column
         if (size(draws%coefficient, 1) < position .or. &
            size(draws%coefficient, 2) /= data%observations .or. &
            size(draws%coefficient, 3) /= draws%retained_draws) then
            out%info = 3
            return
         end if
         do draw = 1, draws%retained_draws
            out%value(draw, :) = draws%coefficient(position, :, draw)
         end do
      case ('covariance', 'vcv')
         if (.not. allocated(draws%covariance)) then
            out%info = 4
            return
         end if
         if (size(draws%covariance, 1) < selected_row .or. &
            size(draws%covariance, 2) < selected_column .or. &
            size(draws%covariance, 3) /= data%observations .or. &
            size(draws%covariance, 4) /= draws%retained_draws) then
            out%info = 4
            return
         end if
         do draw = 1, draws%retained_draws
            do time = 1, data%observations
               out%value(draw, time) = &
                  draws%covariance(selected_row, selected_column, time, draw)
            end do
         end do
      case default
         out%info = 5
      end select
   end function bvarsv_parameter_draws

   function bvarsv_simulate_var1(initial_coefficient, initial_contemporaneous, &
      initial_log_variance, coefficient_state_covariance, &
      contemporaneous_state_covariance, volatility_state_covariance, &
      observations, burnin) result(out)
      !! Simulate bvarsv's VAR(1) with random-walk parameters and volatility.
      real(dp), intent(in) :: initial_coefficient(:, :) !! Initial intercept and lag-one coefficient matrix.
      real(dp), intent(in) :: initial_contemporaneous(:) !! Initial packed strict-lower contemporaneous state.
      real(dp), intent(in) :: initial_log_variance(:) !! Initial structural log variances.
      real(dp), intent(in) :: coefficient_state_covariance(:, :) !! Coefficient random-walk covariance Q.
      real(dp), intent(in) :: contemporaneous_state_covariance(:, :) !! Contemporaneous random-walk covariance S.
      real(dp), intent(in) :: volatility_state_covariance(:, :) !! Log-variance random-walk covariance W.
      integer, intent(in) :: observations !! Number of observations retained after burn-in.
      integer, intent(in), optional :: burnin !! Number of initial observations discarded.
      type(bvarsv_simulation_t) :: out
      real(dp), allocatable :: coefficient(:), contemporaneous(:), log_variance(:)
      real(dp), allocatable :: coefficient_factor(:, :), contemporaneous_factor(:, :)
      real(dp), allocatable :: volatility_factor(:, :), observation_factor(:, :)
      real(dp), allocatable :: normal(:), previous(:), mean(:), covariance(:, :)
      integer :: variables, states, count, discard, time, retained, equation
      integer :: source, position, info

      variables = size(initial_coefficient, 1)
      states = variables + variables**2
      count = variables*(variables - 1)/2
      discard = 1000
      if (present(burnin)) discard = burnin
      if (variables < 2 .or. size(initial_coefficient, 2) /= variables + 1 .or. &
         size(initial_contemporaneous) /= count .or. &
         size(initial_log_variance) /= variables .or. observations < 1 .or. &
         discard < 0 .or. any(shape(coefficient_state_covariance) /= &
         [states, states]) .or. &
         any(shape(contemporaneous_state_covariance) /= [count, count]) .or. &
         any(shape(volatility_state_covariance) /= [variables, variables])) then
         out%info = 1
         return
      end if
      call cholesky_lower_semidefinite(coefficient_state_covariance, &
         coefficient_factor, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      call cholesky_lower_semidefinite(contemporaneous_state_covariance, &
         contemporaneous_factor, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      call cholesky_lower_semidefinite(volatility_state_covariance, &
         volatility_factor, info)
      if (info /= 0) then
         out%info = 4
         return
      end if
      allocate(coefficient(states), contemporaneous(count), log_variance(variables))
      allocate(normal(max(states, max(count, variables))))
      allocate(previous(variables), mean(variables))
      allocate(out%series(variables, observations))
      allocate(out%coefficient(variables, variables + 1, observations))
      allocate(out%covariance(variables, variables, observations))
      allocate(out%contemporaneous(count, observations))
      allocate(out%log_variance(variables, observations))
      coefficient(:variables) = initial_coefficient(:, 1)
      do equation = 1, variables
         position = variables + (equation - 1)*variables
         coefficient(position + 1:position + variables) = &
            initial_coefficient(equation, 2:)
      end do
      contemporaneous = initial_contemporaneous
      log_variance = initial_log_variance
      previous = 0.0_dp
      do time = 1, discard + observations
         call fill_standard_normal(normal(:states))
         coefficient = coefficient + matmul(coefficient_factor, normal(:states))
         if (count > 0) then
            call fill_standard_normal(normal(:count))
            contemporaneous = contemporaneous + &
               matmul(contemporaneous_factor, normal(:count))
         end if
         call fill_standard_normal(normal(:variables))
         log_variance = log_variance + &
            matmul(volatility_factor, normal(:variables))
         mean = coefficient(:variables)
         do equation = 1, variables
            do source = 1, variables
               position = variables + (equation - 1)*variables + source
               mean(equation) = mean(equation) + coefficient(position)*previous(source)
            end do
         end do
         covariance = bvarsv_covariance_from_state(contemporaneous, log_variance)
         call cholesky_lower(covariance, observation_factor, info)
         if (info /= 0) then
            out%info = 5
            return
         end if
         call fill_standard_normal(normal(:variables))
         previous = mean + matmul(observation_factor, normal(:variables))
         if (time > discard) then
            retained = time - discard
            out%series(:, retained) = previous
            out%coefficient(:, 1, retained) = coefficient(:variables)
            do equation = 1, variables
               position = variables + (equation - 1)*variables
               out%coefficient(equation, 2:, retained) = &
                  coefficient(position + 1:position + variables)
            end do
            out%covariance(:, :, retained) = covariance
            out%contemporaneous(:, retained) = contemporaneous
            out%log_variance(:, retained) = log_variance
         end if
      end do
   end function bvarsv_simulate_var1

   pure function bvarsv_prepare(series, lags, training_observations) result(out)
      !! Construct the equation-stacked TVP-VAR design used by bvarsv.
      real(dp), intent(in) :: series(:, :) !! Endogenous variables by chronological time.
      integer, intent(in) :: lags !! Positive autoregressive lag order.
      integer, intent(in), optional :: training_observations !! Initial usable observations reserved for priors.
      type(bvarsv_data_t) :: out
      integer :: training, total, time, output_time, lag, equation, source, position

      training = 0
      if (present(training_observations)) training = training_observations
      out%variables = size(series, 1)
      total = size(series, 2)
      out%lags = lags
      out%training_observations = training
      out%states = out%variables + lags*out%variables*out%variables
      out%observations = total - lags - training
      if (out%variables < 2 .or. lags < 1 .or. training < 0 .or. &
         out%observations < 1) then
         out%info = 1
         return
      end if
      out%series = series
      allocate(out%response(out%variables, out%observations))
      allocate(out%design(out%variables, out%states, out%observations))
      out%design = 0.0_dp
      do output_time = 1, out%observations
         time = lags + training + output_time
         out%response(:, output_time) = series(:, time)
         do equation = 1, out%variables
            out%design(equation, equation, output_time) = 1.0_dp
            do lag = 1, lags
               do source = 1, out%variables
                  position = out%variables + (lag - 1)*out%variables**2 + &
                     (equation - 1)*out%variables + source
                  out%design(equation, position, output_time) = &
                     series(source, time - lag)
               end do
            end do
         end do
      end do
   end function bvarsv_prepare

   function bvarsv_coefficient_update(data, observation_covariance, &
      state_covariance, initial_mean, initial_covariance) result(out)
      !! Draw the complete time-varying VAR coefficient path.
      type(bvarsv_data_t), intent(in) :: data !! Prepared bvarsv response and design arrays.
      real(dp), intent(in) :: observation_covariance(:, :, :) !! Reduced-form covariance by time.
      real(dp), intent(in) :: state_covariance(:, :) !! Coefficient random-walk covariance Q.
      real(dp), intent(in) :: initial_mean(:) !! Prior mean of the initial coefficient state.
      real(dp), intent(in) :: initial_covariance(:, :) !! Prior covariance of the initial coefficient state.
      type(bvarsv_state_draw_t) :: out

      if (data%info /= 0 .or. any(shape(observation_covariance) /= &
         [data%variables, data%variables, data%observations]) .or. &
         any(shape(state_covariance) /= [data%states, data%states]) .or. &
         size(initial_mean) /= data%states .or. &
         any(shape(initial_covariance) /= [data%states, data%states])) then
         out%info = 1
         return
      end if
      out = bvarsv_carter_kohn(data%response, data%design, &
         observation_covariance, state_covariance, initial_mean, initial_covariance)
      if (out%info /= 0) out%info = 10 + out%info
   end function bvarsv_coefficient_update

   function bvarsv_ols_prior(series, lags, covariance_draws) result(out)
      !! Construct bvarsv's OLS and inverse-Wishart training-sample prior moments.
      real(dp), intent(in) :: series(:, :) !! Training series by variable and time.
      integer, intent(in) :: lags !! Positive autoregressive lag order.
      integer, intent(in), optional :: covariance_draws !! Monte Carlo draws for alpha covariance.
      type(bvarsv_prior_t) :: out
      type(bvarsv_data_t) :: data
      real(dp), allocatable :: crossproduct(:, :), crossresponse(:), inverse(:, :)
      real(dp), allocatable :: residual(:, :), residual_covariance(:, :)
      real(dp), allocatable :: residual_precision(:, :), weighted_crossproduct(:, :)
      real(dp), allocatable :: lower(:, :), normalized(:, :), impact(:, :), scale(:)
      real(dp), allocatable :: covariance_draw(:, :), precision_draw(:, :), alpha(:)
      real(dp), allocatable :: alpha_sum(:), alpha_second(:, :)
      integer :: draws, variables, observations, states, count
      integer :: time, item, draw, info

      draws = 2000
      if (present(covariance_draws)) draws = covariance_draws
      data = bvarsv_prepare(series, lags)
      variables = data%variables
      observations = data%observations
      states = data%states
      count = variables*(variables - 1)/2
      if (data%info /= 0 .or. draws < 2 .or. observations <= variables) then
         out%info = 1
         return
      end if
      allocate(crossproduct(states, states), crossresponse(states))
      crossproduct = 0.0_dp
      crossresponse = 0.0_dp
      do time = 1, observations
         crossproduct = crossproduct + matmul(transpose(data%design(:, :, time)), &
            data%design(:, :, time))
         crossresponse = crossresponse + matmul(transpose(data%design(:, :, time)), &
            data%response(:, time))
      end do
      call invert_matrix(crossproduct, inverse, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      out%coefficient_mean = matmul(inverse, crossresponse)
      allocate(residual(variables, observations))
      do time = 1, observations
         residual(:, time) = data%response(:, time) - &
            matmul(data%design(:, :, time), out%coefficient_mean)
      end do
      residual_covariance = matmul(residual, transpose(residual))/ &
         real(observations, dp)
      call invert_matrix(residual_covariance, residual_precision, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      allocate(weighted_crossproduct(states, states))
      weighted_crossproduct = 0.0_dp
      do time = 1, observations
         weighted_crossproduct = weighted_crossproduct + matmul( &
            transpose(data%design(:, :, time)), matmul(residual_precision, &
            data%design(:, :, time)))
      end do
      call invert_matrix(weighted_crossproduct, out%coefficient_covariance, info)
      if (info /= 0) then
         out%info = 4
         return
      end if
      call cholesky_lower(residual_covariance, lower, info)
      if (info /= 0) then
         out%info = 5
         return
      end if
      allocate(scale(variables), out%log_variance_mean(variables))
      do item = 1, variables
         scale(item) = lower(item, item)
      end do
      out%log_variance_mean = log(scale**2)
      normalized = lower
      do item = 1, variables
         normalized(:, item) = normalized(:, item)/scale(item)
      end do
      call invert_matrix(normalized, impact, info)
      if (info /= 0) then
         out%info = 6
         return
      end if
      out%contemporaneous_mean = pack_strict_lower(impact)
      allocate(alpha_sum(count), alpha_second(count, count))
      alpha_sum = 0.0_dp
      alpha_second = 0.0_dp
      do draw = 1, draws
         call inverse_wishart_draw(real(observations, dp)*residual_covariance, &
            real(observations, dp), covariance_draw, precision_draw, info)
         if (info /= 0) then
            out%info = 7
            return
         end if
         call cholesky_lower(covariance_draw, lower, info)
         if (info /= 0) then
            out%info = 8
            return
         end if
         do item = 1, variables
            scale(item) = lower(item, item)
            normalized(:, item) = lower(:, item)/scale(item)
         end do
         call invert_matrix(normalized, impact, info)
         if (info /= 0) then
            out%info = 9
            return
         end if
         alpha = pack_strict_lower(impact)
         alpha_sum = alpha_sum + alpha
         alpha_second = alpha_second + spread(alpha, 2, count)* &
            spread(alpha, 1, count)
      end do
      alpha_sum = alpha_sum/real(draws, dp)
      out%contemporaneous_covariance = alpha_second/real(draws, dp) - &
         spread(alpha_sum, 2, count)*spread(alpha_sum, 1, count)
      out%contemporaneous_covariance = 0.5_dp*( &
         out%contemporaneous_covariance + transpose(out%contemporaneous_covariance))
      out%training_observations = observations
   end function bvarsv_ols_prior

   function bvarsv_gibbs(data, prior, draws, burnin, thin, coefficient_scale, &
      contemporaneous_scale, volatility_scale) result(out)
      !! Run the Primiceri bvarsv Gibbs sampler with random-walk parameters and SV.
      type(bvarsv_data_t), intent(in) :: data !! Prepared estimation-sample data.
      type(bvarsv_prior_t), intent(in) :: prior !! Training-sample OLS prior moments.
      integer, intent(in) :: draws !! Number of retained posterior draws.
      integer, intent(in), optional :: burnin !! Number of discarded Gibbs iterations.
      integer, intent(in), optional :: thin !! Iterations between retained draws.
      real(dp), intent(in), optional :: coefficient_scale !! Innovation prior scale k_Q.
      real(dp), intent(in), optional :: contemporaneous_scale !! Innovation prior scale k_S.
      real(dp), intent(in), optional :: volatility_scale !! Innovation prior scale k_W.
      type(bvarsv_draws_t) :: out
      type(bvarsv_state_draw_t) :: coefficient_draw
      type(bvarsv_contemporaneous_draw_t) :: contemporaneous_draw
      type(bvarsv_volatility_draw_t) :: volatility_draw
      type(bvarsv_covariance_draw_t) :: covariance_draw
      real(dp), allocatable :: coefficient(:, :), contemporaneous(:, :)
      real(dp), allocatable :: log_variance(:, :), covariance(:, :, :)
      real(dp), allocatable :: coefficient_state_covariance(:, :)
      real(dp), allocatable :: contemporaneous_state_covariance(:, :)
      real(dp), allocatable :: volatility_state_covariance(:, :)
      real(dp), allocatable :: coefficient_prior_scale(:, :)
      real(dp), allocatable :: contemporaneous_prior_scale(:, :)
      real(dp), allocatable :: volatility_prior_scale(:, :)
      real(dp), allocatable :: residual(:, :), structural_residual(:, :)
      real(dp), allocatable :: matrix(:, :), coefficient_initial_covariance(:, :)
      real(dp), allocatable :: contemporaneous_initial_covariance(:, :)
      real(dp), allocatable :: volatility_initial_covariance(:, :)
      real(dp), allocatable :: block_degrees(:)
      real(dp) :: coefficient_scale_value, contemporaneous_scale_value
      real(dp) :: volatility_scale_value, coefficient_degrees, volatility_degrees
      integer :: discard, thinning, total, iteration, retained
      integer :: variables, observations, states, count, time, equation
      integer :: block_first, block_size

      discard = 0
      thinning = 1
      coefficient_scale_value = 0.01_dp
      contemporaneous_scale_value = 0.1_dp
      volatility_scale_value = 0.01_dp
      if (present(burnin)) discard = burnin
      if (present(thin)) thinning = thin
      if (present(coefficient_scale)) coefficient_scale_value = coefficient_scale
      if (present(contemporaneous_scale)) contemporaneous_scale_value = &
         contemporaneous_scale
      if (present(volatility_scale)) volatility_scale_value = volatility_scale
      variables = data%variables
      observations = data%observations
      states = data%states
      count = variables*(variables - 1)/2
      if (data%info /= 0 .or. prior%info /= 0 .or. draws < 1 .or. discard < 0 .or. &
         thinning < 1 .or. coefficient_scale_value <= 0.0_dp .or. &
         contemporaneous_scale_value <= 0.0_dp .or. &
         volatility_scale_value <= 0.0_dp .or. &
         size(prior%coefficient_mean) /= states .or. &
         any(shape(prior%coefficient_covariance) /= [states, states]) .or. &
         size(prior%contemporaneous_mean) /= count .or. &
         any(shape(prior%contemporaneous_covariance) /= [count, count]) .or. &
         size(prior%log_variance_mean) /= variables) then
         out%info = 1
         return
      end if
      coefficient_degrees = real(max(prior%training_observations, states + 2), dp)
      volatility_degrees = real(variables + 1, dp)
      coefficient_prior_scale = coefficient_scale_value**2*coefficient_degrees* &
         prior%coefficient_covariance
      contemporaneous_prior_scale = contemporaneous_scale_value**2* &
         prior%contemporaneous_covariance
      volatility_prior_scale = volatility_scale_value**2*volatility_degrees* &
         identity_matrix(variables)
      coefficient_initial_covariance = 4.0_dp*prior%coefficient_covariance
      contemporaneous_initial_covariance = 4.0_dp* &
         prior%contemporaneous_covariance
      volatility_initial_covariance = identity_matrix(variables)
      coefficient = spread(prior%coefficient_mean, 2, observations)
      contemporaneous = spread(prior%contemporaneous_mean, 2, observations)
      log_variance = spread(prior%log_variance_mean, 2, observations)
      coefficient_state_covariance = 0.0001_dp*identity_matrix(states)
      contemporaneous_state_covariance = 0.0001_dp*identity_matrix(count)
      volatility_state_covariance = 0.0001_dp*identity_matrix(variables)
      allocate(covariance(variables, variables, observations))
      do time = 1, observations
         covariance(:, :, time) = bvarsv_covariance_from_state( &
            contemporaneous(:, time), log_variance(:, time))
      end do
      allocate(block_degrees(variables - 1))
      block_first = 1
      do equation = 2, variables
         block_degrees(equation - 1) = real(equation + 1, dp)
         block_size = equation - 1
         contemporaneous_prior_scale(block_first:block_first + block_size - 1, &
            block_first:block_first + block_size - 1) = &
            block_degrees(equation - 1)*contemporaneous_prior_scale( &
            block_first:block_first + block_size - 1, &
            block_first:block_first + block_size - 1)
         block_first = block_first + block_size
      end do
      allocate(out%coefficient(states, observations, draws))
      allocate(out%contemporaneous(count, observations, draws))
      allocate(out%log_variance(variables, observations, draws))
      allocate(out%covariance(variables, variables, observations, draws))
      allocate(out%coefficient_state_covariance(states, states, draws))
      allocate(out%contemporaneous_state_covariance(count, count, draws))
      allocate(out%volatility_state_covariance(variables, variables, draws))
      total = discard + draws*thinning
      retained = 0
      do iteration = 1, total
         coefficient_draw = bvarsv_coefficient_update(data, covariance, &
            coefficient_state_covariance, prior%coefficient_mean, &
            coefficient_initial_covariance)
         if (coefficient_draw%info /= 0) then
            out%info = 10 + coefficient_draw%info
            return
         end if
         coefficient = coefficient_draw%state
         covariance_draw = bvarsv_random_walk_covariance(coefficient, &
            coefficient_prior_scale, coefficient_degrees, degrees_offset=1)
         if (covariance_draw%info /= 0) then
            out%info = 20 + covariance_draw%info
            return
         end if
         coefficient_state_covariance = covariance_draw%covariance
         allocate(residual(variables, observations))
         do time = 1, observations
            residual(:, time) = data%response(:, time) - &
               matmul(data%design(:, :, time), coefficient(:, time))
         end do
         contemporaneous_draw = bvarsv_contemporaneous_update(residual, &
            log_variance, contemporaneous_state_covariance, &
            prior%contemporaneous_mean, contemporaneous_initial_covariance)
         if (contemporaneous_draw%info /= 0) then
            out%info = 30 + contemporaneous_draw%info
            return
         end if
         contemporaneous = contemporaneous_draw%alpha
         covariance_draw = bvarsv_contemporaneous_covariance_update( &
            contemporaneous, contemporaneous_prior_scale, block_degrees)
         if (covariance_draw%info /= 0) then
            out%info = 40 + covariance_draw%info
            return
         end if
         contemporaneous_state_covariance = covariance_draw%covariance
         allocate(structural_residual(variables, observations))
         do time = 1, observations
            matrix = bvarsv_contemporaneous_matrix(contemporaneous(:, time), &
               variables)
            structural_residual(:, time) = matmul(matrix, residual(:, time))
         end do
         volatility_draw = bvarsv_log_volatility_update(structural_residual, &
            log_variance, volatility_state_covariance, prior%log_variance_mean, &
            volatility_initial_covariance)
         if (volatility_draw%info /= 0) then
            out%info = 50 + volatility_draw%info
            return
         end if
         log_variance = volatility_draw%log_variance
         covariance_draw = bvarsv_random_walk_covariance(log_variance, &
            volatility_prior_scale, volatility_degrees, degrees_offset=1)
         if (covariance_draw%info /= 0) then
            out%info = 60 + covariance_draw%info
            return
         end if
         volatility_state_covariance = covariance_draw%covariance
         do time = 1, observations
            covariance(:, :, time) = bvarsv_covariance_from_state( &
               contemporaneous(:, time), log_variance(:, time))
         end do
         if (iteration > discard .and. mod(iteration - discard, thinning) == 0) then
            retained = retained + 1
            out%coefficient(:, :, retained) = coefficient
            out%contemporaneous(:, :, retained) = contemporaneous
            out%log_variance(:, :, retained) = log_variance
            out%covariance(:, :, :, retained) = covariance
            out%coefficient_state_covariance(:, :, retained) = &
               coefficient_state_covariance
            out%contemporaneous_state_covariance(:, :, retained) = &
               contemporaneous_state_covariance
            out%volatility_state_covariance(:, :, retained) = &
               volatility_state_covariance
         end if
         deallocate(residual, structural_residual)
      end do
      out%retained_draws = retained
   end function bvarsv_gibbs

   function bvarsv_forecast(draws, data, horizon, parameter_drift) result(out)
      !! Simulate posterior forecasts with fixed or random-walk parameter states.
      type(bvarsv_draws_t), intent(in) :: draws !! Retained Primiceri posterior draws.
      type(bvarsv_data_t), intent(in) :: data !! Prepared model data and source series.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      logical, intent(in), optional :: parameter_drift !! Propagate parameter random walks.
      type(bvarsv_predictive_t) :: out
      real(dp), allocatable :: coefficient(:), contemporaneous(:), log_variance(:)
      real(dp), allocatable :: history(:, :), mean(:), covariance(:, :), normal(:)
      real(dp), allocatable :: innovation(:), path(:)
      logical :: drift
      integer :: draw, step, lag, equation, source, position, info
      integer :: variables, states, count

      drift = .true.
      if (present(parameter_drift)) drift = parameter_drift
      variables = data%variables
      states = data%states
      count = variables*(variables - 1)/2
      if (draws%info /= 0 .or. data%info /= 0 .or. horizon < 1 .or. &
         draws%retained_draws < 1 .or. .not. allocated(draws%coefficient) .or. &
         size(data%series, 2) < data%lags) then
         out%info = 1
         return
      end if
      allocate(out%mean(variables, horizon, draws%retained_draws))
      allocate(out%covariance(variables, variables, horizon, draws%retained_draws))
      allocate(out%path(variables, horizon, draws%retained_draws))
      allocate(coefficient(states), contemporaneous(count), log_variance(variables))
      allocate(history(variables, data%lags), mean(variables), path(variables))
      allocate(normal(max(states, max(count, variables))), innovation(max(states, &
         max(count, variables))))
      do draw = 1, draws%retained_draws
         coefficient = draws%coefficient(:, data%observations, draw)
         contemporaneous = draws%contemporaneous(:, data%observations, draw)
         log_variance = draws%log_variance(:, data%observations, draw)
         history = data%series(:, size(data%series, 2) - data%lags + 1:)
         do step = 1, horizon
            if (drift) then
               call fill_standard_normal(normal(:states))
               call draw_normal(0.0_dp*coefficient, &
                  draws%coefficient_state_covariance(:, :, draw), normal(:states), &
                  innovation(:states), info)
               if (info /= 0) then
                  out%info = 2
                  return
               end if
               coefficient = coefficient + innovation(:states)
               if (count > 0) then
                  call fill_standard_normal(normal(:count))
                  call draw_normal(0.0_dp*contemporaneous, &
                     draws%contemporaneous_state_covariance(:, :, draw), &
                     normal(:count), innovation(:count), info)
                  if (info /= 0) then
                     out%info = 3
                     return
                  end if
                  contemporaneous = contemporaneous + innovation(:count)
               end if
               call fill_standard_normal(normal(:variables))
               call draw_normal(0.0_dp*log_variance, &
                  draws%volatility_state_covariance(:, :, draw), &
                  normal(:variables), innovation(:variables), info)
               if (info /= 0) then
                  out%info = 4
                  return
               end if
               log_variance = log_variance + innovation(:variables)
            end if
            mean = coefficient(:variables)
            do equation = 1, variables
               do lag = 1, data%lags
                  do source = 1, variables
                     position = variables + (lag - 1)*variables**2 + &
                        (equation - 1)*variables + source
                     mean(equation) = mean(equation) + coefficient(position)* &
                        history(source, data%lags + 1 - lag)
                  end do
               end do
            end do
            covariance = bvarsv_covariance_from_state(contemporaneous, log_variance)
            call fill_standard_normal(normal(:variables))
            call draw_normal(mean, covariance, normal(:variables), path, info)
            if (info /= 0) then
               out%info = 5
               return
            end if
            out%mean(:, step, draw) = mean
            out%covariance(:, :, step, draw) = covariance
            out%path(:, step, draw) = path
            if (data%lags > 1) history(:, :data%lags - 1) = history(:, 2:)
            history(:, data%lags) = path
         end do
      end do
   end function bvarsv_forecast

   function bvarsv_irf(draws, data, horizon, selected_time, scenario) result(out)
      !! Compute bvarsv time-selected impulse responses under three impact scenarios.
      type(bvarsv_draws_t), intent(in) :: draws !! Retained Primiceri posterior draws.
      type(bvarsv_data_t), intent(in) :: data !! Prepared model dimensions.
      integer, intent(in) :: horizon !! Largest nonnegative response horizon.
      integer, intent(in), optional :: selected_time !! State time used for responses.
      integer, intent(in), optional :: scenario !! One: identity, two: Cholesky, three: Primiceri.
      type(bvarsv_irf_t) :: out
      real(dp), allocatable :: companion(:, :), power(:, :), impact(:, :)
      real(dp), allocatable :: contemporaneous(:, :), inverse(:, :), scale(:, :)
      real(dp), allocatable :: mean_standard_deviation(:)
      integer :: time_value, scenario_value, variables, companion_size
      integer :: draw, step, lag, equation, source, position, item, info

      time_value = data%observations
      scenario_value = 2
      if (present(selected_time)) time_value = selected_time
      if (present(scenario)) scenario_value = scenario
      variables = data%variables
      companion_size = variables*data%lags
      if (draws%info /= 0 .or. data%info /= 0 .or. draws%retained_draws < 1 .or. &
         horizon < 0 .or. time_value < 1 .or. time_value > data%observations .or. &
         scenario_value < 1 .or. scenario_value > 3) then
         out%info = 1
         return
      end if
      allocate(out%response(variables, variables, horizon + 1, &
         draws%retained_draws))
      allocate(mean_standard_deviation(variables))
      mean_standard_deviation = 0.0_dp
      if (scenario_value == 3) then
         do draw = 1, draws%retained_draws
            do item = 1, variables
               mean_standard_deviation(item) = mean_standard_deviation(item) + &
                  sum(exp(0.5_dp*draws%log_variance(item, :, draw)))
            end do
         end do
         mean_standard_deviation = mean_standard_deviation/ &
            real(data%observations*draws%retained_draws, dp)
      end if
      do draw = 1, draws%retained_draws
         allocate(companion(companion_size, companion_size))
         companion = 0.0_dp
         do lag = 1, data%lags
            do equation = 1, variables
               do source = 1, variables
                  position = variables + (lag - 1)*variables**2 + &
                     (equation - 1)*variables + source
                  companion(equation, (lag - 1)*variables + source) = &
                     draws%coefficient(position, time_value, draw)
               end do
            end do
         end do
         if (data%lags > 1) then
            companion(variables + 1:, :companion_size - variables) = &
               identity_matrix(companion_size - variables)
         end if
         select case (scenario_value)
         case (1)
            impact = identity_matrix(variables)
         case (2)
            call cholesky_lower(draws%covariance(:, :, time_value, draw), &
               impact, info)
            if (info /= 0) then
               out%info = 2
               return
            end if
         case (3)
            contemporaneous = bvarsv_contemporaneous_matrix( &
               draws%contemporaneous(:, time_value, draw), variables)
            call invert_matrix(contemporaneous, inverse, info)
            if (info /= 0) then
               out%info = 3
               return
            end if
            allocate(scale(variables, variables))
            scale = 0.0_dp
            do item = 1, variables
               scale(item, item) = mean_standard_deviation(item)
            end do
            impact = matmul(inverse, scale)
            deallocate(scale)
         end select
         power = identity_matrix(companion_size)
         do step = 0, horizon
            out%response(:, :, step + 1, draw) = &
               matmul(power(:variables, :variables), impact)
            power = matmul(power, companion)
         end do
         deallocate(companion)
      end do
      out%selected_time = time_value
      out%scenario = scenario_value
   end function bvarsv_irf

   pure function bvarsv_predictive_density(prediction, points, variable, &
      horizon, cumulative) result(value)
      !! Evaluate the draw-averaged Gaussian predictive PDF or CDF.
      type(bvarsv_predictive_t), intent(in) :: prediction !! Draw-wise predictive moments.
      real(dp), intent(in) :: points(:) !! Evaluation points.
      integer, intent(in) :: variable !! One-based response-variable index.
      integer, intent(in) :: horizon !! One-based forecast-horizon index.
      logical, intent(in), optional :: cumulative !! Return CDF rather than PDF values.
      real(dp), allocatable :: value(:)
      real(dp) :: mean, variance, standardized, contribution
      integer :: point, draw, draws
      logical :: use_cdf

      use_cdf = .false.
      if (present(cumulative)) use_cdf = cumulative
      if (prediction%info /= 0 .or. .not. allocated(prediction%mean) .or. &
         .not. allocated(prediction%covariance) .or. variable < 1 .or. &
         variable > size(prediction%mean, 1) .or. horizon < 1 .or. &
         horizon > size(prediction%mean, 2) .or. &
         size(prediction%mean, 3) < 1 .or. &
         variable > size(prediction%covariance, 1) .or. &
         variable > size(prediction%covariance, 2) .or. &
         horizon > size(prediction%covariance, 3) .or. &
         size(prediction%covariance, 4) /= size(prediction%mean, 3)) then
         allocate(value(0))
         return
      end if
      draws = size(prediction%mean, 3)
      allocate(value(size(points)))
      value = 0.0_dp
      do draw = 1, draws
         mean = prediction%mean(variable, horizon, draw)
         variance = prediction%covariance(variable, variable, horizon, draw)
         if (variance <= 0.0_dp) then
            value = 0.0_dp
            return
         end if
         do point = 1, size(points)
            standardized = (points(point) - mean)/sqrt(variance)
            if (use_cdf) then
               contribution = 0.5_dp*(1.0_dp + erf(standardized/sqrt(2.0_dp)))
            else
               contribution = exp(-0.5_dp*standardized**2)/ &
                  sqrt(2.0_dp*acos(-1.0_dp)*variance)
            end if
            value(point) = value(point) + contribution
         end do
      end do
      value = value/real(draws, dp)
   end function bvarsv_predictive_density

   function bvarsv_carter_kohn(response, design, observation_covariance, &
      state_covariance, initial_mean, initial_covariance) result(out)
      !! Draw random-walk states by the Carter-Kohn forward-filter backward-sampler.
      real(dp), intent(in) :: response(:, :) !! Observation vector by time.
      real(dp), intent(in) :: design(:, :, :) !! Observation loading by time.
      real(dp), intent(in) :: observation_covariance(:, :, :) !! Observation covariance by time.
      real(dp), intent(in) :: state_covariance(:, :) !! Random-walk state innovation covariance.
      real(dp), intent(in) :: initial_mean(:) !! Prior mean of the initial state.
      real(dp), intent(in) :: initial_covariance(:, :) !! Prior covariance of the initial state.
      type(bvarsv_state_draw_t) :: out
      real(dp), allocatable :: predicted_mean(:), filtered_mean(:, :), filtered_covariance(:, :, :)
      real(dp), allocatable :: predicted_covariance(:, :), innovation(:), innovation_covariance(:, :)
      real(dp), allocatable :: inverse_innovation(:, :), normal(:), backward_covariance(:, :)
      real(dp), allocatable :: backward_mean(:), factor(:, :), inverse_prediction(:, :)
      integer :: observations, measurement_count, states, time, info

      measurement_count = size(response, 1)
      observations = size(response, 2)
      states = size(initial_mean)
      if (measurement_count < 1 .or. observations < 1 .or. states < 1 .or. &
         any(shape(design) /= [measurement_count, states, observations]) .or. &
         any(shape(observation_covariance) /= &
         [measurement_count, measurement_count, observations]) .or. &
         any(shape(state_covariance) /= [states, states]) .or. &
         any(shape(initial_covariance) /= [states, states])) then
         out%info = 1
         return
      end if
      allocate(filtered_mean(states, observations))
      allocate(filtered_covariance(states, states, observations))
      predicted_mean = initial_mean
      predicted_covariance = initial_covariance
      out%log_likelihood = 0.0_dp
      do time = 1, observations
         innovation = response(:, time) - matmul(design(:, :, time), predicted_mean)
         innovation_covariance = observation_covariance(:, :, time) + matmul( &
            design(:, :, time), matmul(predicted_covariance, &
            transpose(design(:, :, time))))
         call invert_matrix(innovation_covariance, inverse_innovation, info)
         if (info /= 0) then
            out%info = 2
            return
         end if
         call cholesky_lower(innovation_covariance, factor, info)
         if (info /= 0) then
            out%info = 3
            return
         end if
         out%log_likelihood = out%log_likelihood + 2.0_dp*sum(log(diagonal(factor))) + &
            dot_product(innovation, matmul(inverse_innovation, innovation))
         filtered_mean(:, time) = predicted_mean + matmul(matmul(predicted_covariance, &
            transpose(design(:, :, time))), matmul(inverse_innovation, innovation))
         filtered_covariance(:, :, time) = predicted_covariance - matmul( &
            matmul(matmul(predicted_covariance, transpose(design(:, :, time))), &
            inverse_innovation), matmul(design(:, :, time), predicted_covariance))
         filtered_covariance(:, :, time) = 0.5_dp*(filtered_covariance(:, :, time) + &
            transpose(filtered_covariance(:, :, time)))
         if (time < observations) then
            predicted_mean = filtered_mean(:, time)
            predicted_covariance = filtered_covariance(:, :, time) + state_covariance
         end if
      end do
      allocate(out%state(states, observations), normal(states))
      do time = 1, states
         normal(time) = random_standard_normal()
      end do
      call draw_normal(filtered_mean(:, observations), &
         filtered_covariance(:, :, observations), normal, out%state(:, observations), info)
      if (info /= 0) then
         out%info = 4
         return
      end if
      do time = observations - 1, 1, -1
         predicted_covariance = filtered_covariance(:, :, time) + state_covariance
         call invert_matrix(predicted_covariance, inverse_prediction, info)
         if (info /= 0) then
            out%info = 5
            return
         end if
         backward_mean = filtered_mean(:, time) + matmul(matmul(filtered_covariance(:, :, time), &
            inverse_prediction), out%state(:, time + 1) - filtered_mean(:, time))
         backward_covariance = filtered_covariance(:, :, time) - matmul( &
            matmul(filtered_covariance(:, :, time), inverse_prediction), &
            filtered_covariance(:, :, time))
         backward_covariance = 0.5_dp*(backward_covariance + transpose(backward_covariance))
         do info = 1, states
            normal(info) = random_standard_normal()
         end do
         call draw_normal(backward_mean, backward_covariance, normal, out%state(:, time), info)
         if (info /= 0) then
            out%info = 6
            return
         end if
      end do
   end function bvarsv_carter_kohn

   pure function bvarsv_contemporaneous_matrix(alpha, variables) result(matrix)
      !! Reconstruct bvarsv's unit lower-triangular contemporaneous matrix.
      real(dp), intent(in) :: alpha(:) !! Packed strict-lower-triangular coefficients.
      integer, intent(in) :: variables !! Positive system dimension.
      real(dp), allocatable :: matrix(:, :)
      integer :: equation, source, position

      if (variables < 1 .or. size(alpha) /= variables*(variables - 1)/2) then
         allocate(matrix(0, 0))
         return
      end if
      matrix = identity_matrix(variables)
      position = 0
      do equation = 2, variables
         do source = 1, equation - 1
            position = position + 1
            matrix(equation, source) = alpha(position)
         end do
      end do
   end function bvarsv_contemporaneous_matrix

   pure function bvarsv_covariance_from_state(alpha, log_variance) result(covariance)
      !! Form the reduced-form covariance from triangular and log-variance states.
      real(dp), intent(in) :: alpha(:) !! Packed strict-lower-triangular coefficients.
      real(dp), intent(in) :: log_variance(:) !! Orthogonal-shock log variances.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: contemporaneous(:, :), impact(:, :), scale(:, :)
      integer :: variables, item, info

      variables = size(log_variance)
      contemporaneous = bvarsv_contemporaneous_matrix(alpha, variables)
      if (size(contemporaneous, 1) == 0) then
         allocate(covariance(0, 0))
         return
      end if
      call invert_matrix(contemporaneous, impact, info)
      if (info /= 0) then
         allocate(covariance(0, 0))
         return
      end if
      allocate(scale(variables, variables))
      scale = 0.0_dp
      do item = 1, variables
         scale(item, item) = exp(0.5_dp*log_variance(item))
      end do
      impact = matmul(impact, scale)
      covariance = matmul(impact, transpose(impact))
   end function bvarsv_covariance_from_state

   function bvarsv_log_volatility_update(structural_residual, &
      current_log_variance, state_covariance, initial_mean, &
      initial_covariance) result(out)
      !! Draw Primiceri random-walk log volatilities with the seven-component KSC mixture.
      real(dp), intent(in) :: structural_residual(:, :) !! Orthogonal structural residuals by time.
      real(dp), intent(in) :: current_log_variance(:, :) !! Current log-variance states by time.
      real(dp), intent(in) :: state_covariance(:, :) !! Log-volatility random-walk covariance.
      real(dp), intent(in) :: initial_mean(:) !! Prior mean of the initial log variances.
      real(dp), intent(in) :: initial_covariance(:, :) !! Prior covariance of initial log variances.
      type(bvarsv_volatility_draw_t) :: out
      real(dp), parameter :: probability(7) = [0.00730_dp, 0.10556_dp, &
         0.00002_dp, 0.04395_dp, 0.34001_dp, 0.24566_dp, 0.25750_dp]
      real(dp), parameter :: mixture_mean(7) = [-10.12999_dp, -3.97281_dp, &
         -8.56686_dp, 2.77786_dp, 0.61942_dp, 1.79518_dp, -1.08819_dp]
      real(dp), parameter :: mixture_variance(7) = [5.79596_dp, 2.61369_dp, &
         5.17950_dp, 0.16735_dp, 0.64009_dp, 0.34023_dp, 1.26261_dp]
      real(dp), allocatable :: transformed(:, :), adjusted(:, :)
      real(dp), allocatable :: design(:, :, :), observation_covariance(:, :, :)
      real(dp) :: weight(7), total_weight, threshold, cumulative
      type(bvarsv_state_draw_t) :: state_draw
      integer :: variables, observations, variable, time, component

      variables = size(structural_residual, 1)
      observations = size(structural_residual, 2)
      if (variables < 1 .or. observations < 1 .or. &
         any(shape(current_log_variance) /= [variables, observations]) .or. &
         any(shape(state_covariance) /= [variables, variables]) .or. &
         size(initial_mean) /= variables .or. &
         any(shape(initial_covariance) /= [variables, variables])) then
         out%info = 1
         return
      end if
      allocate(out%component(variables, observations))
      allocate(transformed(variables, observations), adjusted(variables, observations))
      transformed = log(0.001_dp + structural_residual**2)
      do time = 1, observations
         do variable = 1, variables
            do component = 1, 7
               weight(component) = probability(component)/ &
                  sqrt(2.0_dp*acos(-1.0_dp)*mixture_variance(component))*exp( &
                  -0.5_dp*(transformed(variable, time) - &
                  current_log_variance(variable, time) - mixture_mean(component) + &
                  1.2704_dp)**2/mixture_variance(component))
            end do
            total_weight = sum(weight)
            if (total_weight <= tiny(1.0_dp)) then
               out%component(variable, time) = maxloc(weight, dim=1)
            else
               threshold = random_uniform()*total_weight
               cumulative = 0.0_dp
               out%component(variable, time) = 7
               do component = 1, 7
                  cumulative = cumulative + weight(component)
                  if (threshold <= cumulative) then
                     out%component(variable, time) = component
                     exit
                  end if
               end do
            end if
            component = out%component(variable, time)
            adjusted(variable, time) = transformed(variable, time) - &
               mixture_mean(component) + 1.2704_dp
         end do
      end do
      allocate(design(variables, variables, observations))
      allocate(observation_covariance(variables, variables, observations))
      design = 0.0_dp
      observation_covariance = 0.0_dp
      do time = 1, observations
         do variable = 1, variables
            design(variable, variable, time) = 1.0_dp
            component = out%component(variable, time)
            observation_covariance(variable, variable, time) = &
               mixture_variance(component)
         end do
      end do
      state_draw = bvarsv_carter_kohn(adjusted, design, &
         observation_covariance, state_covariance, initial_mean, initial_covariance)
      if (state_draw%info /= 0) then
         out%info = 10 + state_draw%info
         return
      end if
      out%log_variance = state_draw%state
      out%standard_deviation = exp(0.5_dp*out%log_variance)
   end function bvarsv_log_volatility_update

   function bvarsv_contemporaneous_update(reduced_residual, log_variance, &
      state_covariance, initial_mean, initial_covariance) result(out)
      !! Draw equation-block random-walk contemporaneous coefficients.
      real(dp), intent(in) :: reduced_residual(:, :) !! Reduced-form residuals by time.
      real(dp), intent(in) :: log_variance(:, :) !! Orthogonal-shock log variances by time.
      real(dp), intent(in) :: state_covariance(:, :) !! Block-diagonal alpha-state covariance.
      real(dp), intent(in) :: initial_mean(:) !! Prior mean of packed initial alpha states.
      real(dp), intent(in) :: initial_covariance(:, :) !! Prior covariance of initial alpha states.
      type(bvarsv_contemporaneous_draw_t) :: out
      real(dp), allocatable :: response(:, :), design(:, :, :)
      real(dp), allocatable :: observation_covariance(:, :, :)
      type(bvarsv_state_draw_t) :: state_draw
      integer :: variables, observations, count, equation, block_size
      integer :: first, last, time

      variables = size(reduced_residual, 1)
      observations = size(reduced_residual, 2)
      count = variables*(variables - 1)/2
      if (variables < 2 .or. observations < 1 .or. &
         any(shape(log_variance) /= [variables, observations]) .or. &
         any(shape(state_covariance) /= [count, count]) .or. &
         size(initial_mean) /= count .or. &
         any(shape(initial_covariance) /= [count, count])) then
         out%info = 1
         return
      end if
      allocate(out%alpha(count, observations))
      first = 1
      do equation = 2, variables
         block_size = equation - 1
         last = first + block_size - 1
         allocate(response(1, observations), design(1, block_size, observations))
         allocate(observation_covariance(1, 1, observations))
         response(1, :) = reduced_residual(equation, :)
         do time = 1, observations
            design(1, :, time) = -reduced_residual(:equation - 1, time)
            observation_covariance(1, 1, time) = exp(log_variance(equation, time))
         end do
         state_draw = bvarsv_carter_kohn(response, design, &
            observation_covariance, state_covariance(first:last, first:last), &
            initial_mean(first:last), initial_covariance(first:last, first:last))
         if (state_draw%info /= 0) then
            out%info = 10 + state_draw%info
            return
         end if
         out%alpha(first:last, :) = state_draw%state
         first = last + 1
         deallocate(response, design, observation_covariance)
      end do
   end function bvarsv_contemporaneous_update

   function bvarsv_random_walk_covariance(state, prior_scale, prior_degrees, &
      degrees_offset) result(out)
      !! Draw a random-walk innovation covariance from its inverse-Wishart posterior.
      real(dp), intent(in) :: state(:, :) !! State vector by time.
      real(dp), intent(in) :: prior_scale(:, :) !! Inverse-Wishart prior scale.
      real(dp), intent(in) :: prior_degrees !! Inverse-Wishart prior degrees of freedom.
      integer, intent(in), optional :: degrees_offset !! Package-specific posterior degree adjustment.
      type(bvarsv_covariance_draw_t) :: out
      real(dp), allocatable :: difference(:, :), posterior_scale(:, :)
      real(dp), allocatable :: precision(:, :)
      integer :: states, observations, offset, info

      states = size(state, 1)
      observations = size(state, 2)
      offset = 0
      if (present(degrees_offset)) offset = degrees_offset
      if (states < 1 .or. observations < 2 .or. &
         any(shape(prior_scale) /= [states, states]) .or. &
         prior_degrees + real(observations - 1 + offset, dp) <= &
         real(states - 1, dp)) then
         out%info = 1
         return
      end if
      difference = state(:, 2:) - state(:, :observations - 1)
      posterior_scale = prior_scale + matmul(difference, transpose(difference))
      posterior_scale = 0.5_dp*(posterior_scale + transpose(posterior_scale))
      call inverse_wishart_draw(posterior_scale, prior_degrees + &
         real(observations - 1 + offset, dp), out%covariance, precision, info)
      if (info /= 0) out%info = 2
   end function bvarsv_random_walk_covariance

   function bvarsv_contemporaneous_covariance_update(alpha, prior_scale, &
      prior_degrees) result(out)
      !! Draw the block-diagonal contemporaneous-state innovation covariance.
      real(dp), intent(in) :: alpha(:, :) !! Packed contemporaneous states by time.
      real(dp), intent(in) :: prior_scale(:, :) !! Block-diagonal inverse-Wishart prior scale.
      real(dp), intent(in) :: prior_degrees(:) !! Equation-block prior degrees of freedom.
      type(bvarsv_covariance_draw_t) :: out
      type(bvarsv_covariance_draw_t) :: block
      integer :: count, variables, equation, block_size, first, last

      count = size(alpha, 1)
      variables = size(prior_degrees) + 1
      if (count /= variables*(variables - 1)/2 .or. size(alpha, 2) < 2 .or. &
         any(shape(prior_scale) /= [count, count])) then
         out%info = 1
         return
      end if
      allocate(out%covariance(count, count))
      out%covariance = 0.0_dp
      first = 1
      do equation = 2, variables
         block_size = equation - 1
         last = first + block_size - 1
         block = bvarsv_random_walk_covariance(alpha(first:last, :), &
            prior_scale(first:last, first:last), prior_degrees(equation - 1))
         if (block%info /= 0) then
            out%info = 10 + block%info
            return
         end if
         out%covariance(first:last, first:last) = block%covariance
         first = last + 1
      end do
   end function bvarsv_contemporaneous_covariance_update

   subroutine draw_normal(mean, covariance, normal, draw, info)
      !! Transform supplied standard normals to a multivariate Gaussian draw.
      real(dp), intent(in) :: mean(:) !! Gaussian mean vector.
      real(dp), intent(in) :: covariance(:, :) !! Positive-definite covariance matrix.
      real(dp), intent(in) :: normal(:) !! Independent standard-normal variates.
      real(dp), intent(out) :: draw(:) !! Generated Gaussian draw.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: factor(:, :)

      call cholesky_lower(covariance, factor, info)
      if (info /= 0) return
      draw = mean + matmul(factor, normal)
   end subroutine draw_normal

   subroutine fill_standard_normal(values)
      !! Fill a vector with independent standard-normal variates.
      real(dp), intent(out) :: values(:) !! Generated standard-normal variates.
      integer :: item

      do item = 1, size(values)
         values(item) = random_standard_normal()
      end do
   end subroutine fill_standard_normal

   subroutine inverse_wishart_draw(scale, degrees_of_freedom, covariance, &
      precision, info)
      !! Draw an inverse-Wishart matrix by Bartlett decomposition.
      real(dp), intent(in) :: scale(:, :) !! Inverse-Wishart scale matrix.
      real(dp), intent(in) :: degrees_of_freedom !! Positive degrees of freedom.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Drawn covariance matrix.
      real(dp), allocatable, intent(out) :: precision(:, :) !! Drawn precision matrix.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: inverse_scale(:, :), lower(:, :)
      real(dp), allocatable :: bartlett(:, :), factor(:, :)
      integer :: variables, row, column

      variables = size(scale, 1)
      call invert_matrix(scale, inverse_scale, info)
      if (info /= 0) return
      call cholesky_lower(inverse_scale, lower, info)
      if (info /= 0) return
      allocate(bartlett(variables, variables))
      bartlett = 0.0_dp
      do row = 1, variables
         bartlett(row, row) = sqrt(2.0_dp*random_gamma(0.5_dp* &
            (degrees_of_freedom - real(row, dp) + 1.0_dp)))
         do column = 1, row - 1
            bartlett(row, column) = random_standard_normal()
         end do
      end do
      factor = matmul(lower, bartlett)
      precision = matmul(factor, transpose(factor))
      call invert_matrix(precision, covariance, info)
   end subroutine inverse_wishart_draw

   pure function diagonal(matrix) result(values)
      !! Extract the diagonal of a square matrix.
      real(dp), intent(in) :: matrix(:, :) !! Square input matrix.
      real(dp) :: values(min(size(matrix, 1), size(matrix, 2)))
      integer :: item

      do item = 1, size(values)
         values(item) = matrix(item, item)
      end do
   end function diagonal

   pure function pack_strict_lower(matrix) result(values)
      !! Pack strict-lower-triangular entries in equation-major order.
      real(dp), intent(in) :: matrix(:, :) !! Square input matrix.
      real(dp), allocatable :: values(:)
      integer :: variables, equation, source, position

      variables = size(matrix, 1)
      if (size(matrix, 2) /= variables) then
         allocate(values(0))
         return
      end if
      allocate(values(variables*(variables - 1)/2))
      position = 0
      do equation = 2, variables
         do source = 1, equation - 1
            position = position + 1
            values(position) = matrix(equation, source)
         end do
      end do
   end function pack_strict_lower

end module bvarsv_mod
