! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Example simulation and estimation of a Gaussian VAR(2).
program var2_fit
   use kind_mod, only: dp
   use multivariate_mod, only: mar_simulation_t, vars_fit_t, mar_simulate, &
      vars_fit, display
   use random_mod, only: set_random_seed
   implicit none

   integer, parameter :: variables = 2
   integer, parameter :: order = 2
   integer, parameter :: observations = 100000
   integer, parameter :: burnin = 500
   real(dp) :: intercept(variables)
   real(dp) :: ar(variables, variables, order)
   real(dp) :: innovation_covariance(variables, variables)
   type(mar_simulation_t) :: simulation
   type(vars_fit_t) :: fit

   intercept = [0.50_dp, -0.25_dp]
   ar(:, :, 1) = reshape([ &
      0.55_dp, -0.15_dp, &
      0.10_dp, 0.40_dp], [variables, variables])
   ar(:, :, 2) = reshape([ &
      -0.20_dp, 0.10_dp, &
      0.05_dp, -0.15_dp], [variables, variables])
   innovation_covariance = reshape([ &
      0.40_dp, 0.12_dp, &
      0.12_dp, 0.30_dp], [variables, variables])

   call set_random_seed(20260717)
   simulation = mar_simulate(intercept, ar, innovation_covariance, &
      observations, burnin)
   call display(simulation)
   fit = vars_fit(simulation%series, order)
   call display(fit)

end program var2_fit
