! SPDX-License-Identifier: MIT
! SPDX-FileComment: Original VAR data-layout infrastructure for this Fortran library.
module time_series_var_utils_mod
   use kind_mod, only: dp
   implicit none
   private

   public :: build_var_data, build_varx_data

contains

   pure subroutine build_var_data(series, lag_order, response, design)
      !! Construct lag-major VAR response and predictor matrices.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: lag_order !! Model lag order.
      real(dp), allocatable, intent(out) :: response(:, :) !! Response observations.
      real(dp), allocatable, intent(out) :: design(:, :) !! Design.
      integer :: observations, variables, lag

      observations = size(series, 1)
      variables = size(series, 2)
      allocate(response(observations - lag_order, variables))
      allocate(design(observations - lag_order, variables*lag_order))
      response = series(lag_order + 1:observations, :)
      do lag = 1, lag_order
         design(:, (lag - 1)*variables + 1:lag*variables) = &
            series(lag_order + 1 - lag:observations - lag, :)
      end do
   end subroutine build_var_data

   pure subroutine build_varx_data(endogenous, exogenous, ar_order, &
      exogenous_order, response, ar_design, x_design, contemporaneous)
      !! Construct aligned VARX data with an optional current exogenous block.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      real(dp), allocatable, intent(out) :: response(:, :) !! Response observations.
      real(dp), allocatable, intent(out) :: ar_design(:, :) !! Autoregressive design.
      real(dp), allocatable, intent(out) :: x_design(:, :) !! X design.
      logical, intent(in), optional :: contemporaneous !! Flag controlling contemporaneous.
      integer :: observations, variables, x_variables, leading, lag, block, x_blocks
      logical :: use_contemporaneous

      observations = size(endogenous, 1)
      variables = size(endogenous, 2)
      x_variables = size(exogenous, 2)
      use_contemporaneous = .false.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      leading = max(ar_order, exogenous_order)
      x_blocks = exogenous_order
      if (use_contemporaneous) x_blocks = x_blocks + 1
      allocate(response(observations - leading, variables))
      allocate(ar_design(observations - leading, variables*ar_order))
      allocate(x_design(observations - leading, x_variables*x_blocks))
      response = endogenous(leading + 1:observations, :)
      do lag = 1, ar_order
         ar_design(:, (lag - 1)*variables + 1:lag*variables) = &
            endogenous(leading + 1 - lag:observations - lag, :)
      end do
      if (use_contemporaneous) x_design(:, 1:x_variables) = exogenous(leading + 1:observations, :)
      do lag = 1, exogenous_order
         block = lag
         if (use_contemporaneous) block = block + 1
         x_design(:, (block - 1)*x_variables + 1:block*x_variables) = &
            exogenous(leading + 1 - lag:observations - lag, :)
      end do
   end subroutine build_varx_data

end module time_series_var_utils_mod
