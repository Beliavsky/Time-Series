! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the R GMDH package translation.
program test_gmdh
   use kind_mod, only: dp
   use gmdh_mod
   implicit none

   type(gmdh_ridge_fit_t) :: ridge
   type(gmdh_model_t) :: classic, revised, invalid_model
   type(gmdh_forecast_t) :: forecast, invalid_forecast
   real(dp) :: predictor(20, 2), response(20), series(100)
   real(dp), allocatable :: prediction(:), values(:)
   real(dp) :: transformed, restored
   integer :: index, transfer, scratch_unit, display_size

   do transfer = gmdh_transfer_polynomial, gmdh_transfer_tangent
      transformed = gmdh_transform(0.37_dp, transfer)
      restored = gmdh_inverse_transform(transformed, transfer)
      call check(abs(restored - 0.37_dp) < 1.0e-12_dp, &
         "transfer-function round trip")
   end do

   do index = 1, size(response)
      predictor(index, :) = [1.0_dp, real(index, dp)/10.0_dp]
      response(index) = 1.5_dp + 2.25_dp*predictor(index, 2)
   end do
   ridge = gmdh_ridge_fit(predictor, response, [0.0_dp, 0.1_dp, 1.0_dp], &
      0.7_dp)
   call check(ridge%info == 0 .and. abs(ridge%lambda) < 1.0e-14_dp .and. &
      maxval(abs(ridge%coefficients - [1.5_dp, 2.25_dp])) < 1.0e-10_dp, &
      "validation-selected ridge regression")

   series(1) = 0.31_dp
   do index = 2, size(series)
      series(index) = 0.2_dp + 2.4_dp*series(index - 1)* &
         (1.0_dp - series(index - 1))
   end do
   classic = gmdh_fit(series, input_count=4, layer_count=1, &
      method=gmdh_method_classic, transfers=[gmdh_transfer_polynomial], &
      lambdas=[0.0_dp])
   call check(classic%info == 0 .and. size(classic%layers) == 1 .and. &
      size(classic%layers(1)%nodes) == 1 .and. &
      size(classic%fitted) == 96 .and. &
      maxval(abs(classic%residuals)) < 2.0e-6_dp, &
      "classical GMDH quadratic-neuron fit")
   prediction = gmdh_predict(classic, series, 5)
   do index = 1, size(prediction)
      if (index == 1) then
         restored = 0.2_dp + 2.4_dp*series(100)*(1.0_dp - series(100))
      else
         restored = 0.2_dp + 2.4_dp*prediction(index - 1)* &
            (1.0_dp - prediction(index - 1))
      end if
      call check(abs(prediction(index) - restored) < 1.0e-7_dp, &
         "recursive classical GMDH prediction")
   end do

   revised = gmdh_fit(series, input_count=4, layer_count=2, &
      method=gmdh_method_revised, transfers=[gmdh_transfer_polynomial, &
      gmdh_transfer_sigmoid], lambdas=[0.0_dp, 0.01_dp])
   call check(revised%info == 0 .and. size(revised%layers) == 2 .and. &
      size(revised%layers(1)%nodes) == 4 .and. &
      size(revised%layers(2)%nodes) == 1 .and. &
      all(revised%layers(1)%nodes%node_type >= gmdh_node_quadratic) .and. &
      all(revised%layers(1)%nodes%node_type <= gmdh_node_feedback) .and. &
      all(abs(revised%residuals) < huge(1.0_dp)), &
      "revised GMDH layered fit")

   forecast = gmdh_forecast(series, horizon=3, input_count=4, layer_count=2, &
      method=gmdh_method_revised, level=90.0_dp, &
      transfers=[gmdh_transfer_polynomial], lambdas=[0.0_dp, 0.01_dp])
   call check(forecast%info == 0 .and. size(forecast%mean) == 3 .and. &
      all(forecast%standard_error >= 0.0_dp) .and. &
      all(forecast%lower <= forecast%mean) .and. &
      all(forecast%upper >= forecast%mean), &
      "rolling-origin GMDH prediction intervals")
   values = gmdh_coefficients(revised)
   call check(size(values) > 0 .and. &
      size(gmdh_fitted_values(revised)) == 96 .and. &
      size(gmdh_residuals(revised)) == 96, "GMDH accessors")

   invalid_model = gmdh_fit(series, input_count=1)
   invalid_forecast = gmdh_forecast(series, horizon=6)
   call check(invalid_model%info /= 0 .and. invalid_forecast%info /= 0, &
      "GMDH input validation")

   open(newunit=scratch_unit, status='scratch', action='write')
   call display(revised, scratch_unit)
   call display(forecast, scratch_unit)
   inquire(unit=scratch_unit, size=display_size)
   close(scratch_unit)
   call check(display_size > 0, "GMDH display methods")

   print '(a)', "GMDH tests passed"

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Test condition.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', "FAILED: "//message
         error stop 1
      end if
   end subroutine check

end program test_gmdh
