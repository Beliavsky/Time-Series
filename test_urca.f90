! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the urca translation.
program test_urca
   use kind_mod, only: dp
   use urca_mod
   implicit none
   real(dp) :: y(40)
   real(dp) :: system(100, 3), exogenous(100, 2), common, u, v
   integer :: i
   type(adf_result_t) :: adf
   type(kpss_result_t) :: kpss
   type(pp_result_t) :: pp
   type(ers_result_t) :: ers
   type(za_result_t) :: za
   type(johansen_result_t) :: johansen

   y(1) = 1.0_dp
   do i = 2, size(y)
      y(i) = 0.72_dp*y(i - 1) + sin(0.31_dp*real(i*i, dp)) &
         + 0.15_dp*cos(1.7_dp*real(i, dp))
   end do
   adf = adf_test(y, 'drift', 2, 'Fixed')
   call check(adf%info == 0, 'ADF status')
   call check(adf%lags == 2, 'ADF lag')
   call check(size(adf%statistic) == 2, 'ADF statistics')
   call check(maxval(abs(adf%statistic - [-2.675507648496013_dp, 3.6557555478186035_dp])) &
      < 1.0e-11_dp, 'ADF R reference')
   call check(maxval(abs(adf%coefficients - [-0.05404722786005757_dp, -0.43843248260017093_dp, &
      -0.037442689619825262_dp, -0.079727993736929709_dp])) < 1.0e-11_dp, &
      'ADF coefficient reference')
   call check(maxval(abs(adf%critical_values(1, :) - [-3.58_dp, -2.93_dp, -2.60_dp])) &
      < epsilon(1.0_dp), &
      'ADF critical values')
   kpss = kpss_test(y, 'mu', use_lag=2)
   call check(kpss%info == 0, 'KPSS status')
   call check(kpss%lags == 2, 'KPSS lag')
   call check(kpss%statistic > 0.0_dp, 'KPSS statistic')
   call check(abs(kpss%statistic - 0.47122363688193619_dp) < 1.0e-12_dp, 'KPSS R reference')
   adf = adf_test(y, 'trend', 4, 'AIC')
   call check(maxval(abs(adf%statistic - [-2.8610097193189885_dp, 2.9049001187072441_dp, &
      4.3358195792238625_dp])) < 1.0e-10_dp, 'ADF AIC R reference')
   pp = pp_test(y, 'Z-tau', 'constant', use_lag=2)
   call check(pp%info == 0, 'PP constant status')
   call check(abs(pp%statistic + 3.0102180942304528_dp) < 1.0e-11_dp, 'PP constant Z-tau reference')
   call check(abs(pp%auxiliary(1) + 0.0985_dp) < 1.0e-12_dp, 'PP constant auxiliary')
   call check(maxval(abs(pp%coefficients - [-0.0099606840043021585_dp, &
      0.60254541146044804_dp])) < 1.0e-12_dp, 'PP constant coefficients')
   call check(maxval(abs(pp%critical_values - [-3.6065512820512819_dp, -2.9378015121630505_dp, &
      -2.606917225509533_dp])) < 1.0e-12_dp, 'PP constant critical values')
   pp = pp_test(y, 'Z-alpha', 'trend', use_lag=2)
   call check(pp%info == 0, 'PP trend status')
   call check(abs(pp%statistic + 16.801468143573363_dp) < 1.0e-10_dp, 'PP trend Z-alpha reference')
   call check(maxval(abs(pp%auxiliary - [-0.0340_dp, -0.8553_dp])) < 1.0e-12_dp, &
      'PP trend auxiliary')
   call check(maxval(abs(pp%coefficients - [-0.0038015406547643447_dp, 0.55269792897163794_dp, &
      -0.0094083278318010458_dp])) < 1.0e-12_dp, 'PP trend coefficients')
   ers = ers_test(y, 'DF-GLS', 'constant', 2)
   call check(ers%info == 0, 'ERS DF-GLS constant status')
   call check(abs(ers%statistic + 2.020721266758525_dp) < 1.0e-11_dp, &
      'ERS DF-GLS constant reference')
   ers = ers_test(y, 'P-test', 'constant', 2)
   call check(ers%info == 0, 'ERS P constant status')
   call check(abs(ers%statistic - 2.9173918055766186_dp) < 1.0e-10_dp, 'ERS P constant reference')
   call check(ers%lags == 1, 'ERS P selected lag')
   ers = ers_test(y, 'DF-GLS', 'trend', 2)
   call check(abs(ers%statistic + 2.4262123909401168_dp) < 1.0e-11_dp, 'ERS DF-GLS trend reference')
   ers = ers_test(y, 'P-test', 'trend', 2)
   call check(abs(ers%statistic - 6.4329373893768942_dp) < 1.0e-10_dp, 'ERS P trend reference')

   za = za_test(y, 'intercept', 2)
   call check(za%info == 0, 'ZA intercept status')
   call check(za%break_point == 18, 'ZA intercept break')
   call check(abs(za%statistic + 4.5118945477153689_dp) < 1.0e-10_dp, 'ZA intercept reference')
   za = za_test(y, 'trend', 2)
   call check(za%break_point == 22, 'ZA trend break')
   call check(abs(za%statistic + 3.7397281039922059_dp) < 1.0e-10_dp, 'ZA trend reference')
   za = za_test(y, 'both', 2)
   call check(za%break_point == 19, 'ZA both break')
   call check(abs(za%statistic + 5.6939359999625703_dp) < 1.0e-10_dp, 'ZA both reference')
   call check(maxval(abs(za%coefficients - [1.0862196037876592_dp, -0.34564149866073196_dp, &
      -0.043174571214751045_dp, 0.36816193205127384_dp, 0.15316280462749859_dp, &
      -2.0384820045345324_dp, 0.15124090554221728_dp])) < 1.0e-9_dp, &
      'ZA both coefficients')

   common = 0.0_dp
   u = 0.0_dp
   v = 0.0_dp
   do i = 1, size(system, 1)
      common = common + 0.15_dp + sin(0.17_dp*real(i*i, dp))
      u = 0.45_dp*u + cos(0.31_dp*real(i*i, dp))
      v = 0.30_dp*v + sin(0.23_dp*real(i*i, dp))
      system(i, :) = [common, common + u, 0.5_dp*common + v]
      exogenous(i, 1) = sin(0.07_dp*real(i*i, dp))
      exogenous(i, 2) = merge(1.0_dp, 0.0_dp, i > 55)
   end do
   johansen = johansen_test(system, 'trace', 'none', 3, 'longrun')
   call check(johansen%info == 0, 'Johansen none status')
   call check(maxval(abs(johansen%eigenvalues - [0.2792715040519849_dp, 0.1247053499623679_dp, &
      0.007042824150148524_dp])) < 1.0e-10_dp, 'Johansen none eigenvalues')
   call check(maxval(abs(johansen%statistic - [0.6855709642314909_dp, 13.60545748595657_dp, &
      45.37225700712754_dp])) < 1.0e-9_dp, 'Johansen trace statistics')
   johansen = johansen_test(system, 'eigen', 'const', 3, 'transitory')
   call check(johansen%info == 0, 'Johansen const status')
   call check(maxval(abs(johansen%eigenvalues - [0.2793407294320964_dp, 0.1257730628006152_dp, &
      0.0895194813653945_dp, 0.0_dp])) < 1.0e-10_dp, 'Johansen const eigenvalues')
   call check(maxval(abs(johansen%statistic - [9.096929308956033_dp, 13.03828249800063_dp, &
      31.77611673896896_dp])) < 1.0e-9_dp, 'Johansen eigen statistics')
   johansen = johansen_test(system, 'trace', 'trend', 3, 'longrun')
   call check(johansen%info == 0, 'Johansen trend status')
   call check(maxval(abs(johansen%eigenvalues - [0.2826069260699755_dp, 0.1339426274729635_dp, &
      0.0596197491929899_dp, 0.0_dp])) < 1.0e-10_dp, 'Johansen trend eigenvalues')
   johansen = johansen_test(system, 'trace', 'const', 3, 'longrun')
   call check(maxval(abs(johansen%pi - reshape([ &
      -0.03760541609614241_dp, 0.4545521595854977_dp, 0.1353985384113341_dp, &
       0.05592298593865674_dp, -0.4968291322701125_dp, 0.2054218631267599_dp, &
      -0.01604815888594499_dp, 0.08509829362434387_dp, -0.647874764206939_dp, &
       0.1356356776321014_dp, 0.2575287803247119_dp, -0.07176132168331739_dp], [3, 4]))) &
      < 1.0e-9_dp, 'Johansen PI reference')
   johansen = johansen_test(system, 'trace', 'const', 3, 'longrun', season=4, exogenous=exogenous)
   call check(johansen%info == 0, 'Johansen seasonal/exogenous status')
   call check(maxval(abs(johansen%eigenvalues - [0.2853745581216316_dp, 0.1504763185970229_dp, &
      0.02939931296478148_dp, 0.0_dp])) < 1.0e-10_dp, &
      'Johansen seasonal/exogenous eigenvalues')
   call check(maxval(abs(johansen%statistic - [2.894493014042241_dp, 18.71320076852811_dp, &
      51.30488368156022_dp])) < 1.0e-9_dp, 'Johansen seasonal/exogenous statistics')
   johansen = johansen_test(system, 'trace', 'none', 3, 'longrun', season=4, exogenous=exogenous)
   call check(maxval(abs(johansen%eigenvalues - [0.2845337602079497_dp, 0.147789569534698_dp, &
      0.0001019879350748765_dp])) < 1.0e-10_dp, 'Johansen dummy none eigenvalues')
   johansen = johansen_test(system, 'trace', 'trend', 3, 'longrun', season=4, exogenous=exogenous)
   call check(maxval(abs(johansen%eigenvalues - [0.2885802588906924_dp, 0.1973289026969577_dp, &
      0.05007037283172314_dp, 0.0_dp])) < 1.0e-10_dp, 'Johansen dummy trend eigenvalues')
   print '(a)', 'All urca_mod tests passed.'

contains

   subroutine check(ok, name)
      ! Stop the test program when a named assertion fails.
      logical, intent(in) :: ok
      character(len=*), intent(in) :: name
      if (.not. ok) then
         print '(a)', 'FAILED: '//name
         error stop 1
      end if
   end subroutine check
end program
