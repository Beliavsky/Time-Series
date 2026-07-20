cmake_minimum_required(VERSION 3.20)
project(time_series LANGUAGES Fortran)

add_library(time_series STATIC
   kind.f90
   linalg.f90
   stats.f90
   spline.f90
   mixture.f90
   utils.f90
   calendar.f90
   fourier.f90
   special_functions.f90
   polynomial.f90
   random.f90
   distribution.f90
   resampling.f90
   optimization.f90
   distribution_fit.f90
   neural_network.f90
   recurrent_network.f90
   time_series_stats.f90
   rolling_forecast.f90
   time_series_var_utils.f90
   forecast.f90
   kfas.f90
   urca.f90
   astsa.f90
   itsmr.f90
   time_series_diagnostics.f90
   arima2.f90
   arfima.f90
   fracdiff.f90
   fracdist_tables.f90
   fracdist.f90
   nsarfima.f90
   garma.f90
   esemifar.f90
   tfarima.f90
   robustarima.f90
   mts.f90
   bayesforecast.f90
   bsts.f90
   bssm.f90
   tsss.f90
   tsissm.f90
   mar.f90
   fcvar.f90
   bigtime.f90
   bigvar.f90
   varshrink.f90
   vars.f90
   var_etp.f90
   bvartools.f90
   bvar.f90
   bvars.f90
   bvarsv.f90
   gmvarkit.f90
   ugmar.f90
   mswm.f90
   nts.f90
   tsdyn.f90
   expar.f90
   exparma.f90
   bentcablear.f90
   baystar.f90
   mixar.f90
   setartree.f90
   tseriestarma_tables.f90
   tseriestarma.f90
   gmdh.f90
   nnfor.f90
   narfima.f90
   nlints.f90
   tslstmplus.f90
   tsann.f90
   echos.f90
   starvars.f90
   rugarch.f90
   rugarch_extensions.f90
   rugarch_diagnostics.f90
   arma.f90
   long_memory.f90
   state_space.f90
   multivariate.f90
   volatility.f90
   spectral.f90
   diagnostics.f90
   forecasting.f90
   regression_time_series.f90
   bayesian_time_series.f90
   markov_switching.f90
   count_time_series.f90
   functional_time_series.f90
)

set_target_properties(time_series PROPERTIES
   Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

target_include_directories(time_series PUBLIC
   "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>"
)

add_executable(example_var2_fit examples/var2_fit.f90)
target_link_libraries(example_var2_fit PRIVATE time_series)

include(CTest)
if(BUILD_TESTING)
   foreach(test_name IN ITEMS
      forecast kfas urca astsa random mixture calendar spline neural_network itsmr arima2 arfima fracdiff fracdist
      nsarfima garma esemifar tfarima robustarima mts bayesforecast bsts bssm
      tsss tsissm mar fcvar
      bigtime bigvar varshrink vars var_etp bvartools bvar bvars bvarsv gmvarkit ugmar
      mswm nts tsdyn expar exparma bentcablear baystar mixar setartree tseriestarma gmdh nnfor narfima nlints tslstmplus tsann echos starvars distribution rugarch facades
   )
      add_executable(test_${test_name} test_${test_name}.f90)
      target_link_libraries(test_${test_name} PRIVATE time_series)
      add_test(NAME ${test_name} COMMAND test_${test_name})
   endforeach()
endif()
