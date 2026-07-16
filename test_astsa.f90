! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the astsa translation.
program test_astsa
   use kind_mod, only: dp
   use astsa_mod
   use time_series_random_mod, only: set_random_seed
   implicit none
   real(dp) :: y(6, 2), observation(2, 2, 1), mu0(2), sigma0(2, 2)
   real(dp) :: transition(2, 2), sq(2, 2), sr(2, 2), matrix(2, 2)
   real(dp) :: correlation(2, 2)
   real(dp), allocatable :: coefficients(:), product(:), powered(:, :)
   type(astsa_filter_t) :: filtered
   type(astsa_smoother_t) :: smoothed
   type(astsa_em_result_t) :: em
   type(astsa_matrix_power_result_t) :: power_result
   type(astsa_ffbs_result_t) :: ffbs
   type(astsa_spectrum_t) :: spectrum
   type(astsa_simulation_t) :: simulation
   type(astsa_arma_diagnostic_t) :: diagnostic
   type(astsa_prewhite_result_t) :: prewhitened
   type(astsa_sarima_likelihood_t) :: likelihood
   type(astsa_sarima_fit_t) :: sarima_estimate
   type(astsa_sarima_fit_t) :: fixed_estimate
   type(astsa_sarima_forecast_t) :: sarima_prediction
   type(astsa_sarima_diagnostics_t) :: sarima_diagnostic
   type(astsa_lag_regression_t) :: lag_fit
   type(astsa_signal_extraction_t) :: extracted
   type(astsa_stochastic_regression_t) :: stochastic_fit
   type(astsa_sv_filter_t) :: volatility_filter
   type(astsa_sv_fit_t) :: volatility_fit
   type(astsa_sv_particle_filter_t) :: particle_filter
   type(astsa_sv_mcmc_t) :: volatility_mcmc
   type(astsa_ssm_result_t) :: scalar_ssm
   type(astsa_ar_bootstrap_t) :: ar_bootstrap
   type(astsa_ar_mcmc_t) :: ar_posterior
   type(astsa_spectrum_ic_t) :: selected_ar_spectrum
   type(astsa_linearity_test_t) :: linearity
   real(dp) :: spectral_data(8, 2)
   real(dp) :: fit_series(20), fit_innovations(20), expected_ar
   real(dp) :: fit_regressors(20, 1), regression_series(20), expected_beta
   real(dp) :: lag_input(8), lag_output(8)
   real(dp) :: impulse_series(101)
   real(dp) :: stochastic_data(8, 3)
   real(dp) :: volatility_returns(12)
   real(dp) :: particle_normals(4, 11), particle_uniforms(4, 11), ancestor_uniforms(11)
   real(dp) :: ssm_data(20)
   integer :: bootstrap_indices(19, 3)
   real(dp) :: ar_normal_draws(2, 5), ar_gamma_draws(5)
   real(dp) :: linear_series(64)
   integer :: i

   y(:, 1) = [1.0_dp, 0.5_dp, 1.2_dp, 0.7_dp, 1.5_dp, 1.1_dp]
   y(:, 2) = [-0.2_dp, 0.1_dp, 0.4_dp, 0.3_dp, 0.6_dp, 0.2_dp]
   observation(:, :, 1) = reshape([1.0_dp, 0.2_dp, 0.3_dp, 1.0_dp], [2, 2])
   mu0 = 0.0_dp
   sigma0 = 0.0_dp
   sigma0(1, 1) = 1.0_dp
   sigma0(2, 2) = 0.5_dp
   transition = reshape([0.8_dp, 0.0_dp, 0.1_dp, 0.6_dp], [2, 2])
   sq = 0.0_dp
   sq(1, 1) = 0.3_dp
   sq(2, 2) = 0.2_dp
   sr = 0.0_dp
   sr(1, 1) = 0.5_dp
   sr(2, 2) = 0.4_dp

   filtered = astsa_kfilter(y, observation, mu0, sigma0, transition, sq, sr)
   call check(filtered%info == 0, 'Kfilter status')
   call check(maxval(abs(filtered%pp(:, :, 1) - reshape([0.735_dp, 0.03_dp, &
      0.03_dp, 0.22_dp], [2, 2]))) < 1.0e-12_dp, 'Kfilter initial covariance')
   call check(maxval(abs(filtered%xf(:, 1) - [0.7287120802343682_dp, &
      -0.1519905849096833_dp])) < 1.0e-12_dp, 'Kfilter state reference')
   call check(abs(filtered%likelihood + 3.157319364822349_dp) < 1.0e-11_dp, &
      'Kfilter likelihood reference')

   smoothed = astsa_ksmooth(y, observation, mu0, sigma0, transition, sq, sr)
   call check(smoothed%info == 0, 'Ksmooth status')
   call check(maxval(abs(smoothed%xs(:, 1) - [0.8634181644817324_dp, &
      -0.1085734366450517_dp])) < 1.0e-11_dp, 'Ksmooth state reference')
   call check(maxval(abs(smoothed%x0n - [0.961239797439466_dp, &
      -0.1371315068178038_dp])) < 1.0e-11_dp, 'Ksmooth initial state reference')

   correlation = reshape([0.15_dp, 0.05_dp, 0.0_dp, -0.1_dp], [2, 2])
   filtered = astsa_kfilter_correlated(y, observation, mu0, sigma0, transition, sq, sr, correlation)
   call check(filtered%info == 0, 'correlated Kfilter status')
   call check(maxval(abs(filtered%xp(:, 2) - [0.5962902642679945_dp, &
      -0.07516905748417976_dp])) < 1.0e-11_dp, 'correlated prediction reference')
   call check(maxval(abs(filtered%pp(:, :, 2) - reshape([0.182121749441793_dp, &
      -0.01515792497851144_dp, -0.01515792497851144_dp, &
      0.08041184606439436_dp], [2, 2]))) < 1.0e-11_dp, &
      'correlated covariance reference')
   call check(abs(filtered%likelihood + 3.222296566579689_dp) < 1.0e-11_dp, &
      'correlated likelihood reference')
   smoothed = astsa_ksmooth_correlated(y, observation, mu0, sigma0, transition, sq, sr, correlation)
   call check(smoothed%info == 0, 'correlated Ksmooth status')
   call check(maxval(abs(smoothed%xs(:, 1) - [0.87432965686772_dp, &
      -0.114180242861945_dp])) < 1.0e-11_dp, 'correlated smoother reference')
   call check(maxval(abs(smoothed%x0n - [0.9740195611778942_dp, &
      -0.1446319270710853_dp])) < 1.0e-11_dp, 'correlated initial state reference')

   coefficients = arma_to_ar([0.5_dp, -0.2_dp], [0.3_dp], 6)
   call check(maxval(abs(coefficients - [-0.8_dp, 0.44_dp, -0.132_dp, &
      0.0396_dp, -0.01188_dp, 0.003564_dp])) < 1.0e-12_dp, &
      'ARMAtoAR reference')
   product = poly_mul([1.0_dp, 2.0_dp, 3.0_dp], [2.0_dp, -1.0_dp])
   call check(maxval(abs(product - [2.0_dp, 3.0_dp, 4.0_dp, -3.0_dp])) &
      < 1.0e-12_dp, 'polyMul reference')
   matrix = 0.0_dp
   matrix(1, 1) = 4.0_dp
   matrix(2, 2) = 9.0_dp
   power_result = symmetric_matrix_power(matrix, 0.5_dp)
   powered = power_result%value
   call check(power_result%info == 0 .and. abs(powered(1, 1) - 2.0_dp) < 1.0e-12_dp &
      .and. abs(powered(2, 2) - 3.0_dp) < 1.0e-12_dp, 'matrix power')
   call check(fdr_cutoff([0.001_dp, 0.04_dp, 0.03_dp, 0.2_dp], &
      0.05_dp) == 1, 'FDR reference')

   em = astsa_em(y, observation, mu0, sigma0, transition, matmul(sq, transpose(sq)), &
      matmul(sr, transpose(sr)), max_iterations=5, tolerance=1.0e-15_dp)
   call check(em%info == 0, 'EM status')
   call check(em%iterations == 5, 'EM iterations')
   call check(maxval(abs(em%transition - reshape([0.9947408524535486_dp, &
      0.06887480658896009_dp, 0.2615727883676138_dp, &
      0.4608235580365132_dp], [2, 2]))) < 1.0e-9_dp, 'EM transition reference')
   call check(maxval(abs(em%q - reshape([0.04149096278771728_dp, &
      0.006708583378070144_dp, 0.006708583378070144_dp, &
      0.01649359688562731_dp], [2, 2]))) < 1.0e-9_dp, 'EM Q reference')
   call check(maxval(abs(em%r - reshape([0.08810427309266976_dp, 0.0_dp, &
      0.0_dp, 0.01757599413820875_dp], [2, 2]))) < 1.0e-9_dp, 'EM R reference')
   call check(maxval(abs(em%mu0 - [1.179443415144531_dp, -0.7468490650937289_dp])) &
      < 1.0e-9_dp, 'EM initial mean reference')
   call check(maxval(abs(em%likelihood - [-3.157319364822349_dp, -6.721038472290012_dp, &
      -8.394807971429008_dp, -9.424340191532089_dp, &
      -10.10828880793605_dp])) < 1.0e-9_dp, 'EM likelihood reference')

   ffbs = astsa_ffbs_draws(y, observation, mu0, sigma0, transition, sq, sr, &
      reshape([(0.0_dp, i=1, 14)], [2, 7]))
   call check(ffbs%info == 0, 'FFBS deterministic status')
   call check(maxval(abs(ffbs%states(:, 1) - [0.8634181644817324_dp, &
      -0.1085734366450517_dp])) < 1.0e-11_dp, 'FFBS first state reference')
   call check(maxval(abs(ffbs%states(:, 6) - [0.8532899049358921_dp, &
      0.08840485491407681_dp])) < 1.0e-11_dp, 'FFBS final state reference')
   call check(maxval(abs(ffbs%initial_state - [0.961239797439466_dp, &
      -0.1371315068178038_dp])) < 1.0e-11_dp, 'FFBS initial state reference')

   spectrum = arma_spectrum([0.7_dp, -0.2_dp], [0.3_dp], &
      noise_variance=1.4_dp, n_frequency=6, sampling_frequency=4.0_dp)
   call check(spectrum%info == 0, 'ARMA spectrum status')
   call check(maxval(abs(spectrum%spectrum(:, 1) - [9.464_dp, 7.490278490234968_dp, &
      2.598163944757945_dp, 0.7339293428687489_dp, 0.2809473810087784_dp, &
      0.1900277008310249_dp])) < 1.0e-11_dp, 'ARMA spectrum reference')
   spectral_data(:, 1) = [1.0_dp, 0.5_dp, 1.2_dp, 0.7_dp, &
      1.5_dp, 1.1_dp, 1.4_dp, 0.9_dp]
   spectral_data(:, 2) = [-0.2_dp, 0.1_dp, 0.4_dp, 0.3_dp, &
      0.6_dp, 0.2_dp, 0.5_dp, 0.1_dp]
   spectrum = mv_periodogram(spectral_data, demean=.true., detrend=.false.)
   call check(spectrum%info == 0, 'periodogram status')
   call check(maxval(abs(spectrum%spectrum(:, 1) - [0.1498896103067893_dp, &
      0.001249999999999997_dp, 0.02261038969321076_dp, 0.45125_dp])) &
      < 1.0e-11_dp, 'periodogram marginal reference')
   call check(maxval(abs(spectrum%phase(1, 2, :) - [-0.7453887094443166_dp, &
      0.1973955598498818_dp, -1.317537928094345_dp, 0.0_dp])) &
      < 1.0e-11_dp, 'periodogram phase reference')
   spectrum = mv_periodogram(spectral_data, demean=.true., detrend=.false., span=3)
   call check(maxval(abs(spectrum%spectrum(:, 1) - [0.1127297077300919_dp, &
      0.04375_dp, 0.1244301948466054_dp, 0.2369301948466054_dp])) &
      < 1.0e-11_dp, 'smoothed spectrum reference')
   call check(abs(spectrum%degrees_freedom - 5.333333333333333_dp) < 1.0e-10_dp, &
      'smoothed spectrum degrees of freedom')
   simulation = sarima_sim_from_innovations(ar=[0.5_dp], n=4, burnin=0, &
      innovations=[1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
   call check(simulation%info == 0 .and. maxval(abs(simulation%series - &
      [1.0_dp, 0.5_dp, 0.25_dp, 0.125_dp])) < 1.0e-12_dp, 'SARIMA AR simulation')
   simulation = sarima_sim_from_innovations(d=1, n=3, burnin=0, &
      innovations=[1.0_dp, 1.0_dp, 1.0_dp])
   call check(simulation%info == 0 .and. maxval(abs(simulation%series - &
      [1.0_dp, 2.0_dp, 3.0_dp])) < 1.0e-12_dp, 'SARIMA integrated simulation')
   simulation = sarima_sim_from_innovations(sar=[0.5_dp], season=2, n=5, burnin=0, &
      innovations=[1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
   call check(simulation%info == 0 .and. simulation%frequency == 2 .and. &
      maxval(abs(simulation%series - [1.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.25_dp])) &
      < 1.0e-12_dp, 'SARIMA seasonal simulation')
   simulation = sarima_sim_from_innovations(ar=[1.0_dp], n=2, burnin=0, &
      innovations=[0.0_dp, 0.0_dp])
   call check(simulation%info == 6, 'SARIMA causality check')
   simulation = sarima_sim_from_innovations(ma=[1.0_dp], n=2, burnin=0, &
      innovations=[0.0_dp, 0.0_dp])
   call check(simulation%info == 8, 'SARIMA invertibility check')
   diagnostic = arma_check(ar=[0.5_dp], ma=[0.4_dp])
   call check(diagnostic%info == 0 .and. diagnostic%causal .and. diagnostic%invertible, &
      'ARMA stability diagnostics')
   call check(abs(diagnostic%ar_roots(1) - cmplx(2.0_dp, 0.0_dp, dp)) < 1.0e-12_dp, &
      'ARMA root value')
   diagnostic = arma_check(ar=[0.5_dp], ma=[-0.5_dp])
   call check(diagnostic%redundancy_checked .and. diagnostic%redundant, &
      'ARMA common factor')
   diagnostic = arma_check(sar=[0.5_dp], sma=[-0.5_dp], season=4)
   call check(diagnostic%seasonal_redundant .and. diagnostic%seasonal_period == 4, &
      'seasonal ARMA common factor')
   diagnostic = arma_check(ar=[1.0_dp])
   call check(.not. diagnostic%causal .and. .not. diagnostic%redundancy_checked, &
      'ARMA noncausal diagnostic')
   prewhitened = pre_white([1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp, 11.0_dp], &
      [2.0_dp, 3.0_dp, 5.0_dp, 8.0_dp, 12.0_dp], differences=1, max_lag=1, order_max=0)
   call check(prewhitened%info == 0 .and. prewhitened%order == 0 .and. &
      prewhitened%differences == 1, 'prewhitening model')
   call check(maxval(abs(prewhitened%first - [-1.5_dp, -0.5_dp, 0.5_dp, 1.5_dp])) &
      < 1.0e-12_dp, 'prewhitening first series')
   call check(maxval(abs(prewhitened%second - [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])) &
      < 1.0e-12_dp, 'prewhitening second series')
   call check(all(prewhitened%lags == [-1, 0, 1]) .and. &
      abs(prewhitened%correlation(2) - 1.0_dp) < 1.0e-12_dp, 'prewhitening CCF')
   prewhitened = pre_white(spectral_data(:, 1), spectral_data(:, 2), order_max=2)
   call check(prewhitened%info == 0 .and. prewhitened%order >= 0 .and. &
      prewhitened%order <= 2 .and. size(prewhitened%first) == 8 - prewhitened%order, &
      'automatic prewhitening order')
   likelihood = sarima_likelihood([1.0_dp, 2.5_dp, 0.25_dp, 0.625_dp], ar=[0.5_dp])
   call check(likelihood%info == 0 .and. likelihood%observations == 3 .and. &
      abs(likelihood%sigma2 - 1.75_dp) < 1.0e-12_dp, 'SARIMA AR likelihood')
   call check(maxval(abs(likelihood%residuals - [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp])) &
      < 1.0e-12_dp, 'SARIMA AR innovations')
   call check(abs(likelihood%log_likelihood + &
      1.5_dp*(log(2.0_dp*acos(-1.0_dp)*1.75_dp) + 1.0_dp)) < 1.0e-12_dp, &
      'SARIMA Gaussian likelihood')
   likelihood = sarima_likelihood([1.0_dp, 2.5_dp, 0.0_dp, 0.0_dp], ma=[0.5_dp])
   call check(likelihood%info == 0 .and. maxval(abs(likelihood%residuals - &
      [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp])) < 1.0e-12_dp, 'SARIMA MA likelihood')
   likelihood = sarima_exact_likelihood([1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp])
   call check(likelihood%info == 0 .and. likelihood%observations == 4 .and. &
      abs(likelihood%sigma2 - 1.5625_dp) < 1.0e-12_dp, 'exact white-noise likelihood')
   fit_innovations = [0.3_dp, -0.7_dp, 0.2_dp, 1.1_dp, -0.4_dp, 0.6_dp, -0.2_dp, &
      0.8_dp, -1.0_dp, 0.5_dp, 0.1_dp, -0.3_dp, 0.9_dp, -0.6_dp, 0.4_dp, &
      -0.1_dp, 0.7_dp, -0.5_dp, 0.2_dp, -0.8_dp]
   fit_series(1) = fit_innovations(1)
   do i = 2, size(fit_series)
      fit_series(i) = 0.6_dp*fit_series(i - 1) + fit_innovations(i)
   end do
   expected_ar = dot_product(fit_series(2:), fit_series(:size(fit_series) - 1))/ &
      sum(fit_series(:size(fit_series) - 1)**2)
   sarima_estimate = sarima_fit(fit_series, 1, 0, 0, 0, 0, 0, 1, &
      max_iterations=100, tolerance=1.0e-7_dp)
   call check(sarima_estimate%converged .and. sarima_estimate%info == 0, &
      'SARIMA BFGS convergence')
   call check(abs(sarima_estimate%coefficients(1) - expected_ar) < 1.0e-5_dp, &
      'SARIMA BFGS AR estimate')
   diagnostic = arma_check(ar=sarima_estimate%coefficients(1:1))
   call check(diagnostic%causal, 'SARIMA transformed AR estimate')
   call check(sarima_estimate%standard_errors(1) > 0.0_dp .and. &
      abs(sarima_estimate%covariance(1, 1) - sarima_estimate%standard_errors(1)**2) &
      < 1.0e-10_dp, 'SARIMA Hessian standard error')
   sarima_estimate = sarima_fit(fit_series, 1, 0, 0, 0, 0, 0, 1, &
      initial=[0.4_dp], estimated=[.false.])
   call check(sarima_estimate%converged .and. sarima_estimate%iterations == 0 .and. &
      abs(sarima_estimate%coefficients(1) - 0.4_dp) < 1.0e-12_dp, 'SARIMA fixed parameter')
   fixed_estimate = sarima_fit(fit_series, 1, 0, 0, 0, 0, 0, 1, &
      initial=[0.6_dp], estimated=[.false.])
   sarima_prediction = sarima_forecast(fixed_estimate, fit_series, 3)
   call check(sarima_prediction%info == 0 .and. maxval(abs(sarima_prediction%mean - &
      [0.6_dp*fit_series(20), 0.6_dp**2*fit_series(20), 0.6_dp**3*fit_series(20)])) &
      < 1.0e-10_dp, 'SARIMA AR forecast means')
   call check(sarima_prediction%standard_error(2) > sarima_prediction%standard_error(1) .and. &
      all(sarima_prediction%lower < sarima_prediction%mean) .and. &
      all(sarima_prediction%upper > sarima_prediction%mean), 'SARIMA forecast intervals')
   fixed_estimate = sarima_fit(fit_series, 1, 0, 0, 0, 0, 0, 1, &
      initial=[0.6_dp], estimated=[.false.], exact_likelihood=.true.)
   call check(fixed_estimate%info == 0 .and. fixed_estimate%uses_exact_likelihood .and. &
      fixed_estimate%likelihood%observations == 20, 'exact SARIMA fit likelihood')
   sarima_diagnostic = sarima_diagnostics(fixed_estimate, max_lag=5)
   call check(sarima_diagnostic%info == 0 .and. &
      size(sarima_diagnostic%standardized_residuals) == 20 .and. &
      size(sarima_diagnostic%lags) == 4, 'SARIMA residual diagnostics')
   call check(abs(sarima_diagnostic%residual_acf(1) - 1.0_dp) < 1.0e-12_dp .and. &
      all(sarima_diagnostic%p_values >= 0.0_dp) .and. &
      all(sarima_diagnostic%p_values <= 1.0_dp), 'SARIMA Ljung-Box diagnostics')
   call check(all(sarima_diagnostic%qq_sample(2:) >= sarima_diagnostic%qq_sample(:19)) .and. &
      abs(sarima_diagnostic%qq_theoretical(1) + sarima_diagnostic%qq_theoretical(20)) &
      < 1.0e-12_dp, 'SARIMA Q-Q diagnostics')
   lag_input = [1.0_dp, -1.0_dp, 2.0_dp, 0.0_dp, -2.0_dp, 1.0_dp, 3.0_dp, -1.0_dp]
   lag_output = 1.0_dp + 2.0_dp*lag_input
   lag_fit = lag_reg(lag_input, lag_output, m=8, span=3, threshold=1.5_dp)
   call check(lag_fit%info == 0 .and. any(lag_fit%selected_lags == 0), &
      'LagReg transfer selection')
   call check(lag_fit%mse < 1.0e-20_dp, 'LagReg fitted output')
   impulse_series = 0.0_dp
   impulse_series(51) = 1.0_dp
   extracted = signal_extract(impulse_series, m=10, max_frequency=0.3_dp, min_frequency=0.05_dp)
   call check(extracted%info == 0 .and. extracted%valid_start == 5 .and. &
      extracted%valid_end == 97, 'SigExtract valid range')
   call check(maxval(abs(extracted%coefficients - extracted%coefficients(9:1:-1))) &
      < 1.0e-12_dp, 'SigExtract symmetric filter')
   call check(abs(extracted%filtered(51) - extracted%coefficients(5)) < 1.0e-12_dp .and. &
      all(abs(extracted%attained_response) < huge(1.0_dp)), 'SigExtract impulse response')
   stochastic_data(:, :2) = spectral_data
   stochastic_data(:, 3) = 0.7_dp*spectral_data(:, 1) - 0.3_dp*spectral_data(:, 2) + &
      [0.1_dp, -0.05_dp, 0.08_dp, -0.02_dp, 0.04_dp, -0.06_dp, 0.03_dp, -0.01_dp]
   stochastic_fit = stochastic_regression(stochastic_data, [1, 2], [1], 3, 3, 8, 0.05_dp)
   call check(stochastic_fit%info == 0 .and. all(stochastic_fit%full_power <= &
      stochastic_fit%reduced_power + 1.0e-10_dp), 'stochastic regression powers')
   call check(all(stochastic_fit%coherence >= 0.0_dp) .and. &
      all(stochastic_fit%coherence <= 1.0_dp) .and. stochastic_fit%critical_f > 0.0_dp, &
      'stochastic regression tests')
   call check(all(shape(stochastic_fit%coefficients) == [7, 2]), &
      'stochastic regression impulse coefficients')
   volatility_returns = [0.02_dp, -0.01_dp, 0.03_dp, -0.025_dp, 0.015_dp, -0.02_dp, &
      0.04_dp, -0.03_dp, 0.01_dp, -0.015_dp, 0.025_dp, -0.02_dp]
   volatility_filter = sv_filter(volatility_returns, 0.0_dp, 0.95_dp, 0.1_dp, &
      -7.0_dp, 1.0_dp, -3.0_dp, 2.0_dp, 0.0_dp)
   call check(volatility_filter%info == 0 .and. &
      all(volatility_filter%prediction_variance > 0.0_dp), 'SV mixture filter')
   call check(maxval(abs(volatility_filter%component_zero_probability + &
      volatility_filter%component_one_probability - 1.0_dp)) < 1.0e-12_dp, &
      'SV mixture probabilities')
   volatility_fit = sv_mle(volatility_returns, max_iterations=40, tolerance=1.0e-5_dp)
   call check(volatility_fit%filter%info == 0 .and. abs(volatility_fit%coefficients(1)) < 1.0_dp .and. &
      volatility_fit%coefficients(2) > 0.0_dp, 'SV maximum likelihood fit')
   particle_normals = 0.0_dp
   particle_uniforms = 0.5_dp
   ancestor_uniforms = 0.5_dp
   particle_filter = sv_particle_filter_draws(volatility_returns, 0.9_dp, 0.04_dp, 0.1_dp, &
      [(0.0_dp, i=1, 12)], particle_normals, particle_uniforms, ancestor_uniforms, 0.5_dp)
   call check(particle_filter%info == 0 .and. &
      maxval(abs(sum(particle_filter%weights, dim=1) - 1.0_dp)) < 1.0e-12_dp, &
      'SV particle weights')
   call check(maxval(abs(particle_filter%trajectory)) < 1.0e-12_dp .and. &
      particle_filter%log_likelihood > -huge(1.0_dp), 'SV particle trajectory')
   call set_random_seed(24680)
   volatility_mcmc = sv_mcmc(volatility_returns, 4, 2, 4)
   call check(volatility_mcmc%info == 0 .and. size(volatility_mcmc%phi) == 4 .and. &
      all(shape(volatility_mcmc%latent) == [4, 12]), 'SV MCMC retained draws')
   call check(all(abs(volatility_mcmc%phi) < 1.0_dp) .and. &
      all(volatility_mcmc%state_sd > 0.0_dp) .and. all(volatility_mcmc%observation_scale > 0.0_dp), &
      'SV MCMC parameter constraints')
   call check(volatility_mcmc%acceptance_rate >= 0.0_dp .and. &
      volatility_mcmc%acceptance_rate <= 1.0_dp .and. all(volatility_mcmc%effective_size > 0.0_dp), &
      'SV MCMC diagnostics')
   ssm_data = [0.2_dp, 0.35_dp, 0.41_dp, 0.52_dp, 0.48_dp, 0.61_dp, 0.73_dp, &
      0.69_dp, 0.82_dp, 0.91_dp, 0.88_dp, 1.02_dp, 1.08_dp, 1.15_dp, &
      1.12_dp, 1.23_dp, 1.31_dp, 1.28_dp, 1.39_dp, 1.46_dp]
   scalar_ssm = ssm_fit(ssm_data, 1.0_dp, 0.8_dp, 0.05_dp, 0.1_dp, 0.15_dp, &
      fix_phi=.true., max_iterations=60)
   call check(size(scalar_ssm%smoothed_state) == 20 .and. scalar_ssm%fixed_phi .and. &
      scalar_ssm%coefficients(2) > 0.0_dp .and. scalar_ssm%coefficients(3) > 0.0_dp, &
      'scalar SSM estimation')
   call check(all(scalar_ssm%predicted_variance > 0.0_dp) .and. &
      all(scalar_ssm%smoothed_variance >= 0.0_dp), 'scalar SSM state variances')
   do i = 1, 19
      bootstrap_indices(i, 1) = i
      bootstrap_indices(i, 2) = 1 + mod(i, 19)
      bootstrap_indices(i, 3) = 1 + mod(i + 5, 19)
   end do
   ar_bootstrap = ar_boot_draws(fit_series, 1, bootstrap_indices, [0.1_dp, 0.5_dp, 0.9_dp])
   call check(ar_bootstrap%info == 0 .and. all(shape(ar_bootstrap%coefficient_draws) == [3, 1]) .and. &
      all(shape(ar_bootstrap%simulated_series) == [20, 3]), 'AR bootstrap shapes')
   call check(all(abs(ar_bootstrap%simulated_series(1, :) - fit_series(1)) < 1.0e-12_dp) .and. &
      size(ar_bootstrap%quantiles, 1) == 3, 'AR bootstrap simulation and quantiles')
   ar_normal_draws = 0.0_dp
   ar_gamma_draws = 10.0_dp
   ar_posterior = ar_mcmc_draws(fit_series, 1, 4, 2, ar_normal_draws, ar_gamma_draws, &
      probabilities=[0.1_dp, 0.5_dp, 0.9_dp])
   call check(ar_posterior%info == 0 .and. all(shape(ar_posterior%coefficient_draws) == [4, 2]) .and. &
      size(ar_posterior%innovation_sd) == 4, 'AR MCMC retained draws')
   call check(all(ar_posterior%innovation_sd > 0.0_dp) .and. &
      all(shape(ar_posterior%quantiles) == [3, 3]) .and. all(ar_posterior%effective_size > 0.0_dp), &
      'AR MCMC posterior summaries')
   selected_ar_spectrum = spectrum_ic(spectral_data(:, 1), order_max=3, use_bic=.true., &
      detrend=.false., frequency_count=16)
   call check(selected_ar_spectrum%info == 0 .and. selected_ar_spectrum%selected_order >= 0 .and. &
      selected_ar_spectrum%selected_order <= 3, 'spectral IC selected order')
   call check(abs(minval(selected_ar_spectrum%aic)) < 1.0e-12_dp .and. &
      abs(minval(selected_ar_spectrum%bic)) < 1.0e-12_dp .and. &
      size(selected_ar_spectrum%spectrum) == 16, 'spectral IC criteria and spectrum')
   do i = 1, size(linear_series)
      linear_series(i) = sin(2.0_dp*acos(-1.0_dp)*real(i, dp)/8.0_dp) + &
         0.2_dp*cos(2.0_dp*acos(-1.0_dp)*real(i, dp)/5.0_dp)
   end do
   linearity = test_linearity(linear_series)
   call check(linearity%info == 0 .and. linearity%block_length == 7 .and. &
      all(shape(linearity%statistic) == [3, 3]), 'linearity bispectrum shape')
   call check(all(linearity%p_values >= 0.0_dp) .and. all(linearity%p_values <= 1.0_dp) .and. &
      linearity%noncentrality >= 0.0_dp, 'linearity p-values')
   fixed_estimate = sarima_fit([1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp], 0, 1, 0, 0, 0, 0, 1, &
      exact_likelihood=.true.)
   call check(fixed_estimate%info == 0 .and. fixed_estimate%uses_exact_likelihood .and. &
      fixed_estimate%likelihood%observations == 3 .and. &
      abs(fixed_estimate%likelihood%sigma2 - 14.0_dp/3.0_dp) < 1.0e-10_dp, &
      'exact diffuse random-walk likelihood')
   fixed_estimate = sarima_fit([1.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, 9.0_dp], &
      0, 0, 0, 0, 1, 0, 2, exact_likelihood=.true.)
   call check(fixed_estimate%info == 0 .and. fixed_estimate%likelihood%observations == 3 .and. &
      abs(fixed_estimate%likelihood%sigma2 - 50.0_dp/3.0_dp) < 1.0e-10_dp, &
      'exact seasonal diffuse likelihood')
   fixed_estimate = sarima_fit([1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp], 0, 1, 0, 0, 0, 0, 1)
   sarima_prediction = sarima_forecast(fixed_estimate, [1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp], 3)
   call check(sarima_prediction%info == 0 .and. &
      maxval(abs(sarima_prediction%mean - 7.0_dp)) < 1.0e-12_dp, &
      'SARIMA integrated forecast means')
   call check(abs(sarima_prediction%standard_error(3)**2 - &
      3.0_dp*fixed_estimate%likelihood%sigma2) < 1.0e-10_dp, &
      'SARIMA integrated forecast variance')
   fit_regressors(:, 1) = [(real(i, dp), i=1, 20)]
   regression_series = 2.0_dp*fit_regressors(:, 1) + fit_innovations
   expected_beta = dot_product(fit_regressors(:, 1), regression_series)/ &
      sum(fit_regressors(:, 1)**2)
   sarima_estimate = sarima_fit(regression_series, 0, 0, 0, 0, 0, 0, 1, &
      regressors=fit_regressors, max_iterations=100)
   call check(sarima_estimate%converged .and. &
      abs(sarima_estimate%coefficients(1) - expected_beta) < 1.0e-5_dp, &
      'SARIMA regression estimate')
   sarima_prediction = sarima_forecast(sarima_estimate, regression_series, 2, &
      regressors=fit_regressors, new_regressors=reshape([21.0_dp, 22.0_dp], [2, 1]))
   call check(sarima_prediction%info == 0 .and. maxval(abs(sarima_prediction%mean - &
      sarima_estimate%coefficients(1)*[21.0_dp, 22.0_dp])) < 1.0e-8_dp, &
      'SARIMA regression forecast')
   print '(a)', 'All astsa_mod tests passed.'

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
