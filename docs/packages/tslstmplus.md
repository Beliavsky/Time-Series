# TSLSTMplus

[Back to the implemented-package index](../../README.md#implemented-packages).

`tslstmplus.f90` translates the
[TSLSTMplus](https://cran.r-project.org/web/packages/TSLSTMplus/index.html) R
package without requiring Keras or TensorFlow. The package adapter prepares
response and exogenous lags as recurrent timesteps or flattened features,
supports independent exogenous lag orders, and retains standard or min-max
input and output scaling. One or more LSTM layers use configurable candidate,
cell, and recurrent-gate activations, optional input dropout, and stateful or
independent sequence execution. Optional hidden dense layers precede the final
linear forecast output. Analytic truncated backpropagation through time trains
the complete network with mini-batch SGD, Adam, or RMSprop and MSE or MAE loss;
terminal validation samples, minimum improvement, and patience provide early
stopping. Rolling fitted values and residuals retain their original alignment,
while recursive forecasts accept known future exogenous observations. The
reusable stand-alone recurrent engine is implemented in
`recurrent_network.f90`. Principal interfaces are available through
`forecasting_mod` and `regression_time_series_mod`. The translation is licensed
under GPL-3.0-only; see `LICENSE-TSLSTMPLUS`.
