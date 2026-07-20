# esemifar

[Back to the implemented-package index](../../README.md#implemented-packages).

`esemifar.f90` translates the
[esemifar](https://cran.r-project.org/web/packages/esemifar/index.html) R
package. It implements boundary-aware local-polynomial trend and derivative
smoothing, ESEMIFAR iterative plug-in bandwidth selection, FARIMA order grids,
finite and infinite filter conversions, and analytic, residual-bootstrap, and
FARIMA-refitted predictive-root forecasts. Parametric fitting and simulation
reuse `fracdiff_mod`; random sampling and numerical utilities use the shared
library modules. The translation is licensed under GPL-3.0-only.
