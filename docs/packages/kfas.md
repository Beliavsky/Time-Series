# KFAS

[Back to the implemented-package index](../../README.md#implemented-packages).

`kfas.f90` translates the
[KFAS](https://cran.r-project.org/web/packages/KFAS/index.html) R package. It
contains typed Gaussian state-space models, Kalman filtering,
smoothing, prediction, exact diffuse filtering for diagonal observation
covariance, standardized innovations, and disturbance means. Covariance
smoothing returns lag-one and forward conditional state covariances using
singular-safe pseudoinverses, while a fast smoother can omit covariance output.
The translation is licensed under GPL-2.0-or-later.
