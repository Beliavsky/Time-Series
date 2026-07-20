# KFAS

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`kfas.f90` translates the
[KFAS](https://cran.r-project.org/web/packages/KFAS/index.html) R package.

## Algorithms and Procedures

It contains typed Gaussian state-space models, Kalman filtering, smoothing,
prediction, exact diffuse filtering for diagonal observation covariance,
standardized innovations, and disturbance means.

Covariance smoothing returns lag-one and forward conditional state covariances
using singular-safe pseudoinverses, while a fast smoother can omit covariance
output.

## Licensing

The translation is licensed under GPL-2.0-or-later.
