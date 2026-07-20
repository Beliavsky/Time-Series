---
project: Fortran Time-Series Library
summary: Time-series analysis algorithms translated from R packages.
author: Beliavsky
src_dir: .
output_dir: build/ford
exclude: test_*.f90
    **/test_*.f90
    fracdist_tables.f90
    tseriestarma_tables.f90
exclude_dir: build*
    .git
    docs
    examples
display: public
sort: alpha
source: false
graph: false
warn: false
---
