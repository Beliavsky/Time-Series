# SPDX-License-Identifier: GPL-3.0-or-later
# Generate a standalone Fortran form of the tseriesTARMA critical-value data.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: generate_tseriestarma_critical_values.R unit_root.rda garch.rda output.f90")
}

load(args[[1L]])
if (!exists("supLMQur") || !identical(dim(supLMQur), c(4L, 9L, 7L, 11L))) {
  stop("supLMQur has an unexpected shape")
}

load(args[[2L]])
if (!exists("ACValues") || !identical(dim(ACValues), c(13L, 62L))) {
  stop("ACValues has an unexpected shape")
}

format_data <- function(values, per_line = 4L) {
  numbers <- sprintf("%.16e_dp", values)
  groups <- split(numbers, ceiling(seq_along(numbers) / per_line))
  vapply(seq_along(groups), function(i) {
    suffix <- if (i < length(groups)) ", &" else "], &"
    paste0("      ", paste(groups[[i]], collapse = ", "), suffix)
  }, character(1L))
}

unit_data_lines <- format_data(as.vector(supLMQur))
garch_array <- array(0, dim = c(3L, 20L, 13L))
for (row in seq_len(13L)) {
  for (parameters in seq_len(20L)) {
    columns <- paste(parameters, c(90, 95, 99), sep = "-")
    garch_array[, parameters, row] <- as.numeric(ACValues[row, columns])
  }
}
garch_data_lines <- format_data(as.vector(garch_array))
garch_probability_lines <- format_data(as.numeric(ACValues[, "pi"]))
garch_probability_lines[[length(garch_probability_lines)]] <- sub(
  "\\], &$", "]", garch_probability_lines[[length(garch_probability_lines)]])

header <- c(
  "! SPDX-License-Identifier: GPL-3.0-or-later",
  "! SPDX-FileComment: Critical-value data translated from the R tseriesTARMA package.",
  "module tseriestarma_tables_mod",
  "   !! Tabulated critical values used by tseriesTARMA hypothesis tests.",
  "   use kind_mod, only: dp",
  "   implicit none",
  "   private",
  "",
  "   real(dp), parameter :: threshold_probabilities(9) = [0.01_dp, 0.05_dp, &",
  "      0.10_dp, 0.15_dp, 0.20_dp, 0.25_dp, 0.30_dp, 0.35_dp, 0.40_dp]",
  "   real(dp), parameter :: ma_coefficients(7) = [-0.9_dp, -0.6_dp, &",
  "      -0.3_dp, 0.0_dp, 0.3_dp, 0.6_dp, 0.9_dp]",
  "   integer, parameter :: sample_sizes(11) = [100, 200, 300, 400, 500, &",
  "      600, 700, 800, 900, 1000, 5000]",
  "   real(dp), parameter :: critical_probabilities(4) = [0.90_dp, 0.95_dp, &",
  "      0.99_dp, 0.999_dp]",
  "   real(dp), parameter :: critical_table(4, 9, 7, 11) = reshape([ &"
)

middle <- c(
  "      [4, 9, 7, 11])",
  "   real(dp), parameter :: garch_threshold_probabilities(13) = [ &",
  garch_probability_lines,
  "   real(dp), parameter :: garch_critical_probabilities(3) = [0.90_dp, &",
  "      0.95_dp, 0.99_dp]",
  "   real(dp), parameter :: garch_critical_table(3, 20, 13) = reshape([ &"
)

footer <- c(
  "      [3, 20, 13])",
  "",
  "   public :: tseriestarma_unit_root_critical_values",
  "   public :: tseriestarma_unit_root_critical_probabilities",
  "   public :: tseriestarma_garch_critical_values",
  "   public :: tseriestarma_garch_critical_probabilities",
  "",
  "contains",
  "",
  "   pure function tseriestarma_unit_root_critical_values(observations, &",
  "      lower_probability, ma_coefficient) result(values)",
  "      !! Return the nearest tabulated unit-root supLM critical values.",
  "      integer, intent(in) :: observations !! Number of time-series observations.",
  "      real(dp), intent(in) :: lower_probability !! Lower threshold-search quantile.",
  "      real(dp), intent(in) :: ma_coefficient !! Fitted moving-average coefficient.",
  "      real(dp) :: values(4)",
  "      integer :: observation_index, probability_index, ma_index",
  "",
  "      observation_index = minloc(abs(real(sample_sizes, dp) - &",
  "         real(observations, dp)), dim=1)",
  "      probability_index = minloc(abs(threshold_probabilities - &",
  "         lower_probability), dim=1)",
  "      ma_index = minloc(abs(ma_coefficients - ma_coefficient), dim=1)",
  "      values = critical_table(:, probability_index, ma_index, observation_index)",
  "   end function tseriestarma_unit_root_critical_values",
  "",
  "   pure function tseriestarma_unit_root_critical_probabilities() result(values)",
  "      !! Return probabilities associated with the unit-root critical values.",
  "      real(dp) :: values(4)",
  "",
  "      values = critical_probabilities",
  "   end function tseriestarma_unit_root_critical_probabilities",
  "",
  "   pure function tseriestarma_garch_critical_values(parameter_count, &",
  "      lower_probability) result(values)",
  "      !! Return the nearest tabulated ARMA-GARCH supLM critical values.",
  "      integer, intent(in) :: parameter_count !! Number of tested threshold parameters.",
  "      real(dp), intent(in) :: lower_probability !! Lower threshold-search quantile.",
  "      real(dp) :: values(3)",
  "      integer :: probability_index",
  "",
  "      values = 0.0_dp",
  "      if (parameter_count < 1 .or. parameter_count > 20) return",
  "      probability_index = minloc(abs(garch_threshold_probabilities - &",
  "         lower_probability), dim=1)",
  "      values = garch_critical_table(:, parameter_count, probability_index)",
  "   end function tseriestarma_garch_critical_values",
  "",
  "   pure function tseriestarma_garch_critical_probabilities() result(values)",
  "      !! Return probabilities associated with the ARMA-GARCH critical values.",
  "      real(dp) :: values(3)",
  "",
  "      values = garch_critical_probabilities",
  "   end function tseriestarma_garch_critical_probabilities",
  "",
  "end module tseriestarma_tables_mod"
)

writeLines(c(header, unit_data_lines, middle, garch_data_lines, footer), args[[3L]], useBytes = TRUE)
