args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: extract_fracdist_tables.R sysdata.rda output.f90")
}

load(args[[1]])
output <- file(args[[2]], open = "wt")
on.exit(close(output))

writeLines(c(
  "! Generated from the GPL-3 fracdist package sysdata.rda.",
  "module fracdist_tables_mod",
  "   use kind_mod, only: dp",
  "   implicit none",
  "   private",
  "   public :: fracdist_table_values",
  "contains",
  "",
  "   pure subroutine fracdist_table_values(iq, iscon, b_values, probabilities, &",
  "      chi_square_quantiles, quantiles, info)",
  "      ! Return one embedded fracdist response-surface table.",
  "      integer, intent(in) :: iq, iscon",
  "      real(dp), intent(out) :: b_values(31), probabilities(221)",
  "      real(dp), intent(out) :: chi_square_quantiles(221), quantiles(221, 31)",
  "      integer, intent(out) :: info",
  "      integer :: i",
  "",
  "      info = 0",
  "      b_values = [ (0.50_dp + 0.05_dp*real(i - 1, dp), i = 1, 31) ]"
), output)

probabilities <- ginv_tab[, "probs"]
write_vector <- function(name, values, indent = "      ", chunk = 100, column = NULL) {
  starts <- seq(1, length(values), by = chunk)
  for (start in starts) {
    finish <- min(length(values), start + chunk - 1)
    lines <- sprintf("%.17g_dp", values[start:finish])
    target <- if (is.null(column)) {
      sprintf("%s(%d:%d)", name, start, finish)
    } else {
      sprintf("%s(%d:%d, %d)", name, start, finish, column)
    }
    prefix <- sprintf("%s%s = [ &", indent, target)
    writeLines(prefix, output)
    groups <- split(lines, ceiling(seq_along(lines) / 4))
    for (group_index in seq_along(groups)) {
      suffix <- if (group_index == length(groups)) " ]" else ", &"
      writeLines(sprintf("%s   %s%s", indent, paste(groups[[group_index]], collapse = ", "), suffix),
                 output)
    }
  }
}

write_vector("probabilities", probabilities)
writeLines(c(
  "      select case (iq)",
  "      case (1:12)",
  "         continue",
  "      case default",
  "         info = 1",
  "         quantiles = 0.0_dp",
  "         chi_square_quantiles = 0.0_dp",
  "         return",
  "      end select"
), output)

for (iq in 1:12) {
  writeLines(sprintf("      if (iq == %d) then", iq), output)
  write_vector("chi_square_quantiles", ginv_tab[, sprintf("iq_%d", iq)], "         ")
  writeLines("      end if", output)
}

writeLines(c(
  "      if (iscon /= 0 .and. iscon /= 1) then",
  "         info = 2",
  "         quantiles = 0.0_dp",
  "         return",
  "      end if",
  "      select case (iscon*12 + iq)"
), output)

for (iscon in 0:1) {
  for (iq in 1:12) {
    object_name <- sprintf("%s%02d", if (iscon == 0) "frmapp" else "frcapp", iq)
    table <- get(object_name)
    matrix_values <- matrix(table$xndf, nrow = 221, ncol = 31)
    writeLines(sprintf("      case (%d)", iscon*12 + iq), output)
    for (column in 1:31) {
      write_vector("quantiles", matrix_values[, column], "         ", column = column)
    }
  }
}

writeLines(c(
  "      end select",
  "   end subroutine fracdist_table_values",
  "",
  "end module fracdist_tables_mod"
), output)
