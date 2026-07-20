"""Tests for the standalone Fortran style auditor."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("fortran_style_audit.py")
SPEC = importlib.util.spec_from_file_location("fortran_style_audit", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDIT
SPEC.loader.exec_module(AUDIT)


class AuditTests(unittest.TestCase):
    """Exercise quote handling and procedure-level rules."""

    def audit(self, source: str):
        """Audit one temporary free-form source file."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.f90"
            path.write_text(source, encoding="utf-8")
            return AUDIT.audit_file(path, "sample.f90", 132)

    def test_clean_function(self):
        """Accept the preferred function form."""
        findings, _ = self.audit(
            "pure elemental function square(x) result(value)\n"
            "   ! Return the square of a scalar.\n"
            "   real(dp), intent(in) :: x\n"
            "   real(dp) :: value\n"
            "   value = x*x\n"
            "end function square\n"
        )
        self.assertEqual([], findings)

    def test_semicolon_ignores_strings_and_comments(self):
        """Only diagnose statement separators in code."""
        findings, _ = self.audit(
            "subroutine show()\n"
            "   ! Display punctuation.\n"
            "   print *, 'a;b' ! another ; here\n"
            "   x = 1; y = 2\n"
            "end subroutine show\n"
        )
        self.assertEqual(["F002"], [item.code for item in findings])

    def test_procedure_contract_diagnostics(self):
        """Diagnose missing docs, unnamed ends, and output function arguments."""
        findings, _ = self.audit(
            "function update(x) result(value)\n"
            "   real(dp), intent(inout) :: x\n"
            "   real(dp) :: value\n"
            "   value = x\n"
            "end function\n"
        )
        self.assertEqual({"F005", "F007", "F010"}, {item.code for item in findings})

    def test_continued_header_doc_comment(self):
        """Locate docs after a continued procedure header."""
        findings, _ = self.audit(
            "pure function add_values(first, &\n"
            "   second) result(value)\n"
            "   ! Add two values.\n"
            "   real(dp), intent(in) :: first, second\n"
            "   real(dp) :: value\n"
            "   value = first + second\n"
            "end function add_values\n"
        )
        self.assertEqual([], findings)

    def test_double_precision_aliases(self):
        """Reject explicit double kinds and package-specific aliases."""
        findings, _ = self.audit(
            "module bad_kinds\n"
            "   integer, parameter :: package_dp = dp\n"
            "   real(real64) :: first\n"
            "   real(package_dp) :: second\n"
            "end module bad_kinds\n"
        )
        self.assertEqual(["F009", "F008", "F008"], [item.code for item in findings])

    def test_cmake_coverage(self):
        """Diagnose missing and stale library and test entries."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sources = [root / "one.f90", root / "test_one.f90"]
            for source in sources:
                source.write_text("", encoding="utf-8")
            (root / "CMakeLists.txt").write_text(
                "add_library(time_series STATIC two.f90)\n"
                "foreach(test_name IN ITEMS two)\n",
                encoding="utf-8",
            )
            findings = AUDIT.cmake_findings(root, sources)
        self.assertEqual({"C001", "C002", "C003", "C004"}, {item.code for item in findings})

    def test_purity_suggestions_are_conservative(self):
        """Suggest a simple candidate but exclude random-number procedures."""
        _, procedures = self.audit(
            "function simple_value(x) result(value)\n"
            "   ! Return a transformed value.\n"
            "   real(dp), intent(in) :: x\n"
            "   real(dp) :: value\n"
            "   value = 2.0_dp*x\n"
            "end function simple_value\n"
            "function random_value() result(value)\n"
            "   ! Draw a random value.\n"
            "   real(dp) :: value\n"
            "   value = random_standard_normal()\n"
            "end function random_value\n"
        )
        suggestions = AUDIT.purity_suggestions("sample.f90", procedures)
        self.assertEqual(1, len(suggestions))
        self.assertIn("simple_value", suggestions[0].message)


if __name__ == "__main__":
    unittest.main()
