#!/usr/bin/env python3
"""Audit the local Fortran style conventions without external dependencies.

The quote-aware scanner design is adapted from
``C:\\python\\fortran\\f2c\\fortran_scan.py``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Sequence


FORTRAN_SUFFIXES = {".f", ".for", ".f90", ".f95", ".f03", ".f08"}
EXCLUDED_DIRS = {
    ".git",
    ".pytest_cache",
    "__pycache__",
    "build",
    "build-debug",
    "build-release",
    "build-verify",
    "build-pure",
    "build-ninja",
    "build-pure-ninja",
}
PROC_START_RE = re.compile(
    r"^\s*(?P<prefix>(?:(?:pure|elemental|impure|recursive|module)\s+)*"
    r"(?:(?:integer|real(?:\s*\([^)]*\))?|logical|complex(?:\s*\([^)]*\))?|"
    r"character(?:\s*\([^)]*\))?|double\s+precision|type\s*\([^)]*\)|"
    r"class\s*\([^)]*\))\s+)?"
    r"(?:(?:pure|elemental|impure|recursive|module)\s+)*)"
    r"(?P<kind>function|subroutine)\s+(?P<name>[a-z][a-z0-9_]*)\s*"
    r"\((?P<args>[^)]*)\)(?:\s*result\s*\(\s*(?P<result>[a-z][a-z0-9_]*)\s*\))?",
    re.IGNORECASE,
)
PROC_END_RE = re.compile(
    r"^\s*end\s+(?P<kind>function|subroutine)(?:\s+(?P<name>[a-z][a-z0-9_]*))?\s*$",
    re.IGNORECASE,
)
DECL_RE = re.compile(r"^\s*(?P<attrs>.+?)\s*::\s*(?P<entities>.+)$", re.IGNORECASE)


@dataclass(order=True)
class Finding:
    """One deterministic audit finding."""

    path: str
    line: int
    severity: str
    code: str
    message: str


@dataclass
class LogicalStatement:
    """One free-form Fortran statement assembled across continuations."""

    start: int
    end: int
    code: str


@dataclass
class Procedure:
    """Procedure metadata needed by the local style checks."""

    name: str
    kind: str
    start: int
    header_end: int
    args: set[str]
    result: str | None
    attrs: set[str]
    statements: list[LogicalStatement] = field(default_factory=list)
    end: int = 0


def split_code_comment(line: str) -> tuple[str, str]:
    """Split a Fortran line at an unquoted exclamation mark."""
    quote: str | None = None
    index = 0
    while index < len(line):
        char = line[index]
        if quote is None:
            if char in {"'", '"'}:
                quote = char
            elif char == "!":
                return line[:index], line[index:]
        elif char == quote:
            if index + 1 < len(line) and line[index + 1] == quote:
                index += 1
            else:
                quote = None
        index += 1
    return line, ""


def unquoted_semicolon_column(line: str) -> int | None:
    """Return the first statement-separating semicolon column, if present."""
    code, _ = split_code_comment(line)
    quote: str | None = None
    index = 0
    while index < len(code):
        char = code[index]
        if quote is None:
            if char in {"'", '"'}:
                quote = char
            elif char == ";":
                return index + 1
        elif char == quote:
            if index + 1 < len(code) and code[index + 1] == quote:
                index += 1
            else:
                quote = None
        index += 1
    return None


def logical_statements(lines: Sequence[str]) -> list[LogicalStatement]:
    """Assemble free-form continuation lines into logical statements."""
    statements: list[LogicalStatement] = []
    pieces: list[str] = []
    start = 0
    continuing = False
    for number, raw in enumerate(lines, 1):
        code, _ = split_code_comment(raw)
        stripped = code.strip()
        if not stripped and not continuing:
            continue
        if not continuing:
            start = number
            pieces = []
        if continuing:
            stripped = stripped.lstrip()
            if stripped.startswith("&"):
                stripped = stripped[1:].lstrip()
        has_continuation = stripped.endswith("&")
        if has_continuation:
            stripped = stripped[:-1].rstrip()
        pieces.append(stripped)
        continuing = has_continuation
        if not continuing:
            statements.append(LogicalStatement(start, number, " ".join(pieces)))
            pieces = []
    if pieces:
        statements.append(LogicalStatement(start, len(lines), " ".join(pieces)))
    return statements


def split_top_level(text: str) -> list[str]:
    """Split a comma list while respecting parentheses and quoted strings."""
    items: list[str] = []
    start = 0
    depth = 0
    quote: str | None = None
    index = 0
    while index < len(text):
        char = text[index]
        if quote is None:
            if char in {"'", '"'}:
                quote = char
            elif char == "(":
                depth += 1
            elif char == ")":
                depth = max(0, depth - 1)
            elif char == "," and depth == 0:
                items.append(text[start:index].strip())
                start = index + 1
        elif char == quote:
            if index + 1 < len(text) and text[index + 1] == quote:
                index += 1
            else:
                quote = None
        index += 1
    items.append(text[start:].strip())
    return [item for item in items if item]


def parse_procedures(statements: Sequence[LogicalStatement]) -> list[Procedure]:
    """Parse named function and subroutine blocks from logical statements."""
    procedures: list[Procedure] = []
    stack: list[Procedure] = []
    for statement in statements:
        code = statement.code.strip()
        end_match = PROC_END_RE.match(code)
        if end_match:
            if stack:
                procedure = stack.pop()
                procedure.end = statement.end
            continue
        if code.lower().startswith("module procedure"):
            continue
        start_match = PROC_START_RE.match(code)
        if start_match:
            prefix = start_match.group("prefix").lower()
            attrs = {
                word
                for word in ("pure", "elemental", "impure", "recursive", "module")
                if re.search(rf"\b{word}\b", prefix)
            }
            args = {
                item.strip().lower()
                for item in start_match.group("args").split(",")
                if item.strip()
            }
            procedure = Procedure(
                name=start_match.group("name").lower(),
                kind=start_match.group("kind").lower(),
                start=statement.start,
                header_end=statement.end,
                args=args,
                result=(start_match.group("result") or "").lower() or None,
                attrs=attrs,
            )
            procedures.append(procedure)
            stack.append(procedure)
            continue
        for procedure in stack:
            procedure.statements.append(statement)
    return procedures


def entity_name(entity: str) -> str:
    """Return the declared identifier at the start of an entity declaration."""
    match = re.match(r"\s*([a-z][a-z0-9_]*)", entity, re.IGNORECASE)
    return match.group(1).lower() if match else ""


def dummy_declarations(procedure: Procedure) -> dict[str, tuple[str, str]]:
    """Map dummy names to their declaration attributes and entity text."""
    declarations: dict[str, tuple[str, str]] = {}
    for statement in procedure.statements:
        match = DECL_RE.match(statement.code)
        if not match:
            continue
        attrs = match.group("attrs").lower()
        for entity in split_top_level(match.group("entities")):
            name = entity_name(entity)
            if name in procedure.args or name == procedure.result:
                declarations[name] = (attrs, entity.lower())
    return declarations


def has_short_doc(lines: Sequence[str], header_end: int) -> bool:
    """Check for a comment immediately after a procedure header."""
    index = header_end
    while index < len(lines) and not lines[index].strip():
        index += 1
    return index < len(lines) and lines[index].lstrip().startswith("!")


def audit_file(path: Path, display: str, max_length: int) -> tuple[list[Finding], list[Procedure]]:
    """Run deterministic style checks on one Fortran source file."""
    findings: list[Finding] = []
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = path.read_text(encoding="cp1252")
    lines = text.splitlines()
    if not lines or not re.fullmatch(r"!\s*SPDX-License-Identifier:\s*\S(?:.*\S)?", lines[0]):
        findings.append(Finding(display, 1, "error", "F011", "line 1 must contain an SPDX license identifier"))
    if len(lines) < 2 or not re.fullmatch(r"!\s*SPDX-FileComment:\s*\S(?:.*\S)?", lines[1]):
        findings.append(Finding(display, 2, "error", "F012", "line 2 must contain an SPDX provenance comment"))
    for number, line in enumerate(lines, 1):
        if len(line) > max_length:
            findings.append(Finding(display, number, "error", "F001", f"line exceeds {max_length} columns"))
        column = unquoted_semicolon_column(line)
        if column is not None:
            findings.append(Finding(display, number, "error", "F002", f"semicolon separates statements at column {column}"))
        if "\t" in line:
            findings.append(Finding(display, number, "error", "F003", "tab character"))
        if line.rstrip() != line:
            findings.append(Finding(display, number, "error", "F004", "trailing whitespace"))
        code, _ = split_code_comment(line)
        if re.search(
            r"\bdouble\s+precision\b|\breal\s*\*\s*8\b|"
            r"\breal\s*\(\s*(?:kind\s*=\s*)?(?:8|real64|[a-z][a-z0-9_]*_dp)\s*\)",
            code,
            re.IGNORECASE,
        ):
            findings.append(Finding(display, number, "error", "F008", "use real(dp) for double precision"))
        alias = re.search(
            r"\binteger\b[^!]*\bparameter\b[^!]*::[^!]*\b([a-z][a-z0-9_]*_dp)\s*=",
            code,
            re.IGNORECASE,
        )
        if alias and alias.group(1).lower() != "dp":
            findings.append(Finding(display, number, "error", "F009", f"use shared dp instead of alias {alias.group(1)}"))

    statements = logical_statements(lines)
    procedures = parse_procedures(statements)
    stack: list[tuple[str, str, int]] = []
    for statement in statements:
        code = statement.code.strip()
        end_match = PROC_END_RE.match(code)
        if end_match:
            end_name = (end_match.group("name") or "").lower()
            if not end_name:
                findings.append(Finding(display, statement.start, "error", "F005", f"unnamed end {end_match.group('kind').lower()}"))
            if stack:
                kind, name, _ = stack.pop()
                if end_match.group("kind").lower() != kind or (end_name and end_name != name):
                    findings.append(Finding(display, statement.start, "error", "F006", f"procedure end does not match {kind} {name}"))
            continue
        if code.lower().startswith("module procedure"):
            continue
        start_match = PROC_START_RE.match(code)
        if start_match:
            stack.append((start_match.group("kind").lower(), start_match.group("name").lower(), statement.start))

    for kind, name, start in stack:
        findings.append(Finding(display, start, "error", "F006", f"{kind} {name} has no matching end"))

    for procedure in procedures:
        if not has_short_doc(lines, procedure.header_end):
            findings.append(Finding(display, procedure.start, "error", "F007", f"{procedure.kind} {procedure.name} lacks a following comment"))
        if procedure.kind == "function":
            for name, (attrs, _) in dummy_declarations(procedure).items():
                if name in procedure.args and re.search(r"\bintent\s*\(\s*(?:out|inout)\s*\)", attrs):
                    findings.append(Finding(display, procedure.start, "error", "F010", f"function argument {name} is not intent(in)"))
    return findings, procedures


def purity_suggestions(path: str, procedures: Sequence[Procedure]) -> list[Finding]:
    """Return conservative, non-failing suggestions for purity attributes."""
    findings: list[Finding] = []
    declared_pure = {proc.name for proc in procedures if proc.attrs & {"pure", "elemental"}}
    forbidden = re.compile(
        r"^\s*(?:print\b|read\b|write\b|open\b|close\b|rewind\b|backspace\b|"
        r"stop\b|error\s+stop\b)|\brandom_[a-z0-9_]*\b|\bsave\b",
        re.IGNORECASE,
    )
    call_re = re.compile(r"\bcall\s+([a-z][a-z0-9_]*)\b", re.IGNORECASE)
    intrinsic_pure_calls = {"move_alloc"}
    for procedure in procedures:
        if procedure.attrs & {"pure", "elemental", "impure"} or procedure.end == 0:
            continue
        body = [statement.code for statement in procedure.statements]
        if any(code.strip().lower() == "contains" for code in body):
            continue
        if procedure.name.startswith("random_") or any(forbidden.search(code) for code in body):
            continue
        calls = {match.group(1).lower() for code in body for match in call_re.finditer(code)}
        if calls - declared_pure - intrinsic_pure_calls:
            continue
        if not any(re.match(r"^\s*(?:if\b|do\b|select\b|call\b|[a-z][a-z0-9_%]*(?:\([^=]*\))?\s*=)", code, re.IGNORECASE) for code in body):
            continue
        findings.append(Finding(path, procedure.start, "suggestion", "P001", f"review whether {procedure.kind} {procedure.name} can be pure"))

    for procedure in procedures:
        if "pure" not in procedure.attrs or "elemental" in procedure.attrs or procedure.end == 0:
            continue
        declarations = dummy_declarations(procedure)
        if not procedure.args or not all(name in declarations for name in procedure.args):
            continue
        scalar = True
        for name in procedure.args:
            attrs, entity = declarations[name]
            if "dimension" in attrs or "allocatable" in attrs or "pointer" in attrs or "(" in entity:
                scalar = False
        result_name = procedure.result or (procedure.name if procedure.kind == "function" else None)
        if result_name and result_name in declarations:
            attrs, entity = declarations[result_name]
            if "dimension" in attrs or "allocatable" in attrs or "pointer" in attrs or "(" in entity:
                scalar = False
        elif procedure.kind == "function":
            scalar = False
        if scalar:
            findings.append(Finding(path, procedure.start, "suggestion", "P002", f"review whether pure {procedure.kind} {procedure.name} can be elemental"))
    return findings


def source_paths(inputs: Sequence[Path]) -> list[Path]:
    """Expand input files and directories to a sorted unique source list."""
    paths: set[Path] = set()
    for item in inputs:
        resolved = item.resolve()
        if resolved.is_file() and resolved.suffix.lower() in FORTRAN_SUFFIXES:
            paths.add(resolved)
        elif resolved.is_dir():
            for candidate in resolved.rglob("*"):
                if candidate.is_file() and candidate.suffix.lower() in FORTRAN_SUFFIXES:
                    if not any(
                        part.lower() in EXCLUDED_DIRS
                        or part.lower().startswith("build")
                        or part.startswith(".")
                        for part in candidate.relative_to(resolved).parts[:-1]
                    ):
                        paths.add(candidate.resolve())
    return sorted(paths, key=lambda value: str(value).lower())


def cmake_findings(root: Path, sources: Sequence[Path]) -> list[Finding]:
    """Check that root library and test sources appear in CMakeLists.txt."""
    cmake = root / "CMakeLists.txt"
    if not cmake.is_file():
        return []
    text = cmake.read_text(encoding="utf-8")
    findings: list[Finding] = []
    library_match = re.search(r"add_library\s*\(\s*time_series\s+STATIC(?P<body>.*?)\)", text, re.IGNORECASE | re.DOTALL)
    listed_sources = set(re.findall(r"\b([a-z0-9_]+\.f(?:90|95|03|08|or)?)\b", library_match.group("body"), re.IGNORECASE)) if library_match else set()
    test_match = re.search(r"foreach\s*\(\s*test_name\s+IN\s+ITEMS(?P<body>.*?)\)", text, re.IGNORECASE | re.DOTALL)
    listed_tests = set(re.findall(r"\b[a-z][a-z0-9_]*\b", test_match.group("body"), re.IGNORECASE)) if test_match else set()
    listed_executable_sources = set(
        re.findall(
            r"add_executable\s*\(\s*[a-z0-9_]+\s+([a-z0-9_./-]+\.f(?:90|95|03|08|or)?)",
            text,
            re.IGNORECASE,
        )
    )
    for source in sources:
        if source.parent != root:
            continue
        name = source.name.lower()
        if name.startswith("test_"):
            test_name = source.stem[5:].lower()
            if test_name not in {item.lower() for item in listed_tests}:
                findings.append(Finding("CMakeLists.txt", 1, "error", "C001", f"test {name} is not in the test list"))
        elif name in {Path(item).name.lower() for item in listed_executable_sources}:
            continue
        elif name not in {item.lower() for item in listed_sources}:
            findings.append(Finding("CMakeLists.txt", 1, "error", "C002", f"source {name} is not in time_series"))
    root_sources = [
        source
        for source in root.iterdir()
        if source.is_file() and source.suffix.lower() in FORTRAN_SUFFIXES
    ]
    available_sources = {
        source.name.lower()
        for source in root_sources
        if not source.name.lower().startswith("test_")
    }
    available_tests = {
        source.stem[5:].lower()
        for source in root_sources
        if source.name.lower().startswith("test_")
    }
    for name in sorted({item.lower() for item in listed_sources} - available_sources):
        findings.append(Finding("CMakeLists.txt", 1, "error", "C003", f"listed source {name} does not exist"))
    for name in sorted({item.lower() for item in listed_tests} - available_tests):
        findings.append(Finding("CMakeLists.txt", 1, "error", "C004", f"listed test test_{name}.f90 does not exist"))
    return findings


def display_path(path: Path, root: Path) -> str:
    """Return a stable path relative to the audit root when possible."""
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def print_report(findings: Sequence[Finding], files_checked: int) -> None:
    """Print findings and a compact summary for command-line use."""
    for finding in findings:
        print(f"{finding.path}:{finding.line}: {finding.severity} {finding.code}: {finding.message}")
    errors = sum(finding.severity == "error" for finding in findings)
    suggestions = sum(finding.severity == "suggestion" for finding in findings)
    print(f"Checked {files_checked} Fortran files: {errors} errors, {suggestions} suggestions")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line style audit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path, default=[Path(".")], help="source files or directories")
    parser.add_argument("--max-line-length", type=int, default=132, help="maximum physical line length")
    parser.add_argument("--suggest-purity", action="store_true", help="emit heuristic pure/elemental suggestions")
    parser.add_argument("--no-cmake", action="store_true", help="skip CMake source-list checks")
    parser.add_argument("--json", action="store_true", help="write findings as JSON")
    args = parser.parse_args(argv)
    inputs = args.paths or [Path(".")]
    root = Path.cwd().resolve()
    sources = source_paths(inputs)
    findings: list[Finding] = []
    for source in sources:
        shown = display_path(source, root)
        file_findings, procedures = audit_file(source, shown, args.max_line_length)
        findings.extend(file_findings)
        if args.suggest_purity:
            findings.extend(purity_suggestions(shown, procedures))
    if not args.no_cmake:
        findings.extend(cmake_findings(root, sources))
    findings.sort()
    if args.json:
        print(json.dumps({"files_checked": len(sources), "findings": [asdict(item) for item in findings]}, indent=2))
    else:
        print_report(findings, len(sources))
    return 1 if any(item.severity == "error" for item in findings) else 0


if __name__ == "__main__":
    raise SystemExit(main())
