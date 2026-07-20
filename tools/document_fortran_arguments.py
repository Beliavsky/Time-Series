#!/usr/bin/env python3
"""Add consistent FORD comments to Fortran procedure arguments."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Sequence

from fortran_style_audit import (
    DECL_RE,
    Procedure,
    entity_name,
    logical_statements,
    parse_procedures,
    split_code_comment,
    split_top_level,
)


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
    "build-varshrink",
}

EXACT_DESCRIPTIONS = {
    "alpha": "Significance, smoothing, or model coefficient.",
    "ar": "Autoregressive coefficients.",
    "ar_order": "Autoregressive order.",
    "ar_polynomial": "Autoregressive polynomial coefficients.",
    "beta": "Regression or model coefficients.",
    "burn": "Number of initial simulation draws to discard.",
    "burnin": "Number of initial simulation draws to discard.",
    "coefficients": "Model coefficients.",
    "components": "Model components.",
    "convergence_tolerance": "Convergence tolerance.",
    "covariance": "Covariance matrix.",
    "d": "Fractional-differencing parameter or differencing order.",
    "degrees_of_freedom": "Degrees of freedom.",
    "delta": "Model increment or differencing parameter.",
    "distribution": "Probability-distribution specification.",
    "endogenous": "Endogenous time-series observations.",
    "exogenous": "Exogenous predictor observations.",
    "fit": "Previously fitted model.",
    "first": "First operand.",
    "forecast": "Forecast values.",
    "frequency": "Number of observations per seasonal cycle.",
    "gamma": "Model coefficient or scale parameter.",
    "horizon": "Number of periods to forecast.",
    "include_intercept": "Whether to include an intercept.",
    "include_mean": "Whether to include a mean term.",
    "index": "Element or observation index.",
    "info": "Status code; zero indicates success.",
    "initial": "Initial value.",
    "initial_covariance": "Initial state covariance matrix.",
    "initial_mean": "Initial state mean.",
    "initial_parameters": "Initial parameter values.",
    "initial_state": "Initial state vector.",
    "initial_variance": "Initial variance.",
    "innovations": "Model innovations.",
    "intercept": "Model intercept.",
    "iterations": "Number of algorithm iterations.",
    "lag": "Lag index or number of lags.",
    "lag_order": "Model lag order.",
    "lambda": "Penalty or shrinkage parameter.",
    "lambdas": "Candidate penalty or shrinkage parameters.",
    "level": "Model level or confidence level.",
    "log_density": "Log-density value.",
    "log_likelihood": "Log-likelihood value.",
    "ma": "Moving-average coefficients.",
    "ma_order": "Moving-average order.",
    "ma_polynomial": "Moving-average polynomial coefficients.",
    "matrix": "Input matrix.",
    "max_iterations": "Maximum number of algorithm iterations.",
    "max_lag": "Maximum lag to consider.",
    "mean": "Mean value or vector.",
    "method": "Algorithm or estimation method.",
    "model": "Model specification.",
    "n": "Number of observations or elements.",
    "normal_draws": "Independent standard-normal draws.",
    "normals": "Independent standard-normal draws.",
    "observation": "Observed value or vector.",
    "observation_loading": "Observation loading matrix.",
    "observation_variance": "Observation-error variance.",
    "observations": "Observed time-series values.",
    "offset": "Known additive offset.",
    "order": "Model or polynomial order.",
    "options": "Algorithm options.",
    "out": "Procedure result.",
    "p": "Autoregressive order or model dimension.",
    "parameters": "Model parameter values.",
    "particles": "Number of particles.",
    "period": "Seasonal period.",
    "phi": "Autoregressive or model coefficient.",
    "predictors": "Predictor matrix.",
    "predicted": "Predicted values.",
    "prior": "Prior-distribution specification.",
    "probabilities": "Probability values.",
    "probability": "Probability value.",
    "proposal_normals": "Standard-normal proposal draws.",
    "q": "Model order, dimension, or parameter.",
    "rank": "Matrix or cointegration rank.",
    "regressors": "Regression design matrix.",
    "residuals": "Model residuals.",
    "response": "Response observations.",
    "seed": "Random-number seed.",
    "second": "Second operand.",
    "series": "Time-series observations.",
    "sigma": "Scale parameter or standard deviation.",
    "simulations": "Number of simulation draws.",
    "state": "State vector or state sequence.",
    "state_normal_draws": "Independent standard-normal state draws.",
    "structure": "Model-structure specification.",
    "threshold": "Decision or truncation threshold.",
    "training": "Training observations.",
    "time": "Observation times.",
    "tolerance": "Numerical convergence tolerance.",
    "transition": "State transition matrix.",
    "u": "Input vector or random variate.",
    "value": "Input value.",
    "values": "Input values.",
    "variables": "Number or indices of variables.",
    "variance": "Variance value or matrix.",
    "weights": "Observation or objective weights.",
    "x": "Input data or predictor values.",
    "y": "Response or time-series observations.",
    "actual": "Observed values used for evaluation.",
    "bandwidth": "Smoothing or spectral bandwidth.",
    "denominator": "Denominator polynomial coefficients.",
    "numerator": "Numerator polynomial coefficients.",
    "pgram": "Periodogram values and frequencies.",
    "transforms": "Transformation specifications.",
}

WORD_REPLACEMENTS = {
    "acf": "autocorrelation",
    "ar": "autoregressive",
    "arma": "ARMA",
    "arima": "ARIMA",
    "cdf": "cumulative distribution",
    "df": "degrees of freedom",
    "garch": "GARCH",
    "hessian": "Hessian",
    "id": "identifier",
    "iekf": "iterated extended Kalman filter",
    "irf": "impulse-response",
    "ll": "log-likelihood",
    "ma": "moving-average",
    "mcmc": "MCMC",
    "pdf": "probability density",
    "rng": "random-number generator",
    "sd": "standard deviation",
    "sv": "stochastic-volatility",
    "var": "VAR",
}

DECLARATION_PREFIX_RE = re.compile(
    r"^(?:integer|real|logical|complex|character|double\s+precision|"
    r"type\s*\(|class\s*\(|procedure\s*\()",
    re.IGNORECASE,
)


def humanize(name: str) -> str:
    """Convert an identifier into a short documentation phrase."""
    words = name.strip("_").split("_")
    rendered = [WORD_REPLACEMENTS.get(word, word) for word in words]
    return " ".join(rendered)


def attributes_intent(attributes: str) -> str:
    """Return the INTENT value from a declaration attribute list."""
    match = re.search(r"\bintent\s*\(\s*(inout|in|out)\s*\)", attributes, re.I)
    return match.group(1).lower() if match else ""


def describe_argument(name: str, attributes: str) -> str:
    """Create a conservative argument description from its declaration."""
    if name in EXACT_DESCRIPTIONS:
        description = EXACT_DESCRIPTIONS[name]
    elif name.startswith("include_"):
        description = f"Whether to include the {humanize(name[8:])}."
    elif name.startswith("use_"):
        description = f"Whether to use the {humanize(name[4:])}."
    elif name.startswith("return_"):
        description = f"Whether to return the {humanize(name[7:])}."
    elif name.startswith("estimate_"):
        description = f"Whether to estimate the {humanize(name[9:])}."
    elif name.startswith("max_") or name.startswith("maximum_"):
        description = f"Maximum {humanize(name.removeprefix('max_').removeprefix('maximum_'))}."
    elif name.startswith("min_") or name.startswith("minimum_"):
        description = f"Minimum {humanize(name.removeprefix('min_').removeprefix('minimum_'))}."
    elif name.startswith("initial_"):
        description = f"Initial {humanize(name[8:])}."
    elif name.startswith("lower_"):
        description = f"Lower bound for {humanize(name[6:])}."
    elif name.startswith("upper_"):
        description = f"Upper bound for {humanize(name[6:])}."
    elif name.startswith("first_"):
        description = f"First {humanize(name[6:])}."
    elif name.startswith("last_"):
        description = f"Last {humanize(name[5:])}."
    elif name.startswith("prior_"):
        description = f"Prior {humanize(name[6:])}."
    elif name.startswith("observation_"):
        description = f"Observation {humanize(name[12:])}."
    elif name.startswith("state_"):
        description = f"State {humanize(name[6:])}."
    elif name.startswith("proposal_"):
        description = f"Proposal {humanize(name[9:])}."
    elif name.startswith("future_"):
        description = f"Future {humanize(name[7:])}."
    elif name.endswith("_order"):
        description = f"{humanize(name[:-6]).capitalize()} order."
    elif name.endswith("_variance"):
        description = f"{humanize(name[:-9]).capitalize()} variance."
    elif name.endswith("_covariance"):
        description = f"{humanize(name[:-11]).capitalize()} covariance matrix."
    elif name.endswith("_matrix"):
        description = f"{humanize(name[:-7]).capitalize()} matrix."
    elif name.endswith("_polynomial"):
        description = f"{humanize(name[:-11]).capitalize()} polynomial coefficients."
    elif name.endswith("_draws"):
        description = f"{humanize(name[:-6]).capitalize()} simulation draws."
    elif name.endswith("_probability"):
        description = f"{humanize(name[:-12]).capitalize()} probability."
    elif name.endswith("_probabilities"):
        description = f"{humanize(name[:-14]).capitalize()} probabilities."
    elif name.endswith("_df"):
        description = f"{humanize(name[:-3]).capitalize()} degrees of freedom."
    elif name.endswith("_tolerance"):
        description = f"{humanize(name[:-10]).capitalize()} tolerance."
    elif name.endswith("_iterations"):
        description = f"Number of {humanize(name[:-11])} iterations."
    elif name.endswith("_count"):
        description = f"Number of {humanize(name[:-6])}."
    elif name.endswith("_index"):
        description = f"Index of {humanize(name[:-6])}."
    elif name.endswith("_lag"):
        description = f"{humanize(name[:-4]).capitalize()} lag."
    elif name.startswith("n_"):
        description = f"Number of {humanize(name[2:])}."
    elif "procedure" in attributes.lower():
        description = f"{humanize(name).capitalize()} callback procedure."
    elif re.search(r"\blogical\b", attributes, re.I):
        description = f"Flag controlling {humanize(name)}."
    else:
        description = f"{humanize(name).capitalize()}."

    intent = attributes_intent(attributes)
    if intent == "inout":
        description = description[:-1] + ", updated in place."
    return description


def immediate_doc_block(lines: Sequence[str], header_end: int) -> range:
    """Return the comment block immediately following a procedure header."""
    index = header_end
    while index < len(lines) and not lines[index].strip():
        index += 1
    start = index
    while index < len(lines) and lines[index].lstrip().startswith("!"):
        index += 1
    return range(start, index)


def declaration_replacements(
    lines: Sequence[str], procedures: Sequence[Procedure]
) -> dict[int, tuple[int, list[str]]]:
    """Build source-line replacements for dummy-argument declarations."""
    replacements: dict[int, tuple[int, list[str]]] = {}
    for procedure in procedures:
        nested_ranges = [
            (other.start, other.end)
            for other in procedures
            if procedure.start < other.start and other.end <= procedure.end
        ]
        for statement in procedure.statements:
            if any(start <= statement.start <= end for start, end in nested_ranges):
                continue
            match = DECL_RE.match(statement.code)
            if not match:
                continue
            if not DECLARATION_PREFIX_RE.match(match.group("attrs").strip()):
                continue
            entities = split_top_level(match.group("entities"))
            dummy_entities = [
                entity for entity in entities if entity_name(entity) in procedure.args
            ]
            if not dummy_entities:
                continue

            first_line = lines[statement.start - 1]
            indent = first_line[: len(first_line) - len(first_line.lstrip())]
            attributes = match.group("attrs").strip()
            local_entities = [
                entity for entity in entities if entity_name(entity) not in procedure.args
            ]
            new_lines = []
            for entity in dummy_entities:
                name = entity_name(entity)
                comment = describe_argument(name, attributes)
                new_lines.append(
                    f"{indent}{attributes} :: {entity.strip()} !! {comment}"
                )
            if local_entities:
                new_lines.append(
                    f"{indent}{attributes} :: {', '.join(local_entities)}"
                )
            replacements[statement.start - 1] = (statement.end, new_lines)
    return replacements


def document_file(path: Path) -> bool:
    """Document procedure arguments in one Fortran source file."""
    original = path.read_text(encoding="utf-8")
    trailing_newline = original.endswith("\n")
    lines = original.splitlines()
    procedures = parse_procedures(logical_statements(lines))
    replacements = declaration_replacements(lines, procedures)

    for procedure in procedures:
        for index in immediate_doc_block(lines, procedure.header_end):
            stripped = lines[index].lstrip()
            if stripped.startswith("!!"):
                continue
            indent = lines[index][: len(lines[index]) - len(stripped)]
            replacements[index] = (index + 1, [f"{indent}!!{stripped[1:]}"])

    output: list[str] = []
    index = 0
    while index < len(lines):
        replacement = replacements.get(index)
        if replacement is None:
            output.append(lines[index])
            index += 1
        else:
            end, new_lines = replacement
            output.extend(new_lines)
            index = end

    updated = "\n".join(output) + ("\n" if trailing_newline else "")
    if updated == original:
        return False
    path.write_text(updated, encoding="utf-8", newline="\n")
    return True


def source_paths(arguments: Sequence[str]) -> list[Path]:
    """Expand command-line files and directories into Fortran sources."""
    candidates = [Path(argument) for argument in arguments] if arguments else [Path(".")]
    paths: set[Path] = set()
    for candidate in candidates:
        if candidate.is_file() and candidate.suffix.lower() == ".f90":
            paths.add(candidate)
        elif candidate.is_dir():
            for path in candidate.rglob("*.f90"):
                if not any(part in EXCLUDED_DIRS for part in path.parts):
                    paths.add(path)
    return sorted(paths)


def main() -> int:
    """Document selected Fortran sources."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", help="Fortran files or source directories")
    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        help="file name to exclude; may be specified more than once",
    )
    args = parser.parse_args()
    excluded = set(args.exclude)
    changed = 0
    for path in source_paths(args.paths):
        if path.name in excluded:
            continue
        if document_file(path):
            print(path)
            changed += 1
    print(f"Updated {changed} Fortran files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
