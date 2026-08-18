#!/usr/bin/env python3
"""Run loaded increments from every hex case using its production dictionaries."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import tempfile
from pathlib import Path

from validate_solids4foam import ignore_generated


ROOT = Path(__file__).resolve().parents[1]


def run_case(executable: str, case: Path, timeout: int) -> tuple[int, str, bool]:
    process = subprocess.Popen(
        [executable],
        cwd=case,
        env=os.environ.copy(),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        start_new_session=True,
    )
    try:
        output, _ = process.communicate(timeout=timeout)
        return process.returncode, output, False
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
        return process.returncode, output, True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--canonical-only", action="store_true")
    args = parser.parse_args()
    if args.steps < 1:
        parser.error("--steps must be positive")
    executable = shutil.which("solids4Foam")
    if executable is None:
        parser.error("solids4Foam is not available in PATH")

    case_sets = [("canonical_mesh3", ROOT / "cases")]
    if not args.canonical_only:
        case_sets.extend(
            (f"mesh{level}", ROOT / "refinements" / f"mesh{level}")
            for level in range(1, 7)
        )

    report = {}
    with tempfile.TemporaryDirectory(prefix="arostica_loaded_hex_suite_") as temporary:
        temporary = Path(temporary)
        for label, case_root in case_sets:
            report[label] = {}
            for geometry in ("monoventricle", "biventricle"):
                source = case_root / geometry
                target = temporary / label / geometry
                shutil.copytree(source, target, ignore=ignore_generated)
                control = target / "system" / "controlDict"
                control_text = control.read_text()
                match = re.search(r"(?m)^deltaT\s+([^;]+);", control_text)
                if match is None:
                    raise RuntimeError(f"deltaT missing in {control}")
                delta_t = float(match.group(1))
                control_text = re.sub(
                    r"(?m)^endTime\s+[^;]+;",
                    f"endTime {delta_t * args.steps:.16g};",
                    control_text,
                )
                control.write_text(control_text)
                code, log, timed_out = run_case(executable, target, args.timeout)
                (source / "log.loaded.hexSetup").write_text(log)
                norms = [
                    float(value)
                    for value in re.findall(r"SNES Function norm ([0-9.eE+-]+)", log)
                ]
                time_values = re.findall(r"(?m)^Time = ([0-9.eE+-]+)$", log)
                completed = code == 0 and "End" in log and "FOAM FATAL" not in log
                mesh_failure = any(
                    token in log.lower()
                    for token in (
                        "negative cell",
                        "negative jacobian",
                        "invalid supplied aróstica",
                    )
                ) or re.search(r"(^|[^a-z])nan([^a-z]|$)", log.lower()) is not None
                report[label][geometry] = {
                    "case": str(source.relative_to(ROOT)),
                    "cell_count": json.loads((source / "meshMetadata.json").read_text())[
                        "cell_count"
                    ],
                    "delta_t": delta_t,
                    "requested_steps": args.steps,
                    "return_code": code,
                    "timed_out": timed_out,
                    "completed": completed,
                    "last_time": float(time_values[-1]) if time_values else None,
                    "nonlinear_converged_count": log.count("Nonlinear solve converged"),
                    "initial_snes_norm": norms[0] if norms else None,
                    "final_reported_snes_norm": norms[-1] if norms else None,
                    "mesh_or_nan_signature": mesh_failure,
                }
                print(
                    f"{label}/{geometry}: completed={completed}, "
                    f"finalNorm={report[label][geometry]['final_reported_snes_norm']}",
                    flush=True,
                )

    report["case_count"] = 2 * len(case_sets)
    report["all_passed"] = all(
        report[label][geometry]["completed"]
        and report[label][geometry]["nonlinear_converged_count"] == args.steps
        and not report[label][geometry]["mesh_or_nan_signature"]
        for label, _ in case_sets
        for geometry in ("monoventricle", "biventricle")
    )
    destination = ROOT / "reports" / "loaded_hex_cases.json"
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"loaded hex validation: {report['case_count']} cases, "
        f"all_passed={report['all_passed']}",
        flush=True,
    )
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
