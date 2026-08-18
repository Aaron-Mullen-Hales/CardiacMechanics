#!/usr/bin/env python3
"""Run one pressure-loaded increment and classify solver versus setup failure."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
from pathlib import Path

from validate_solids4foam import ignore_generated


ROOT = Path(__file__).resolve().parents[1]


def run_with_timeout(executable, target, timeout=120):
    process = subprocess.Popen(
        [executable],
        cwd=target,
        env=os.environ.copy(),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        start_new_session=True,
    )
    timed_out = False
    try:
        output, _ = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
    return process.returncode, output, timed_out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--classify-existing", action="store_true")
    args = parser.parse_args()
    executable = shutil.which("solids4Foam")
    if executable is None:
        print("ERROR: solids4Foam is not available in PATH", file=sys.stderr)
        return 2
    results = {}
    with tempfile.TemporaryDirectory(prefix="arostica_hex_loaded_") as temporary:
        temporary = Path(temporary)
        for name in ("monoventricle", "biventricle"):
            source = ROOT / "cases" / name
            target = temporary / name
            shutil.copytree(source, target, ignore=ignore_generated)
            # Bound diagnostic cost.  The production file remains at 1000; a
            # 40-iteration cap is enough to expose the approximate-Jacobian
            # preconditioner limitation without launching a long campaign.
            options = target / "petscOptions.hex"
            text = re.sub(r"-ksp_max_it\s+\d+", "-ksp_max_it 40", options.read_text())
            options.write_text(text)
            if args.classify_existing:
                log = (source / "log.smoke.loaded").read_text()
                return_code = 0 if "End" in log and "FOAM FATAL" not in log else 1
                timed_out = False
            else:
                return_code, log, timed_out = run_with_timeout(executable, target)
                (source / "log.smoke.loaded").write_text(log)
            norms = [float(value) for value in re.findall(r"SNES Function norm ([0-9.eE+-]+)", log)]
            entered = "Solving the momentum equation for D using PETSc SNES" in log
            setup_ok = all(token in log for token in (
                "Selecting mechanical law electroMechanicalLaw",
                "Selecting mechanical law ArosticaHolzapfelOgdenViscoelastic",
                "solvePressure = true",
            ))
            lower_log = log.lower()
            mesh_failure = any(token in lower_log for token in (
                "negative cell", "negative jacobian", "invalid supplied aróstica",
            )) or re.search(r"(^|[^a-z])nan([^a-z]|$)", lower_log) is not None
            linear_limit = any(token in log for token in (
                "DIVERGED_LINEAR_SOLVE", "Linear solve did not converge due to DIVERGED_ITS",
            ))
            completed = return_code == 0 and "End" in log
            classification = (
                "completed" if completed else
                "solver_preconditioner_limited" if entered and setup_ok and linear_limit and not mesh_failure else
                "timeout_after_snes_entry" if timed_out and entered and setup_ok and not mesh_failure else
                "mesh_or_field_failure" if mesh_failure else
                "setup_or_unclassified_failure"
            )
            results[name] = {
                "return_code": return_code,
                "timed_out": timed_out,
                "entered_snes": entered,
                "setup_and_fields_instantiated": setup_ok,
                "initial_snes_residual": norms[0] if norms else None,
                "lowest_reported_snes_residual": min(norms) if norms else None,
                "reported_snes_iterations": max(0, len(norms) - 1),
                "mesh_or_nan_signature": mesh_failure,
                "linear_iteration_limit_signature": linear_limit,
                "classification": classification,
                "diagnostic_ksp_max_it": 40,
            }
            print(f"{name}: {classification}, norms={norms}")
    acceptable = {"completed", "solver_preconditioner_limited", "timeout_after_snes_entry"}
    results["no_mesh_or_setup_failure"] = all(
        results[name]["classification"] in acceptable
        for name in ("monoventricle", "biventricle")
    )
    (ROOT / "reports" / "loaded_smoke.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0 if results["no_mesh_or_setup_failure"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
