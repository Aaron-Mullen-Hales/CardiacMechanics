#!/usr/bin/env python3
"""Run one zero-load solids4foam step on temporary copies of both cases."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def ignore_generated(directory, names):
    ignored = []
    for name in names:
        if name == "VTK" or name.startswith("log."):
            ignored.append(name)
            continue
        try:
            float(name)
        except ValueError:
            continue
        if name != "0":
            ignored.append(name)
    return ignored


def validate_case_root(cases_root: Path, log_name: str = "log.smoke.clean") -> dict:
    """Instantiate both cases below *cases_root* without changing their loads."""

    executable = shutil.which("solids4Foam")
    if executable is None:
        raise RuntimeError("solids4Foam is not available in PATH")
    cases_root = cases_root.resolve()
    results = {}
    with tempfile.TemporaryDirectory(prefix="hex_meshes_solids4foam_") as temporary:
        temporary = Path(temporary)
        for name in ("monoventricle", "biventricle"):
            source = cases_root / name
            if not source.is_dir():
                raise FileNotFoundError(f"case directory does not exist: {source}")
            target = temporary / name
            shutil.copytree(
                source,
                target,
                ignore=ignore_generated,
            )
            # Registration/compatibility smoke: remove chamber load in the
            # temporary copy only.  The checked-in case retains the benchmark
            # pressure series; loaded behaviour is diagnosed separately.
            d_path = target / "0" / "D"
            d_text = d_path.read_text()
            d_text = re.sub(
                r'pressureSeries\s*\{[^}]*\}',
                'pressure uniform 0',
                d_text,
            )
            d_path.write_text(d_text)
            process = subprocess.run(
                [executable],
                cwd=target,
                env=os.environ.copy(),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            log = process.stdout
            (source / log_name).write_text(log)
            results[name] = {
                "return_code": process.returncode,
                "completed": process.returncode == 0 and "FOAM FATAL" not in log and "End" in log,
                "dynamic_mesh_selected": "Selecting dynamicFvMesh staticFvMesh" in log,
                "mechanical_model_created": "Creating the mechanicalModel" in log,
                "solid_model_selected":
                    "Selecting solidModel nonLinearGeometryTotalLagrangianTotalDisplacement" in log,
                "arostica_viscoelastic_law_selected":
                    "Selecting mechanical law ArosticaHolzapfelOgdenViscoelastic" in log,
                "clean_activation_wrapper_selected":
                    "Selecting mechanical law electroMechanicalLaw" in log,
                "mixed_pressure_enabled": "solvePressure = true" in log,
                "cell_and_face_triads_validated": log.count("determinant range =") == 2,
                "spring_dashpot_conditions_created": log.count("Creating solidTraction on") == 2,
                "snes_entered": "Solving the momentum equation for D using PETSc SNES" in log,
                "zero_load_step_converged": "Nonlinear solve converged" in log,
                "temporary_test_only": True,
            }
            print(
                f"{cases_root.name}/{name}: "
                f"solids4foam completed={results[name]['completed']}"
            )
    results["all_passed"] = all(
        all(value for key, value in results[name].items() if key not in ("return_code",))
        for name in ("monoventricle", "biventricle")
    )
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases-root", type=Path, default=ROOT / "cases")
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "reports" / "solids4foam_compatibility.json",
    )
    parser.add_argument("--log-name", default="log.smoke.clean")
    args = parser.parse_args()
    try:
        results = validate_case_root(args.cases_root, args.log_name)
    except (FileNotFoundError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(results, indent=2) + "\n")
    return 0 if results["all_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
