#!/usr/bin/env python3
"""Validate retained canonical loaded results and deformation quality."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

from fibres import read_internal_scalars, read_internal_vectors


ROOT = Path(__file__).resolve().parents[1]


def tensor_determinants(path: Path) -> list[float]:
    text = path.read_text()
    match = re.search(
        r"internalField\s+nonuniform\s+List<tensor>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if match is None:
        raise ValueError(f"cannot parse tensor field {path}")
    values = []
    for row in re.findall(r"\(([^()]*)\)", match.group(2)):
        a = [float(value) for value in row.split()]
        if len(a) != 9:
            continue
        values.append(
            a[0] * (a[4] * a[8] - a[5] * a[7])
            - a[1] * (a[3] * a[8] - a[5] * a[6])
            + a[2] * (a[3] * a[7] - a[4] * a[6])
        )
    if len(values) != int(match.group(1)):
        raise ValueError(f"tensor count mismatch in {path}")
    return values


def main() -> int:
    report = {}
    for geometry in ("monoventricle", "biventricle"):
        case = ROOT / "cases" / geometry
        metadata = json.loads((case / "meshMetadata.json").read_text())
        expected_cells = metadata["cell_count"]
        times = sorted(
            float(path.name)
            for path in case.iterdir()
            if path.is_dir() and path.name != "0" and re.fullmatch(r"[0-9.eE+-]+", path.name)
        )
        log = (case / "log.loaded.production").read_text()
        vtk_log = (case / "log.foamToVTK.loaded").read_text()
        final = case / "0.001"
        displacement = read_internal_vectors(final / "D")
        pressure = read_internal_scalars(final / "p")
        jacobian = tensor_determinants(final / "F")
        max_displacement = max(
            math.sqrt(sum(component * component for component in value))
            for value in displacement
        )
        passed = all(
            (
                len(times) == 20,
                abs(times[-1] - 0.001) < 1e-14,
                log.count("Nonlinear solve converged") == 20,
                "End\n" in log,
                "FOAM FATAL" not in log,
                len(displacement) == expected_cells,
                len(pressure) == expected_cells,
                len(jacobian) == expected_cells,
                min(jacobian) > 0.0,
                "Number of cells/points in mesh and field do not match" not in vtk_log,
                "FOAM FATAL" not in vtk_log,
                "End:" in vtk_log,
            )
        )
        norms = [float(value) for value in re.findall(r"SNES Function norm ([0-9.eE+-]+)", log)]
        report[geometry] = {
            "passed": passed,
            "cell_count": expected_cells,
            "time_step_count": len(times),
            "final_time": times[-1],
            "nonlinear_convergence_count": log.count("Nonlinear solve converged"),
            "final_snes_norm": norms[-1],
            "maximum_displacement_m": max_displacement,
            "pressure_range_pa": [min(pressure), max(pressure)],
            "deformation_jacobian_range": [min(jacobian), max(jacobian)],
            "vtk_export_passed": "End:" in vtk_log and "FOAM FATAL" not in vtk_log,
        }
        print(f"{geometry}: retained loaded run passed={passed}")
    report["all_passed"] = all(report[name]["passed"] for name in report)
    (ROOT / "reports" / "canonical_loaded_run.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
