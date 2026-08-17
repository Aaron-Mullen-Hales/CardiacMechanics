#!/usr/bin/env python3
"""Validate exact doubling, all-hex quality, and fibre fields at every level."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from fibres import load_mesh_data, read_internal_scalars, read_internal_vectors
from validate import checkmesh_diagnostics, continuity, vector_diagnostics


ROOT = Path(__file__).resolve().parents[1]


def validate_case(case: Path, name: str) -> dict:
    meta = load_mesh_data(case)
    quality = checkmesh_diagnostics((case / "log.checkMesh").read_text())
    scalar = read_internal_scalars(case / "0" / "t")
    fibres = read_internal_vectors(case / "0" / "f0")
    sheets = read_internal_vectors(case / "0" / "s0")
    normals = read_internal_vectors(case / "0" / "n0")
    face_fibres = read_internal_vectors(case / "0" / "f0f")
    face_sheets = read_internal_vectors(case / "0" / "s0f")
    face_normals = read_internal_vectors(case / "0" / "n0f")
    triad = vector_diagnostics(fibres, sheets, normals)
    face_triad = vector_diagnostics(face_fibres, face_sheets, face_normals)
    fibre_continuity = continuity(fibres, meta["internal_owner_neighbour"])
    failures = []
    if not quality["mesh_ok"] or quality["failed_checks"]:
        failures.append("checkMesh did not report an unqualified Mesh OK")
    if quality["cells"] != meta["cell_count"] or quality["hexahedra"] != meta["cell_count"]:
        failures.append("cell/hexahedron count does not match metadata")
    if quality["regions"] != 1 or meta["non_manifold_boundary_edge_count"]:
        failures.append("mesh is not one manifold region")
    if quality["volume_min"] is None or quality["volume_min"] <= 0.0:
        failures.append("non-positive cell volume")
    if quality["determinant_min"] is None or quality["determinant_min"] < 1.0e-3:
        failures.append("minimum cell determinant is below 1e-3")
    if quality["non_orthogonality_max"] is None or quality["non_orthogonality_max"] > 85.0:
        failures.append("maximum non-orthogonality exceeds 85 degrees")
    # The revised mapped charts are deliberately held below 2 at every level.
    if quality["skewness_max"] is None or quality["skewness_max"] >= 2.0:
        failures.append("maximum skewness is not below 2")
    if any(len(values) != meta["cell_count"] for values in (scalar, fibres, sheets, normals)):
        failures.append("cell field length does not match the mesh")
    if any(
        len(values) != meta["internal_face_count"]
        for values in (face_fibres, face_sheets, face_normals)
    ):
        failures.append("face field length does not match the mesh")
    if min(scalar) < -1.0e-10 or max(scalar) > 1.0 + 1.0e-10:
        failures.append("transmural coordinate is outside [0,1]")
    for label, diagnostics in (("cell", triad), ("face", face_triad)):
        if (
            not diagnostics["finite"]
            or diagnostics["max_norm_error"] > 2.0e-12
            or diagnostics["max_abs_pairwise_dot"] > 2.0e-12
            or diagnostics["max_abs_determinant_error"] > 2.0e-12
        ):
            failures.append(f"{label} fibre triad is not finite and orthonormal")
    if fibre_continuity["sign_invariant_p95_degrees"] > 60.0:
        failures.append("p95 neighbour fibre jump exceeds 60 degrees")
    if name == "monoventricle":
        log = (case / "log.setFibreField").read_text()
        if "Done" not in log or "FOAM FATAL" in log:
            failures.append("setFibreFieldArostica did not complete")
    return {
        "cell_count": meta["cell_count"],
        "quality": quality,
        "triad": triad,
        "face_triad": face_triad,
        "fibre_continuity": fibre_continuity,
        "passed": not failures,
        "failures": failures,
    }


def main() -> int:
    report = {}
    previous = {"monoventricle": None, "biventricle": None}
    for level in range(1, 7):
        level_report = {}
        for name in ("monoventricle", "biventricle"):
            result = validate_case(ROOT / "refinements" / f"mesh{level}" / name, name)
            expected = previous[name] * 2 if previous[name] is not None else result["cell_count"]
            if result["cell_count"] != expected:
                result["failures"].append(
                    f"cell count {result['cell_count']} is not expected exact-2x value {expected}"
                )
                result["passed"] = False
            previous[name] = result["cell_count"]
            level_report[name] = result
            q = result["quality"]
            print(
                f"mesh{level} {name}: passed={result['passed']} cells={result['cell_count']} "
                f"skew={q['skewness_max']:.4f} nonOrth={q['non_orthogonality_max']:.3f}"
            )
            for failure in result["failures"]:
                print(f"  FAIL: {failure}")
        report[f"mesh{level}"] = level_report
    report["all_passed"] = all(
        report[f"mesh{level}"][name]["passed"]
        for level in range(1, 7)
        for name in ("monoventricle", "biventricle")
    )
    (ROOT / "reports" / "refinement_quality.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
