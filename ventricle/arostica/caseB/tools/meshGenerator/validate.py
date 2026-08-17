#!/usr/bin/env python3
"""Acceptance checks for mesh topology, geometry, and fibre fields."""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

from fibres import cross, dot, load_mesh_data, norm, read_internal_scalars, read_internal_vectors, utility_basis


ROOT = Path(__file__).resolve().parents[1]


def percentile(values, fraction):
    values = sorted(values)
    if not values:
        return None
    position = fraction * (len(values) - 1)
    lower = int(math.floor(position)); upper = int(math.ceil(position))
    if lower == upper:
        return values[lower]
    weight = position - lower
    return values[lower] * (1.0 - weight) + values[upper] * weight


def wrap_degrees(value):
    return (value + 180.0) % 360.0 - 180.0


def vector_diagnostics(fibres, sheets, normals):
    norm_errors = []
    orthogonality = []
    determinants = []
    finite = True
    for f, s, n in zip(fibres, sheets, normals):
        finite = finite and all(math.isfinite(value) for vector in (f, s, n) for value in vector)
        norm_errors.extend(abs(math.sqrt(dot(value, value)) - 1.0) for value in (f, s, n))
        orthogonality.extend(abs(value) for value in (dot(f, s), dot(f, n), dot(s, n)))
        determinants.append(dot(f, cross(s, n)))
    determinant_mean = sum(determinants) / len(determinants)
    return {
        "finite": finite,
        "max_norm_error": max(norm_errors),
        "max_abs_pairwise_dot": max(orthogonality),
        "determinant_min": min(determinants),
        "determinant_max": max(determinants),
        "max_abs_determinant_error": max(abs(abs(value) - 1.0) for value in determinants),
        "handedness": "right" if determinant_mean > 0.0 else "left",
    }


def regex_number(log, pattern, group=1, cast=float):
    match = re.search(pattern, log)
    return cast(match.group(group).rstrip(".,;")) if match else None


def checkmesh_diagnostics(log):
    return {
        "points": regex_number(log, r"points:\s+(\d+)", cast=int),
        "faces": regex_number(log, r"faces:\s+(\d+)", cast=int),
        "cells": regex_number(log, r"cells:\s+(\d+)", cast=int),
        "hexahedra": regex_number(log, r"hexahedra:\s+(\d+)", cast=int),
        "regions": regex_number(log, r"Number of regions:\s+(\d+)", cast=int),
        "aspect_ratio_max": regex_number(log, r"Max aspect ratio = ([^ ]+)"),
        "volume_min": regex_number(log, r"Min volume = ([^ ]+)"),
        "volume_max": regex_number(log, r"Max volume = ([^ ]+)"),
        "volume_total": regex_number(log, r"Total volume = ([^\.\s]+(?:\.[^\.\s]+)?)"),
        "non_orthogonality_max": regex_number(log, r"non-orthogonality Max: ([^ ]+)"),
        "non_orthogonality_average": regex_number(log, r"non-orthogonality Max: [^ ]+ average: ([^\s]+)"),
        "severe_non_orthogonal_faces": regex_number(
            log, r"severely non-orthogonal \(> 70 degrees\) faces: (\d+)", cast=int
        ) or 0,
        "skewness_max": regex_number(log, r"Max skewness = ([^ ]+)"),
        "determinant_min": regex_number(log, r"Cell determinant .*?minimum: ([^ ]+)"),
        "determinant_average": regex_number(log, r"Cell determinant .*?average: ([^\s]+)"),
        "face_weight_min": regex_number(log, r"Face interpolation weight : minimum: ([^ ]+)"),
        "volume_ratio_min": regex_number(log, r"Face volume ratio : minimum: ([^ ]+)"),
        "mesh_ok": "Mesh OK." in log,
        "failed_checks": "Failed " in log,
        "concavity_ok": "Concave cell check OK." in log,
    }


def continuity(values, owner_neighbour):
    angles = []
    for owner, neighbour in owner_neighbour:
        cosine = max(0.0, min(1.0, abs(dot(values[owner], values[neighbour]))))
        angles.append(math.degrees(math.acos(cosine)))
    return {
        "sign_invariant_mean_degrees": sum(angles) / len(angles),
        "sign_invariant_p95_degrees": percentile(angles, 0.95),
        "sign_invariant_max_degrees": max(angles),
        "faces_over_60_degrees": sum(value > 60.0 for value in angles),
        "fraction_over_60_degrees": sum(value > 60.0 for value in angles) / len(angles),
    }


def validate_case(name):
    case = ROOT / "cases" / name
    meta = load_mesh_data(case)
    log = (case / "log.checkMesh").read_text()
    quality = checkmesh_diagnostics(log)
    mesh = {
        "cell_count": meta["cell_count"],
        "checkMesh": quality,
        "non_manifold_boundary_edge_count": meta["non_manifold_boundary_edge_count"],
        "patch_face_counts": meta["patch_face_counts"],
        "surface_vertex_error": meta.get("surface_vertex_error"),
        "patch_reconstruction": meta.get("patch_reconstruction"),
        "layer_alignment": meta.get("layer_alignment"),
    }
    t_values = read_internal_scalars(case / "0" / "t")
    fibres = read_internal_vectors(case / "0" / "f0")
    sheets = read_internal_vectors(case / "0" / "s0")
    normals = read_internal_vectors(case / "0" / "n0")
    face_fibres = read_internal_vectors(case / "0" / "f0f")
    face_sheets = read_internal_vectors(case / "0" / "s0f")
    face_normals = read_internal_vectors(case / "0" / "n0f")
    counts = {
        "t": len(t_values), "f0": len(fibres), "s0": len(sheets), "n0": len(normals),
        "f0f": len(face_fibres), "s0f": len(face_sheets), "n0f": len(face_normals),
    }
    triad = vector_diagnostics(fibres, sheets, normals)
    face_triad = vector_diagnostics(face_fibres, face_sheets, face_normals)
    fibre_continuity = continuity(fibres, meta["internal_owner_neighbour"])
    result = {
        "mesh": mesh,
        "field_counts": counts,
        "transmural": {"min": min(t_values), "mean": sum(t_values) / len(t_values), "max": max(t_values)},
        "triad": triad,
        "internal_face_triad": face_triad,
        "fibre_continuity": fibre_continuity,
    }
    boundary_report_path = ROOT / "reports" / "boundary_audit.json"
    boundary_report = json.loads(boundary_report_path.read_text()) if boundary_report_path.exists() else {}
    result["boundary_orientation"] = boundary_report.get(name)
    if name == "monoventricle":
        angle_errors = []
        for point, t, fibre in zip(meta["cell_centres"], t_values, fibres):
            e_mu, e_theta = utility_basis(point, t)
            measured = math.degrees(math.atan2(dot(fibre, e_mu), dot(fibre, e_theta)))
            expected = -60.0 + 120.0 * t
            angle_errors.append(abs(wrap_degrees(measured - expected)))
        result["helix_angle_error_degrees"] = {
            "mean": sum(angle_errors) / len(angle_errors),
            "p95": percentile(angle_errors, 0.95),
            "max": max(angle_errors),
        }
        solve_log = (case / "log.setFibreField").read_text()
        result["setFibreField"] = {
            "completed": "Done" in solve_log and "FOAM FATAL" not in solve_log,
            "laplace_final_residual": float(re.search(r"Final residual = ([^,]+)", solve_log).group(1)),
            "laplace_iterations": int(re.search(r"No Iterations (\d+)", solve_log).group(1)),
        }

    smoke_log_path = case / "log.smoke.clean"
    smoke_log = smoke_log_path.read_text() if smoke_log_path.exists() else ""
    result["clean_build_smoke"] = {
        "completed": "End" in smoke_log and "FOAM FATAL" not in smoke_log,
        "selected_standard_mixed_model":
            "Selecting solidModel nonLinearGeometryTotalLagrangianTotalDisplacement" in smoke_log,
        "selected_arostica_viscoelastic_law":
            "Selecting mechanical law ArosticaHolzapfelOgdenViscoelastic" in smoke_log,
        "selected_clean_activation_wrapper":
            "Selecting mechanical law electroMechanicalLaw" in smoke_log,
        "solve_pressure": "solvePressure = true" in smoke_log,
        "entered_snes": "Solving the momentum equation for D using PETSc SNES" in smoke_log,
        "snes_converged": "Nonlinear solve converged" in smoke_log,
        "support_conditions_created": smoke_log.count("Creating solidTraction on") == 2,
    }

    failures = []
    if not quality["mesh_ok"] or quality["failed_checks"]:
        failures.append("checkMesh did not finish with an unqualified Mesh OK")
    if quality["hexahedra"] != meta["cell_count"]:
        failures.append("not every cell is reported as a hexahedron")
    if quality["regions"] != 1 or mesh["non_manifold_boundary_edge_count"]:
        failures.append("mesh connectivity/manifold check failed")
    if quality["volume_min"] is None or quality["volume_min"] <= 0.0:
        failures.append("mesh has a non-positive cell volume")
    if quality["determinant_min"] is None or quality["determinant_min"] < 1e-3:
        failures.append("minimum checkMesh cell determinant is below 1e-3")
    if quality["non_orthogonality_max"] is None or quality["non_orthogonality_max"] > 85.0:
        failures.append("maximum non-orthogonality exceeds 85 degrees")
    if quality["skewness_max"] is None or quality["skewness_max"] > 4.0:
        failures.append("maximum skewness exceeds 4")
    if any(counts[name] != meta["cell_count"] for name in ("t", "f0", "s0", "n0")):
        failures.append("volume field length does not match cell count")
    if any(counts[name] != meta["internal_face_count"] for name in ("f0f", "s0f", "n0f")):
        failures.append("face field length does not match internal face count")
    if min(t_values) < -1e-10 or max(t_values) > 1.0 + 1e-10:
        failures.append("transmural field is outside [0,1]")
    if not triad["finite"] or triad["max_norm_error"] > 2e-12:
        failures.append("triad contains non-finite or non-unit vectors")
    if triad["max_abs_pairwise_dot"] > 2e-12 or triad["max_abs_determinant_error"] > 2e-12:
        failures.append("triad is not orthonormal")
    if not face_triad["finite"] or face_triad["max_norm_error"] > 2e-12:
        failures.append("face triad contains non-finite or non-unit vectors")
    if face_triad["max_abs_pairwise_dot"] > 2e-12 or face_triad["max_abs_determinant_error"] > 2e-12:
        failures.append("face triad is not orthonormal")
    if fibre_continuity["sign_invariant_p95_degrees"] > 60.0:
        failures.append("fibre field has excessive cell-to-cell directional jumps")
    if name == "monoventricle":
        if result["helix_angle_error_degrees"]["max"] > 0.05:
            failures.append("monoventricle helix angle law was not reproduced")
        if not result["setFibreField"]["completed"]:
            failures.append("setFibreField did not complete")
    if not all(result["clean_build_smoke"].values()):
        failures.append("clean-build zero-load mixed/viscoelastic smoke test did not pass")
    if not result["boundary_orientation"] or any(
        patch["reversed_face_count"] for patch in result["boundary_orientation"].values()
    ):
        failures.append("one or more physical patch faces has a reversed normal")
    warnings = []
    if quality["severe_non_orthogonal_faces"]:
        warnings.append(
            f"{quality['severe_non_orthogonal_faces']} faces exceed 70 degrees non-orthogonality"
        )
    if fibre_continuity["faces_over_60_degrees"]:
        warnings.append(
            f"{fibre_continuity['faces_over_60_degrees']} internal faces exceed a 60 degree sign-invariant fibre jump"
        )
    result["acceptance"] = {"passed": not failures, "failures": failures, "warnings": warnings}
    return result


def main():
    reports = {name: validate_case(name) for name in ("monoventricle", "biventricle")}
    compatibility_path = ROOT / "reports" / "solids4foam_compatibility.json"
    compatibility = json.loads(compatibility_path.read_text()) if compatibility_path.exists() else {"all_passed": False}
    reports["solids4foam_compatibility"] = compatibility
    loaded_path = ROOT / "reports" / "loaded_smoke.json"
    loaded = json.loads(loaded_path.read_text()) if loaded_path.exists() else {"no_mesh_or_setup_failure": False}
    reports["loaded_smoke"] = loaded
    reports["all_passed"] = (
        all(reports[name]["acceptance"]["passed"] for name in ("monoventricle", "biventricle"))
        and compatibility.get("all_passed", False)
        and loaded.get("no_mesh_or_setup_failure", False)
    )
    output = ROOT / "reports" / "acceptance.json"
    output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(reports, indent=2) + "\n")
    for name in ("monoventricle", "biventricle"):
        report = reports[name]
        print(
            f"{name}: passed={report['acceptance']['passed']} cells={report['mesh']['cell_count']} "
            f"p95_fibre_jump={report['fibre_continuity']['sign_invariant_p95_degrees']:.3f} deg"
        )
        for failure in report["acceptance"]["failures"]:
            print(f"  FAIL: {failure}")
    if not compatibility.get("all_passed", False):
        print("solids4foam compatibility: FAIL")
    if not loaded.get("no_mesh_or_setup_failure", False):
        print("loaded smoke classification: FAIL")
    return 0 if reports["all_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
