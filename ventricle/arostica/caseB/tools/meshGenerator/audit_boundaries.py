#!/usr/bin/env python3
"""Check patch face orientation against the analytical benchmark surfaces."""

from __future__ import annotations

import json
import math
from pathlib import Path

from conforming_hexmesh import BIV, BIV_MAPPED, MONO, ellipsoid_q
from mesh_common import face_normal, mean_point, read_boundary, read_faces, read_points


ROOT = Path(__file__).resolve().parents[1]


def unit(value):
    magnitude = math.sqrt(sum(x * x for x in value))
    return tuple(x / magnitude for x in value)


def ellipsoid_gradient(point, definition):
    centre, axes = definition
    return tuple((point[d] - centre[d]) / (axes[d] * axes[d]) for d in range(3))


def expected_normal(case, patch, point):
    if case == "monoventricle":
        if patch == "base":
            return (1.0, 0.0, 0.0)
        definition = (
            ((0.0, 0.0, 0.0), (MONO["r_long_endo"], MONO["r_short_endo"], MONO["r_short_endo"]))
            if patch == "endocardium" else
            ((0.0, 0.0, 0.0), (MONO["r_long_epi"], MONO["r_short_epi"], MONO["r_short_epi"]))
        )
        value = ellipsoid_gradient(point, definition)
        return tuple(-x for x in value) if patch == "endocardium" else value
    if patch == "base":
        return (-1.0, 0.0, 0.0)
    if patch == "endocardiumLV":
        definition = BIV_MAPPED["lv"]
        centre = definition["endo_centre"]
        z_axis = definition["endo_z_free"] if point[2] < centre[2] else definition["endo_z_facing"]
        axes = (*definition["endo_axes"], z_axis)
        return tuple(-x for x in ellipsoid_gradient(point, (centre, axes)))
    if patch == "endocardiumRV":
        definition = BIV_MAPPED["rv"]
        centre = definition["endo_centre"]
        z_axis = definition["endo_z_free"] if point[2] > centre[2] else definition["endo_z_facing"]
        axes = (*definition["endo_axes"], z_axis)
        return tuple(-x for x in ellipsoid_gradient(point, (centre, axes)))
    # Free-wall epicardium follows the mapped LV/RV ellipsoids.  The narrow
    # connecting strip follows their common disk-diameter curve; its outward
    # direction is captured by the same global radial expression.
    centre = (0.0, 0.0, 0.010)
    return (
        point[0] / (0.0775 * 0.0775),
        point[1] / (0.0385 * 0.0385),
        (point[2] - centre[2]) / (0.069 * 0.069),
    )


def audit(case_name):
    poly = ROOT / "cases" / case_name / "constant" / "polyMesh"
    points = read_points(poly / "points")
    faces = read_faces(poly / "faces")
    result = {}
    for patch, start, count in read_boundary(poly / "boundary"):
        cosines = []
        for face_i in range(start, start + count):
            face = faces[face_i]
            centre = mean_point(points, face)
            actual = unit(face_normal(points, face))
            expected = unit(expected_normal(case_name, patch, centre))
            cosines.append(sum(actual[d] * expected[d] for d in range(3)))
        result[patch] = {
            "face_count": count,
            "minimum_expected_normal_cosine": min(cosines),
            "mean_expected_normal_cosine": sum(cosines) / len(cosines),
            "reversed_face_count": sum(value <= 0.0 for value in cosines),
        }
    return result


def main():
    report = {case: audit(case) for case in ("monoventricle", "biventricle")}
    report["all_faces_outward"] = all(
        patch["reversed_face_count"] == 0
        for case in ("monoventricle", "biventricle")
        for patch in report[case].values()
    )
    (ROOT / "reports" / "boundary_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["all_faces_outward"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
