#!/usr/bin/env python3
"""Generate and complete fibre triads on the independent hex meshes."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from conforming_hexmesh import BIV, MONO, bivent_piecewise_point, ellipsoid_q


ROOT = Path(__file__).resolve().parents[1]


def load_mesh_data(case):
    values = json.loads((case / "meshMetadata.json").read_text())
    values.update(json.loads((case / "meshData.json").read_text()))
    return values


def dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a):
    magnitude = math.sqrt(dot(a, a))
    if magnitude < 1e-14:
        raise ValueError(f"cannot normalise near-zero vector {a}")
    return tuple(value / magnitude for value in a)


def add_scaled(a, scale, b):
    return tuple(a[i] + scale * b[i] for i in range(3))


def mono_t(point):
    x, y, z = point
    lo, hi = 0.0, 1.0
    for _ in range(60):
        t = 0.5 * (lo + hi)
        rs = MONO["r_short_endo"] + (MONO["r_short_epi"] - MONO["r_short_endo"]) * t
        rl = MONO["r_long_endo"] + (MONO["r_long_epi"] - MONO["r_long_endo"]) * t
        q = math.sqrt((x / rl) ** 2 + (y * y + z * z) / (rs * rs))
        if q > 1.0:
            lo = t
        else:
            hi = t
    return 0.5 * (lo + hi)


def utility_basis(point, t):
    x, y, z = point
    rs = MONO["r_short_endo"] + (MONO["r_short_epi"] - MONO["r_short_endo"]) * t
    rl = MONO["r_long_endo"] + (MONO["r_long_epi"] - MONO["r_long_endo"]) * t
    u = math.atan2(math.sqrt(y * y + z * z) / rs, x / rl)
    v = 0.0 if abs(u) < 1e-7 else math.pi - math.atan2(z, -y)
    e_mu = norm((-rl * math.sin(u), rs * math.cos(u) * math.cos(v), rs * math.cos(u) * math.sin(v)))
    raw_theta = (0.0, -rs * math.sin(u) * math.sin(v), rs * math.sin(u) * math.cos(v))
    e_theta = (0.0, 0.0, 1.0) if math.sqrt(dot(raw_theta, raw_theta)) < 1e-14 else norm(raw_theta)
    return e_mu, e_theta


def mono_expected_triad(point, t, fibre=None):
    e_mu, e_theta = utility_basis(point, t)
    if fibre is None:
        alpha = math.radians(-60.0 + 120.0 * t)
        fibre = norm(tuple(math.sin(alpha) * e_mu[i] + math.cos(alpha) * e_theta[i] for i in range(3)))
    else:
        fibre = norm(fibre)
    normal = cross(e_mu, e_theta)
    # Gram-Schmidt protects the triad from round-off and a cell-centre t that
    # is not exactly the analytical radial coordinate.
    normal = norm(add_scaled(normal, -dot(normal, fibre), fibre))
    sheet = norm(cross(fibre, normal))
    return fibre, sheet, normal


def gradient(point, definition):
    centre, axes = definition
    return norm(tuple((point[i] - centre[i]) / (axes[i] * axes[i]) for i in range(3)))


def distance_like(point, definition):
    return abs(ellipsoid_q(point, definition) - 1.0) * min(definition[1])


def biv_triad(point):
    inner_defs = (BIV["lv_endo"], BIV["rv_endo"])
    outer_defs = (BIV["lv_epi"], BIV["rv_epi"])
    inner = min(inner_defs, key=lambda definition: distance_like(point, definition))
    outer = min(outer_defs, key=lambda definition: distance_like(point, definition))
    d_endo = distance_like(point, inner)
    d_epi = distance_like(point, outer)
    t = d_endo / max(d_endo + d_epi, 1e-14)

    normal = gradient(point, inner)
    # On the overlap (septal) side, tissue is inward from an inner ellipsoid.
    if ellipsoid_q(point, inner) < 1.0:
        normal = tuple(-value for value in normal)
    longitudinal = add_scaled((1.0, 0.0, 0.0), -dot((1.0, 0.0, 0.0), normal), normal)
    if math.sqrt(dot(longitudinal, longitudinal)) < 1e-10:
        fallback = (0.0, 1.0, 0.0)
        longitudinal = add_scaled(fallback, -dot(fallback, normal), normal)
    longitudinal = norm(longitudinal)
    circumferential = norm(cross(normal, longitudinal))
    alpha = math.radians(-60.0 + 120.0 * t)
    fibre = norm(tuple(math.sin(alpha) * longitudinal[i] + math.cos(alpha) * circumferential[i] for i in range(3)))
    sheet = norm(cross(fibre, normal))
    return t, fibre, sheet, normal


def field_header(class_name, name, dimensions="[0 0 0 0 0 0 0]"):
    return (
        "FoamFile\n{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        f"    class {class_name};\n"
        "    location \"0\";\n"
        f"    object {name};\n"
        "}\n"
        f"dimensions {dimensions};\n"
    )


def write_scalar(path, name, values, patches, boundary_types=None, patch_owners=None):
    boundary_types = boundary_types or {}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write(field_header("volScalarField", name))
        stream.write(f"internalField nonuniform List<scalar> {len(values)}\n(\n")
        stream.write("\n".join(f"{value:.16g}" for value in values))
        stream.write("\n);\n\nboundaryField\n{\n")
        for patch in patches:
            kind, value = boundary_types.get(patch, ("calculated", None))
            stream.write(f"    {patch}\n    {{\n        type {kind};\n")
            if value is not None:
                stream.write(f"        value uniform {value};\n")
            elif kind == "calculated" and patch_owners is not None:
                boundary_values = [values[owner] for owner in patch_owners[patch]]
                stream.write(f"        value nonuniform List<scalar> {len(boundary_values)}\n        (\n")
                for boundary_value in boundary_values:
                    stream.write(f"            {boundary_value:.16g}\n")
                stream.write("        );\n")
            elif kind == "calculated":
                stream.write("        value uniform 0;\n")
            stream.write("    }\n")
        stream.write("}\n")


def write_vector(path, name, values, patches, patch_owners=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write(field_header("volVectorField", name))
        stream.write(f"internalField nonuniform List<vector> {len(values)}\n(\n")
        for value in values:
            stream.write(f"({value[0]:.16g} {value[1]:.16g} {value[2]:.16g})\n")
        stream.write(");\n\nboundaryField\n{\n")
        for patch in patches:
            stream.write(f"    {patch}\n    {{\n        type calculated;\n")
            if patch_owners is None:
                stream.write("        value uniform (0 0 0);\n")
            else:
                boundary_values = [values[owner] for owner in patch_owners[patch]]
                stream.write(f"        value nonuniform List<vector> {len(boundary_values)}\n        (\n")
                for value in boundary_values:
                    stream.write(f"            ({value[0]:.16g} {value[1]:.16g} {value[2]:.16g})\n")
                stream.write("        );\n")
            stream.write("    }\n")
        stream.write("}\n")


def write_surface_vector(path, name, internal_values, patch_values):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write(field_header("surfaceVectorField", name))
        stream.write(f"internalField nonuniform List<vector> {len(internal_values)}\n(\n")
        for value in internal_values:
            stream.write(f"({value[0]:.16g} {value[1]:.16g} {value[2]:.16g})\n")
        stream.write(");\n\nboundaryField\n{\n")
        for patch, values in patch_values.items():
            stream.write(f"    {patch}\n    {{\n        type calculated;\n")
            stream.write(f"        value nonuniform List<vector> {len(values)}\n        (\n")
            for value in values:
                stream.write(f"            ({value[0]:.16g} {value[1]:.16g} {value[2]:.16g})\n")
            stream.write("        );\n    }\n")
        stream.write("}\n")


def right_handed_triad(fibre, sheet, normal):
    fibre = norm(fibre)
    sheet = norm(add_scaled(sheet, -dot(sheet, fibre), fibre))
    result_normal = norm(cross(fibre, sheet))
    if dot(result_normal, normal) < 0.0:
        sheet = tuple(-value for value in sheet)
        result_normal = tuple(-value for value in result_normal)
    return fibre, sheet, result_normal


def mapped_biv_shell_triad(chamber, parameter, transmural):
    """Rule-based triad on one mapped chamber column."""
    disk_y, disk_z = parameter
    inner = bivent_piecewise_point(chamber, "endo", disk_y, disk_z)
    outer = bivent_piecewise_point(chamber, "outer", disk_y, disk_z)
    sheet_normal = norm(tuple(outer[d] - inner[d] for d in range(3)))
    longitudinal = add_scaled((1.0, 0.0, 0.0), -dot((1.0, 0.0, 0.0), sheet_normal), sheet_normal)
    if math.sqrt(dot(longitudinal, longitudinal)) < 1.0e-10:
        longitudinal = add_scaled((0.0, 1.0, 0.0), -dot((0.0, 1.0, 0.0), sheet_normal), sheet_normal)
    longitudinal = norm(longitudinal)
    circumferential = norm(cross(sheet_normal, longitudinal))
    alpha = math.radians(-60.0 + 120.0 * transmural)
    fibre = norm(tuple(
        math.sin(alpha) * longitudinal[d] + math.cos(alpha) * circumferential[d]
        for d in range(3)
    ))
    sheet = norm(cross(sheet_normal, fibre))
    return right_handed_triad(fibre, sheet, sheet_normal)


def mapped_biv_triads(meta):
    """Triads for both chamber walls and the interpolating septal block."""
    triads = []
    scalar_t = []
    for chamber, parameter, fraction in zip(
        meta["cell_chamber"], meta["cell_parametric_centre"], meta["cell_transmural_fraction"]
    ):
        if chamber in ("LV", "RV"):
            triads.append(mapped_biv_shell_triad(chamber, parameter, fraction))
            scalar_t.append(fraction)
            continue
        lv = mapped_biv_shell_triad("LV", parameter, 1.0)
        rv = mapped_biv_shell_triad("RV", (parameter[0], -parameter[1]), 1.0)
        rv_f = rv[0] if dot(lv[0], rv[0]) >= 0.0 else tuple(-x for x in rv[0])
        rv_s = rv[1] if dot(lv[1], rv[1]) >= 0.0 else tuple(-x for x in rv[1])
        fibre = tuple((1.0 - fraction) * lv[0][d] + fraction * rv_f[d] for d in range(3))
        sheet = tuple((1.0 - fraction) * lv[1][d] + fraction * rv_s[d] for d in range(3))
        normal_hint = tuple((1.0 - fraction) * lv[2][d] + fraction * rv[2][d] for d in range(3))
        if math.sqrt(dot(normal_hint, normal_hint)) < 1.0e-10:
            normal_hint = cross(fibre, sheet)
        triads.append(right_handed_triad(fibre, sheet, normal_hint))
        scalar_t.append(1.0)
    return triads, scalar_t


def build_face_triads(meta, fibres, sheets, normals, handedness="left"):
    internal = []
    for owner, neighbour in meta["internal_owner_neighbour"]:
        f_owner, f_neighbour = fibres[owner], fibres[neighbour]
        if dot(f_owner, f_neighbour) < 0.0:
            f_neighbour = tuple(-value for value in f_neighbour)
        face_f = norm(tuple(f_owner[i] + f_neighbour[i] for i in range(3)))
        if handedness == "right":
            s_owner, s_neighbour = sheets[owner], sheets[neighbour]
            if dot(s_owner, s_neighbour) < 0.0:
                s_neighbour = tuple(-value for value in s_neighbour)
            face_f, face_s, face_n = right_handed_triad(
                face_f,
                tuple(s_owner[i] + s_neighbour[i] for i in range(3)),
                tuple(normals[owner][i] + normals[neighbour][i] for i in range(3)),
            )
        else:
            n_owner, n_neighbour = normals[owner], normals[neighbour]
            if dot(n_owner, n_neighbour) < 0.0:
                n_neighbour = tuple(-value for value in n_neighbour)
            face_n = tuple(n_owner[i] + n_neighbour[i] for i in range(3))
            face_n = norm(add_scaled(face_n, -dot(face_n, face_f), face_f))
            face_s = norm(cross(face_f, face_n))
        internal.append((face_f, face_s, face_n))
    boundary = {}
    for patch, owners in meta["patch_owner_cells"].items():
        boundary[patch] = [(fibres[owner], sheets[owner], normals[owner]) for owner in owners]
    return internal, boundary


def write_face_triad(case, meta, fibres, sheets, normals, handedness="left"):
    internal, boundary = build_face_triads(meta, fibres, sheets, normals, handedness)
    for index, name in ((0, "f0f"), (1, "s0f"), (2, "n0f")):
        write_surface_vector(
            case / "0" / name,
            name,
            [triad[index] for triad in internal],
            {patch: [triad[index] for triad in triads] for patch, triads in boundary.items()},
        )


def read_internal_vectors(path):
    text = path.read_text()
    match = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"cannot parse internal vector list in {path}")
    values = [tuple(float(x) for x in item.split()) for item in re.findall(r"\(([^()]*)\)", match.group(2))]
    if len(values) != int(match.group(1)):
        raise ValueError(f"vector count mismatch in {path}")
    return values


def read_internal_scalars(path):
    text = path.read_text()
    match = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"cannot parse internal scalar list in {path}")
    values = [float(x) for x in match.group(2).split()]
    if len(values) != int(match.group(1)):
        raise ValueError(f"scalar count mismatch in {path}")
    return values


def write_solids4foam_fields(case, patches):
    basal_spring = 1e6 if case.name == "biventricle" else 1e5
    with (case / "0" / "D").open("w") as stream:
        stream.write(field_header("volVectorField", "D", "[0 1 0 0 0 0 0]"))
        stream.write("internalField uniform (0 0 0);\n\nboundaryField\n{\n")
        for patch in patches:
            if patch == "base":
                stream.write(
                    f"    {patch}\n    {{\n"
                    "        type arosticaVectorSpringDashpotTraction;\n"
                    f"        springCoefficient springCoefficientValue [1 -2 -2 0 0 0 0] {basal_spring:.16g};\n"
                    "        dashpotCoefficient dashpotCoefficientValue [1 -2 -1 0 0 0 0] 5e3;\n"
                    "        useUndeformedArea true;\n        writeDiagnostics false;\n"
                    "        value uniform (0 0 0);\n    }\n"
                )
            elif patch == "epicardium":
                stream.write(
                    f"    {patch}\n    {{\n"
                    "        type arosticaNormalSpringDashpotTraction;\n"
                    "        springCoefficient springCoefficientValue [1 -2 -2 0 0 0 0] 1e8;\n"
                    "        dashpotCoefficient dashpotCoefficientValue [1 -2 -1 0 0 0 0] 5e3;\n"
                    "        useUndeformedArea true;\n        writeDiagnostics false;\n"
                    "        tangentialTraction uniform (0 0 0);\n"
                    "        value uniform (0 0 0);\n    }\n"
                )
            else:
                if case.name == "monoventricle":
                    pressure_file = "pressure_caseB.dat"
                elif patch == "endocardiumLV":
                    pressure_file = "pressure_LV.dat"
                else:
                    pressure_file = "pressure_RV.dat"
                stream.write(
                    f"    {patch}\n    {{\n        type solidTraction;\n"
                    "        traction uniform (0 0 0);\n"
                    f"        pressureSeries {{ file \"$FOAM_CASE/constant/loadCurves/{pressure_file}\"; outOfBounds clamp; }}\n"
                    "        useUndeformedArea false;\n"
                    "        value uniform (0 0 0);\n    }\n"
                )
        stream.write("}\n")
    with (case / "0" / "pointD").open("w") as stream:
        stream.write(field_header("pointVectorField", "pointD", "[0 1 0 0 0 0 0]"))
        stream.write("internalField uniform (0 0 0);\n\nboundaryField\n{\n")
        for patch in patches:
            stream.write(f"    {patch} {{ type calculated; value uniform (0 0 0); }}\n")
        stream.write("}\n")
    with (case / "0" / "p").open("w") as stream:
        stream.write(
            "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            "    class volScalarField;\n    location \"0\";\n    object p;\n}\n"
            "dimensions [1 -1 -2 0 0 0 0];\ninternalField uniform 0;\n\nboundaryField\n{\n"
        )
        for patch in patches:
            stream.write(f"    {patch} {{ type zeroGradient; }}\n")
        stream.write("}\n")


def prepare(cases_root=None):
    cases_root = Path(cases_root) if cases_root is not None else ROOT / "cases"
    mono_case = cases_root / "monoventricle"
    mono_meta = load_mesh_data(mono_case)
    mono_values = [mono_t(point) for point in mono_meta["cell_centres"]]
    write_scalar(
        mono_case / "0" / "t", "t", mono_values,
        ("endocardium", "epicardium", "base"),
        {
            "endocardium": ("fixedValue", "0"),
            "epicardium": ("fixedValue", "1"),
            "base": ("zeroGradient", None),
        },
    )
    write_solids4foam_fields(mono_case, ("endocardium", "epicardium", "base"))

    biv_case = cases_root / "biventricle"
    biv_meta = load_mesh_data(biv_case)
    triads, geometry_values = mapped_biv_triads(biv_meta)
    patches = ("endocardiumLV", "endocardiumRV", "epicardium", "base")
    write_solids4foam_fields(biv_case, patches)
    write_scalar(biv_case / "0" / "t", "t", geometry_values, patches, patch_owners=biv_meta["patch_owner_cells"])
    for index, name in ((0, "f0"), (1, "s0"), (2, "n0")):
        write_vector(biv_case / "0" / name, name, [value[index] for value in triads], patches, biv_meta["patch_owner_cells"])
    write_face_triad(
        biv_case,
        biv_meta,
        [value[0] for value in triads],
        [value[1] for value in triads],
        [value[2] for value in triads],
        handedness="right",
    )
    (biv_case / "fibreMethod.json").write_text(json.dumps({
        "method": "mapped LV/RV rule-based fibre charts with an explicitly interpolated septal chart",
        "source": "new structured chamber and septal coordinates in meshData.json",
        "angles_degrees": {"endocardial": -60.0, "epicardial": 60.0},
        "interpolation": "column-local helix law followed by right-handed Gram-Schmidt orthonormalisation",
        "face_interpolation": "sign-aligned owner/neighbour cell triads followed by right-handed Gram-Schmidt orthonormalisation",
        "septum": "sign-aligned interpolation between the LV and RV epicardial-side chart values",
        "transmural_scalar_note": "exact logical shell fraction; septal cells use t=1 while vector triads are interpolated across the septal block",
    }, indent=2) + "\n")


def complete_mono(cases_root=None):
    cases_root = Path(cases_root) if cases_root is not None else ROOT / "cases"
    case = cases_root / "monoventricle"
    meta = load_mesh_data(case)
    t_values = read_internal_scalars(case / "0" / "t")
    # The utility's Laplace field supplies the exact mesh-dependent
    # transmural coordinate.  Reconstruct the requested helix angle from that
    # scalar instead of retaining its rounded intermediate vector output.
    triads = [mono_expected_triad(point, t) for point, t in zip(meta["cell_centres"], t_values)]
    patches = ("endocardium", "epicardium", "base")
    # Rewrite all three volume fields together.  The utility computes useful
    # boundary values for f0 independently, but the clean Aróstica law checks
    # each boundary triad as well as the internal field.  Owner-cell values
    # keep every boundary tuple exactly consistent with the completed triad.
    write_vector(case / "0" / "f0", "f0", [value[0] for value in triads], patches, meta["patch_owner_cells"])
    write_vector(case / "0" / "s0", "s0", [value[1] for value in triads], patches, meta["patch_owner_cells"])
    write_vector(case / "0" / "n0", "n0", [value[2] for value in triads], patches, meta["patch_owner_cells"])
    write_face_triad(
        case, meta,
        [value[0] for value in triads],
        [value[1] for value in triads],
        [value[2] for value in triads],
    )
    (case / "fibreMethod.json").write_text(json.dumps({
        "method": "setFibreFieldArostica Arostica2025 FVM Laplace solve plus orthonormal triad completion",
        "angles_degrees": {"endocardial": -60.0, "epicardial": 60.0},
        "transmural_range": [min(t_values), max(t_values)],
    }, indent=2) + "\n")
    # The utility's scalar diagnostic `a` aliases solids4foam's tensor `A` on
    # case-insensitive filesystems. None of these diagnostics is required by
    # the material law, so keep the production case namespace unambiguous.
    for auxiliary in ("a", "alphaRadians", "b", "q", "rl", "rs", "uu", "vv"):
        path = case / "0" / auxiliary
        if path.exists():
            path.unlink()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("prepare", "complete-mono"))
    parser.add_argument("--cases-root", type=Path, default=ROOT / "cases")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.cases_root)
    else:
        complete_mono(args.cases_root)


if __name__ == "__main__":
    main()
