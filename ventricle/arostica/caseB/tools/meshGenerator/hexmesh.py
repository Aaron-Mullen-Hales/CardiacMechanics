#!/usr/bin/env python3
"""Independent all-hexahedral meshes for the Arostica benchmark geometries.

Only the analytical dimensions published in Sections 3.1 and 4.2 are used.
The distributed tetrahedral meshes and their Gmsh generator are not inputs.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import deque
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

MONO = {
    "r_short_endo": 0.025,
    "r_short_epi": 0.035,
    "r_long_endo": 0.090,
    "r_long_epi": 0.097,
    "base_x": 0.025,
}

BIV = {
    "lv_epi": ((0.0, 0.0, 0.0), (0.080, 0.039, 0.039)),
    "rv_epi": ((0.0, 0.0, 0.020), (0.075, 0.038, 0.059)),
    "lv_endo": ((0.0, 0.0, 0.0), (0.069, 0.025, 0.025)),
    "rv_endo": ((0.0, 0.0, 0.020), (0.070, 0.033, 0.054)),
    "base_x": 0.0,
}


def ellipsoid_q(point, definition):
    centre, axes = definition
    return math.sqrt(sum(((point[i] - centre[i]) / axes[i]) ** 2 for i in range(3)))


def mono_inside(point):
    x, y, z = point
    if x >= MONO["base_x"]:
        return False
    q_epi = math.sqrt((x / MONO["r_long_epi"]) ** 2 + (y * y + z * z) / MONO["r_short_epi"] ** 2)
    q_endo = math.sqrt((x / MONO["r_long_endo"]) ** 2 + (y * y + z * z) / MONO["r_short_endo"] ** 2)
    return q_epi < 1.0 and q_endo >= 1.0


def biv_inside(point):
    if point[0] >= BIV["base_x"]:
        return False
    outer = ellipsoid_q(point, BIV["lv_epi"]) < 1.0 or ellipsoid_q(point, BIV["rv_epi"]) < 1.0
    in_lv = ellipsoid_q(point, BIV["lv_endo"]) < 1.0
    in_rv = ellipsoid_q(point, BIV["rv_endo"]) < 1.0
    # The exclusive parts of the intersecting inner ellipsoids are the two
    # cavities. Their overlap is the interventricular septum.
    cavity = in_lv != in_rv
    return outer and not cavity


def largest_face_connected_component(selected):
    remaining = set(selected)
    components = []
    directions = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    while remaining:
        seed = remaining.pop()
        component = {seed}
        queue = deque([seed])
        while queue:
            cell = queue.popleft()
            for direction in directions:
                other = tuple(cell[d] + direction[d] for d in range(3))
                if other in remaining:
                    remaining.remove(other)
                    component.add(other)
                    queue.append(other)
        components.append(component)
    components.sort(key=len, reverse=True)
    return components[0], [len(c) for c in components[1:]]


def classification_confidence(case_name, point):
    """Distance-like cost for changing one voxel's inside/outside state."""
    if case_name == "monoventricle":
        x, y, z = point
        q_endo = math.sqrt((x / MONO["r_long_endo"]) ** 2 + (y * y + z * z) / MONO["r_short_endo"] ** 2)
        q_epi = math.sqrt((x / MONO["r_long_epi"]) ** 2 + (y * y + z * z) / MONO["r_short_epi"] ** 2)
        return min(
            abs(q_endo - 1.0) * MONO["r_short_endo"],
            abs(q_epi - 1.0) * MONO["r_short_epi"],
            abs(x - MONO["base_x"]),
        )
    values = [abs(point[0] - BIV["base_x"])]
    for name in ("lv_epi", "rv_epi", "lv_endo", "rv_endo"):
        values.append(abs(ellipsoid_q(point, BIV[name]) - 1.0) * min(BIV[name][1]))
    return min(values)


def repair_edge_diagonals(case_name, selected, counts, origins, spacing):
    """Remove voxel checkerboards that create non-manifold boundary edges.

    Around a Cartesian grid edge there are four possible incident cells. Two
    diagonally opposite material cells make the union touch only along that
    edge. Flip the least-confident cell classification to obtain a manifold
    one- or three-cell pattern while perturbing the analytical geometry by
    the smallest available amount.
    """
    selected = set(selected)
    flips = []
    axes = ((0, 1, 2), (1, 0, 2), (2, 0, 1))
    while True:
        changed = False
        for along, across_a, across_b in axes:
            for line in range(counts[along]):
                for va in range(1, counts[across_a]):
                    for vb in range(1, counts[across_b]):
                        square = []
                        for da, db in ((-1, -1), (0, -1), (0, 0), (-1, 0)):
                            cell = [0, 0, 0]
                            cell[along] = line
                            cell[across_a] = va + da
                            cell[across_b] = vb + db
                            square.append(tuple(cell))
                        pattern = [cell in selected for cell in square]
                        if pattern not in ([True, False, True, False], [False, True, False, True]):
                            continue
                        def confidence(cell):
                            centre = tuple(origins[d] + spacing * (cell[d] + 0.5) for d in range(3))
                            return classification_confidence(case_name, centre)
                        # Monotone removal cannot oscillate: shave the less
                        # securely classified of the two edge-touching voxels.
                        target = min((cell for cell in square if cell in selected), key=confidence)
                        selected.remove(target)
                        flips.append({"cell": target, "action": "removed", "confidence_m": confidence(target)})
                        changed = True
                        break
                    if changed:
                        break
                if changed:
                    break
            if changed:
                break
        if not changed:
            return selected, flips


def prune_rank_deficient_cells(selected):
    """Remove cells whose face-neighbour directions do not span 3-D.

    OpenFOAM reports these tip or one-voxel ligament cells as having a zero
    cell determinant. Removing them is a conservative erosion of less than
    one grid spacing and makes the finite-volume reconstruction well posed.
    """
    selected = set(selected)
    removed = []
    while True:
        deficient = []
        for cell in selected:
            represented_axes = 0
            for axis in range(3):
                minus = list(cell); minus[axis] -= 1
                plus = list(cell); plus[axis] += 1
                if tuple(minus) in selected or tuple(plus) in selected:
                    represented_axes += 1
            if represented_axes < 3:
                deficient.append(cell)
        if not deficient:
            return selected, removed
        # Simultaneous removal avoids ordering-dependent survival at tips.
        for cell in deficient:
            selected.remove(cell)
        removed.extend(deficient)


def local_face_components(cells):
    remaining = set(cells)
    components = []
    while remaining:
        seed = remaining.pop()
        component = {seed}
        queue = [seed]
        while queue:
            cell = queue.pop()
            for axis in range(3):
                for step in (-1, 1):
                    other = list(cell); other[axis] += step
                    other = tuple(other)
                    if other in remaining:
                        remaining.remove(other)
                        component.add(other)
                        queue.append(other)
        components.append(component)
    return components


def repair_point_contacts(case_name, selected, counts, origins, spacing):
    """Remove vertex-only solid/void contacts in 2x2x2 neighbourhoods."""
    selected = set(selected)
    removed = []

    def confidence(cell):
        centre = tuple(origins[d] + spacing * (cell[d] + 0.5) for d in range(3))
        return classification_confidence(case_name, centre)

    queue = deque(
        (i, j, k)
        for i in range(counts[0] + 1)
        for j in range(counts[1] + 1)
        for k in range(counts[2] + 1)
    )
    queued = set(queue)
    while queue:
        i, j, k = queue.popleft()
        queued.discard((i, j, k))
        local = {
            (i + di, j + dj, k + dk)
            for di in (-1, 0) for dj in (-1, 0) for dk in (-1, 0)
        }
        occupied = local & selected
        empty = local - occupied
        occupied_components = local_face_components(occupied)
        empty_components = local_face_components(empty)
        if len(occupied_components) <= 1 and len(empty_components) <= 1:
            continue
        candidates = []
        if len(occupied_components) > 1:
            largest = max(occupied_components, key=len)
            candidates.extend(occupied - largest)
        if len(empty_components) > 1:
            old_count = len(empty_components)
            empty_bridges = []
            for cell in occupied:
                if len(local_face_components(empty | {cell})) < old_count:
                    empty_bridges.append(cell)
            # Opposite empty cube corners need a two-cell path, so the first
            # removal cannot reduce the component count by itself.
            candidates.extend(empty_bridges if empty_bridges else occupied)
        if not candidates:
            continue
        target = min(set(candidates), key=confidence)
        selected.remove(target)
        removed.append(target)
        # A voxel state affects only its eight corner neighbourhoods.
        for di in (0, 1):
            for dj in (0, 1):
                for dk in (0, 1):
                    vertex = (target[0] + di, target[1] + dj, target[2] + dk)
                    if vertex not in queued:
                        queue.append(vertex)
                        queued.add(vertex)
    return selected, removed


def fmt(value):
    return f"{value:.16g}"


def foam_header(class_name, object_name, location="constant/polyMesh"):
    return (
        "FoamFile\n{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        f"    class {class_name};\n"
        f"    location \"{location}\";\n"
        f"    object {object_name};\n"
        "}\n\n"
    )


def write_list(path, class_name, object_name, values, formatter):
    with path.open("w") as stream:
        stream.write(foam_header(class_name, object_name))
        stream.write(f"{len(values)}\n(\n")
        for value in values:
            stream.write(formatter(value) + "\n")
        stream.write(")\n")


def classify_boundary(case_name, centre, direction, spacing):
    face = tuple(centre[i] + 0.5 * spacing * direction[i] for i in range(3))
    if direction == (1, 0, 0):
        base_x = MONO["base_x"] if case_name == "monoventricle" else BIV["base_x"]
        if abs(face[0] - base_x) < 0.51 * spacing:
            return "base"
    if case_name == "monoventricle":
        x, y, z = face
        q_endo = math.sqrt((x / MONO["r_long_endo"]) ** 2 + (y * y + z * z) / MONO["r_short_endo"] ** 2)
        q_epi = math.sqrt((x / MONO["r_long_epi"]) ** 2 + (y * y + z * z) / MONO["r_short_epi"] ** 2)
        endo_error = abs(q_endo - 1.0) * min(MONO["r_long_endo"], MONO["r_short_endo"])
        epi_error = abs(q_epi - 1.0) * min(MONO["r_long_epi"], MONO["r_short_epi"])
        return "endocardium" if endo_error < epi_error else "epicardium"
    errors = {
        "endocardiumLV": abs(ellipsoid_q(face, BIV["lv_endo"]) - 1.0) * min(BIV["lv_endo"][1]),
        "endocardiumRV": abs(ellipsoid_q(face, BIV["rv_endo"]) - 1.0) * min(BIV["rv_endo"][1]),
        "epicardium": min(
            abs(ellipsoid_q(face, BIV["lv_epi"]) - 1.0) * min(BIV["lv_epi"][1]),
            abs(ellipsoid_q(face, BIV["rv_epi"]) - 1.0) * min(BIV["rv_epi"][1]),
        ),
    }
    return min(errors, key=errors.get)


def surface_error(case_name, patch, point):
    if patch == "base":
        base_x = MONO["base_x"] if case_name == "monoventricle" else BIV["base_x"]
        return abs(point[0] - base_x)
    if case_name == "monoventricle":
        if patch == "endocardium":
            q = math.sqrt((point[0] / MONO["r_long_endo"]) ** 2 + (point[1] ** 2 + point[2] ** 2) / MONO["r_short_endo"] ** 2)
            return abs(q - 1.0) * MONO["r_short_endo"]
        q = math.sqrt((point[0] / MONO["r_long_epi"]) ** 2 + (point[1] ** 2 + point[2] ** 2) / MONO["r_short_epi"] ** 2)
        return abs(q - 1.0) * MONO["r_short_epi"]
    definition = {
        "endocardiumLV": BIV["lv_endo"],
        "endocardiumRV": BIV["rv_endo"],
    }.get(patch)
    if definition:
        return abs(ellipsoid_q(point, definition) - 1.0) * min(definition[1])
    return min(
        abs(ellipsoid_q(point, BIV["lv_epi"]) - 1.0) * min(BIV["lv_epi"][1]),
        abs(ellipsoid_q(point, BIV["rv_epi"]) - 1.0) * min(BIV["rv_epi"][1]),
    )


def build_mesh(case_name, spacing, output):
    if case_name == "monoventricle":
        bounds = ((-0.101, MONO["base_x"]), (-0.042, 0.042), (-0.042, 0.042))
        inside = mono_inside
        patch_order = ("endocardium", "epicardium", "base")
    else:
        bounds = ((-0.084, BIV["base_x"]), (-0.045, 0.045), (-0.045, 0.084))
        inside = biv_inside
        patch_order = ("endocardiumLV", "endocardiumRV", "epicardium", "base")

    counts = tuple(int(math.ceil((hi - lo) / spacing)) for lo, hi in bounds)
    origins = tuple(bounds[d][1] - counts[d] * spacing for d in range(3))
    selected = set()
    for i in range(counts[0]):
        for j in range(counts[1]):
            for k in range(counts[2]):
                centre = tuple(origins[d] + spacing * ((i, j, k)[d] + 0.5) for d in range(3))
                if inside(centre):
                    selected.add((i, j, k))
    topology_repairs = []
    point_repairs = []
    determinant_prunes = []
    while True:
        selected, repairs = repair_edge_diagonals(case_name, selected, counts, origins, spacing)
        selected, corner_repairs = repair_point_contacts(case_name, selected, counts, origins, spacing)
        selected, prunes = prune_rank_deficient_cells(selected)
        topology_repairs.extend(repairs)
        point_repairs.extend(corner_repairs)
        determinant_prunes.extend(prunes)
        if not repairs and not corner_repairs and not prunes:
            break
    selected, discarded_components = largest_face_connected_component(selected)
    cells = sorted(selected)
    cell_id = {cell: i for i, cell in enumerate(cells)}
    centres = [tuple(origins[d] + spacing * (cell[d] + 0.5) for d in range(3)) for cell in cells]

    point_ids = {}
    points = []

    def vertex_id(key):
        if key not in point_ids:
            point_ids[key] = len(points)
            points.append(tuple(origins[d] + spacing * key[d] for d in range(3)))
        return point_ids[key]

    local_vertices = (
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
    )
    # Each face is outward-oriented from its first (and therefore owner) cell.
    local_faces = (
        ((0, 3, 2, 1), (0, 0, -1)),
        ((4, 5, 6, 7), (0, 0, 1)),
        ((0, 1, 5, 4), (0, -1, 0)),
        ((1, 2, 6, 5), (1, 0, 0)),
        ((2, 3, 7, 6), (0, 1, 0)),
        ((3, 0, 4, 7), (-1, 0, 0)),
    )
    face_map = {}
    for cell in cells:
        owner = cell_id[cell]
        vertex_keys = [tuple(cell[d] + offset[d] for d in range(3)) for offset in local_vertices]
        vertices = [vertex_id(key) for key in vertex_keys]
        for local, direction in local_faces:
            oriented = tuple(vertices[index] for index in local)
            key = tuple(sorted(oriented))
            if key in face_map:
                face_map[key]["neighbour"] = owner
            else:
                face_map[key] = {"face": oriented, "owner": owner, "neighbour": None, "direction": direction}

    internal = [record for record in face_map.values() if record["neighbour"] is not None]
    internal.sort(key=lambda record: (record["owner"], record["neighbour"]))
    boundary = {patch: [] for patch in patch_order}
    error_values = {patch: [] for patch in patch_order}
    for record in face_map.values():
        if record["neighbour"] is not None:
            continue
        centre = centres[record["owner"]]
        patch = classify_boundary(case_name, centre, record["direction"], spacing)
        face_centre = tuple(sum(points[index][d] for index in record["face"]) / 4.0 for d in range(3))
        record["faceCentre"] = face_centre
        boundary[patch].append(record)
        error_values[patch].append(surface_error(case_name, patch, face_centre))
    ordered = internal + [record for patch in patch_order for record in boundary[patch]]

    poly = output / "constant" / "polyMesh"
    poly.mkdir(parents=True, exist_ok=True)
    write_list(poly / "points", "vectorField", "points", points, lambda p: f"({fmt(p[0])} {fmt(p[1])} {fmt(p[2])})")
    write_list(poly / "faces", "faceList", "faces", ordered, lambda r: f"4({' '.join(str(x) for x in r['face'])})")
    write_list(poly / "owner", "labelList", "owner", ordered, lambda r: str(r["owner"]))
    write_list(poly / "neighbour", "labelList", "neighbour", internal, lambda r: str(r["neighbour"]))
    with (poly / "boundary").open("w") as stream:
        stream.write(foam_header("polyBoundaryMesh", "boundary"))
        stream.write(f"{len(patch_order)}\n(\n")
        start = len(internal)
        for patch in patch_order:
            stream.write(f"    {patch}\n    {{\n        type wall;\n        nFaces {len(boundary[patch])};\n        startFace {start};\n    }}\n")
            start += len(boundary[patch])
        stream.write(")\n")

    all_boundary_edges = {}
    for patch in patch_order:
        for record in boundary[patch]:
            face = record["face"]
            for a, b in zip(face, face[1:] + face[:1]):
                edge = tuple(sorted((a, b)))
                all_boundary_edges[edge] = all_boundary_edges.get(edge, 0) + 1
    bad_edges = sum(count != 2 for count in all_boundary_edges.values())
    error_report = {}
    for patch, values in error_values.items():
        error_report[patch] = {
            "count": len(values),
            "mean_m": sum(values) / len(values) if values else None,
            "max_m": max(values) if values else None,
        }
    metadata = {
        "case": case_name,
        "method": "independent cell-centred implicit-domain Cartesian voxel mesh",
        "spacing_m": spacing,
        "grid_origin": origins,
        "grid_counts": counts,
        "cell_count": len(cells),
        "point_count": len(points),
        "face_count": len(ordered),
        "internal_face_count": len(internal),
        "patch_face_counts": {patch: len(boundary[patch]) for patch in patch_order},
        "patch_owner_cells": {patch: [record["owner"] for record in boundary[patch]] for patch in patch_order},
        "discarded_component_sizes": discarded_components,
        "topology_repair_flip_count": len(topology_repairs),
        "topology_repairs": topology_repairs,
        "point_contact_repair_count": len(point_repairs),
        "point_contact_repaired_cells": point_repairs,
        "determinant_prune_count": len(determinant_prunes),
        "determinant_pruned_cells": determinant_prunes,
        "non_manifold_boundary_edge_count": bad_edges,
        "surface_face_centre_error": error_report,
        "cell_grid_indices": cells,
        "cell_centres": centres,
        "internal_owner_neighbour": [[r["owner"], r["neighbour"]] for r in internal],
    }
    output.mkdir(parents=True, exist_ok=True)
    data_keys = {
        "cell_grid_indices", "cell_centres", "internal_owner_neighbour",
        "patch_owner_cells", "topology_repairs", "point_contact_repaired_cells",
        "determinant_pruned_cells",
    }
    mesh_data = {key: metadata[key] for key in data_keys}
    mesh_summary = {key: value for key, value in metadata.items() if key not in data_keys}
    (output / "meshData.json").write_text(json.dumps(mesh_data, sort_keys=True, separators=(",", ":")) + "\n")
    (output / "meshMetadata.json").write_text(json.dumps(mesh_summary, indent=2) + "\n")
    (output / "case.foam").touch()
    return metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=("monoventricle", "biventricle", "all"), default="all")
    parser.add_argument("--spacing", type=float, default=0.003)
    args = parser.parse_args()
    names = ("monoventricle", "biventricle") if args.case == "all" else (args.case,)
    reports = {}
    for name in names:
        reports[name] = build_mesh(name, args.spacing, ROOT / "cases" / name)
        print(f"{name}: {reports[name]['cell_count']} all-hex cells")
    (ROOT / "reports").mkdir(exist_ok=True)
    bulky = {
        "cell_grid_indices", "cell_centres", "internal_owner_neighbour",
        "patch_owner_cells", "topology_repairs", "point_contact_repaired_cells",
        "determinant_pruned_cells",
    }
    compact = {name: {key: value for key, value in report.items() if key not in bulky} for name, report in reports.items()}
    (ROOT / "reports" / "mesh_generation.json").write_text(json.dumps(compact, indent=2) + "\n")


if __name__ == "__main__":
    main()
