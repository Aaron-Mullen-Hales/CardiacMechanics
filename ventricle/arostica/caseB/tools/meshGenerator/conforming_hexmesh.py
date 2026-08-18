#!/usr/bin/env python3
"""Mapped multi-block all-hex meshes for the Aróstica benchmark.

Both meshes are built from quadrilateral surface charts and explicit
transmural columns.  No Cartesian cell selection and no tetrahedron-to-hex
subdivision is used.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from mesh_common import HEX_FACES, write_hex_poly_mesh


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

# Smooth, non-intersecting two-cavity approximation used by the mapped bivent
# topology.  It retains the published chamber long axes and z bounds; the
# common outer x/y scales differ from the extrema by at most 3.2%.  Piecewise
# facing-side radii replace the analytical XOR intersection with a finite
# septum and valid endocardium-to-epicardium columns.
BIV_MAPPED = {
    "lv": {
        "endo_centre": (0.0, 0.0, 0.0),
        "outer_centre": (0.0, 0.0, 0.007),
        "endo_axes": (0.069, 0.025),
        "endo_z_free": 0.025,
        "endo_z_facing": 0.006,
        "outer_axes": (0.0775, 0.0385),
        "outer_facing_axis_delta": (-0.002, -0.002),
        "outer_z_free": 0.046,
        "outer_z_facing": 0.005,
    },
    "rv": {
        "endo_centre": (0.0, 0.0, 0.020),
        "outer_centre": (0.0, 0.0, 0.016),
        "endo_axes": (0.070, 0.033),
        "endo_z_free": 0.055,
        "endo_z_facing": 0.002,
        "outer_axes": (0.0775, 0.0385),
        "outer_facing_axis_delta": (-0.007, -0.005),
        "outer_z_free": 0.063,
        "outer_z_facing": 0.002,
    },
    "base_x": 0.0,
}


def ellipsoid_q(point, definition):
    centre, axes = definition
    return math.sqrt(sum(((point[i] - centre[i]) / axes[i]) ** 2 for i in range(3)))


def mono_point(disk_y, disk_z, transmural):
    radial_parameter = math.hypot(disk_y, disk_z)
    a = MONO["r_long_endo"] + transmural * (MONO["r_long_epi"] - MONO["r_long_endo"])
    b = MONO["r_short_endo"] + transmural * (MONO["r_short_epi"] - MONO["r_short_endo"])
    x = -a + (MONO["base_x"] + a) * radial_parameter * radial_parameter
    cross_radius = b * math.sqrt(max(0.0, 1.0 - (x / a) ** 2))
    if radial_parameter < 1.0e-14:
        return (-a, 0.0, 0.0)
    return (
        x,
        cross_radius * disk_y / radial_parameter,
        cross_radius * disk_z / radial_parameter,
    )


def cap_quad_mesh(divisions, longitudinal_divisions=20, cap_fraction=0.08):
    """Central cap block plus a regular annular quad grid.

    The former large square core made four chevron seams visible over most of
    the ventricle.  Here the extraordinary topology is confined to the distal
    ``cap_fraction`` of the long axis; the remainder consists of smooth,
    near-longitudinal rings.
    """
    if divisions < 4 or divisions % 2:
        raise ValueError("surface divisions must be an even integer >= 4")
    if longitudinal_divisions < 4:
        raise ValueError("longitudinal divisions must be >= 4")
    inner_half_width = math.sqrt(0.5 * cap_fraction)
    points = []
    point_ids = {}

    def point_id(y, z):
        key = (round(y, 13), round(z, 13))
        if key not in point_ids:
            point_ids[key] = len(points)
            points.append((y, z))
        return point_ids[key]

    quads = []
    outer_quads = set()
    # Small central square.  Keeping its boundary slightly non-circular avoids
    # the three-collinear-corner quads produced when every square-boundary
    # node is forced onto one ellipsoidal ring.
    core = {}
    for j in range(divisions + 1):
        square_z = -1.0 + 2.0 * j / divisions
        for i in range(divisions + 1):
            square_y = -1.0 + 2.0 * i / divisions
            y, z = inner_half_width * square_y, inner_half_width * square_z
            core[(i, j)] = point_id(y, z)
    for j in range(divisions):
        for i in range(divisions):
            quads.append((core[(i, j)], core[(i + 1, j)], core[(i + 1, j + 1)], core[(i, j + 1)]))

    # Four annular blocks.  Their inner edge is the circular image of a core
    # edge; constant-angle rays continue to the base.  Squared radius is
    # distributed uniformly, which is close to uniform long-axis spacing for
    # the ellipsoidal map x = x_apex + L*r^2.
    blocks = (
        ((1.0, -1.0), (1.0, 1.0)),
        ((1.0, 1.0), (-1.0, 1.0)),
        ((-1.0, 1.0), (-1.0, -1.0)),
        ((-1.0, -1.0), (1.0, -1.0)),
    )
    for inner0, inner1 in blocks:
        grid = {}
        for j in range(longitudinal_divisions + 1):
            axial = j / longitudinal_divisions
            radius = math.sqrt(cap_fraction + (1.0 - cap_fraction) * axial)
            for i in range(divisions + 1):
                q = i / divisions
                square = tuple((1.0 - q) * inner0[d] + q * inner1[d] for d in range(2))
                square_radius = math.hypot(*square)
                direction = (square[0] / square_radius, square[1] / square_radius)
                inner_radius = inner_half_width * square_radius
                radius = math.sqrt(inner_radius * inner_radius + (1.0 - inner_radius * inner_radius) * axial)
                y, z = radius * direction[0], radius * direction[1]
                grid[(i, j)] = point_id(y, z)
        for j in range(longitudinal_divisions):
            for i in range(divisions):
                quad = (grid[(i, j)], grid[(i + 1, j)], grid[(i + 1, j + 1)], grid[(i, j + 1)])
                quad_i = len(quads)
                quads.append(quad)
                if j == longitudinal_divisions - 1:
                    outer_quads.add(quad_i)
    return points, quads, outer_quads, longitudinal_divisions


def build_monoventricle(output: Path, surface_divisions: int, longitudinal_divisions: int, transmural_layers: int):
    if surface_divisions < 6:
        raise ValueError("surface-divisions must be >= 6")
    if transmural_layers < 3:
        raise ValueError("at least three transmural layers are required")

    cap_points, cap_quads, outer_quads, radial_divisions = cap_quad_mesh(
        surface_divisions, longitudinal_divisions
    )
    points = []
    point_id = {}
    for k in range(transmural_layers + 1):
        t = k / transmural_layers
        for cap_i, (disk_y, disk_z) in enumerate(cap_points):
            point_id[(cap_i, k)] = len(points)
            points.append(mono_point(disk_y, disk_z, t))

    cells = []
    boundary_keys = {}
    cell_logical_indices = []
    for k in range(transmural_layers):
        for quad_i, quad in enumerate(cap_quads):
            cell = tuple(point_id[(cap_i, layer)] for layer in (k, k + 1) for cap_i in quad)
            cells.append(cell)
            cell_logical_indices.append([quad_i, k])
            if k == 0:
                boundary_keys[tuple(sorted(cell[n] for n in HEX_FACES[0]))] = "endocardium"
            if k == transmural_layers - 1:
                boundary_keys[tuple(sorted(cell[n] for n in HEX_FACES[1]))] = "epicardium"
            if quad_i in outer_quads:
                boundary_keys[tuple(sorted(cell[n] for n in HEX_FACES[4]))] = "base"

    connectivity = write_hex_poly_mesh(
        output, points, cells, ("endocardium", "epicardium", "base"), boundary_keys
    )
    surface_errors = {}
    for patch, layer in (("endocardium", 0), ("epicardium", transmural_layers)):
        a = MONO["r_long_endo"] if layer == 0 else MONO["r_long_epi"]
        b = MONO["r_short_endo"] if layer == 0 else MONO["r_short_epi"]
        errors = []
        for cap_i in range(len(cap_points)):
            x, y, z = points[point_id[(cap_i, layer)]]
            errors.append(abs(math.sqrt((x / a) ** 2 + (y * y + z * z) / (b * b)) - 1.0))
        surface_errors[patch] = {"max_implicit_error": max(errors), "mean_implicit_error": sum(errors) / len(errors)}
    base_errors = [abs(points[index][0] - MONO["base_x"]) for key, patch in boundary_keys.items() if patch == "base" for index in key]
    surface_errors["base"] = {"max_plane_error_m": max(base_errors), "mean_plane_error_m": sum(base_errors) / len(base_errors)}

    summary = {
        "case": "monoventricle",
        "method": "mapped five-block ellipsoidal shell; compact concentric cap and regular longitudinal rings",
        "geometry_source": "Aróstica paper nested truncated ellipsoids",
        "surface_divisions_per_square_axis": surface_divisions,
        "longitudinal_divisions": radial_divisions,
        "surface_quads_per_layer": len(cap_quads),
        "transmural_layers": transmural_layers,
        "surface_vertex_error": surface_errors,
        **{key: value for key, value in connectivity.items() if key not in {
            "cell_centres", "internal_owner_neighbour", "patch_owner_cells"
        }},
    }
    data = {
        "cell_centres": connectivity["cell_centres"],
        "internal_owner_neighbour": connectivity["internal_owner_neighbour"],
        "patch_owner_cells": connectivity["patch_owner_cells"],
        "cell_logical_indices": cell_logical_indices,
    }
    return summary, data


def ellipsoid_disk_point(disk_y, disk_z, definition):
    """Map a unit quad disk to the positive-x half of an ellipsoid."""
    centre, axes = definition
    radius = math.hypot(disk_y, disk_z)
    x = centre[0] + axes[0] * (1.0 - radius * radius)
    if radius < 1.0e-14:
        return x, centre[1], centre[2]
    cross_radius = math.sqrt(max(0.0, 1.0 - ((x - centre[0]) / axes[0]) ** 2))
    return (
        x,
        centre[1] + axes[1] * cross_radius * disk_y / radius,
        centre[2] + axes[2] * cross_radius * disk_z / radius,
    )


def bivent_piecewise_point(chamber, surface, disk_y, disk_z):
    """Nested piecewise ellipsoid for a free wall and septal-side wall.

    Each surface keeps common long/cross-fibre axes on both sides.  The z
    semi-axis changes at the equator and the outer centre may be mildly
    shifted, while the chart remains continuous and columns do not cross.
    """
    definition = BIV_MAPPED[chamber.lower()]
    facing = disk_z >= 0.0 if chamber == "LV" else disk_z <= 0.0
    axes_xy = definition[f"{surface}_axes"]
    if surface == "outer" and facing:
        delta = definition["outer_facing_axis_delta"]
        weight = abs(disk_z)
        axes_xy = tuple(axes_xy[d] + weight * delta[d] for d in range(2))
    z_axis = definition[f"{surface}_z_{'facing' if facing else 'free'}"]
    return ellipsoid_disk_point(
        disk_y, disk_z,
        (definition[f"{surface}_centre"], (axes_xy[0], axes_xy[1], z_axis)),
    )


def layer_transition_cross_section(free_layers, facing_layers, variant=8, smoothing_iterations=20):
    """All-quad two-row transition from a fine to a coarse wall stack.

    A polygon with half the requested boundary divisions is triangulated and
    each triangle is split into three quads.  The split doubles every boundary
    edge, giving exactly ``free_layers`` edges at b=0, ``facing_layers`` at
    b=1, and two surface-chart rows at t=0 and t=1.
    """
    if free_layers < 4 or facing_layers < 2 or free_layers % 2 or facing_layers % 2:
        raise ValueError("bivent wall-transition layer counts must be even (free >= 4, facing >= 2)")
    left_macro = free_layers // 2
    right_macro = facing_layers // 2
    polygon = [(0.0, 0.0), (1.0, 0.0)]
    polygon.extend((1.0, j / right_macro) for j in range(1, right_macro + 1))
    polygon.extend((0.0, j / left_macro) for j in range(left_macro, 0, -1))
    if free_layers == 6 and facing_layers == 2:
        variants = (
            ((0, 1, 2), (0, 2, 4), (0, 4, 5), (2, 3, 4)),
            ((0, 1, 2), (0, 2, 5), (2, 3, 4), (2, 4, 5)),
            ((0, 1, 2), (0, 2, 5), (2, 3, 5), (3, 4, 5)),
            ((0, 1, 4), (0, 4, 5), (1, 2, 3), (1, 3, 4)),
            ((0, 1, 4), (0, 4, 5), (1, 2, 4), (2, 3, 4)),
            ((0, 1, 5), (1, 2, 3), (1, 3, 4), (1, 4, 5)),
            ((0, 1, 5), (1, 2, 3), (1, 3, 5), (3, 4, 5)),
            ((0, 1, 5), (1, 2, 4), (1, 4, 5), (2, 3, 4)),
            ((0, 1, 5), (1, 2, 5), (2, 3, 4), (2, 4, 5)),
            ((0, 1, 5), (1, 2, 5), (2, 3, 5), (3, 4, 5)),
        )
        if not 0 <= variant < len(variants):
            raise ValueError(f"transition variant must be in [0, {len(variants) - 1}]")
        triangles = variants[variant]
    else:
        triangles = [(1, k, k + 1) for k in range(2, len(polygon) - 1)]
        triangles.append((1, len(polygon) - 1, 0))
    points = []
    point_ids = {}

    def point_id(value):
        key = tuple(round(component, 13) for component in value)
        if key not in point_ids:
            point_ids[key] = len(points)
            points.append(value)
        return point_ids[key]

    quads = []
    for triangle in triangles:
        a, b, c = (polygon[index] for index in triangle)
        signed_area = sum(
            value[0] * following[1] - following[0] * value[1]
            for value, following in zip((a, b, c), (b, c, a))
        )
        if signed_area < 0.0:
            b, c = c, b
        ia, ib, ic = (point_id(value) for value in (a, b, c))
        iab = point_id(tuple((a[d] + b[d]) / 2.0 for d in range(2)))
        ibc = point_id(tuple((b[d] + c[d]) / 2.0 for d in range(2)))
        ica = point_id(tuple((c[d] + a[d]) / 2.0 for d in range(2)))
        centre = point_id(tuple((a[d] + b[d] + c[d]) / 3.0 for d in range(2)))
        quads.extend(((ia, iab, centre, ica), (ib, ibc, centre, iab), (ic, ica, centre, ibc)))

    # Relax only interior reference points. Boundary points encode the exact
    # six/two layer and two-row interface counts and must remain fixed.
    edge_counts = {}
    neighbours = {index: set() for index in range(len(points))}
    for quad in quads:
        for a, b in zip(quad, quad[1:] + quad[:1]):
            edge = tuple(sorted((a, b)))
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
            neighbours[a].add(b)
            neighbours[b].add(a)
    boundary_points = {
        index for edge, count in edge_counts.items() if count == 1 for index in edge
    }
    for _ in range(smoothing_iterations):
        old = list(points)
        points = [
            value if index in boundary_points else tuple(
                sum(old[other][d] for other in neighbours[index]) / len(neighbours[index])
                for d in range(2)
            )
            for index, value in enumerate(old)
        ]
    return points, quads


def build_biventricle(
    output: Path,
    surface_divisions: int,
    longitudinal_divisions: int,
    apex_cap_fraction: float,
    transmural_layers: int,
    facing_layers: int,
    septal_layers: int,
    transition_variant: int,
    transition_smoothing: int,
    transition_physical_smoothing: int,
):
    """Build mapped chamber free walls and a resolution-balanced septum."""
    if transmural_layers < 2:
        raise ValueError("bivent free-wall layers must be >= 2")
    if facing_layers < 2:
        raise ValueError("bivent facing-wall layers must be >= 2")
    if facing_layers > transmural_layers:
        raise ValueError("facing-wall layers cannot exceed free-wall layers")
    if facing_layers < transmural_layers and (
        transmural_layers < 4 or transmural_layers % 2 or facing_layers % 2
    ):
        raise ValueError("transitioned bivent layer counts must be even (free >= 4, facing >= 2)")
    if septal_layers < 1:
        raise ValueError("septal connector layers must be >= 1")
    if not 0.05 <= apex_cap_fraction <= 0.35:
        raise ValueError("bivent apex-cap fraction must be in [0.05, 0.35]")
    cap_points, cap_quads, outer_quads, _ = cap_quad_mesh(
        surface_divisions, longitudinal_divisions, cap_fraction=apex_cap_fraction
    )
    required = ("endocardiumLV", "endocardiumRV", "epicardium", "base")
    use_layer_transition = facing_layers < transmural_layers
    points = []
    point_ids = {}

    def global_point_id(chamber, cap_i, fraction):
        fraction = round(fraction, 13)
        key = (chamber, cap_i, fraction)
        if key in point_ids:
            return point_ids[key]
        disk_y, disk_z = cap_points[cap_i]
        inner = bivent_piecewise_point(chamber, "endo", disk_y, disk_z)
        outer = bivent_piecewise_point(chamber, "outer", disk_y, disk_z)
        point = tuple((1.0 - fraction) * inner[d] + fraction * outer[d] for d in range(3))
        point_ids[key] = len(points)
        points.append(point)
        return point_ids[key]

    cells = []
    cell_chambers = []
    cell_transmural_layers = []
    cell_parametric_centres = []
    cell_transmural_fractions = []
    boundary_keys = {}

    # The two chart rows touching disk_z=0 carry an all-quad 6-to-2 layer
    # transition.  Everything farther into the free wall keeps six layers;
    # everything farther into the chamber-facing half uses two.
    seam_records = {}
    transition_quads = set()
    for quad_i, quad in enumerate(cap_quads):
        if not use_layer_transition:
            break
        z_values = [cap_points[index][1] for index in quad]
        zero_edges = []
        for edge_i in range(4):
            if (
                abs(z_values[edge_i]) < 1.0e-12
                and abs(z_values[(edge_i + 1) % 4]) < 1.0e-12
            ):
                zero_edges.append(edge_i)
        if not zero_edges or not (min(z_values) < -1.0e-12 or max(z_values) > 1.0e-12):
            continue
        if len(zero_edges) != 1:
            raise ValueError(f"unexpected transition chart topology in quad {quad_i}")
        edge_i = zero_edges[0]
        s0, s1 = quad[edge_i], quad[(edge_i + 1) % 4]
        side = "negative" if sum(z_values) < 0.0 else "positive"
        record = {
            "quad_i": quad_i,
            "off": {s0: quad[(edge_i - 1) % 4], s1: quad[(edge_i + 2) % 4]},
        }
        seam_records.setdefault(tuple(sorted((s0, s1))), {})[side] = record
        transition_quads.add(quad_i)

    if use_layer_transition and (
        not seam_records
        or any(set(record) != {"negative", "positive"} for record in seam_records.values())
    ):
        raise ValueError("incomplete bivent free-to-facing transition strip")

    for chamber in ("LV", "RV"):
        for quad_i, quad in enumerate(cap_quads):
            if quad_i in transition_quads:
                continue
            mean_z = sum(cap_points[index][1] for index in quad) / 4.0
            free_wall = mean_z < 0.0 if chamber == "LV" else mean_z > 0.0
            layer_count = transmural_layers if free_wall or not use_layer_transition else facing_layers
            for layer in range(layer_count):
                cell = tuple(
                    global_point_id(chamber, cap_i, shell_layer / layer_count)
                    for shell_layer in (layer, layer + 1)
                    for cap_i in quad
                )
                cells.append(cell)
                cell_chambers.append(chamber)
                cell_transmural_layers.append(layer)
                cell_parametric_centres.append([
                    sum(cap_points[index][d] for index in quad) / 4.0 for d in range(2)
                ])
                cell_transmural_fractions.append((layer + 0.5) / layer_count)
                if layer == 0:
                    patch = "endocardiumLV" if chamber == "LV" else "endocardiumRV"
                    boundary_keys[tuple(sorted(cell[n] for n in HEX_FACES[0]))] = patch
                if layer == layer_count - 1 and free_wall:
                    boundary_keys[tuple(sorted(cell[n] for n in HEX_FACES[1]))] = "epicardium"
                if quad_i in outer_quads:
                    boundary_keys[tuple(sorted(cell[n] for n in HEX_FACES[4]))] = "base"

    # Replace the skipped chart rows with conformal all-hex transition prisms.
    if use_layer_transition:
        transition_points, transition_section_quads = layer_transition_cross_section(
            transmural_layers, facing_layers, transition_variant, transition_smoothing
        )
    else:
        transition_points, transition_section_quads = [], []
    transition_ids = {}
    transition_cell_start = len(cells)

    def transition_point_id(chamber, seam_i, free_i, facing_i, b_value, t_value):
        b_value, t_value = round(b_value, 13), round(t_value, 13)
        if abs(b_value) < 1.0e-12:
            return global_point_id(chamber, free_i, t_value)
        if abs(b_value - 1.0) < 1.0e-12:
            return global_point_id(chamber, facing_i, t_value)
        if abs(b_value - 0.5) < 1.0e-12 and (abs(t_value) < 1.0e-12 or abs(t_value - 1.0) < 1.0e-12):
            return global_point_id(chamber, seam_i, t_value)
        key = (chamber, seam_i, free_i, facing_i, b_value, t_value)
        if key in transition_ids:
            return transition_ids[key]
        if b_value < 0.5:
            left_i, right_i, weight = free_i, seam_i, 2.0 * b_value
        else:
            left_i, right_i, weight = seam_i, facing_i, 2.0 * (b_value - 0.5)
        inner_ends = [bivent_piecewise_point(chamber, "endo", *cap_points[index]) for index in (left_i, right_i)]
        outer_ends = [bivent_piecewise_point(chamber, "outer", *cap_points[index]) for index in (left_i, right_i)]
        inner = tuple((1.0 - weight) * inner_ends[0][d] + weight * inner_ends[1][d] for d in range(3))
        outer = tuple((1.0 - weight) * outer_ends[0][d] + weight * outer_ends[1][d] for d in range(3))
        transition_ids[key] = len(points)
        points.append(tuple((1.0 - t_value) * inner[d] + t_value * outer[d] for d in range(3)))
        return transition_ids[key]

    edge_to_side_face = {
        (0, 1): HEX_FACES[2],
        (1, 2): HEX_FACES[3],
        (2, 3): HEX_FACES[4],
        (0, 3): HEX_FACES[5],
    }
    for chamber in ("LV", "RV"):
        free_side = "negative" if chamber == "LV" else "positive"
        facing_side = "positive" if chamber == "LV" else "negative"
        endocardium_patch = "endocardiumLV" if chamber == "LV" else "endocardiumRV"
        for seam_edge, records in seam_records.items():
            seam_nodes = tuple(sorted(seam_edge, key=lambda index: cap_points[index][0]))
            free_off = records[free_side]["off"]
            facing_off = records[facing_side]["off"]
            for section_quad in transition_section_quads:
                cell = tuple(
                    transition_point_id(
                        chamber,
                        seam_i,
                        free_off[seam_i],
                        facing_off[seam_i],
                        transition_points[section_i][0],
                        transition_points[section_i][1],
                    )
                    for seam_i in seam_nodes
                    for section_i in section_quad
                )
                cells.append(cell)
                cell_chambers.append(chamber)
                cell_transmural_layers.append(-1)
                parameters = []
                for seam_i in seam_nodes:
                    for section_i in section_quad:
                        b_value = transition_points[section_i][0]
                        if b_value < 0.5:
                            left, right, weight = cap_points[free_off[seam_i]], cap_points[seam_i], 2.0 * b_value
                        else:
                            left, right, weight = cap_points[seam_i], cap_points[facing_off[seam_i]], 2.0 * (b_value - 0.5)
                        parameters.append(tuple((1.0 - weight) * left[d] + weight * right[d] for d in range(2)))
                cell_parametric_centres.append([
                    sum(parameter[d] for parameter in parameters) / len(parameters) for d in range(2)
                ])
                cell_transmural_fractions.append(
                    sum(transition_points[index][1] for index in section_quad) / 4.0
                )
                if math.isclose(sum(value * value for value in cap_points[seam_nodes[0]]), 1.0, abs_tol=1.0e-10):
                    boundary_keys[tuple(sorted(cell[n] for n in HEX_FACES[0]))] = "base"
                if math.isclose(sum(value * value for value in cap_points[seam_nodes[1]]), 1.0, abs_tol=1.0e-10):
                    boundary_keys[tuple(sorted(cell[n] for n in HEX_FACES[1]))] = "base"
                for edge in ((0, 1), (1, 2), (2, 3), (0, 3)):
                    values = [transition_points[section_quad[position]] for position in edge]
                    local_face = edge_to_side_face[tuple(sorted(edge))]
                    if all(abs(value[1]) < 1.0e-12 for value in values):
                        boundary_keys[tuple(sorted(cell[n] for n in local_face))] = endocardium_patch
                    elif all(abs(value[1] - 1.0) < 1.0e-12 for value in values):
                        if sum(value[0] for value in values) / 2.0 < 0.5:
                            boundary_keys[tuple(sorted(cell[n] for n in local_face))] = "epicardium"

    # Smooth only the newly introduced interior transition points in physical
    # space. Interface, endocardial, epicardial, and basal vertices are all
    # global chart points and remain fixed, so geometry and conformity are
    # preserved while warped interior transition faces are relaxed.
    transition_edges = (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
    )
    movable_transition_points = set(transition_ids.values())
    transition_neighbours = {index: set() for index in movable_transition_points}
    for cell in cells[transition_cell_start:]:
        for a, b in transition_edges:
            pa, pb = cell[a], cell[b]
            if pa in transition_neighbours:
                transition_neighbours[pa].add(pb)
            if pb in transition_neighbours:
                transition_neighbours[pb].add(pa)
    physical_transition_smoothing = transition_physical_smoothing
    for _ in range(physical_transition_smoothing):
        old = list(points)
        for index, neighbours in transition_neighbours.items():
            average = tuple(sum(old[other][d] for other in neighbours) / len(neighbours) for d in range(3))
            points[index] = tuple(0.5 * old[index][d] + 0.5 * average[d] for d in range(3))

    # Fill the septum with its own mapped block.  It connects the upper half
    # of the LV outer chart to the reflected lower half of the RV chart.  The
    # two chamber interfaces become internal faces; only the disk diameter is
    # epicardial and the outer semicircle is basal.
    coordinate_to_cap = {
        (round(value[0], 13), round(value[1], 13)): index
        for index, value in enumerate(cap_points)
    }

    def mirrored_cap_id(cap_i):
        disk_y, disk_z = cap_points[cap_i]
        key = (round(disk_y, 13), round(-disk_z, 13))
        if key in coordinate_to_cap:
            return coordinate_to_cap[key]
        # Independent annular blocks can differ at their shared reflected
        # edge by a few final floating-point bits at some division counts.
        # Resolve that harmless representation difference geometrically.
        nearest = min(
            range(len(cap_points)),
            key=lambda index: (cap_points[index][0] - disk_y) ** 2
            + (cap_points[index][1] + disk_z) ** 2,
        )
        error = math.hypot(cap_points[nearest][0] - disk_y, cap_points[nearest][1] + disk_z)
        if error > 1.0e-10:
            raise ValueError(f"cap chart is not reflection symmetric at point {cap_i}: error {error}")
        return nearest

    # Two connector layers are sufficient at the checked coarse level.  The
    # earlier four-layer connector visibly over-resolved the tissue between
    # the two already layered chamber walls.
    septal_ids = {}

    def septal_point_id(cap_i, layer):
        key = (cap_i, layer)
        if key in septal_ids:
            return septal_ids[key]
        if layer == 0:
            value = global_point_id("LV", cap_i, 1.0)
        elif layer == septal_layers:
            value = global_point_id("RV", mirrored_cap_id(cap_i), 1.0)
        else:
            left = points[global_point_id("LV", cap_i, 1.0)]
            right = points[global_point_id("RV", mirrored_cap_id(cap_i), 1.0)]
            fraction = layer / septal_layers
            value = len(points)
            points.append(tuple((1.0 - fraction) * left[d] + fraction * right[d] for d in range(3)))
        septal_ids[key] = value
        return value

    edge_to_side_face = {
        (0, 1): HEX_FACES[2],
        (1, 2): HEX_FACES[3],
        (2, 3): HEX_FACES[4],
        (0, 3): HEX_FACES[5],
    }
    for quad_i, quad in enumerate(cap_quads):
        mean_z = sum(cap_points[index][1] for index in quad) / 4.0
        if mean_z <= 1.0e-12:
            continue
        zero_edges = []
        for edge in ((0, 1), (1, 2), (2, 3), (0, 3)):
            if all(abs(cap_points[quad[position]][1]) < 1.0e-12 for position in edge):
                zero_edges.append(edge)
        for layer in range(septal_layers):
            cell = tuple(
                septal_point_id(cap_i, slab_layer)
                for slab_layer in (layer, layer + 1)
                for cap_i in quad
            )
            cells.append(cell)
            cell_chambers.append("septum")
            cell_transmural_layers.append(layer)
            cell_parametric_centres.append([
                sum(cap_points[index][d] for index in quad) / 4.0 for d in range(2)
            ])
            cell_transmural_fractions.append((layer + 0.5) / septal_layers)
            if quad_i in outer_quads:
                boundary_keys[tuple(sorted(cell[n] for n in HEX_FACES[4]))] = "base"
            for edge in zero_edges:
                local_face = edge_to_side_face[tuple(sorted(edge))]
                boundary_keys[tuple(sorted(cell[n] for n in local_face))] = "epicardium"

    connectivity = write_hex_poly_mesh(output, points, cells, required, boundary_keys)
    summary = {
        "case": "biventricle",
        "method": "mapped nested LV and RV shell blocks connected by a separately layered conformal septal block",
        "geometry_source": "Aróstica analytical dimensions with documented non-intersecting two-cavity approximation",
        "geometry_parameters": BIV_MAPPED,
        "geometry_approximation": (
            "retains benchmark chamber long axes, z extent, free-wall scales, and distinct LV/RV endocardial "
            "patches; common outer x/y scales are within 3.2 percent of the published extrema; replaces the "
            "intersecting XOR chart by non-intersecting piecewise ellipsoidal cavities and an explicit "
            "finite-thickness septal block"
        ),
        "surface_divisions_per_cap_block": surface_divisions,
        "longitudinal_divisions": longitudinal_divisions,
        "apex_cap_fraction": apex_cap_fraction,
        "surface_quads_per_chamber": len(cap_quads),
        "transmural_layers": transmural_layers,
        "facing_transmural_layers": facing_layers,
        "septal_layers": septal_layers,
        "interventricular_wall_layers": 2 * facing_layers + septal_layers,
        "layer_transition_variant": transition_variant,
        "layer_transition_smoothing_iterations": transition_smoothing,
        "layer_transition_physical_smoothing_iterations": physical_transition_smoothing,
        "layer_alignment": (
            f"{transmural_layers} free-wall layers transition conformally to "
            f"{facing_layers} chamber-facing layers on each side of the "
            f"{septal_layers}-layer septal connector"
            if use_layer_transition
            else f"uniform {transmural_layers}-layer chamber walls with a {septal_layers}-layer septal connector"
        ),
        **{key: value for key, value in connectivity.items() if key not in {
            "cell_centres", "internal_owner_neighbour", "patch_owner_cells"
        }},
    }
    data = {
        "cell_centres": connectivity["cell_centres"],
        "internal_owner_neighbour": connectivity["internal_owner_neighbour"],
        "patch_owner_cells": connectivity["patch_owner_cells"],
        "cell_chamber": cell_chambers,
        "cell_transmural_layer": cell_transmural_layers,
        "cell_parametric_centre": cell_parametric_centres,
        "cell_transmural_fraction": cell_transmural_fractions,
    }
    return summary, data


def write_metadata(output: Path, summary, data):
    output.mkdir(parents=True, exist_ok=True)
    (output / "meshMetadata.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "meshData.json").write_text(json.dumps(data, separators=(",", ":")) + "\n")
    (output / "case.foam").touch()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=("monoventricle", "biventricle", "all"), default="all")
    parser.add_argument("--surface-divisions", type=int, default=12)
    parser.add_argument("--longitudinal-divisions", type=int, default=20)
    parser.add_argument("--bivent-apex-cap-fraction", type=float, default=0.187)
    parser.add_argument("--transmural-layers", type=int, default=6)
    parser.add_argument("--facing-layers", type=int, default=2)
    parser.add_argument("--septal-layers", type=int, default=2)
    parser.add_argument("--transition-variant", type=int, default=1)
    parser.add_argument("--transition-smoothing", type=int, default=20)
    parser.add_argument("--transition-physical-smoothing", type=int, default=3)
    args = parser.parse_args()
    names = ("monoventricle", "biventricle") if args.case == "all" else (args.case,)
    reports = {}
    for name in names:
        output = ROOT / "cases" / name
        if name == "monoventricle":
            summary, data = build_monoventricle(
                output, args.surface_divisions, args.longitudinal_divisions, args.transmural_layers
            )
        else:
            summary, data = build_biventricle(
                output,
                args.surface_divisions,
                args.longitudinal_divisions,
                args.bivent_apex_cap_fraction,
                args.transmural_layers,
                args.facing_layers,
                args.septal_layers,
                args.transition_variant,
                args.transition_smoothing,
                args.transition_physical_smoothing,
            )
        write_metadata(output, summary, data)
        reports[name] = summary
        print(f"{name}: {summary['cell_count']} geometry-conforming hexahedra")
    (ROOT / "reports").mkdir(exist_ok=True)
    (ROOT / "reports" / "mesh_generation.json").write_text(json.dumps(reports, indent=2) + "\n")


if __name__ == "__main__":
    main()
