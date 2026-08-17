#!/usr/bin/env python3
"""Generate reproducible exterior and true cell-intersection audit drawings."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

from fibres import read_internal_scalars
from mesh_common import read_boundary, read_faces, read_labels, read_points


ROOT = Path(__file__).resolve().parents[1]


def project(point, axes):
    return point[axes[0]], point[axes[1]]


def bounds(polygons):
    values = [value for polygon in polygons for value in polygon]
    xs = [p[0] for p in values]
    ys = [p[1] for p in values]
    return min(xs), max(xs), min(ys), max(ys)


def svg(path, polygons, fills, strokes, title, width=1200, height=900):
    x0, x1, y0, y1 = bounds(polygons)
    padding = 55
    scale = min((width - 2 * padding) / (x1 - x0), (height - 2 * padding) / (y1 - y0))

    def screen(p):
        x = padding + (p[0] - x0) * scale
        y = height - padding - (p[1] - y0) * scale
        return x, y

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="30" y="35" font-family="sans-serif" font-size="24">{title}</text>',
    ]
    for polygon, fill, stroke in zip(polygons, fills, strokes):
        coords = " ".join(f"{x:.2f},{y:.2f}" for x, y in map(screen, polygon))
        lines.append(
            f'<polygon points="{coords}" fill="{fill}" fill-opacity="0.72" '
            f'stroke="{stroke}" stroke-width="0.55" stroke-linejoin="round"/>'
        )
    lines.append(
        f'<text x="30" y="{height - 18}" font-family="sans-serif" font-size="16" '
        f'fill="#333">extent: {x0:.4g}..{x1:.4g} m, {y0:.4g}..{y1:.4g} m</text>'
    )
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n")


def colour(value):
    value = max(0.0, min(1.0, value))
    # Endocardium blue -> epicardium red.
    r = int(45 + 200 * value)
    g = int(105 + 75 * (1.0 - abs(2.0 * value - 1.0)))
    b = int(225 - 175 * value)
    return f"#{r:02x}{g:02x}{b:02x}"


def intersect_edge(a, b, axis, level, tolerance=1e-13):
    da = a[axis] - level
    db = b[axis] - level
    if abs(da) <= tolerance and abs(db) <= tolerance:
        return None
    if da * db > 0.0 or abs(da - db) <= tolerance:
        return None
    fraction = da / (da - db)
    if fraction < -tolerance or fraction > 1.0 + tolerance:
        return None
    return tuple(a[d] + fraction * (b[d] - a[d]) for d in range(3))


def cell_intersections(points, faces, owners, neighbours, axis, level):
    cell_count = max(owners + neighbours) + 1
    cell_edges = [set() for _ in range(cell_count)]
    for face_i, face in enumerate(faces):
        edges = {tuple(sorted(edge)) for edge in zip(face, face[1:] + face[:1])}
        cell_edges[owners[face_i]].update(edges)
        if face_i < len(neighbours):
            cell_edges[neighbours[face_i]].update(edges)
    axes = tuple(value for value in range(3) if value != axis)
    result = []
    cell_ids = []
    for cell_i, edges in enumerate(cell_edges):
        hits = []
        for p0, p1 in edges:
            hit = intersect_edge(points[p0], points[p1], axis, level)
            if hit is not None and not any(
                sum((hit[d] - old[d]) ** 2 for d in range(3)) < 1e-22 for old in hits
            ):
                hits.append(hit)
        if len(hits) < 3:
            continue
        values = [project(hit, axes) for hit in hits]
        centre = (
            sum(p[0] for p in values) / len(values),
            sum(p[1] for p in values) / len(values),
        )
        values.sort(key=lambda p: math.atan2(p[1] - centre[1], p[0] - centre[0]))
        result.append(values)
        cell_ids.append(cell_i)
    return result, cell_ids


def read_set(path):
    if not path.exists():
        return []
    text = path.read_text()
    match = re.search(r"\n\s*\d+\s*\n\s*\((.*?)\)", text, re.S)
    return [int(value) for value in re.findall(r"\d+", match.group(1))] if match else []


def render_case(name, cuts):
    case = ROOT / "cases" / name
    poly = case / "constant" / "polyMesh"
    points = read_points(poly / "points")
    faces = read_faces(poly / "faces")
    owners = read_labels(poly / "owner")
    neighbours = read_labels(poly / "neighbour")
    patches = read_boundary(poly / "boundary")
    out = ROOT / "reports" / "visual" / name
    out.mkdir(parents=True, exist_ok=True)

    patch_faces = {}
    for patch, start, count in patches:
        patch_faces[patch] = list(range(start, start + count))

    exterior_axes = (0, 2) if name == "biventricle" else (0, 1)
    exterior = [project(points[p], exterior_axes) for i in patch_faces["epicardium"] for p in faces[i]]
    exterior = [exterior[i:i + 4] for i in range(0, len(exterior), 4)]
    svg(
        out / "exterior.svg", exterior,
        ["#d9e8f5"] * len(exterior), ["#224b6e"] * len(exterior),
        f"{name}: epicardial surface mesh",
    )

    endo_names = [patch for patch in patch_faces if patch.startswith("endocardium")]
    endo_polygons = []
    endo_colours = []
    for patch_i, patch in enumerate(endo_names):
        for face_i in patch_faces[patch]:
            endo_polygons.append([project(points[p], exterior_axes) for p in faces[face_i]])
            endo_colours.append("#ef8a62" if patch_i == 0 else "#fdd49e")
    svg(
        out / "endocardium.svg", endo_polygons, endo_colours,
        ["#7f0000"] * len(endo_polygons), f"{name}: endocardial surface mesh",
    )

    # This is the ParaView +x-axis view: the planar base projected into y-z.
    # It makes the number of cells crossing the interventricular wall explicit.
    positive_x = [
        [project(points[p], (1, 2)) for p in faces[face_i]]
        for face_i in patch_faces["base"]
    ]
    svg(
        out / "positive_x_base.svg",
        positive_x,
        ["#d9e8f5"] * len(positive_x),
        ["#224b6e"] * len(positive_x),
        f"{name}: +x-axis view of base mesh",
    )

    t_values = read_internal_scalars(case / "0" / "t")
    cut_report = {}
    for label, axis, level in cuts:
        polygons, cell_ids = cell_intersections(points, faces, owners, neighbours, axis, level)
        svg(
            out / f"cut_{label}.svg", polygons,
            [colour(t_values[cell_i]) for cell_i in cell_ids],
            ["#202020"] * len(polygons),
            f"{name}: true cell cut {label}, axis {axis} = {level:g} m",
        )
        cut_report[label] = {"axis": axis, "level_m": level, "cut_cell_count": len(cell_ids)}

    severe_faces = read_set(poly / "sets" / "nonOrthoFaces")
    severe_centres = []
    for face_i in severe_faces:
        face = faces[face_i]
        severe_centres.append([
            sum(points[p][d] for p in face) / len(face) for d in range(3)
        ])
    report = {
        "case": name,
        "views": ["exterior.svg", "endocardium.svg", "positive_x_base.svg"]
        + [f"cut_{c[0]}.svg" for c in cuts],
        "cuts": cut_report,
        "severe_non_orthogonal_face_count": len(severe_faces),
        "severe_non_orthogonal_face_centres_m": severe_centres,
    }
    (out / "visualAudit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    reports = {
        "monoventricle": render_case(
            "monoventricle", (("longitudinal", 2, 0.00031), ("transverse", 0, -0.04))
        ),
        "biventricle": render_case(
            "biventricle", (("longitudinal_septal", 1, 0.00023), ("transverse", 0, 0.04))
        ),
    }
    (ROOT / "reports" / "visual_audit.json").write_text(json.dumps(reports, indent=2) + "\n")
    for name, report in reports.items():
        print(name, report["cuts"])


if __name__ == "__main__":
    main()
