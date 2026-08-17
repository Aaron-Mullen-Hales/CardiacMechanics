#!/usr/bin/env python3
"""Small dependency-free OpenFOAM polyMesh helpers."""

from __future__ import annotations

import math
import re
import shutil
from pathlib import Path


HEX_FACES = (
    (0, 3, 2, 1),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)


def fmt(value: float) -> str:
    return f"{value:.16g}"


def foam_header(class_name: str, object_name: str, location: str = "constant/polyMesh") -> str:
    return (
        "FoamFile\n{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        f"    class {class_name};\n"
        f"    location \"{location}\";\n"
        f"    object {object_name};\n"
        "}\n\n"
    )


def write_list(path: Path, class_name: str, object_name: str, values, formatter) -> None:
    with path.open("w") as stream:
        stream.write(foam_header(class_name, object_name))
        stream.write(f"{len(values)}\n(\n")
        for value in values:
            stream.write(formatter(value) + "\n")
        stream.write(")\n")


def mean_point(points, ids):
    return tuple(sum(points[index][d] for index in ids) / len(ids) for d in range(3))


def sub(a, b):
    return tuple(a[d] - b[d] for d in range(3))


def dot(a, b):
    return sum(a[d] * b[d] for d in range(3))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def face_normal(points, face):
    """Area-weighted polygon normal for an ordered face."""
    origin = points[face[0]]
    result = (0.0, 0.0, 0.0)
    for i in range(1, len(face) - 1):
        value = cross(sub(points[face[i]], origin), sub(points[face[i + 1]], origin))
        result = tuple(result[d] + value[d] for d in range(3))
    return result


def write_hex_poly_mesh(output: Path, points, cells, patch_order, boundary_patch_by_key):
    """Write an all-hex polyMesh and return connectivity metadata.

    ``boundary_patch_by_key`` maps the sorted point IDs of every boundary
    quadrilateral to its physical patch. Faces are oriented out of the owner
    cell using geometry, so the supplied cell vertex order need only follow
    the topological HEX_FACES convention.
    """
    cell_centres = [mean_point(points, cell) for cell in cells]
    face_map = {}
    for cell_i, cell in enumerate(cells):
        for local in HEX_FACES:
            face = tuple(cell[index] for index in local)
            key = tuple(sorted(face))
            if key in face_map:
                face_map[key]["neighbour"] = cell_i
                continue
            fc = mean_point(points, face)
            normal = face_normal(points, face)
            if dot(normal, sub(fc, cell_centres[cell_i])) < 0.0:
                face = tuple(reversed(face))
            face_map[key] = {"face": face, "owner": cell_i, "neighbour": None}

    internal = [record for record in face_map.values() if record["neighbour"] is not None]
    internal.sort(key=lambda record: (record["owner"], record["neighbour"]))
    boundary = {patch: [] for patch in patch_order}
    for key, record in face_map.items():
        if record["neighbour"] is not None:
            continue
        if key not in boundary_patch_by_key:
            raise ValueError(
                f"unclassified boundary quadrilateral {key} "
                f"(owner cell {record['owner']}, ordered face {record['face']})"
            )
        boundary[boundary_patch_by_key[key]].append(record)
    for records in boundary.values():
        records.sort(key=lambda record: (record["owner"], record["face"]))
    ordered = internal + [record for patch in patch_order for record in boundary[patch]]

    poly = output / "constant" / "polyMesh"
    poly.mkdir(parents=True, exist_ok=True)
    if (poly / "sets").exists():
        shutil.rmtree(poly / "sets")
    write_list(
        poly / "points", "vectorField", "points", points,
        lambda p: f"({fmt(p[0])} {fmt(p[1])} {fmt(p[2])})",
    )
    write_list(
        poly / "faces", "faceList", "faces", ordered,
        lambda r: f"4({' '.join(str(x) for x in r['face'])})",
    )
    write_list(poly / "owner", "labelList", "owner", ordered, lambda r: str(r["owner"]))
    write_list(poly / "neighbour", "labelList", "neighbour", internal, lambda r: str(r["neighbour"]))
    with (poly / "boundary").open("w") as stream:
        stream.write(foam_header("polyBoundaryMesh", "boundary"))
        stream.write(f"{len(patch_order)}\n(\n")
        start = len(internal)
        for patch in patch_order:
            stream.write(
                f"    {patch}\n    {{\n        type wall;\n"
                f"        nFaces {len(boundary[patch])};\n        startFace {start};\n    }}\n"
            )
            start += len(boundary[patch])
        stream.write(")\n")

    edge_counts = {}
    for records in boundary.values():
        for record in records:
            face = record["face"]
            for a, b in zip(face, face[1:] + face[:1]):
                edge = tuple(sorted((a, b)))
                edge_counts[edge] = edge_counts.get(edge, 0) + 1

    return {
        "cell_count": len(cells),
        "point_count": len(points),
        "face_count": len(ordered),
        "internal_face_count": len(internal),
        "cell_centres": cell_centres,
        "internal_owner_neighbour": [[r["owner"], r["neighbour"]] for r in internal],
        "patch_face_counts": {patch: len(boundary[patch]) for patch in patch_order},
        "patch_owner_cells": {patch: [r["owner"] for r in boundary[patch]] for patch in patch_order},
        "non_manifold_boundary_edge_count": sum(count != 2 for count in edge_counts.values()),
    }


def _foam_list_body(path: Path):
    text = path.read_text()
    match = re.search(r"\n\s*(\d+)\s*\n\s*\(\s*\n(.*?)\n\s*\)\s*$", text, re.S)
    if not match:
        raise ValueError(f"cannot parse OpenFOAM list {path}")
    return int(match.group(1)), match.group(2)


def read_points(path: Path):
    count, body = _foam_list_body(path)
    values = [tuple(float(x) for x in match.groups()) for match in re.finditer(
        r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)", body
    )]
    if len(values) != count:
        raise ValueError(f"point count mismatch in {path}: {len(values)} != {count}")
    return values


def read_faces(path: Path):
    count, body = _foam_list_body(path)
    values = []
    for match in re.finditer(r"(\d+)\s*\(([^)]*)\)", body):
        n = int(match.group(1))
        face = tuple(int(x) for x in match.group(2).split())
        if len(face) != n:
            raise ValueError(f"face size mismatch in {path}")
        values.append(face)
    if len(values) != count:
        raise ValueError(f"face count mismatch in {path}: {len(values)} != {count}")
    return values


def read_labels(path: Path):
    count, body = _foam_list_body(path)
    values = [int(x) for x in re.findall(r"[-+]?\d+", body)]
    if len(values) != count:
        raise ValueError(f"label count mismatch in {path}: {len(values)} != {count}")
    return values


def read_boundary(path: Path):
    text = path.read_text()
    result = []
    for match in re.finditer(
        r"\n\s*([A-Za-z][A-Za-z0-9_]*)\s*\n\s*\{(.*?)\}", text, re.S
    ):
        name, body = match.groups()
        n_faces = re.search(r"\bnFaces\s+(\d+)\s*;", body)
        start = re.search(r"\bstartFace\s+(\d+)\s*;", body)
        if n_faces and start:
            result.append((name, int(start.group(1)), int(n_faces.group(1))))
    return result


def read_tet_poly_mesh(poly: Path):
    points = read_points(poly / "points")
    faces = read_faces(poly / "faces")
    owners = read_labels(poly / "owner")
    neighbours = read_labels(poly / "neighbour")
    cell_count = max(owners + neighbours) + 1
    cell_vertices = [set() for _ in range(cell_count)]
    for face_i, face in enumerate(faces):
        cell_vertices[owners[face_i]].update(face)
        if face_i < len(neighbours):
            cell_vertices[neighbours[face_i]].update(face)
    cells = [tuple(sorted(vertices)) for vertices in cell_vertices]
    if any(len(cell) != 4 for cell in cells):
        raise ValueError("source polyMesh is not purely tetrahedral")
    patch_by_face = {}
    patch_order = []
    for patch, start, count in read_boundary(poly / "boundary"):
        patch_order.append(patch)
        for face_i in range(start, start + count):
            patch_by_face[face_i] = patch
    return points, faces, owners, neighbours, cells, patch_order, patch_by_face
