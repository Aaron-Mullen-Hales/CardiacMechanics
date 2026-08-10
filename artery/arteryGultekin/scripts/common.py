#!/usr/bin/env python3
"""Small, dependency-free OpenFOAM mesh and field helpers."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

def dot(a, b):
    return sum(x*y for x,y in zip(a,b))

def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def mesh_hash(case: Path) -> str:
    digest = hashlib.sha256()
    for name in ("points", "faces", "owner", "neighbour", "boundary"):
        path = case / "constant/polyMesh" / name
        digest.update(name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()

def strip_comments(text: str) -> str:
    text = re.sub(r"//[^\n]*", "", text)
    return re.sub(r"/\*.*?\*/", "", text, flags=re.S)

def _list_body(text: str, name: str) -> str:
    match = re.search(rf"\bobject\s+{name}\s*;.*?\n\s*\d+\s*\(", text, re.S)
    if not match:
        raise ValueError(f"Cannot find {name} list")
    opening = text.find("(", match.start())
    depth = 1
    end = opening + 1
    while end < len(text) and depth:
        depth += (text[end] == "(") - (text[end] == ")")
        end += 1
    return text[opening + 1:end - 1]

def _tuples(body: str):
    return [tuple(float(x) for x in item.split()) for item in re.findall(r"\(([^()]*)\)", body)]

def _ints(body: str):
    return [int(x) for x in re.findall(r"[-+]?\d+", body)]

def load_mesh(case: Path) -> dict:
    root = case / "constant/polyMesh"
    points = _tuples(strip_comments((root / "points").read_text()))
    face_body = _list_body(strip_comments((root / "faces").read_text()), "faces")
    faces = [tuple(int(x) for x in item.split()) for item in re.findall(r"\d+\s*\(([^()]*)\)", face_body)]
    owner = _ints(_list_body(strip_comments((root / "owner").read_text()), "owner"))
    neighbour_path = root / "neighbour"
    neighbour = _ints(_list_body(strip_comments(neighbour_path.read_text()), "neighbour")) if neighbour_path.exists() else []
    boundary_text = strip_comments((root / "boundary").read_text())
    boundary = {}
    for match in re.finditer(r"(\w+)\s*\{([^{}]*)\}", boundary_text, re.S):
        block = match.group(2)
        nf = re.search(r"\bnFaces\s+(\d+)\s*;", block)
        sf = re.search(r"\bstartFace\s+(\d+)\s*;", block)
        if nf and sf:
            boundary[match.group(1)] = {"nFaces": int(nf.group(1)), "startFace": int(sf.group(1))}
    n_internal = len(neighbour)
    face_centres = []
    face_area_vectors = []
    for face in faces:
        vertices = [points[i] for i in face]
        face_centres.append(tuple(sum(v[j] for v in vertices) / len(vertices) for j in range(3)))
        area = [0.0, 0.0, 0.0]
        for a, b in zip(vertices, vertices[1:] + vertices[:1]):
            area[0] += a[1]*b[2] - a[2]*b[1]
            area[1] += a[2]*b[0] - a[0]*b[2]
            area[2] += a[0]*b[1] - a[1]*b[0]
        face_area_vectors.append(tuple(0.5*x for x in area))
    n_cells = max(owner) + 1 if owner else 0
    cells = [[] for _ in range(n_cells)]
    for face_i, cell_i in enumerate(owner):
        cells[cell_i].extend(faces[face_i])
    for face_i, cell_i in enumerate(neighbour):
        cells[cell_i].extend(faces[face_i])
    cell_centres = []
    for ids in cells:
        unique = list(dict.fromkeys(ids))
        cell_centres.append(tuple(sum(points[i][j] for i in unique)/len(unique) for j in range(3)))
    return {"points": points, "faces": faces, "owner": owner, "neighbour": neighbour,
            "boundary": boundary, "n_internal_faces": n_internal, "n_cells": n_cells,
            "face_centres": face_centres, "face_area_vectors": face_area_vectors,
            "cell_centres": cell_centres}

def foam_header(field_class: str, location: str, object_name: str) -> str:
    return ("/*--------------------------------*- C++ -*----------------------------------*\\\n"
            "| solids4foam clean artery benchmark                                  |\n"
            "\\*---------------------------------------------------------------------------*/\n"
            "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            f"    class {field_class};\n    location \"{location}\";\n    object {object_name};\n}}\n")

def _vector(v):
    return "(" + " ".join(f"{x:.16g}" for x in v) + ")"

def write_vol_vector_field(path: Path, name: str, values, boundary_values: dict):
    lines = [foam_header("volVectorField", path.parent.name, name),
             "dimensions [0 0 0 0 0 0 0];\n",
             f"internalField nonuniform List<vector>\n{len(values)}\n(\n"]
    lines += [f"    {_vector(v)}\n" for v in values]
    lines += [");\n\nboundaryField\n{\n"]
    for patch, vals in boundary_values.items():
        lines += [f"    {patch}\n    {{\n        type calculated;\n        value nonuniform List<vector>\n        {len(vals)}\n        (\n"]
        lines += [f"            {_vector(v)}\n" for v in vals]
        lines += ["        );\n    }\n"]
    lines += ["}\n"]
    path.write_text("".join(lines), encoding="utf-8")

def write_surface_vector_field(path: Path, name: str, values, boundary_values: dict):
    lines = [foam_header("surfaceVectorField", path.parent.name, name),
             "dimensions [0 0 0 0 0 0 0];\n",
             f"internalField nonuniform List<vector>\n{len(values)}\n(\n"]
    lines += [f"    {_vector(v)}\n" for v in values]
    lines += [");\n\nboundaryField\n{\n"]
    for patch, vals in boundary_values.items():
        lines += [f"    {patch}\n    {{\n        type calculated;\n        value nonuniform List<vector>\n        {len(vals)}\n        (\n"]
        lines += [f"            {_vector(v)}\n" for v in vals]
        lines += ["        );\n    }\n"]
    path.write_text("".join(lines + ["}\n"]), encoding="utf-8")
