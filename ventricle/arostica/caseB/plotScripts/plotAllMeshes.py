#!/usr/bin/env python3
"""
Plot Aróstica Case B probe displacements for all mesh subdirectories together.

Run this from the caseB directory, e.g.

    cd .../run/CardiacMechanics/ventricle/arostica/caseB
    python3 plotAllMeshes.py

The script looks for mesh*/probeCase/postProcessing/probesDict/.../D
and overlays the curves from each mesh on the same six-panel figure.
"""

from pathlib import Path
import re
import sys

try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError as exc:
    sys.exit(
        f"Missing Python package: {exc.name}\n"
        "Install/activate numpy and matplotlib, then run this script again."
    )


def find_probe_file(mesh_dir: Path) -> Path:
    preferred = mesh_dir / "probeCase" / "postProcessing" / "probesDict" / "0" / "D"
    if preferred.is_file():
        return preferred

    candidates = list(
        (mesh_dir / "probeCase" / "postProcessing" / "probesDict").glob("*/D")
    )
    candidates += list(
        (mesh_dir / "postProcessing" / "probesDict").glob("*/D")
    )

    if not candidates:
        candidates = [
            p for p in mesh_dir.rglob("D")
            if p.is_file()
            and "postProcessing" in p.parts
            and ("probesDict" in p.parts or "probes" in p.parts)
        ]

    if not candidates:
        raise FileNotFoundError(
            f"Could not find a probe displacement file in {mesh_dir}"
        )

    def sort_key(path: Path):
        parent = path.parent.name
        try:
            t = float(parent)
        except ValueError:
            t = float("inf")
        return (0 if "probesDict" in path.parts else 1, t, str(path))

    return sorted(candidates, key=sort_key)[0]


def read_probe_file(path: Path):
    rows = []
    number_pattern = re.compile(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    )

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            values = [float(x) for x in number_pattern.findall(stripped)]
            if len(values) >= 7:
                rows.append(values[:7])

    if not rows:
        raise ValueError(f"No usable probe rows were found in {path}")

    data = np.asarray(rows, dtype=float)
    order = np.argsort(data[:, 0], kind="stable")
    data = data[order]

    unique = {}
    for row in data:
        unique[float(row[0])] = row

    return np.asarray([unique[t] for t in sorted(unique)], dtype=float)


def main():
    case_dir = Path.cwd().resolve()

    mesh_dirs = sorted(
        [d for d in case_dir.iterdir() if d.is_dir() and d.name.lower().startswith("mesh")],
        key=lambda p: p.name
    )

    if not mesh_dirs:
        sys.exit(
            f"No mesh directories found inside:\n  {case_dir}\n\n"
            "Run this from the caseB directory containing mesh1, mesh2, mesh3, ..."
        )

    mesh_data = {}
    for mesh_dir in mesh_dirs:
        try:
            probe_file = find_probe_file(mesh_dir)
            data = read_probe_file(probe_file)
            mesh_data[mesh_dir.name] = (probe_file, data)
        except Exception as exc:
            print(f"Skipping {mesh_dir.name}: {exc}")

    if not mesh_data:
        sys.exit(
            "No usable probe datasets were found in any mesh directory.\n"
            "If needed, run ./extractProbes.sh inside each mesh first."
        )

    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)

    component_map = [
        (1, "p0 - ux"),
        (4, "p1 - ux"),
        (2, "p0 - uy"),
        (5, "p1 - uy"),
        (3, "p0 - uz"),
        (6, "p1 - uz"),
    ]

    for mesh_name, (probe_file, data) in mesh_data.items():
        t = data[:, 0]
        for ax, (col, title) in zip(axes.flat, component_map):
            ax.plot(t, data[:, col], linewidth=1.8, label=mesh_name)
            ax.set_title(title)
            ax.set_ylabel("Displacement [m]")
            ax.grid(True, alpha=0.3)

    for ax in axes.flat:
        ax.legend()

    axes[2, 0].set_xlabel("Time [s]")
    axes[2, 1].set_xlabel("Time [s]")

    fig.suptitle("Aróstica Case B - all meshes")
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    output = case_dir / "caseB_allMeshes_probes.png"
    fig.savefig(output, dpi=200, bbox_inches="tight")

    print(f"Case directory : {case_dir}")
    print(f"Meshes plotted : {', '.join(mesh_data.keys())}")
    for mesh_name, (probe_file, data) in mesh_data.items():
        print(
            f"  {mesh_name:>8} : {probe_file.relative_to(case_dir)}"
            f"   [{len(data)} samples, t={data[:,0].min():.6g}->{data[:,0].max():.6g} s]"
        )
    print(f"Saved plot     : {output.name}")

    plt.show()


if __name__ == "__main__":
    main()
