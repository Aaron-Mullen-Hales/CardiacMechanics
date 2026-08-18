#!/usr/bin/env python3
"""
Plot Aróstica Case B probe displacements for whichever mesh directory
this script is run from.

Usage:
    cd .../caseB/mesh1
    python3 plotProbes.py

The script searches the CURRENT DIRECTORY for the probe output rather than
hard-coding mesh1/mesh2/mesh3.
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


def find_probe_file(case_dir: Path) -> Path:
    preferred = case_dir / "probeCase" / "postProcessing" / "probesDict" / "0" / "D"
    if preferred.is_file():
        return preferred

    candidates = list(
        (case_dir / "probeCase" / "postProcessing" / "probesDict").glob("*/D")
    )
    candidates += list(
        (case_dir / "postProcessing" / "probesDict").glob("*/D")
    )

    if not candidates:
        candidates = [
            p for p in case_dir.rglob("D")
            if p.is_file()
            and "postProcessing" in p.parts
            and ("probesDict" in p.parts or "probes" in p.parts)
        ]

    if not candidates:
        raise FileNotFoundError(
            "Could not find a probe displacement file D.\n\n"
            f"Searched inside:\n  {case_dir}\n\n"
            "Expected something like:\n"
            "  probeCase/postProcessing/probesDict/0/D\n\n"
            "If probes have not been extracted, run:\n"
            "  ./extractProbes.sh"
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
        raise ValueError(
            f"No usable probe rows were found in:\n  {path}\n\n"
            "Expected rows containing:\n"
            "  time (p0x p0y p0z) (p1x p1y p1z)"
        )

    data = np.asarray(rows, dtype=float)
    order = np.argsort(data[:, 0], kind="stable")
    data = data[order]

    unique = {}
    for row in data:
        unique[float(row[0])] = row

    return np.asarray([unique[t] for t in sorted(unique)], dtype=float)


def main():
    case_dir = Path.cwd().resolve()
    mesh_name = case_dir.name

    try:
        probe_file = find_probe_file(case_dir)
        data = read_probe_file(probe_file)
    except (FileNotFoundError, ValueError) as exc:
        sys.exit(str(exc))

    t = data[:, 0]

    p0x, p0y, p0z = data[:, 1], data[:, 2], data[:, 3]
    p1x, p1y, p1z = data[:, 4], data[:, 5], data[:, 6]

    curves = [
        (p0x, "p0 - ux"),
        (p1x, "p1 - ux"),
        (p0y, "p0 - uy"),
        (p1y, "p1 - uy"),
        (p0z, "p0 - uz"),
        (p1z, "p1 - uz"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)

    for ax, (values, title) in zip(axes.flat, curves):
        ax.plot(t, values, linewidth=1.8)
        ax.set_title(title)
        ax.set_ylabel("Displacement [m]")
        ax.grid(True, alpha=0.3)

    axes[2, 0].set_xlabel("Time [s]")
    axes[2, 1].set_xlabel("Time [s]")

    fig.suptitle(f"Aróstica Case B - {mesh_name}")
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    output = case_dir / f"caseB_{mesh_name}_probes.png"
    fig.savefig(output, dpi=200, bbox_inches="tight")

    print(f"Mesh directory : {case_dir}")
    print(f"Probe file     : {probe_file.relative_to(case_dir)}")
    print(f"Samples        : {len(t)}")
    print(f"Time range     : {t.min():.6g} -> {t.max():.6g} s")
    print(f"Saved plot     : {output.name}")

    plt.show()


if __name__ == "__main__":
    main()
