#!/usr/bin/env python3
"""Generate the six exact-2x Aróstica mesh levels as complete case copies."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from conforming_hexmesh import build_biventricle, build_monoventricle, write_metadata


ROOT = Path(__file__).resolve().parents[1]
REFINEMENTS = ROOT / "refinements"

# mesh3 is the checked default mesh. Every adjacent cell count is exactly 2x.
# All six biventricle levels use uniform chamber-wall stacks (free == facing)
# because this is substantially better conditioned than routing any layer
# count through the abrupt 6-to-2 transition template: forcing that template
# even at the coarse/canonical level produced a near-degenerate cell
# (Jacobian determinant 0.0051) and 168 severely non-orthogonal faces, while
# the uniform construction below is checkMesh-clean end to end. Free/facing
# and septum layers are held at 6/2 across levels 2-6 (3/1 at level 1, the
# smallest mesh that stays non-degenerate at that layer count); only surface
# and longitudinal divisions vary, chosen by exact search so every level is
# precisely 2x the previous cell count.
LEVELS = (
    {
        "level": 1,
        "monoventricle": {"surface": 8, "longitudinal": 9, "wall": 3},
        "biventricle": {"surface": 8, "longitudinal": 9, "free": 3, "facing": 3, "septum": 1},
    },
    {
        "level": 2,
        "monoventricle": {"surface": 8, "longitudinal": 20, "wall": 3},
        "biventricle": {"surface": 8, "longitudinal": 9, "free": 6, "facing": 6, "septum": 2},
    },
    {
        "level": 3,
        "monoventricle": {"surface": 8, "longitudinal": 20, "wall": 6},
        "biventricle": {"surface": 8, "longitudinal": 20, "free": 6, "facing": 6, "septum": 2},
    },
    {
        "level": 4,
        "monoventricle": {"surface": 8, "longitudinal": 42, "wall": 6},
        "biventricle": {"surface": 8, "longitudinal": 42, "free": 6, "facing": 6, "septum": 2},
    },
    {
        "level": 5,
        "monoventricle": {"surface": 8, "longitudinal": 42, "wall": 12},
        "biventricle": {
            "surface": 16, "longitudinal": 40, "free": 6, "facing": 6, "septum": 2,
            "apex_cap_fraction": 0.155,
        },
    },
    {
        "level": 6,
        "monoventricle": {"surface": 8, "longitudinal": 86, "wall": 12},
        "biventricle": {
            "surface": 16, "longitudinal": 84, "free": 6, "facing": 6, "septum": 2,
            "apex_cap_fraction": 0.155,
        },
    },
)


def copy_case_template(name: str, destination: Path) -> None:
    source = ROOT / "cases" / name
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(
        source,
        destination,
        ignore=shutil.ignore_patterns("VTK", "processor*", "log.*", "sets"),
    )
    for path in destination.glob("[1-9]*"):
        if path.is_dir():
            shutil.rmtree(path)
    (destination / "case.foam").touch()


def main() -> None:
    reports = {}
    previous_counts = {"monoventricle": None, "biventricle": None}
    for definition in LEVELS:
        level = definition["level"]
        level_root = REFINEMENTS / f"mesh{level}"
        mono_case = level_root / "monoventricle"
        biv_case = level_root / "biventricle"
        copy_case_template("monoventricle", mono_case)
        copy_case_template("biventricle", biv_case)

        mono = definition["monoventricle"]
        mono_summary, mono_data = build_monoventricle(
            mono_case, mono["surface"], mono["longitudinal"], mono["wall"]
        )
        write_metadata(mono_case, mono_summary, mono_data)

        biv = definition["biventricle"]
        biv_summary, biv_data = build_biventricle(
            biv_case,
            biv["surface"],
            biv["longitudinal"],
            biv.get("apex_cap_fraction", 0.187),
            biv["free"],
            biv["facing"],
            biv["septum"],
            1,
            20,
            3,
        )
        write_metadata(biv_case, biv_summary, biv_data)

        summaries = {"monoventricle": mono_summary, "biventricle": biv_summary}
        for name, summary in summaries.items():
            count = summary["cell_count"]
            previous = previous_counts[name]
            if previous is not None and count != 2 * previous:
                raise RuntimeError(
                    f"{name} mesh{level} has {count} cells; expected exactly {2 * previous}"
                )
            previous_counts[name] = count
        reports[f"mesh{level}"] = {
            "parameters": definition,
            "monoventricle": mono_summary,
            "biventricle": biv_summary,
        }
        print(
            f"mesh{level}: mono={mono_summary['cell_count']} hex, "
            f"biv={biv_summary['cell_count']} hex"
        )

    (ROOT / "reports").mkdir(exist_ok=True)
    (ROOT / "reports" / "refinement_generation.json").write_text(
        json.dumps(reports, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
