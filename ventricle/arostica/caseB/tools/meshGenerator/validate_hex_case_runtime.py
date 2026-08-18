#!/usr/bin/env python3
"""Validate clean-build runtime selection on canonical and refinement cases."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from validate_solids4foam import validate_case_root


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--canonical-only",
        action="store_true",
        help="check cases/* only, omitting refinements/mesh1..mesh6",
    )
    args = parser.parse_args()

    case_sets = [("canonical_mesh3", ROOT / "cases")]
    if not args.canonical_only:
        case_sets.extend(
            (f"mesh{level}", ROOT / "refinements" / f"mesh{level}")
            for level in range(1, 7)
        )

    report = {}
    try:
        for label, cases_root in case_sets:
            report[label] = validate_case_root(
                cases_root,
                log_name="log.runtime.cleanBuild",
            )
    except (FileNotFoundError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    report["all_passed"] = all(
        report[label]["all_passed"] for label, _ in case_sets
    )
    report["case_set_count"] = len(case_sets)
    report["case_count"] = 2 * len(case_sets)
    destination = ROOT / "reports" / "hex_case_runtime.json"
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"clean-build runtime validation: {report['case_count']} cases, "
        f"all_passed={report['all_passed']}"
    )
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
