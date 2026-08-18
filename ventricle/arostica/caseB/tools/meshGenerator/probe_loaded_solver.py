#!/usr/bin/env python3
"""Run a pressure-loaded increment on a temporary case with PETSc variants."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import tempfile
from pathlib import Path

from validate_solids4foam import ignore_generated


ROOT = Path(__file__).resolve().parents[1]

COMMON = """\
-snes_type newtonls
-snes_linesearch_type bt
-snes_rtol 1e-7
-snes_atol 1e-12
-snes_stol 1e-8
-snes_max_it 20
-snes_monitor
-snes_converged_reason
-snes_max_funcs 1000
-snes_mf
-snes_mf_operator
-ksp_type lgmres
-ksp_gmres_restart 200
-ksp_rtol 1e-2
-ksp_max_it 500
-ksp_converged_reason
"""

ASSEMBLED_COMMON = """\
-snes_type newtonls
-snes_linesearch_type bt
-snes_rtol 1e-7
-snes_atol 1e-12
-snes_stol 1e-10
-snes_max_it 200
-snes_monitor
-snes_converged_reason
-snes_max_funcs 2000
-ksp_type lgmres
-ksp_gmres_restart 200
-ksp_rtol 1e-8
-ksp_max_it 1000
-ksp_converged_reason
"""

VARIANTS = {
    "fieldsplit": COMMON + """\
-pc_type fieldsplit
-pc_fieldsplit_block_size 4
-pc_fieldsplit_0_fields 0,1,2
-pc_fieldsplit_1_fields 3
-pc_fieldsplit_type schur
-pc_fieldsplit_schur_factorization_type full
-pc_fieldsplit_schur_precondition selfp
-fieldsplit_0_ksp_type preonly
-fieldsplit_0_pc_type hypre
-fieldsplit_0_pc_hypre_type boomeramg
-fieldsplit_1_ksp_type preonly
-fieldsplit_1_pc_type hypre
-fieldsplit_1_pc_hypre_type boomeramg
""",
    "hypre": COMMON + """\
-pc_type hypre
-pc_hypre_type boomeramg
-pc_hypre_boomeramg_max_iter 1
-pc_hypre_boomeramg_strong_threshold 0.7
-pc_hypre_boomeramg_grid_sweeps_up 1
-pc_hypre_boomeramg_grid_sweeps_down 1
-pc_hypre_boomeramg_agg_nl 1
-pc_hypre_boomeramg_agg_num_paths 1
-pc_hypre_boomeramg_max_levels 25
-pc_hypre_boomeramg_coarsen_type HMIS
-pc_hypre_boomeramg_interp_type ext+i
-pc_hypre_boomeramg_P_max 1
-pc_hypre_boomeramg_truncfactor 0.3
""",
    "asm_ilu": COMMON + """\
-pc_type asm
-sub_pc_type ilu
-sub_pc_factor_levels 2
""",
    "bjacobi_lu": COMMON + """\
-pc_type bjacobi
-sub_pc_type lu
""",
    "mumps_lu": COMMON + """\
-pc_type lu
-pc_factor_mat_solver_type mumps
""",
    "assembled_hypre": ASSEMBLED_COMMON + """\
-pc_type hypre
-pc_hypre_type boomeramg
-pc_hypre_boomeramg_max_iter 1
-pc_hypre_boomeramg_strong_threshold 0.7
""",
    "assembled_mumps": ASSEMBLED_COMMON + """\
-pc_type lu
-pc_factor_mat_solver_type mumps
""",
    "assembled_mumps_basic": ASSEMBLED_COMMON.replace(
        "-snes_linesearch_type bt", "-snes_linesearch_type basic"
    ).replace("-snes_max_it 200", "-snes_max_it 1000") + """\
-pc_type lu
-pc_factor_mat_solver_type mumps
""",
    "jfnk_fieldsplit_inexact": COMMON.replace(
        "-ksp_rtol 1e-2", "-ksp_rtol 1e-1\n-snes_ksp_ew"
    ).replace("-ksp_max_it 500", "-ksp_max_it 2000") + """\
-pc_type fieldsplit
-pc_fieldsplit_block_size 4
-pc_fieldsplit_0_fields 0,1,2
-pc_fieldsplit_1_fields 3
-pc_fieldsplit_type schur
-pc_fieldsplit_schur_factorization_type full
-pc_fieldsplit_schur_precondition selfp
-fieldsplit_0_ksp_type preonly
-fieldsplit_0_pc_type hypre
-fieldsplit_1_ksp_type preonly
-fieldsplit_1_pc_type hypre
""",
    "fieldsplit_production": COMMON.replace(
        "-snes_rtol 1e-7", "-snes_rtol 1e-6"
    ).replace("-ksp_max_it 500", "-ksp_max_it 1000") + """\
-pc_type fieldsplit
-pc_fieldsplit_block_size 4
-pc_fieldsplit_0_fields 0,1,2
-pc_fieldsplit_1_fields 3
-pc_fieldsplit_type schur
-pc_fieldsplit_schur_factorization_type full
-pc_fieldsplit_schur_precondition selfp
-fieldsplit_0_ksp_type preonly
-fieldsplit_0_pc_type hypre
-fieldsplit_0_pc_hypre_type boomeramg
-fieldsplit_1_ksp_type preonly
-fieldsplit_1_pc_type hypre
-fieldsplit_1_pc_hypre_type boomeramg
""",
}


def run(executable: str, case: Path, timeout: int) -> tuple[int, str, bool]:
    process = subprocess.Popen(
        [executable],
        cwd=case,
        env=os.environ.copy(),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        start_new_session=True,
    )
    try:
        output, _ = process.communicate(timeout=timeout)
        return process.returncode, output, False
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
        return process.returncode, output, True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("variant", choices=tuple(VARIANTS))
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--delta-t", type=float)
    parser.add_argument("--steps", type=int, default=1)
    args = parser.parse_args()

    executable = shutil.which("solids4Foam")
    if executable is None:
        parser.error("solids4Foam is not available in PATH")
    source = args.case.resolve()
    if not (source / "system" / "controlDict").is_file():
        parser.error(f"not an OpenFOAM case: {source}")

    label = "_".join(source.relative_to(ROOT).parts)
    report_dir = ROOT / "reports" / "solver_probes"
    report_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="arostica_loaded_probe_") as temporary:
        target = Path(temporary) / source.name
        shutil.copytree(source, target, ignore=ignore_generated)
        (target / "petscOptions.hex").write_text(VARIANTS[args.variant])
        if args.delta_t is not None:
            control = target / "system" / "controlDict"
            control_text = control.read_text()
            control_text = re.sub(
                r"(?m)^deltaT\s+[^;]+;",
                f"deltaT {args.delta_t:.16g};",
                control_text,
            )
            control_text = re.sub(
                r"(?m)^endTime\s+[^;]+;",
                f"endTime {args.delta_t * args.steps:.16g};",
                control_text,
            )
            control.write_text(control_text)
        code, log, timed_out = run(executable, target, args.timeout)

    norms = [float(v) for v in re.findall(r"SNES Function norm ([0-9.eE+-]+)", log)]
    linear_iterations = [
        int(v) for v in re.findall(r"Linear solve .* iterations (\d+)", log)
    ]
    completed = code == 0 and "End" in log and "FOAM FATAL" not in log
    converged = "Nonlinear solve converged" in log
    mesh_failure = any(
        token in log.lower()
        for token in ("negative cell", "negative jacobian", "invalid supplied aróstica")
    ) or re.search(r"(^|[^a-z])nan([^a-z]|$)", log.lower()) is not None
    result = {
        "case": str(source.relative_to(ROOT)),
        "variant": args.variant,
        "delta_t": args.delta_t,
        "steps": args.steps,
        "return_code": code,
        "timed_out": timed_out,
        "completed": completed,
        "nonlinear_converged": converged,
        "snes_norms": norms,
        "linear_iterations": linear_iterations,
        "mesh_or_nan_signature": mesh_failure,
    }
    time_tag = "defaultDt" if args.delta_t is None else f"dt{args.delta_t:g}"
    stem = f"{label}__{args.variant}__{time_tag}"
    (report_dir / f"{stem}.log").write_text(log)
    (report_dir / f"{stem}.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if completed else 1


if __name__ == "__main__":
    raise SystemExit(main())
