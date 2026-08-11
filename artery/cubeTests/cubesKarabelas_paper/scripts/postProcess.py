#!/usr/bin/env python3
"""Post-process one Gultekin cube case.

Run from inside a case directory (or pass the case directory as an argument).
Reads the material parameters straight out of constant/mechanicalProperties
and the loading direction out of 0/D, so there is nothing to keep in sync by
hand, then for every written time writes

    postProcessing/paperQuantities.csv

containing the quantities the paper plots in Figures 3 and 4 --

    lambda      applied stretch, = 1 + time
    Javg        volume-average Jacobian                        [-]
    Uavg        volume-average volumetric energy   U(J)        [kPa]
    PsiAniAvg   volume-average anisotropic energy  Psi_ani     [kPa]
    PsiIsoAvg   volume-average isotropic energy    Psi_iso     [kPa]

-- each beside the closed-form value from cubeReference.py, plus the error.
A PASS/FAIL summary is printed.

The deformation here is spatially homogeneous, so the exact solution carries
no discretisation error and the simulation should reproduce it to solver
tolerance. That is what makes this a meaningful check rather than a plot.
"""
from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cubeReference import (  # noqa: E402
    solve_fibre_direction,
    solve_isotropic_direction,
)

NUMBER = r"[-+0-9.eE]+"

# Tolerances. These are homogeneous deformations, so achievable accuracy is
# set by the nonlinear solve alone (-snes_rtol 1e-5). p is recovered from
# p = 0.5*K*(1-J^2)/J, so an error in J is amplified by the bulk modulus
# (dp ~ K*dJ); its tolerance therefore tracks TOL_J rather than being an
# independent requirement.
TOL_J = 1.0e-4
TOL_I4 = 1.0e-5


def read_field(path: Path, width: int) -> list[list[float]]:
    text = path.read_text(encoding="utf-8")
    marker = text.find("internalField")
    section = text[marker: text.find("boundaryField", marker)]
    # "uniform" is a substring of "nonuniform"; match the keyword, not text
    if re.search(r"internalField\s+uniform\b", section):
        vals = [float(x) for x in
                re.findall(NUMBER, section.split("uniform", 1)[1].split(";", 1)[0])]
        return [vals[:width]]
    count = int(re.search(r"nonuniform\s+List<[^>]+>\s+(\d+)", section).group(1))
    body = section[section.index("("):]
    if width == 1:
        return [[v] for v in [float(x) for x in re.findall(NUMBER, body)][:count]]
    return [[float(x) for x in re.findall(NUMBER, t)][:width]
            for t in re.findall(r"\(([^()]*)\)", body)[:count]]


def read_material(case: Path) -> dict:
    text = (case / "constant" / "mechanicalProperties").read_text(encoding="utf-8")
    out = {}
    for key in ("mu", "k1", "k2", "bulkModulus"):
        m = re.search(rf"^\s*{key}\s+\w+\s*\[[^\]]*\]\s*({NUMBER})\s*;", text, re.M)
        if not m:
            raise SystemExit(f"could not read '{key}' from {case}/constant/mechanicalProperties")
        out[key] = float(m.group(1))
    return out


def read_direction(case: Path) -> str:
    """Fibre is always M=(1,0,0) i.e. x. The loaded face is whichever patch
    carries the displacementSeries ramp, so read that rather than assume."""
    text = (case / "0" / "D").read_text(encoding="utf-8")
    idx = text.find("displacementSeries")
    if idx < 0:
        raise SystemExit(f"no displacementSeries ramp found in {case}/0/D")
    patch = None
    for m in re.finditer(r"^\s{4}(\w+)\s*$", text[:idx], re.M):
        patch = m.group(1)
    if patch is None:
        raise SystemExit(f"could not identify the loaded patch in {case}/0/D")
    return "fibre" if patch.startswith("x") else "isotropic"


def main() -> int:
    case = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    mat = read_material(case)
    mu, k1, k2, K = mat["mu"], mat["k1"], mat["k2"], mat["bulkModulus"]
    direction = read_direction(case)
    solver = solve_fibre_direction if direction == "fibre" else solve_isotropic_direction

    times = sorted(
        (p for p in case.iterdir()
         if p.is_dir() and re.fullmatch(r"[0-9]+(\.[0-9]+)?", p.name) and float(p.name) > 0),
        key=lambda p: float(p.name),
    )
    if not times:
        print(f"No solved time directories in {case}. Run ./Allrun first.", file=sys.stderr)
        return 1

    rows, worst_j, worst_i4 = [], 0.0, 0.0
    for td in times:
        lam = 1.0 + float(td.name)
        J = [v[0] for v in read_field(td / "GTF_J", 1)]
        I4 = [v[0] for v in read_field(td / "GTF_I4", 1)]
        F = read_field(td / "F", 9)
        n = len(J)

        # I1bar = J^(-2/3) tr(b), with tr(b) = F : F
        I1bar = [J[i] ** (-2.0 / 3.0) * sum(f * f for f in F[i]) for i in range(n)]

        Javg = sum(J) / n
        I4avg = sum(I4) / n
        Uavg = sum(K * (j - math.log(j) - 1.0) for j in J) / n
        PsiIso = sum(0.5 * mu * (b - 3.0) for b in I1bar) / n
        PsiAni = sum(k1 / (2.0 * k2) * (math.exp(k2 * (i - 1.0) ** 2) - 1.0) for i in I4) / n

        ref = solver(lam, mu, k1, k2, K)
        j_err, i4_err = abs(Javg - ref["J"]), abs(I4avg - ref["I4"])
        worst_j, worst_i4 = max(worst_j, j_err), max(worst_i4, i4_err)

        rows.append({
            "lambda": f"{lam:.4f}",
            "Javg_sim": f"{Javg:.8f}", "Javg_ref": f"{ref['J']:.8f}", "Javg_err": f"{j_err:.3e}",
            "I4avg_sim": f"{I4avg:.8f}", "I4avg_ref": f"{ref['I4']:.8f}", "I4avg_err": f"{i4_err:.3e}",
            "Uavg_kPa_sim": f"{Uavg/1e3:.6f}", "Uavg_kPa_ref": f"{ref['U_Pa']/1e3:.6f}",
            "PsiAniAvg_kPa_sim": f"{PsiAni/1e3:.6f}", "PsiAniAvg_kPa_ref": f"{ref['PsiAni_Pa']/1e3:.6f}",
            "PsiIsoAvg_kPa_sim": f"{PsiIso/1e3:.6f}", "PsiIsoAvg_kPa_ref": f"{ref['PsiIso_Pa']/1e3:.6f}",
        })

    out_dir = case / "postProcessing"
    out_dir.mkdir(exist_ok=True)
    out = out_dir / "paperQuantities.csv"
    with out.open("w", newline="", encoding="utf-8") as s:
        w = csv.DictWriter(s, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # Did every increment converge, and did we reach the end of the ramp?
    log = case / "log.solids4Foam"
    n_bad = log.read_text(errors="ignore").count("did not converge") if log.exists() else -1
    lam_end = float(rows[-1]["lambda"])

    ok = worst_j <= TOL_J and worst_i4 <= TOL_I4 and n_bad == 0 and abs(lam_end - 1.5) < 1e-9
    print(f"  loading            : {direction} direction")
    print(f"  increments written : {len(rows)}   (lambda 1.0 -> {lam_end:.3f})")
    print(f"  not converged      : {'unknown (no log)' if n_bad < 0 else n_bad}")
    print(f"  max |Javg  error|  : {worst_j:.3e}   (tol {TOL_J:g})")
    print(f"  max |I4avg error|  : {worst_i4:.3e}   (tol {TOL_I4:g})")
    print(f"  csv                : {out.relative_to(case) if out.is_relative_to(case) else out}")
    print(f"  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
