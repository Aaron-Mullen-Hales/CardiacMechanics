#!/usr/bin/env python3
"""Independent homogeneous-deformation reference for the Gultekin cube1/cube2
uniaxial-extension cases.

This is NOT a reproduction of the paper's Q1P0/Q1P0+AL/Q1P0+WAS finite
element results. It is the exact closed-form solution of the boundary value
problem that solids4foam's own mixed finite-volume model actually solves for
this geometry/loading, derived from:

  * the constitutive stress read directly from
    GultekinTwoFibreElastic.C::evaluateConstitutive (unsplit fibre invariant,
    since anisotropicSplit=false):
        sigmaIso   = (mu/J) * dev(bBar),  bBar = J^(-2/3) * F F^T
        sigmaFibre = 2*k1*(I4-1)*exp(k2*(I4-1)^2)/J * (m . m^T),  m = F.M
        sigmaPassive = sigmaIso + sigmaFibre

  * the volumetric/assembly equations read directly from the current (as
    built) nonLinGeomTotalLagTotalDispSolid.C:
        pressure residual (uniform field, stabilisation term vanishes):
            -p/K - 0.5*(J^2 - 1)/J = 0  =>  p = 0.5*K*(1 - J^2)/J
        total Cauchy stress:
            sigma = dev(sigmaPassive) - p*I

For a unit cube with a single fibre family M=(1,0,0), loaded by a prescribed
stretch along one Cartesian axis with the other two faces traction-free, the
true deformation is homogeneous and diagonal in the cube axes (by material
and geometric symmetry), F = diag(lambda1, lambda2, lambda3). Axis 1 is
always the fibre axis. Two loading cases are supported:

  fibre:     lambda1 = lambda_x prescribed (cube1, Sec. 4.1); lambda2=lambda3
             unknown and equal by transverse-isotropy symmetry (1 unknown).
  isotropic: lambda2 = lambda_y prescribed (cube2, Sec. 4.2); lambda1
             (fibre axis) and lambda3 unknown and generally different
             (2 unknowns).

The unknown stretch(es) are found by requiring zero Cauchy traction on the
free faces, using plain-Python bisection/Newton (no scipy dependency, to
match the rest of this repo's analytical checkers).

Also reports the paper's own energy quantities (U, Psi_ani, Psi_iso) computed
from the resulting fields, purely as descriptive diagnostics for comparing
trends against Fig. 3/4 of Gultekin, Dal & Holzapfel (2019) -- these are not
what solids4foam's mixed formulation actually minimises (it has no stored
volumetric potential; p is a Lagrange-multiplier-like unknown), so they are
not a pass/fail tolerance.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path


def parse_material_case(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    values = {}
    for key in ("kappa", "mu", "k1", "k2"):
        match = re.search(rf"\b{key}\s+{key}\s+\[[^\]]*\]\s+([-+0-9.eE]+)\s*;", text)
        if not match:
            # k2 is dimensionless: "k2 k2 [0 0 0 0 0 0 0] 2;"
            match = re.search(rf"\b{key}\s+{key}\s+\[[^\]]*\]\s+([-+0-9.eE]+)", text)
        if not match:
            raise ValueError(f"could not find {key} in {path}")
        values[key] = float(match.group(1))
    return values


def stresses(l1: float, l2: float, l3: float, mu: float, k1: float, k2: float, K: float):
    """Return (J, I4, p, sigmaTotal[3]) for principal stretches (l1,l2,l3),
    axis 1 = fibre direction, using the exact code equations above."""
    J = l1 * l2 * l3
    Jm23 = J ** (-2.0 / 3.0)
    b = (l1 * l1, l2 * l2, l3 * l3)
    bBar = tuple(Jm23 * x for x in b)
    trBBar = sum(bBar)
    sigmaIso = tuple((mu / J) * (bBar[i] - trBBar / 3.0) for i in range(3))

    I4 = l1 * l1
    fibreStrain = I4 - 1.0
    coeff = 2.0 * k1 * fibreStrain * math.exp(k2 * fibreStrain ** 2) / J
    sigmaFibre = (coeff * I4, 0.0, 0.0)

    sigmaPassive = tuple(sigmaIso[i] + sigmaFibre[i] for i in range(3))

    # sigma = sigmaPassive - p*I, i.e. the FULL passive stress including the
    # hydrostatic part of the unsplit fibre term. The solid model assembles
    # dev(sigmaToProject) + sigmaPreserved - p*I, and GultekinTwoFibreElastic
    # now returns sigmaPreserved = sph(sigmaPassive), so the two are
    # identical. Before that fix the fibre's hydrostatic part was discarded,
    # which removed the fibres' resistance to volume change and produced
    # large spurious volume growth (J -> 1.36 at lambda_x = 1.5 instead of
    # ~1.001) -- the Q1P0 pathology the unsplit form exists to avoid.
    p = 0.5 * K * (1.0 - J * J) / J
    sigmaTotal = tuple(sigmaPassive[i] - p for i in range(3))

    return J, I4, p, sigmaIso, sigmaFibre, sigmaPassive, sigmaTotal


def bisect(f, lo: float, hi: float, tol: float = 1e-14, maxit: int = 200) -> float:
    flo, fhi = f(lo), f(hi)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0.0:
        # scan for a bracket
        n = 400
        xs = [lo + (hi - lo) * i / n for i in range(n + 1)]
        fs = [f(x) for x in xs]
        bracket = None
        for i in range(n):
            if fs[i] == 0.0:
                return xs[i]
            if fs[i] * fs[i + 1] < 0.0:
                bracket = (xs[i], xs[i + 1])
                break
        if bracket is None:
            raise RuntimeError(f"no sign change found for f on [{lo},{hi}]")
        lo, hi = bracket
        flo, fhi = f(lo), f(hi)
    for _ in range(maxit):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if fmid == 0.0 or (hi - lo) < tol:
            return mid
        if flo * fmid <= 0.0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return 0.5 * (lo + hi)


def newton2d(F, x0, tol: float = 1e-13, maxit: int = 200):
    x = list(x0)
    for _ in range(maxit):
        fx = F(x)
        if max(abs(v) for v in fx) < tol:
            return x
        h = 1e-6
        jac = [[0.0, 0.0], [0.0, 0.0]]
        for j in range(2):
            xp = list(x)
            xp[j] += h
            fp = F(xp)
            jac[0][j] = (fp[0] - fx[0]) / h
            jac[1][j] = (fp[1] - fx[1]) / h
        det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]
        dx0 = (jac[1][1] * (-fx[0]) - jac[0][1] * (-fx[1])) / det
        dx1 = (jac[0][0] * (-fx[1]) - jac[1][0] * (-fx[0])) / det
        step = 1.0
        while step > 1e-6:
            trial = [x[0] + step * dx0, x[1] + step * dx1]
            if trial[0] > 1e-3 and trial[1] > 1e-3:
                break
            step *= 0.5
        x = [x[0] + step * dx0, x[1] + step * dx1]
    return x


def solve_fibre_direction(lam_x: float, mu: float, k1: float, k2: float, K: float) -> dict:
    """cube1: lambda1=lambda_x prescribed along the fibre; lambda2=lambda3
    unknown and equal by transverse-isotropy symmetry."""

    def residual(lt: float) -> float:
        _, _, _, _, _, _, sigmaTotal = stresses(lam_x, lt, lt, mu, k1, k2, K)
        return sigmaTotal[1]

    guess = 1.0 / math.sqrt(max(lam_x, 1e-6))
    lam_t = bisect(residual, max(1e-3, 0.2 * guess), 3.0 * guess + 1.0)

    J, I4, p, sigmaIso, sigmaFibre, sigmaPassive, sigmaTotal = stresses(
        lam_x, lam_t, lam_t, mu, k1, k2, K
    )
    Jm23 = J ** (-2.0 / 3.0)
    I1bar = Jm23 * (lam_x * lam_x + 2.0 * lam_t * lam_t)
    U = K * (J - math.log(J) - 1.0)
    psiIso = 0.5 * mu * (I1bar - 3.0)
    psiAni = k1 / (2.0 * k2) * (math.exp(k2 * (I4 - 1.0) ** 2) - 1.0)

    return {
        "lambda_x": lam_x,
        "lambda_t": lam_t,
        "J": J,
        "I4": I4,
        "p_Pa": p,
        "U_Pa": U,
        "PsiIso_Pa": psiIso,
        "PsiAni_Pa": psiAni,
        "sigmaTotal_free_face_Pa": sigmaTotal[1],
        "sigmaTotal_loaded_face_Pa": sigmaTotal[0],
    }


def solve_isotropic_direction(lam_y: float, mu: float, k1: float, k2: float, K: float) -> dict:
    """cube2: lambda2=lambda_y prescribed, orthogonal to the fibre;
    lambda1 (fibre axis) and lambda3 unknown and generally different."""

    def residuals(x):
        l1, l3 = x
        _, _, _, _, _, _, sigmaTotal = stresses(l1, lam_y, l3, mu, k1, k2, K)
        return (sigmaTotal[0], sigmaTotal[2])

    guess = 1.0 / math.sqrt(max(lam_y, 1e-6))
    lam_x, lam_z = newton2d(residuals, (guess, guess))

    J, I4, p, sigmaIso, sigmaFibre, sigmaPassive, sigmaTotal = stresses(
        lam_x, lam_y, lam_z, mu, k1, k2, K
    )
    Jm23 = J ** (-2.0 / 3.0)
    I1bar = Jm23 * (lam_x * lam_x + lam_y * lam_y + lam_z * lam_z)
    U = K * (J - math.log(J) - 1.0)
    psiIso = 0.5 * mu * (I1bar - 3.0)
    psiAni = k1 / (2.0 * k2) * (math.exp(k2 * (I4 - 1.0) ** 2) - 1.0)

    return {
        "lambda_y": lam_y,
        "lambda_x": lam_x,
        "lambda_z": lam_z,
        "J": J,
        "I4": I4,
        "p_Pa": p,
        "U_Pa": U,
        "PsiIso_Pa": psiIso,
        "PsiAni_Pa": psiAni,
        "sigmaTotal_free_face_x_Pa": sigmaTotal[0],
        "sigmaTotal_free_face_z_Pa": sigmaTotal[2],
        "sigmaTotal_loaded_face_Pa": sigmaTotal[1],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("direction", choices=("fibre", "isotropic"))
    parser.add_argument("material_case", type=Path, help="materialCases/caseX.dict")
    parser.add_argument("--lambda-min", type=float, default=1.0)
    parser.add_argument("--lambda-max", type=float, default=1.5)
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    props = parse_material_case(args.material_case)
    mu, k1, k2, K = props["mu"], props["k1"], props["k2"], props["kappa"]

    rows = []
    for i in range(args.steps + 1):
        lam = args.lambda_min + (args.lambda_max - args.lambda_min) * i / args.steps
        if args.direction == "fibre":
            rows.append(solve_fibre_direction(lam, mu, k1, k2, K))
        else:
            rows.append(solve_isotropic_direction(lam, mu, k1, k2, K))

    fieldnames = list(rows[0].keys())
    out = args.out or sys.stdout
    if isinstance(out, Path):
        stream = out.open("w", newline="", encoding="utf-8")
    else:
        stream = out
    writer = csv.DictWriter(stream, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    if isinstance(out, Path):
        stream.close()
        print(f"wrote {len(rows)} rows to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
