# cube1 — Section 4.1 / Figure 3: uniaxial extension ALONG the fibre

Stretch is applied along `x`, which is also the fibre direction, so the
exponential fibre term is engaged directly from the first increment. This is
the harder of the two cubes and the one where the paper shows plain Q1P0
swelling most severely (`Javg` climbing to ~1.14 for case a and ~1.37 for
case b by `lambda = 1.5`).

## Cases

| directory | Table 1 case | mu | k1 | differs from case a by |
| --- | --- | --- | --- | --- |
| `caseA` | a | 10 kPa | 50 kPa | — (baseline) |
| `caseB` | b | 10 kPa | **500 kPa** | 10x stiffer fibres |

Both share the same geometry, mesh, fibre field, boundary conditions and
loading; only the material parameters differ.

## Setup

- Unit cube `(0,0,0)`–`(1,1,1)`, meshed `2 x 2 x 2` = 8 hexahedra (the paper's
  "8 unstructured hexahedral elements").
- One fibre family, `M = (1,0,0)` along `x`, supplied by the `0/f0` (cell) and
  `0/f0f` (face) fields.
- **Loading: stretch along `x`, the fibre axis.** `xMin` is fixed in its normal direction,
  `xMax` is driven by the displacement ramp in
  `constant/load/timeVsDisplacement` (0 → 0.5 over t = 0 → 0.5, so
  `lambda = 1 + time`, 1.0 → 1.5). `yMin`/`zMin` are fixed in their normal directions and `yMax`/`zMax` are traction-free.
- Material law `GultekinTwoFibreElastic` with `anisotropicSplit false`
  (unsplit invariant `q = I4`, the analogue of the paper's WAS treatment) and
  `useSecondFibreFamily false`.
- Solid model `nonLinearGeometryTotalLagrangianTotalDisplacement`, mixed
  displacement–pressure, PETSc SNES with a matrix-free (JFNK) Newton step.
- Load increment `deltaT = 0.05` (10 increments).

## Running

```bash
./Allrun          # both cases
./caseA/Allrun    # one case
./Allclean
```

Each case writes `postProcessing/paperQuantities.csv` and prints a PASS/FAIL
line comparing against the closed-form reference in `../scripts/`.

## Expected result

Both cases reach the full range with every increment genuinely converged
and `Javg` staying at ~1.001 — i.e. on the paper's Q1P0+AL / Q1P0+WAS
branch, not the plain-Q1P0 branch.

| case | increments | λ reached | max abs `Javg` error |
| --- | --- | --- | --- |
| `caseA` | 10 | 1.500 | 1.1e-06 |
| `caseB` | 10 | 1.500 | 4.6e-06 |

See `../README.md` for the physics being tested, why `Javg` is the most
diagnostic quantity, and the tuned-but-non-physical numerical settings.
