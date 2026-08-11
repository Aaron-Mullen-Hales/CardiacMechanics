# cube2 — Section 4.2 / Figure 4: uniaxial extension TRANSVERSE to the fibre

Stretch is applied along `y`, an isotropic direction, while the fibre remains
along `x`. The fibre is therefore only loaded indirectly — it shortens
slightly (`I4 < 1`) — so the exponential term stays mild and this is the
better-conditioned of the two cubes. The paper uses it to show that all
three formulations agree closely in an isotropic direction.

## Cases

| directory | Table 1 case | mu | k1 | differs from case a by |
| --- | --- | --- | --- | --- |
| `caseA` | a | 10 kPa | 50 kPa | — (baseline) |
| `caseC` | c | **100 kPa** | 50 kPa | 10x stiffer matrix |

Both share the same geometry, mesh, fibre field, boundary conditions and
loading; only the material parameters differ.

## Setup

- Unit cube `(0,0,0)`–`(1,1,1)`, meshed `2 x 2 x 2` = 8 hexahedra (the paper's
  "8 unstructured hexahedral elements").
- One fibre family, `M = (1,0,0)` along `x`, supplied by the `0/f0` (cell) and
  `0/f0f` (face) fields.
- **Loading: stretch along `y`, transverse to the fibre.** `yMin` is fixed in its normal direction,
  `yMax` is driven by the displacement ramp in
  `constant/load/timeVsDisplacement` (0 → 0.5 over t = 0 → 0.5, so
  `lambda = 1 + time`, 1.0 → 1.5). `xMin`/`zMin` are fixed in their normal directions and `xMax`/`zMax` are traction-free.
- Material law `GultekinTwoFibreElastic` with `anisotropicSplit false`
  (unsplit invariant `q = I4`, the analogue of the paper's WAS treatment) and
  `useSecondFibreFamily false`.
- Solid model `nonLinearGeometryTotalLagrangianTotalDisplacement`, mixed
  displacement–pressure, PETSc SNES with a matrix-free (JFNK) Newton step.
- Load increment `deltaT = 0.025` (20 increments).

## Running

```bash
./Allrun          # both cases
./caseA/Allrun    # one case
./Allclean
```

Each case writes `postProcessing/paperQuantities.csv` and prints a PASS/FAIL
line comparing against the closed-form reference in `../scripts/`.

## Expected result

Both cases reach the full range with every increment genuinely converged.
Accuracy here is the best of the four cases (errors ~1e-07), reflecting the
milder fibre engagement.

| case | increments | λ reached | max abs `Javg` error |
| --- | --- | --- | --- |
| `caseA` | 20 | 1.500 | 1.7e-07 |
| `caseC` | 20 | 1.500 | 2.3e-07 |

See `../README.md` for the physics being tested, why `Javg` is the most
diagnostic quantity, and the tuned-but-non-physical numerical settings.
