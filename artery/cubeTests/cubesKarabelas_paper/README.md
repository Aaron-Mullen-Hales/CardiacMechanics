# Gültekin–Dal–Holzapfel unit-cube benchmarks

solids4foam cases reproducing the two unit-cube benchmarks in Sections 4.1
and 4.2 of:

> O. Gültekin, H. Dal and G. A. Holzapfel, *On the quasi-incompressible finite
> element analysis of anisotropic hyperelastic materials*, Computational
> Mechanics **63** (2019), 443–453. DOI 10.1007/s00466-018-1602-9

## Layout

```
cubes/
├── Allrun / Allclean          run or clean everything
├── scripts/
│   ├── cubeReference.py       closed-form homogeneous-deformation solution
│   └── postProcess.py         per-case post-processing and PASS/FAIL check
├── cube1/                     Sec. 4.1 / Fig. 3 — stretch ALONG the fibre
│   ├── Allrun / Allclean      run or clean both material cases
│   ├── caseA/                 Table 1 case a
│   └── caseB/                 Table 1 case b  (10x stiffer fibres)
└── cube2/                     Sec. 4.2 / Fig. 4 — stretch TRANSVERSE to the fibre
    ├── Allrun / Allclean
    ├── caseA/                 Table 1 case a
    └── caseC/                 Table 1 case c  (10x stiffer matrix)
```

Every `caseX/` is a complete, standalone solids4foam case: `cd` into it and
run `./Allrun`. The `Allrun` at each level above simply runs the ones below
and prints a combined summary.

## Running

```bash
source /path/to/OpenFOAM-v2312/etc/bashrc   # must also provide PETSC_DIR
cd cubes
./Allrun                 # all four cases
./cube1/Allrun           # just Sec. 4.1
./cube1/caseA/Allrun     # just one case
./Allclean
```

Each case takes a few seconds. `Allrun` runs `blockMesh`, `solids4Foam`, then
post-processing, which writes `postProcessing/paperQuantities.csv` and prints
a PASS/FAIL line.

## Requirements

- OpenFOAM v2312 with PETSc (`PETSC_DIR` set); the cases use
  `solutionAlgorithm PETScSNES`.
- solids4foam providing the `GultekinTwoFibreElastic` mechanical law and the
  `nonLinearGeometryTotalLagrangianTotalDisplacement` solid model, on a build
  that includes the `correctStressComponents` override in
  `GultekinTwoFibreElastic` (see "Required solids4foam version" below).
- Python 3, standard library only — no numpy or scipy. It is used for the
  post-processing check; if it is unavailable the solve still runs and only
  the automatic PASS/FAIL comparison is skipped.

### Required solids4foam version

`GultekinTwoFibreElastic` must override `correctStressComponents` to return
`sigmaPreserved = sph(sigmaToProject)`. Without it the solid model's `dev()`
projection discards the hydrostatic part of the unsplit fibre stress, the
fibres lose their resistance to volume change, and `J` grows to ~1.36 at
`lambda = 1.5` instead of ~1.001 — the cases will still run, and will still
report converged, but the physics is wrong. See "The physics being tested".

To check the build you are about to use actually has it:

```bash
grep -c correctStressComponents \
  "$SOLIDS4FOAM_SRC"/materialModels/mechanicalModel/mechanicalLaws/\
nonLinearGeometryLaws/GultekinTwoFibreElastic/GultekinTwoFibreElastic.C
# 2 = present (cell and face overrides); 0 = missing, do not trust results
```

If you have just pulled that change, rebuild before running:

```bash
cd "$SOLIDS4FOAM_SRC" && wmake libso
```

Equivalently, `Javg` in the generated `postProcessing/paperQuantities.csv`
should stay at ~1.001 for every increment. If it climbs past ~1.1, the build
is missing the override.

## The benchmark

A unit cube of transversely isotropic material, meshed with 8 hexahedral
elements (2×2×2), reinforced by a **single** fibre family `M = (1,0,0)` along
`x`. One face is pulled by a displacement ramp while the opposite transverse
faces are left traction-free, so the applied stretch is `lambda = 1 + time`,
running 1.0 → 1.5.

The two cubes differ only in the loading direction relative to the fibre:

| | loaded axis | fibre axis | what it probes |
| --- | --- | --- | --- |
| **cube1** (§4.1, Fig. 3) | `x` | `x` | stretch *along* the fibre — engages the exponential fibre term directly |
| **cube2** (§4.2, Fig. 4) | `y` | `x` | stretch *transverse* to the fibre — an isotropic direction |

Material parameters are Table 1 of the paper (quoted there in kPa; the
dictionaries use Pa):

| case | κ | μ | k₁ | k₂ | used by | varies |
| --- | --- | --- | --- | --- | --- | --- |
| a | 5000 kPa | 10 kPa | 50 kPa | 2.0 | cube1 **and** cube2 | baseline |
| b | 5000 kPa | 10 kPa | **500 kPa** | 2.0 | cube1 | stiffer fibres |
| c | 5000 kPa | **100 kPa** | 50 kPa | 2.0 | cube2 | stiffer matrix |

This pairing follows the paper: case b probes anisotropic stiffening so it is
paired with the fibre-direction test, case c probes isotropic stiffening so it
is paired with the transverse test, and case a appears in both as the common
baseline.

## The physics being tested

The paper's subject is what happens to the **volume** when an anisotropic
quasi-incompressible material is stretched. Its central result is that the
plain `Q1P0` element produces a "tremendous increase" in the average Jacobian
`Javg` — the cube swells unphysically — because splitting the anisotropic
term makes the volumetric and fibre energies compete during minimisation.
Its two remedies, `Q1P0+AL` (augmented Lagrangian) and `Q1P0+WAS` (no
anisotropic split), both keep `Javg` at ~1.

These cases use `anisotropicSplit false`, the unsplit invariant `q = I4`,
which is the analogue of the paper's **WAS** treatment. They should therefore
land on the AL/WAS branch, and they do: `Javg` stays at ~1.001 and `Uavg` at
~0 across the whole range, against roughly 1.14 (case a) and 1.37 (case b)
for plain Q1P0 in Figs. 3a/3b.

That makes `Javg` the single most diagnostic number here. If it climbs toward
1.1+ with stretch, the fibre hydrostatic stress is being discarded somewhere
(see Requirements above) and the case is reproducing the paper's *broken*
formulation rather than its recommended one.

## What is and is not reproduced

solids4foam has no Q1P0 finite element and no augmented-Lagrangian option; it
is a mixed finite-volume displacement–pressure solver. So the paper's
three-way `Q1P0` / `Q1P0+AL` / `Q1P0+WAS` comparison cannot be reproduced —
there is one formulation here, closest in spirit to WAS.

Accordingly the pass/fail check is **not** a comparison against the paper's
numbers. It is a comparison against `scripts/cubeReference.py`, an independent
closed-form solution of the same boundary value problem. Because these
deformations are spatially homogeneous, that exact solution carries no
discretisation error, and a correct solver must reproduce it to solver
tolerance on any mesh. The paper's figures are used as a qualitative
cross-check (correct branch, right energy magnitudes), not as a tolerance.

## Verified results

All four cases reach the full range with every increment genuinely converged
(PETSc reporting `CONVERGED_FNORM_RELATIVE` against the true nonlinear
residual at `-snes_atol 1e-6 -snes_rtol 1e-5`):

| case | increments | λ reached | max abs `Javg` error | max abs `I4avg` error |
| --- | --- | --- | --- | --- |
| cube1 / caseA | 10 | 1.500 | 1.1e-06 | 2.1e-08 |
| cube1 / caseB | 10 | 1.500 | 4.6e-06 | 5.3e-07 |
| cube2 / caseA | 20 | 1.500 | 1.7e-07 | 5.2e-07 |
| cube2 / caseC | 20 | 1.500 | 2.3e-07 | 7.9e-07 |

Values at `lambda = 1.5`, beside the paper's plotted AL/WAS curves (read off
the figures by eye, so treat as shape agreement, not a numerical match):

| case | Javg | Uavg (kPa) | PsiAniAvg (kPa) | PsiIsoAvg (kPa) | paper (AL/WAS) |
| --- | --- | --- | --- | --- | --- |
| cube1 a | 1.0011 | 0.003 | **272.0** | 2.91 | Ψ_ani ≈ 270 kPa (Fig. 3a) |
| cube1 b | 1.0011 | 0.003 | 2720.0 | 2.91 | Ψ_ani ≈ 2600 kPa (Fig. 3b) |
| cube2 a | 1.0015 | 0.006 | 0.1 | **3.33** | Ψ_iso ≈ 2.9 kPa (Fig. 4a) |
| cube2 c | 1.0122 | 0.369 | 1.3 | **29.21** | Ψ_iso ≈ 29 kPa (Fig. 4c) |

## Numerical settings that are tuned, not physical

Two settings are per-case and were chosen by experiment. Neither changes the
converged answer:

- **`implicitShearModulus`** (`impK`) in `constant/mechanicalProperties`. It is
  added implicitly and subtracted explicitly, and divided back out of the
  `solidTraction` update, so it cancels at convergence. Its default estimate
  `mu + 2*k1` linearises the fibre response about the *undeformed* state,
  while the true tangent grows exponentially with fibre stretch (~1000× by
  `lambda = 1.5` for case a). Since the `solidTraction` update is a
  fixed-point iteration with contraction factor ≈ `|1 − trueTangent/impK|`,
  the default makes it divergent well before the end of the ramp. Note the
  optimum is bounded on *both* sides: too small diverges, too large
  over-damps.
- **`deltaT`** in `system/controlDict`, i.e. how the quasi-static ramp is
  walked. cube1 uses 0.05 (10 increments); cube2 uses 0.025 (20), because at
  0.05 its case c was erratically sensitive to `impK`.

Both are recorded per case in the dictionaries with comments explaining the
value.

## Known limitation

The cases run the paper's own 2×2×2 mesh. A 4×4×4 refinement was also tested
during development: cube1 caseA and both cube2 cases complete the full range
on it, but **cube1 caseB stops at `lambda = 1.45`** on 4×4×4. Case b at
`lambda = 1.5` develops ~64 MPa of fibre stress against a 5 MPa bulk modulus,
i.e. the material is an order of magnitude stiffer along the fibre than in
bulk, which inverts the assumption a mixed displacement–pressure scheme rests
on. `impK`, load increment, line-search type and mesh level were all swept
without recovering that last increment. The failure is a deliberate guard
(`refusing to accept a diverged solution as the answer`), not a crash, so no
unconverged field is ever written.
