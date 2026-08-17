# Aróstica Case B — monoventricle, Step 0, Case B

Aróstica et al., *A software benchmark for cardiac elastodynamics*, CMAME **435** (2025)
117485. Benchmark 1, monoventricle, **Step 0, Case B**: passive response, `τ(t) = 0`,
endocardial pressure only.

This directory is a **self-contained, runnable case**. Mesh and fibres are pre-generated
and committed here, so there is no meshing step.

---

## Quick start

```bash
./Allrun                  # run ALL meshes: mesh1, mesh2, mesh3
./Allrun mesh1            # run just one
./Allrun mesh1 mesh3      # run a subset
./Allclean                # remove results, KEEP the meshes and the case runnable
```

Each mesh runs in its own directory (`mesh1/`, `mesh2/`, `mesh3/`) assembled from the
shared case setup plus that mesh, so the three resolutions never overwrite each other.
Results, `log.solids4Foam` and probe output all stay inside the per-mesh directory.

`Allrun` sources the OpenFOAM environment itself if you have not already, runs
`solids4Foam` into `log.solids4Foam`, and then calls `./extractProbes.sh`.

To run a shorter test, edit `endTime` in `system/controlDict`. The scripts
deliberately never rewrite that file: `foamDictionary -set` would strip the header
and every comment out of it.

`Allclean` deliberately does **not** use OpenFOAM's `cleanCase`, because `cleanCase`
deletes `constant/polyMesh`. That is right for a tutorial that rebuilds its mesh with
`blockMesh`, but this case ships a pre-generated mesh, so `cleanCase` leaves it
unrunnable (`Cannot find file "points" in directory "polyMesh"`). A pristine copy is
also kept at `constant/polyMesh.orig` in case the mesh is ever lost.

Probe output (the benchmark's two observation points) ends up in:

```
probeCase/postProcessing/probesDict/0/D
```

Columns are `time`, then `u(p0)`, then `u(p1)` as `(ux uy uz)`, with
`p0 = (0.025, 0.03, 0)` and `p1 = (0, 0.03, 0)`.

---

## What is in here

```
0/                      D, p, pointD                 initial/boundary fields
                        f0, s0, n0                   cell fibre / sheet / sheet-normal
                        f0f, s0f, n0f                face fibre triad
                        t                            transmural coordinate (0 endo, 1 epi)
meshes/
    mesh1/              1056 cells  (8 x  9 x 3)   polyMesh + its fibre fields
    mesh2/              2048 cells  (8 x 14 x 4)
    mesh3/              4224 cells  (8 x 20 x 6)
constant/
    polyMesh/           default mesh, overwritten per run by Allrun
    mechanicalProperties  material law + all Table 1/2 coefficients
    solidProperties     solid model, mixed p-D settings, stabilisation
    physicsProperties, g, dynamicMeshDict
    loadCurves/         tabulated chamber pressure (reference only; the BC
                        integrates the ODE at runtime)
system/
    controlDict         dt, endTime, write interval, and the `libs` entry
    fvSchemes           BDF2 in time, leastSquaresS4f gradients
    fvSolution
    probesDict          the p0/p1 probe definition
lib/                    case-local boundary-condition library.  NOT tracked in
                        git (platform-specific build artefact) -- build it with
                        the command under "Case-local boundary conditions".
petscOptions.snes       PETSc/SNES options
tools/neutraliseD.py    helper used by extractProbes.sh
```

### Case-local boundary conditions

Three custom BCs live in `lib/libArosticaHexBoundaryConditions.dylib`, loaded through
`libs` in `system/controlDict`:

| patch | BC type | what it does |
|---|---|---|
| `endocardium` | `arosticaActivationPressureTraction` | integrates the Bestel pressure ODE (Eqs. 7–8) at runtime and applies `−p n` on the **current** area (follower load) |
| `epicardium` | `arosticaNormalSpringDashpotTraction` | `−[α(D·N) + β(Ḋ·N)]N` on the **reference** area, normal only, zero shear |
| `base` | `arosticaVectorSpringDashpotTraction` | `−[αD + βḊ]` on the **reference** area, full vector |

Source: `../HexMeshesTests/src/arosticaBoundaryConditions/`.  The library is a
build artefact and is deliberately not committed; build it with:

```bash
export S4F_CARDIAC_CLEAN=/Volumes/OpenFoam/aaronmullen-hales-v2312/s4f-cardiac-clean
export AROSTICA_HEX_LIBBIN=$PWD/lib
cd ../HexMeshesTests/src/arosticaBoundaryConditions && wmake libso
cp $AROSTICA_HEX_LIBBIN/libArosticaHexBoundaryConditions.dylib <caseB>/lib/
```

Without those two environment variables `wmake` fails on `/etc/wmake-options` and
**silently leaves the previous library in place**, which looks like a successful
build. Check the timestamp.

### Physical parameters (all verified against the paper and against Simula's implementation)

| | value | | value |
|---|---|---|---|
| `rho` | 1000 kg/m³ | `a`, `b` | 59 Pa, 8.023 |
| `eta` | 100 Pa·s | `a_f`, `b_f` | 18472 Pa, 16.026 |
| `kappa` | 1e6 Pa | `a_s`, `b_s` | 2481 Pa, 11.120 |
| `k` (switch) | 100 | `a_fs`, `b_fs` | 216 Pa, 11.436 |
| `alpha_epi` | 1e8 Pa/m | `alpha_top` | 1e5 Pa/m |
| `beta_epi` | 5e3 Pa·s/m | `beta_top` | 5e3 Pa·s/m |
| fibre angles | −60° endo → +60° epi | `tau` | **0** (Case B) |

### Meshes

Three resolutions, generated by the case-local generator in `tools/meshGenerator/`.
Cell count is `(d² + 4·d·l)·m` for `--surface-divisions d --longitudinal-divisions l
--transmural-layers m`:

| | d × l × m | cells | through wall | epi cell size min–max [mm] | max/min |
|---|---|---|---|---|---|
| `mesh1` | 6 × 15 × 3 | 1188 | 3 | 3.72 – 9.31 | 2.50 |
| `mesh2` | 8 × 20 × 4 | 2816 | 4 | 2.78 – 7.05 | 2.54 |
| `mesh3` | 10 × 26 × 5 | 5700 | 5 | 2.22 – 5.55 | 2.50 |
| *(paper's own tet mesh, for reference)* | — | 4052 | — | 2.67 – 6.24 | 2.34 |

Regenerate or add a resolution with:

```bash
./tools/makeMesh.sh 8 20 4 meshes/mesh2
```

which runs the whole validated pipeline (geometry generator → analytic transmural seed →
`setFibreFieldArostica` Laplace solve → orthonormal triad completion).

**These parameters were chosen deliberately, and the earlier ones were wrong.** The
previous ladder (8×9×3, 8×14×4, 8×20×6) used far too few longitudinal divisions relative
to surface divisions (`l/d ≈ 1.1`), which produced two defects:

* surface cells varying by up to **4.3×** in size on one patch, and
* a **minimum** cell size that was *frozen* under refinement — 2.78 mm on both the 1056-
  and 4224-cell meshes, because `d` was held at 8. Refining added cells everywhere except
  around the apex, so the size contrast there got *worse* with refinement.

Measurement showed the apex cell size is set by `d` and the equatorial size by `l`, so the
fix is to scale **both** together at `l/d ≈ 2.5`. The ladder above now has:

* **uniformity constant at ~2.50** across all three resolutions (and better than the
  paper's tet mesh at the coarsest level), and
* a minimum cell size that **shrinks properly** with refinement: 3.72 → 2.78 → 2.22 mm.

That second point matters: it removes a mesh-side candidate for the otherwise odd
observation that refining made `J_min` *worse*.

All three: `checkMesh` **Mesh OK**, max skewness ≈ 0.69, endo/epi surface deviation from
the analytic ellipsoids **0.00000 mm**, apex exactly at `x = −0.097000`, basal rims exactly
at `+0.0264706` (endo) and `+0.0242500` (epi), fibre triads orthonormal to 2.2e-16 with
handedness `(f×s)·n = −1.000000`.

All three share the same **butterfly / O-grid apex** (no polar singularity — the apex
vertex is shared by exactly 4 cells), and the paper-faithful **sloped** basal truncation
(`cos μ = 5/17` endo, `5/20` epi).

---

## The one knob worth knowing about

In `0/D`, the `epicardium` patch has:

```
supportDisplacementMode cellReconstructedFace;   // or: ownerCell
```

This selects **where the epicardial spring/dashpot reads the displacement**:

* `ownerCell` — the adjacent cell-centre value `D_c`. First order in `h`.
* `cellReconstructedFace` — `D_c + (Cf − Cc)·grad(D)`, i.e. the reconstructed **surface**
  value. Second order.

It matters because `alpha_epi = 1e8 Pa/m` is enormous: the two measures differ by ~15× at
`t = 0.26`, which is ~28 kPa of traction against a 16.1 kPa peak chamber pressure, in the
term that supplies ~81% of the axial reaction.

`cellReconstructedFace` is the correct choice — the reference FEM implementation evaluates
the Robin condition on the finite-element **trace** `u|_Γ`, not a cell value — and it is
also more robust here. The gradient used is the production `leastSquaresS4f` scheme, which
sets `useBoundaryFaceValues = false` on every patch, so the reconstruction is built from
cell values only and never feeds the patch's own value back into itself.

---

## Current status — what works and what does not

**Works.** The case runs, converges (typically ~5–6 SNES iterations per step, zero domain
rejections), and produces a physically sensible early response.

**Does not work.** *No configuration yet completes the full `t = 1.0` cycle.* Runs stop
around **`t ≈ 0.21–0.27`**. Two distinct failure modes have been seen:

| configuration | dt | reached | failure reason |
|---|---|---|---|
| coarse HEX, `ownerCell` | 1e-3 | 0.215 | `DIVERGED_LINEAR_SOLVE` |
| coarse HEX, `cellReconstructedFace` | 1e-3 | 0.264 | `DIVERGED_LINEAR_SOLVE` |
| coarse HEX, `cellReconstructedFace` | 5e-4 | 0.2585 | `DIVERGED_LINEAR_SOLVE` |
| coarse HEX, `cellReconstructedFace` | 2e-4 | **0.262** | `DIVERGED_LINE_SEARCH` |
| coarse HEX, `cellReconstructedFace` | 1e-4 | 0.2303 | `DIVERGED_DTOL` |
| mesh3 HEX (4224 cells) | 5e-5 | 0.268 | `DIVERGED_MAX_IT` + element inversion |
| mesh3 HEX (4224 cells) | 2.5e-5 | 0.269 | `DIVERGED_MAX_IT` |

Note that the failure time is **t ≈ 0.215–0.269 in every case, across a 40× range of
timestep (1e-3 down to 2.5e-5) and a 4× range of mesh resolution, with FOUR different
PETSc failure reasons**. That is the single most important clue in this case: a solver
setting would not reproduce the same stopping time through four different failure
mechanisms. It points at a *state the solution reaches*, not at the solver.

**What happens physically.** `J = det F` collapses in the **apical** cells:

| t | J_min (coarse, dt=1e-3) | cells with J < 0.95 |
|---|---|---|
| 0.21 | 0.990214 | 0 |
| 0.23 | 0.971525 | 0 |
| 0.24 | 0.904046 | 20 |
| 0.25 | 0.837251 | 44 |
| 0.26 | 0.817992 | 56 |

At `J ≈ 0.8` with `κ = 1e6 Pa` the local pressure reaches ~220 kPa — about 14× the peak
chamber pressure. The mixed pressure `p` tracks `−κ/2(J − 1/J)` to 1% everywhere, so the
u–p coupling is sound; the compression is a genuine solution of the discrete equations, not
a formulation artefact.

**Why the apex compresses is still open.** It is *not*:

* apex mesh topology — the apex is a butterfly O-grid, tip aspect ratio 2.83 (better than
  the mesh median 3.73), no vanishing edges/areas, no coincident vertices;
* mesh resolution — refining makes `J_min` **worse**, not better (the opposite of an
  under-resolution artefact); see the cell-size table above for a likely reason;
* the fibre field — triads exact to 2e-16 including in the collapsing cells;
* timestep — halving `dt` changes the peak by 0.0006% and `J_min` by 0.1%.

The current best explanation is that it is **downstream of a displacement trajectory that
runs ~20–30 ms ahead of the benchmark**: by `t = 0.26` our `ux` is at ≈ −0.023 while the
published ensemble median is only ≈ −0.005. We arrive early at a deformation state the
participants do not reach until much later, and crush the apex when we get there.

---

## How our result compares to the benchmark

Against the ten official participant datasets (in
`../reference/cardiac_benchmark_toolkit/results/data`):

| quantity | ours | official range | |
|---|---|---|---|
| `p0 ux` first peak | 7.89e-03 | 1.019e-02 … 1.181e-02 | **−22.6%** |
| `p1 ux` first peak | 5.35e-03 | 7.495e-03 … 8.422e-03 | **−28.6%** |
| peak **time** | 0.210 | 0.230 … 0.240 | **20–30 ms early** |

The timing error is the more diagnostic of the two, and it is insensitive to mesh, base
geometry, stabilisation, timestep and the support discretisation.

## What has already been checked and excluded

An external cross-check was done against the **Simula** reference implementation
(`github.com/finsberg/cardiac_benchmark`) and, for the viscous law, **Ambit**. Written up in
`../HexMeshesTests/debug_collaboration/AROSTICA_EXTERNAL_IMPLEMENTATION_CROSSCHECK.md`.
Summary: **no material physical mismatch was found.**

* passive strain energy — **term-by-term identical** (all 8 coefficients, isochoric `I1`
  with raw `I4f/I4s/I8fs`, same logistic switch as a multiplicative prefactor)
* chamber pressure waveform — agrees to **13 significant figures** at 10 sampled times
* geometry — **exact** (their `mu_base = −acos(5/17)`, `−acos(5/20)` are our basal rims)
* viscous law — `S_v = η Ė`, confirmed against **two** independent implementations
* Robin BCs, fibre construction and handedness, units, every parameter — all match

The two remaining differences from the reference are **formulation-level, not bugs**:
they use P2 displacement-only with a `κ/4(J²−1−2lnJ)` penalty and generalized-α time
integration; we use mixed p–D cell-centred FV with BDF2.

---

## Suggested things to try

1. **The timestep has already been swept** — `../caseB_dtStudy/` holds the dt = 5e-4,
   2e-4, 1e-4 runs; combined with 1e-3 here and 5e-5/2.5e-5 on mesh3 the barrier does
   not move. Treat this as closed unless you want to re-check it.
2. **Compare the three resolutions** — `./Allrun` runs all of them; the interesting
   quantity is `J_min(t)` and where it collapses, not just the probe curves.
3. **Flip `supportDisplacementMode`** in `0/D` between `ownerCell` and
   `cellReconstructedFace` and compare — one line, isolates the Robin discretisation.
4. **Watch `J_min` rather than the residual.** Add `writeInterval 10` and look at where
   and when `det(F)` starts to fall; that is the real early-warning signal.
5. **Compare against the official data** with
   `../HexMeshesTests/debug_collaboration/ensemble/full_compare.py <probeDir> <label>`,
   which plots all six curves against the ten-participant envelope.

## Known rough edges

* **Restart does not reproduce the state.** Starting from a written time gives an initial
  residual ~400× larger than the continuous run. Prefer clean runs from `t = 0`.
* **`postProcess` cannot re-read `0/D` directly** because of the case-local pressure BC;
  that is why `extractProbes.sh` copies fields into `probeCase/` and neutralises the
  boundary types first (values are preserved exactly).
* The linear solver is expensive at `dt = 1e-3` (~800 KSP iterations/step). This is a
  preconditioner limitation, deliberately left alone for now.
