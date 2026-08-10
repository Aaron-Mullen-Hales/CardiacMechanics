# Clean Gültekin artery benchmark

This maintained case implements the Gültekin–Dal–Holzapfel extension–inflation–torsion artery benchmark from Karabelas et al. (2019), using the mixed finite-volume formulation: solid model `nonLinearGeometryTotalLagrangianTotalDisplacement` (PETSc SNES, `solvePressure true`) with material law `GultekinTwoFibreElastic`, two direct ±40° fibre families, finite bulk modulus, and Rhie–Chow pressure stabilisation. This case is built and run against a solids4foam checkout, where the historical benchmark-specific solid model `nonLinGeomTotalLagTotalDispGultekinSolid` has been merged into the general-purpose `nonLinGeomTotalLagTotalDispSolid` class (same mixed-pressure options: `scaleMixedPetScFields`, `pressureUnknownScale`, `pressureScaleFactor`, `stabilisation`). This mirrors the solid model used by the Land ventricle cases one directory up (`ventricle/land/problem1` etc.).

Nothing in this case (`scripts/`, `base/`, `Allrun`, `Allclean`) hardcodes a path — it only depends on `solids4Foam`/`checkMesh`/`blockMesh` being resolvable via `PATH`, `DYLD_LIBRARY_PATH`/`LD_LIBRARY_PATH`, and `FOAM_USER_LIBBIN`/`FOAM_LIBBIN`/`FOAM_EXT_LIBBIN`.

The pressure is the prescribed conversion of 500 mmHg: 66,661.184 Pa. The tube has inner radius 8 mm, outer radius 10 mm, and height 10 mm. The final top motion is 2 mm axial extension and 60° rotation about `z`.

## Building solids4foam

This case requires the `feature/active-passive-stress-split` branch of solids4foam, from `git@github.com:Aaron-Mullen-Hales/solids4foam_dev.git`, which is where `GultekinTwoFibreElastic` and the merged mixed-pressure `nonLinGeomTotalLagTotalDispSolid` live. Verified against commit `bd6a2311` (2026-08-10) on that branch; `master` does not have these changes. Clone/checkout that repo and branch anywhere you like:

```bash
source /path/to/OpenFOAM-v2312/etc/bashrc
git clone -b feature/active-passive-stress-split git@github.com:Aaron-Mullen-Hales/solids4foam_dev.git
cd solids4foam_dev
export FOAM_MODULE_APPBIN="$PWD/platforms/$WM_OPTIONS/bin"
export FOAM_MODULE_LIBBIN="$PWD/platforms/$WM_OPTIONS/lib"
./Allwmake -j -s
```

`FOAM_MODULE_APPBIN`/`FOAM_MODULE_LIBBIN` keep the build local to that checkout's own `platforms/` directory instead of the shared `FOAM_USER_APPBIN`/`FOAM_USER_LIBBIN` location, which matters if you have more than one solids4foam checkout on the same machine. Before `./Allrun`, put that `bin`/`lib` pair on `PATH`/`DYLD_LIBRARY_PATH` (and `FOAM_USER_LIBBIN`, which is what the runtime gate searches):

```bash
S4F=/path/to/solids4foam_dev
export PATH="$S4F/platforms/$WM_OPTIONS/bin:$PATH"
export DYLD_LIBRARY_PATH="$S4F/platforms/$WM_OPTIONS/lib:$DYLD_LIBRARY_PATH"
export FOAM_USER_LIBBIN="$S4F/platforms/$WM_OPTIONS/lib"
```

## Meshes

The seven meshes follow Table 1 of Karabelas et al., *An accurate, robust, and efficient finite element framework for anisotropic, nearly and fully incompressible elasticity* (the mesh-convergence study for this same artery benchmark): a uniform geometric refinement radial × circumferential × axial = `4s × 24s × 10s` for level `s = 1..7`, giving hexahedron counts 960, 7680, 25920, 61440, 120000, 207360, and 329280 — an exact match to the paper's element and node counts.

| mesh | level | radial | circumferential | axial | cells |
| --- | --- | --- | --- | --- | --- |
| mesh1 | 1 | 4 | 24 | 10 | 960 |
| mesh2 | 2 | 8 | 48 | 20 | 7680 |
| mesh3 | 3 | 12 | 72 | 30 | 25920 |
| mesh4 | 4 | 16 | 96 | 40 | 61440 |
| mesh5 | 5 | 20 | 120 | 50 | 120000 |
| mesh6 | 6 | 24 | 144 | 60 | 207360 |
| mesh7 | 7 | 28 | 168 | 70 | 329280 |

`mesh1` is the quick regression mesh and matches the accepted `fullRigid_0` topology. `mesh1 → mesh7` is a uniform full-geometry refinement sequence (unlike the previous ad-hoc `mesh1`-`mesh4` scheme, where only `mesh2 → mesh3 → mesh4` refined radially). `mesh7` is the finest paper level, used there for the published stress-distribution figure; at 329,280 cells it is expensive to solve on a single machine (see cost note below).

Generated cases are exactly `runs/mesh1` through `runs/mesh7`.

## Commands

```text
./Allrun
./Allrun --mesh mesh2
./Allrun --meshes mesh1,mesh2,mesh3,mesh4,mesh5,mesh6,mesh7
./Allrun --generate-only
./Allrun --validate-only
./Allrun --run-only
./Allrun --postprocess-only
./Allrun --plot-only
./Allrun --runtime-dir /path/to/staged/runtime
```

`./Allrun` selects `mesh2`. Unknown mesh names are rejected. Generation and validation are separate from solving; the runtime gate checks the active or staged executable/library provenance before a solve begins. A staged runtime must contain `solids4Foam` and `libsolids4FoamModels.dylib` with the hashes recorded in `referenceData/runtime_manifest.json`; neither is copied into this tree. Against an actively sourced (non-staged) runtime, the gate only requires the required model/solid-model markers to be present in the library, not a hash match, so a freshly built solids4foam runtime is accepted without updating the archived hashes.

Mesh cost grows with cell count: `mesh1` (960 cells) solves in well under a minute on this machine. `mesh5`-`mesh7` (120k-330k cells) are substantially more expensive — plan for a long single-machine run, or run them unattended/in the background, before requesting all seven levels at once.

## Fields and stress convention

The default postprocessor writes exactly these final-time ParaView fields:

`sigmaRR_referenceBasis_kPa`, `sigmaTT_referenceBasis_kPa`, and `sigmaZZ_referenceBasis_kPa`.

These three fields are chosen specifically to match the quantities the paper reports for this benchmark: the radial `σrr`, circumferential `σθθ`, and axial/longitudinal `σzz` components of the Cauchy stress tensor, resolved in the cylindrical `(r, θ, z)` basis. Both Karabelas et al. (Table 1 / Figure 2, the mesh-convergence paper this case's 7 mesh levels come from) and the original Gültekin–Dal–Holzapfel (2019) benchmark report exactly these three stress components — `σrr`, `σθθ`, `σzz` — as the comparison quantities for this problem, so writing them under matching names is what makes a direct numeric/visual comparison against the published figures possible. `referenceData/published_comparison_data.csv` carries the Gültekin–Dal–Holzapfel Figure 6 colour-bar ranges for these same three fields (in kPa) as descriptive comparison data, not pass/fail tolerances.

They are final spatial Cauchy-stress components resolved in the cylindrical basis attached to the undeformed/reference position (hence `_referenceBasis`), so the `r`/`θ`/`z` directions used for the projection do not rotate with the deformation. They are kPa visualization copies. The default output does not restore current/deformed-basis fields. Use `--diagnostics` with `--postprocess-only` only when those extra diagnostic fields are explicitly wanted.

The validation report checks the reference-basis identities

```text
sigmaRR_referenceBasis + sigmaTT_referenceBasis = sigmaXX + sigmaYY
sigmaRR_referenceBasis + sigmaTT_referenceBasis + sigmaZZ_referenceBasis = trace(sigma)
```

Open `runs/mesh2/case.foam` in ParaView after a completed run to view the one comparable final-time session.

## Provenance

`referenceData/runtime_manifest.json` records the accepted executable and library hashes, model/law names, and archived provenance labels. Those hashes are from the historical archived runtime and are not required for an ordinary run against a freshly built solids4foam. On another machine, install/source a compatible solids4Foam environment (the `feature/active-passive-stress-split` branch, see above) containing `nonLinearGeometryTotalLagrangianTotalDisplacement` and `GultekinTwoFibreElastic`, then source it before `./Allrun`. The custom executable and library are intentionally not included in this repository. `--runtime-dir` is reserved for an exact archived replay and requires the recorded hashes.

`referenceData/accepted_results.csv` lists the paper's mesh sizes for all seven levels; every row is marked pending until replayed against this maintained topology and that solids4foam build. Published colour-bar ranges are descriptive comparison data in `referenceData/published_comparison_data.csv`, not pass/fail tolerances.

## Known open item: mesh1 smoke test against the new build

A `mesh1` end-to-end run (`./Allrun --mesh mesh1`) against the freshly built solids4foam runtime converges numerically (SNES reaches `CONVERGED_FNORM_RELATIVE` every step, momentum equation converges in all time-steps) and passes every workflow gate. However the Jacobian field `GTF_J` at full load spans roughly `[0.91, 1.70]` (mean ≈ 1.32), a materially larger deviation from `J = 1` than the historical accepted `mesh1` replay (`min_J = 0.879`, `max_J = 1.055`) despite the same material parameters (`K/mu = 500`). The mean stress fields also differ from the historical accepted values by 30-40%, not by a uniform scale factor. `result_gate` in `scripts/validate_case.py` only checks that `J` stays positive, so this passes validation without being physically equivalent to the old accepted case. This is worth checking against the `nonLinGeomTotalLagTotalDispSolid` mixed-pressure formulation (`pressureUnknownScale`, `pressureScaleFactor`, `pressureScaleByTwoMu`, `stabilisation`) before treating any new mesh1-mesh7 results as validated; it may simply need re-tuned stabilisation/scale settings for the merged solid model rather than the values carried over from the retired `nonLinGeomTotalLagTotalDispGultekinSolid`.

Historical research directories and their binaries remain outside this maintained tree and are not modified by `Allrun` or `Allclean`.
