#!/bin/bash
#
# Generate one Aróstica monoventricle HEX mesh + its fibre fields.
#
#   ./tools/makeMesh.sh <surfaceDiv> <longitudinalDiv> <transmuralLayers> <outDir>
#
# e.g.  ./tools/makeMesh.sh 8  9 3 meshes/mesh1     ->  1056 cells
#       ./tools/makeMesh.sh 8 14 4 meshes/mesh2     ->  2048 cells
#       ./tools/makeMesh.sh 8 20 6 meshes/mesh3     ->  4224 cells
#
# Cell count = (d^2 + 4*d*l) * m.
#
# Writes <outDir>/polyMesh and <outDir>/0/{f0,s0,n0,f0f,s0f,n0f,t}.
# Fibres go through the full validated pipeline:
#   analytic transmural seed -> setFibreFieldArostica FVM Laplace solve
#   -> orthonormal triad completion (alpha = -60 + 120 t).
d=$1; l=$2; m=$3; OUT=$4
HERE="$(cd "$(dirname "$0")/.." && pwd)"
GEN="$HERE/tools/meshGenerator"   # case-local generator (radially graded core)
STAGE=$(mktemp -d)

# NOTE: source the OpenFOAM bashrc BEFORE enabling `set -e`.  The bashrc
# returns non-zero from some of its helper functions, which aborts the script
# if errexit is already active.
if [ -z "${WM_PROJECT_DIR:-}" ]; then . /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc; fi
set -e
export DYLD_LIBRARY_PATH="$HERE/lib:$FOAM_USER_LIBBIN:$FOAM_LIBBIN:${DYLD_LIBRARY_PATH:-}"

cp -R "$GEN" "$STAGE/scripts"
cd "$STAGE"
python3 scripts/makeHexMesh.py --case monoventricle \
    --surface-divisions "$d" --longitudinal-divisions "$l" --transmural-layers "$m"

python3 - "$STAGE" <<'PY'
import sys, pathlib
sys.path.insert(0, sys.argv[1] + "/scripts")
import fibres
case = pathlib.Path(sys.argv[1]) / "cases" / "monoventricle"
meta = fibres.load_mesh_data(case)
vals = [fibres.mono_t(p) for p in meta["cell_centres"]]
fibres.write_scalar(case/"0"/"t", "t", vals, ("endocardium","epicardium","base"),
    {"endocardium":("fixedValue","0"),"epicardium":("fixedValue","1"),"base":("zeroGradient",None)})
fibres.write_solids4foam_fields(case, ("endocardium","epicardium","base"))
PY

# setFibreFieldArostica needs a system/ directory (controlDict, fvSchemes,
# fvSolution, setFibreFieldDict); the fibre "prepare" step only writes 0/.
cp -R "$HERE/system" "$STAGE/cases/monoventricle/system"

cd "$STAGE/cases/monoventricle"
setFibreFieldArostica > log.setFibreField 2>&1
python3 - "$STAGE" <<'PY'
import sys
sys.path.insert(0, sys.argv[1] + "/scripts")
import fibres
fibres.complete_mono(sys.argv[1] + "/cases")
PY

mkdir -p "$HERE/$OUT/0"
cp -R "$STAGE/cases/monoventricle/constant/polyMesh" "$HERE/$OUT/"
for f in f0 s0 n0 f0f s0f n0f t; do cp "$STAGE/cases/monoventricle/0/$f" "$HERE/$OUT/0/"; done
rm -rf "$STAGE"
echo "$OUT: $(( (d*d + 4*d*l) * m )) cells"
