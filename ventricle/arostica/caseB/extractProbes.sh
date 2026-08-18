#!/bin/bash
#
# Extract u_h at the two benchmark points
#     p0 = (0.025, 0.03, 0)      p1 = (0, 0.03, 0)
# using OpenFOAM's own cellPoint interpolation (the benchmark requires an
# interpolation algorithm; neither point is a mesh point).
#
# The written D files are copied into ./probeCase with their boundary
# conditions neutralised first, because the case-local
# arosticaActivationPressureTraction cannot be re-read from a written file.
# Boundary VALUES are preserved exactly, so the probe result is unaffected.
#
# Usage:  ./extractProbes.sh        (run from inside a case or a mesh run dir)

# Source the OpenFOAM environment if it is not already active, BEFORE errexit:
# the bashrc returns non-zero from some of its helper functions.
if [ -z "${WM_PROJECT_DIR:-}" ]; then
    . /Volumes/OpenFoam/OpenFOAM-v2312/etc/bashrc
fi
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

export DYLD_LIBRARY_PATH="$HERE/lib:$FOAM_USER_LIBBIN:$FOAM_LIBBIN:${DYLD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$HERE/lib:$FOAM_USER_LIBBIN:$FOAM_LIBBIN:${LD_LIBRARY_PATH:-}"

rm -rf probeCase
mkdir -p probeCase
cp -R constant system probeCase/

n=0
for d in [0-9]*; do
    [ -d "$d" ] || continue
    [ -f "$d/D" ] || continue
    mkdir -p "probeCase/$d"
    cp "$d/D" "probeCase/$d/D"
    python3 tools/neutraliseD.py "probeCase/$d/D" > /dev/null
    n=$((n+1))
done

if [ "$n" -eq 0 ]; then
    echo "No written time directories with a D field found in $HERE"
    echo "Run the solver first."
    exit 1
fi

echo "Neutralised $n time directories; running postProcess..."
cd probeCase
postProcess -func probesDict > log.postProcess 2>&1 || true
cd "$HERE"

F=$(ls probeCase/postProcessing/probesDict/*/D 2>/dev/null | head -1)
if [ -n "$F" ]; then
    echo "OK  ->  $F   ($(grep -vc '^#' "$F") samples)"
else
    echo "FAILED - see probeCase/log.postProcess"
    exit 1
fi
