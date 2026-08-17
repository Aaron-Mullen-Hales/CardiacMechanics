#!/bin/bash
#
# Extract u_h at the two benchmark points
#     p0 = (0.025, 0.03, 0)      p1 = (0, 0.03, 0)
# using OpenFOAM's own cellPoint interpolation (the benchmark requires an
# interpolation algorithm; neither point is a mesh point).
#
# The written D files are copied into ./probeCase with the boundary conditions
# neutralised first, because the case-local arosticaActivationPressureTraction
# cannot currently be re-read from a written file.  Boundary VALUES are
# preserved exactly, so the probe result is unaffected.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
rm -rf "$HERE/probeCase"; mkdir -p "$HERE/probeCase"
cp -R "$HERE/constant" "$HERE/system" "$HERE/probeCase/"
for d in "$HERE"/[0-9]*; do
    b=$(basename "$d"); [ -f "$d/D" ] || continue
    mkdir -p "$HERE/probeCase/$b"; cp "$d/D" "$HERE/probeCase/$b/D"
    python3 "$HERE/tools/neutraliseD.py" "$HERE/probeCase/$b/D" > /dev/null
done
cd "$HERE/probeCase" && postProcess -func probesDict > /dev/null 2>&1
echo "probe data: probeCase/postProcessing/probesDict/0/D"
