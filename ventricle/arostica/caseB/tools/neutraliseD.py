#!/usr/bin/env python3
"""Make a written D file readable by a plain utility.

The case-local arosticaActivationPressureTraction BC writes 'alphaMin [dims] v'
but reads a named dimensionedScalar, so 0.02/D cannot be re-read (pre-existing
defect, out of scope). This rewrites each boundaryField entry to 'calculated'
keeping its written 'value'. gradD is unaffected: leastSquaresS4fGrad has
useBoundaryFaceValues = false, so boundary D values do not enter the gradient.
Verified by the F gate in the audit utility.
"""
import re, sys, pathlib
p = pathlib.Path(sys.argv[1])
s = p.read_text()
i = s.index("boundaryField")
head, body = s[:i], s[i:]

def repl(m):
    name, blk = m.group(1), m.group(2)
    v = re.search(r"(value\s+(?:nonuniform\s+List<vector>\s*\n\d+\s*\n\(.*?\n\)|uniform\s+\([^)]*\))\s*;)", blk, re.S)
    val = v.group(1) if v else "value uniform (0 0 0);"
    return f"\n    {name}\n    {{\n        type            calculated;\n        {val}\n    }}"

body = re.sub(r"\n    (\w+)\s*\n    \{(.*?)\n    \}", repl, body, flags=re.S)
p.write_text(head + body)
print(f"rewrote {p} boundaryField -> calculated")
