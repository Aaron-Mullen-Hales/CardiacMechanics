#!/usr/bin/env python3
"""Generate the four maintained annular meshes with blockMesh."""
from __future__ import annotations
import argparse, json, math, shutil, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object blockMeshDict;
}
"""

def point(radius, theta, z):
    return (radius*math.cos(theta), radius*math.sin(theta), z)

def fmt(v):
    return "(" + " ".join(f"{x:.16g}" for x in v) + ")"

def dictionary(geometry, counts):
    nr, nt, nz = counts["radial"], counts["circumferential"], counts["axial"]
    if nt % 8: raise ValueError("circumferential count must be divisible by 8")
    ri, ro, height = geometry["inner_radius"], geometry["outer_radius"], geometry["height"]
    vertices = [point(r, i*math.pi/4, z) for z in (0, height) for r in (ri, ro) for i in range(8)]
    def vid(zi, ri_, ti): return zi*16 + ri_*8 + ti%8
    blocks, edges, patches = [], [], {"innerWall": [], "outerWall": [], "bottom": [], "top": []}
    for i in range(8):
        j=(i+1)%8
        v=[vid(0,0,i),vid(0,1,i),vid(0,1,j),vid(0,0,j),vid(1,0,i),vid(1,1,i),vid(1,1,j),vid(1,0,j)]
        blocks.append(f"    hex ({' '.join(map(str,v))}) ({nr} {nt//8} {nz}) simpleGrading (1 1 1)")
        for zi,z in enumerate((0,height)):
            edges += [f"    arc {vid(zi,0,i)} {vid(zi,0,j)} {fmt(point(ri, (i+.5)*math.pi/4,z))}",
                      f"    arc {vid(zi,1,i)} {vid(zi,1,j)} {fmt(point(ro, (i+.5)*math.pi/4,z))}"]
        patches["innerWall"].append(f"            ({v[0]} {v[4]} {v[7]} {v[3]})")
        patches["outerWall"].append(f"            ({v[1]} {v[2]} {v[6]} {v[5]})")
        patches["bottom"].append(f"            ({v[0]} {v[3]} {v[2]} {v[1]})")
        patches["top"].append(f"            ({v[4]} {v[5]} {v[6]} {v[7]})")
    def patch(name):
        return f"    {name}\n    {{\n        type patch;\n        faces\n        (\n" + "\n".join(patches[name]) + "\n        );\n    }\n"
    return HEADER + "scale 1;\n\nvertices\n(\n" + "\n".join("    "+fmt(v) for v in vertices) + \
        "\n);\n\nblocks\n(\n" + "\n".join(blocks) + "\n);\n\nedges\n(\n" + \
        "\n".join(edges) + "\n);\n\nboundary\n(\n" + "".join(patch(name) for name in patches) + \
        ");\n\nmergePatchPairs ();\n"

def generate(mesh_name, force=False):
    definitions=json.loads((ROOT/"meshes/mesh_definitions.json").read_text())
    if mesh_name not in definitions["meshes"]: raise ValueError(f"Unknown mesh name: {mesh_name}")
    case=ROOT/"runs"/mesh_name
    if case.exists() and force: shutil.rmtree(case)
    (case/"system").mkdir(parents=True, exist_ok=True)
    (case/"constant").mkdir(exist_ok=True)
    (case/"system/blockMeshDict").write_text(dictionary(definitions["geometry_m"],definitions["meshes"][mesh_name]),encoding="utf-8")
    template=(ROOT/"base/system/controlDict").read_text(encoding="utf-8")
    template="/*--------------------------------*- C++ -*----------------------------------*\\\n| solids4foam: clean artery benchmark                                          |\n\\*---------------------------------------------------------------------------*/\n"+template.replace("}\napplication", "}\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\napplication")
    (case/"system/controlDict").write_text(template.replace("application solids4Foam;","application blockMesh;"),encoding="utf-8")
    for name in ("fvSchemes", "fvSolution"):
        shutil.copy2(ROOT/"base/system"/name, case/"system"/name)
    if not (case/"constant/polyMesh").is_dir():
        if not shutil.which("blockMesh"):
            raise RuntimeError("blockMesh is not available; source the OpenFOAM v2312 environment")
        result=subprocess.run(["blockMesh","-case",str(case)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (case/"blockMesh.log").write_text(result.stdout,encoding="utf-8")
        if result.returncode: raise RuntimeError(f"blockMesh failed for {mesh_name}; see {case/'blockMesh.log'}")
    if not shutil.which("checkMesh"):
        raise RuntimeError("checkMesh is not available; source the OpenFOAM v2312 environment")
    quality=subprocess.run(["checkMesh","-case",str(case),"-allGeometry","-allTopology"],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (case/"checkMesh.log").write_text(quality.stdout,encoding="utf-8")
    (case/"meshQuality.json").write_text(json.dumps({"return_code":quality.returncode,"passed":quality.returncode==0,"command":"checkMesh -allGeometry -allTopology"},indent=2)+"\n",encoding="utf-8")
    if quality.returncode: raise RuntimeError(f"checkMesh failed for {mesh_name}; see {case/'checkMesh.log'}")
    return case

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--meshes",default="mesh1,mesh2,mesh3,mesh4,mesh5,mesh6,mesh7"); ap.add_argument("--force",action="store_true")
    args=ap.parse_args(); definitions=json.loads((ROOT/"meshes/mesh_definitions.json").read_text())
    selected=[x.strip() for x in args.meshes.split(",") if x.strip()]
    unknown=sorted(set(selected)-set(definitions["meshes"]))
    if unknown: ap.error("Unknown mesh names: "+",".join(unknown))
    for name in selected: generate(name,args.force)
    print(json.dumps({name: definitions["meshes"][name] for name in selected},indent=2))
if __name__ == "__main__": main()
