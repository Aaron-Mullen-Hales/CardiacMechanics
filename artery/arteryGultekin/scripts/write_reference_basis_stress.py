#!/usr/bin/env python3
"""Write final spatial Cauchy stress in the reference-attached cylindrical basis."""
from __future__ import annotations
import argparse, json, math, re
from pathlib import Path
from common import load_mesh, foam_header

ROOT=Path(__file__).resolve().parents[1]

def values(path):
    text=re.sub(r"//[^\n]*|/\*.*?\*/","",path.read_text(errors="replace"),flags=re.S)
    match=re.search(r"internalField\s+nonuniform\s+List<\w+>\s+(\d+)\s*\((.*?)\)\s*;",text,re.S)
    if not match:
        m=re.search(r"internalField\s+uniform\s+([^;]+);",text); token=m.group(1).strip()
        return [tuple(map(float,token[1:-1].split())) if token.startswith("(") else float(token)]
    tuples=re.findall(r"\(([^()]*)\)",match.group(2))
    return [tuple(map(float,item.split())) for item in tuples] if tuples else [float(x) for x in match.group(2).split()]

def patch_values(path,patch):
    text=re.sub(r"//[^\n]*|/\*.*?\*/","",path.read_text(errors="replace"),flags=re.S); m=re.search(rf"\b{patch}\s*\{{(.*?)\}}",text,re.S)
    if not m:return None
    b=m.group(1); u=re.search(r"value\s+uniform\s+([^;]+);",b)
    if u:
        t=u.group(1).strip(); return [tuple(map(float,t[1:-1].split())) if t.startswith("(") else float(t)]
    n=re.search(r"value\s+nonuniform\s+List<\w+>\s+\d+\s*\((.*?)\)\s*;",b,re.S)
    if not n:return None
    tuples=re.findall(r"\(([^()]*)\)",n.group(1)); return [tuple(map(float,item.split())) for item in tuples]

def tensor(v):
    if len(v)==6:
        xx,xy,xz,yy,yz,zz=v; return ((xx,xy,xz),(xy,yy,yz),(xz,yz,zz))
    return (tuple(v[0:3]),tuple(v[3:6]),tuple(v[6:9]))
def comp(a,x,y): return sum(x[i]*a[i][j]*y[j] for i in range(3) for j in range(3))

def write_scalar(path,name,internal,patches):
    text=foam_header("volScalarField",path.parent.name,name)+"dimensions [0 0 0 0 0 0 0];\ninternalField nonuniform List<scalar>\n%d\n(\n"%len(internal)
    text+="".join(f"    {x:.16g}\n" for x in internal)+");\n\nboundaryField\n{\n"
    for patch,vals in patches.items():
        text+=f"    {patch}\n    {{\n        type calculated;\n        value nonuniform List<scalar>\n        {len(vals)}\n        (\n"+"".join(f"            {x:.16g}\n" for x in vals)+"        );\n    }\n"
    path.write_text(text+"}\n",encoding="utf-8")

def process(mesh_name,diagnostics=False):
    case=ROOT/"runs"/mesh_name; mesh=load_mesh(case); times=sorted([p for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?",p.name)],key=lambda p:float(p.name))
    if not times: raise RuntimeError(f"{mesh_name}: no solver time directories")
    final=times[-1]; sigma_path=final/"sigma"
    if not sigma_path.is_file(): raise RuntimeError(f"{mesh_name}: final sigma field missing")
    sigmas=[tensor(x) for x in values(sigma_path)]; positions=mesh["cell_centres"]
    if len(sigmas)==1: sigmas=sigmas*len(positions)
    names=("sigmaRR_referenceBasis_kPa","sigmaTT_referenceBasis_kPa","sigmaZZ_referenceBasis_kPa")
    data={name:[] for name in names}; traces=[]
    for pos,s in zip(positions,sigmas):
        r=math.hypot(pos[0],pos[1]); er=(pos[0]/r,pos[1]/r,0); et=(-pos[1]/r,pos[0]/r,0); ez=(0,0,1)
        rr,tt,zz=[comp(s,v,v)/1000 for v in (er,et,ez)]
        data[names[0]].append(rr); data[names[1]].append(tt); data[names[2]].append(zz)
        traces.append((rr+tt-(s[0][0]+s[1][1])/1000,rr+tt+zz-(s[0][0]+s[1][1]+s[2][2])/1000))
    patches={}
    for patch,d in mesh["boundary"].items():
        sv=patch_values(sigma_path,patch); start=d["startFace"]
        if sv is None:
            vals=[sigmas[mesh["owner"][face]] for face in range(start,start+d["nFaces"])]
        else:
            vals=[tensor(x) for x in sv]
            if len(vals)==1: vals=vals*d["nFaces"]
        pvals={name:[] for name in names}
        for face,s in zip(range(start,start+d["nFaces"]),vals):
            pos=mesh["face_centres"][face]; r=math.hypot(pos[0],pos[1]); er=(pos[0]/r,pos[1]/r,0); et=(-pos[1]/r,pos[0]/r,0); ez=(0,0,1)
            for name,v in zip(names,(er,et,ez)): pvals[name].append(comp(s,v,v)/1000)
        for name in names: patches.setdefault(name,{})[patch]=pvals[name]
    for name in names: write_scalar(final/name,name,data[name],patches[name])
    result={"mesh":mesh_name,"final_time":final.name,"fields_written":list(names),"reference_basis":"er/etheta/ez from undeformed cell or face centres","maximum_inplane_trace_identity_error_kPa":max(abs(x[0]) for x in traces),"maximum_full_trace_identity_error_kPa":max(abs(x[1]) for x in traces)}
    (case/"referenceBasisStress.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--meshes",default="mesh2"); ap.add_argument("--diagnostics",action="store_true"); args=ap.parse_args(); print(json.dumps([process(x.strip(),args.diagnostics) for x in args.meshes.split(",") if x.strip()],indent=2))
