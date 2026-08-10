#!/usr/bin/env python3
"""Write directly evaluated ±40-degree reference fibre fields."""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
from common import dot, file_hash, load_mesh, mesh_hash, write_surface_vector_field, write_vol_vector_field

ROOT=Path(__file__).resolve().parents[1]

def direction(position, alpha, sign):
    x,y,_=position; radius=math.hypot(x,y)
    if radius<=1e-14: raise ValueError("fibre point lies on cylinder axis")
    et=(-y/radius,x/radius,0.0)
    return (math.cos(alpha)*et[0],math.cos(alpha)*et[1],sign*math.sin(alpha))

def generate(mesh_name):
    case=ROOT/"runs"/mesh_name; mesh=load_mesh(case); alpha=math.radians(40.0)
    f0=[direction(c,alpha,1) for c in mesh["cell_centres"]]; f1=[direction(c,alpha,-1) for c in mesh["cell_centres"]]
    f0f=[direction(c,alpha,1) for c in mesh["face_centres"]]; f1f=[direction(c,alpha,-1) for c in mesh["face_centres"]]
    b0={}; b1={}
    for patch,data in mesh["boundary"].items():
        indices=range(data["startFace"],data["startFace"]+data["nFaces"])
        b0[patch]=[f0f[i] for i in indices]; b1[patch]=[f1f[i] for i in indices]
    write_vol_vector_field(case/"0/f0","f0",f0,b0); write_vol_vector_field(case/"0/f1","f1",f1,b1)
    write_surface_vector_field(case/"0/f0f","f0f",f0f[:mesh["n_internal_faces"]],b0)
    write_surface_vector_field(case/"0/f1f","f1f",f1f[:mesh["n_internal_faces"]],b1)
    errors=[abs(math.sqrt(dot(v,v))-1) for v in f0+f1+f0f+f1f]
    radial=[]; angles=[]
    for c,v,sign in [(c,v,1) for c,v in zip(mesh["cell_centres"],f0)]+[(c,v,-1) for c,v in zip(mesh["cell_centres"],f1)]:
        r=math.hypot(c[0],c[1]); er=(c[0]/r,c[1]/r,0)
        radial.append(abs(dot(v,er))); et=(-c[1]/r,c[0]/r,0)
        angles.append(math.degrees(math.atan2(sign*v[2],dot(v,et))))
    result={"mesh":mesh_name,"mesh_hash":mesh_hash(case),"fibre_angle_degrees":40.0,
            "maximum_unit_error":max(errors),"maximum_radial_component":max(radial),
            "angle_min_degrees":min(angles),"angle_max_degrees":max(angles),
            "field_sha256":{n:file_hash(case/"0"/n) for n in ("f0","f1","f0f","f1f")}}
    (case/"fibreMetadata.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--meshes",default="mesh1,mesh2,mesh3,mesh4,mesh5,mesh6,mesh7")
    print(json.dumps([generate(x.strip()) for x in ap.parse_args().meshes.split(",") if x.strip()],indent=2))
