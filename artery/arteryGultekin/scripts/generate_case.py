#!/usr/bin/env python3
"""Instantiate a named maintained case from the accepted mechanics template."""
from __future__ import annotations
import argparse, json, shutil
from pathlib import Path
from common import file_hash, mesh_hash
from generate_mesh import generate as generate_mesh
from generate_fibres import generate as generate_fibres

ROOT=Path(__file__).resolve().parents[1]
FINAL_LOADS={"pressure_Pa":66661.184,"extension_m":0.002,"rotation_degrees":60.0}

def generate(mesh_name,force=False):
    case=generate_mesh(mesh_name,force=force)
    base=ROOT/"base"
    for source in (base/"0",base/"constant",base/"system"):
        destination=case/source.relative_to(base)
        destination.mkdir(parents=True,exist_ok=True)
        for item in source.iterdir():
            if item.is_file(): shutil.copy2(item,destination/item.name)
    shutil.copy2(base/"petscOptions.split",case/"petscOptions.split")
    load=case/"constant/load"; load.mkdir(parents=True,exist_ok=True)
    (load/"timeVsPressure").write_text("(\n    (0 0)\n    (1 66661.184)\n)\n",encoding="utf-8")
    (load/"timeVsTranslation").write_text("(\n    (0 (0 0 0))\n    (1 (0 0 0.002))\n)\n",encoding="utf-8")
    (load/"timeVsRotationDegrees").write_text("(\n    (0 0)\n    (1 60)\n)\n",encoding="utf-8")
    generate_fibres(mesh_name)
    (case/"case.foam").write_text("",encoding="utf-8")
    definitions=json.loads((ROOT/"meshes/mesh_definitions.json").read_text())
    counts=definitions["meshes"][mesh_name]
    manifest={
        "schema_version":1,"benchmark":"Gultekin-Dal-Holzapfel extension-inflation-torsion artery",
        "mesh":mesh_name,"mesh_counts":{k:counts[k] for k in ("radial","circumferential","axial")},
        "geometry_m":definitions["geometry_m"],"fibre_angle_degrees":40.0,
        "final_loads":FINAL_LOADS,
        "material":{"rho_kg_m3":1000.0,"mu_Pa":10000.0,"k1_Pa":500000.0,"k2":2.0,"bulk_modulus_Pa":5000000.0,"implicit_shear_modulus_Pa":2010000.0,"impKcoeff":1.0,"runtime_type":"GultekinTwoFibreElastic","use_second_fibre_family":True,"fibres_tension_only":False,"clip_exponent":False,"exponent_limit":650.0,"constitutive_stress_contract":"passive_cauchy_stress_only"},
        "solid_model":"nonLinearGeometryTotalLagrangianTotalDisplacement",
        "stabilisation":{"momentum_type":"diffStencilLaplacian","momentum_scale":1.0,"pressure_type":"RhieChow","pressure_scale":10.0,"pressure_jacobian_scale":10.0},
        "numerics":{"delta_t":0.005,"end_time":1.0,"dry_run_end_time":0.02,"snes_rtol":1e-6,"snes_stol":1e-12,"ksp_rtol":1e-3,"preconditioner":"fieldsplit/amg","accepted_snes_reasons":["CONVERGED_FNORM_RELATIVE","CONVERGED_FNORM_ABS"]},
        "boundary_conditions":{"innerWall":"follower internal pressure; useUndeformedArea false","outerWall":"traction free","bottom":"fixed displacement","top":"fixed rotation about z plus z extension"},
        "mesh_hash":mesh_hash(case)
    }
    (case/"caseManifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    required=["0/D","0/p","0/pointD","0/f0","0/f1","0/f0f","0/f1f","constant/physicsProperties","constant/solidProperties","constant/mechanicalProperties","constant/dynamicMeshDict","constant/g","system/controlDict","system/controlDict.dryRun","system/fvSchemes","system/fvSolution","system/decomposeParDict","petscOptions.split","system/blockMeshDict","constant/load/timeVsPressure","constant/load/timeVsTranslation","constant/load/timeVsRotationDegrees"]
    lock={rel:file_hash(case/rel) for rel in required}
    (case/"inputLock.json").write_text(json.dumps(lock,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return manifest

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--meshes",default="mesh1,mesh2,mesh3,mesh4,mesh5,mesh6,mesh7"); ap.add_argument("--force",action="store_true")
    print(json.dumps([generate(x.strip(),ap.parse_args().force) for x in ap.parse_args().meshes.split(",") if x.strip()],indent=2))
