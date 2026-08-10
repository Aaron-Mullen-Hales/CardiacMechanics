#!/usr/bin/env python3
"""Static, geometric, runtime-provenance, and result acceptance gates."""
from __future__ import annotations
import argparse, json, math, os, re, shutil, subprocess
from pathlib import Path
from common import file_hash, load_mesh, mesh_hash

ROOT=Path(__file__).resolve().parents[1]
RUNTIME=json.loads((ROOT/"referenceData/runtime_manifest.json").read_text())

def runtime_audit(runtime_dir=None):
    if runtime_dir:
        root=Path(runtime_dir).expanduser().resolve()
        exes=list(root.rglob("solids4Foam")); libs=list(root.rglob("libsolids4FoamModels.*"))
        exe=exes[0] if exes else None; lib=libs[0] if libs else None
        source="staged runtime"
    else:
        found=shutil.which("solids4Foam"); exe=Path(found).resolve() if found else None
        lib=None
        search_dirs=[]
        if exe: search_dirs.append(exe.parent.parent/"lib")
        for variable in ("FOAM_USER_LIBBIN","FOAM_LIBBIN","FOAM_EXT_LIBBIN"):
            if os.environ.get(variable): search_dirs.append(Path(os.environ[variable]))
        for variable in ("DYLD_LIBRARY_PATH","LD_LIBRARY_PATH"):
            search_dirs.extend(Path(item) for item in os.environ.get(variable," ").split(":") if item)
        candidates=[]
        for directory in search_dirs:
            if not directory.is_dir(): continue
            candidates.extend(sorted(directory.glob("libsolids4FoamModels.*")))
        # Prefer the candidate that actually registers the required law; this
        # matters when a machine has both a core and a user solids4Foam lib.
        for candidate in dict.fromkeys(candidates):
            marker_text=subprocess.run(["strings",str(candidate)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True).stdout
            if "GultekinTwoFibreElastic" in marker_text and "nonLinGeomTotalLagTotalDispGultekinSolid" in marker_text:
                lib=candidate.resolve(); break
        if lib is None and candidates: lib=candidates[0].resolve()
        source="active sourced environment"
    result={"source":source,"executable":str(exe) if exe else None,"library":str(lib) if lib else None,"passed":False,"preflight":"not attempted"}
    if not exe or not lib:
        result["error"]="solids4Foam executable/library not found"; return result
    result["executable_sha256"]=file_hash(exe); result["library_sha256"]=file_hash(lib)
    result["expected_executable_sha256"]=RUNTIME["accepted_executable_sha256"]
    result["expected_library_sha256"]=RUNTIME["accepted_library_sha256"]
    strings=subprocess.run(["strings",str(lib)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True).stdout
    markers=["GultekinTwoFibreElastic","nonLinearGeometryTotalLagrangianTotalDisplacement","implicitShearModulus","fibreUnitTolerance"]
    result["model_markers"]={m:(m in strings) for m in markers}
    if shutil.which("otool"):
        loader=subprocess.run(["otool","-L",str(exe)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        result["loader_dependencies"]=loader.stdout
    elif shutil.which("ldd"):
        loader=subprocess.run(["ldd",str(exe)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        result["loader_dependencies"]=loader.stdout
    else:
        result["loader_dependencies"]="otool/ldd unavailable"
    result["loader_diagnostics"]={"DYLD_PRINT_LIBRARIES":"available for exact replay","library_search_directory":str(lib.parent)}
    result["accepted_executable_hash_match"]=(result["executable_sha256"]==RUNTIME["accepted_executable_sha256"])
    result["accepted_library_hash_match"]=(result["library_sha256"]==RUNTIME["accepted_library_sha256"])
    result["hashes_match"]=(result["accepted_executable_hash_match"] and result["accepted_library_hash_match"])
    result["preflight"]="static model/material symbol preflight"
    if runtime_dir:
        result["hash_policy"]="strict accepted archived runtime hashes"
        result["passed"]=result["hashes_match"] and all(result["model_markers"].values())
        if not result["passed"]: result["error"]="staged runtime differs from accepted provenance or lacks required model markers"
    else:
        result["hash_policy"]="active sourced runtime accepted when required model/material markers are present"
        result["passed"]=all(result["model_markers"].values())
        if not result["passed"]: result["error"]="active runtime lacks required model/material markers"
    return result

def parse_internal(path):
    text=re.sub(r"//[^\n]*|/\*.*?\*/","",path.read_text(errors="replace"),flags=re.S)
    match=re.search(r"internalField\s+nonuniform\s+List<\w+>\s+(\d+)\s*\((.*?)\)\s*;",text,re.S)
    if not match:
        uniform=re.search(r"internalField\s+uniform\s+([^;]+);",text)
        if not uniform:return None
        token=uniform.group(1).strip(); return [tuple(map(float,token[1:-1].split())) if token.startswith("(") else float(token)]
    tuples=re.findall(r"\(([^()]*)\)",match.group(2))
    values=[tuple(map(float,x.split())) for x in tuples] if tuples else [float(x) for x in match.group(2).split()]
    return values

def result_gate(case):
    errors=[]; times=sorted([p for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?",p.name)],key=lambda p:float(p.name))
    if not times or float(times[-1].name)<1-1e-12: errors.append("full-load time 1 is absent")
    final=times[-1] if times else None
    if final:
        for name in ("sigmaRR_referenceBasis_kPa","sigmaTT_referenceBasis_kPa","sigmaZZ_referenceBasis_kPa"):
            vals=parse_internal(final/name) if (final/name).is_file() else None
            if not vals: errors.append(f"missing default stress field {name}")
            elif any(not math.isfinite(float(v)) for v in vals): errors.append(f"non-finite values in {name}")
        jpath=final/"GTF_J" if (final/"GTF_J").is_file() else final/"J"
        vals=parse_internal(jpath) if jpath.is_file() else None
        if vals and min(float(v) for v in vals)<=0: errors.append("J is not positive at full load")
        elif vals is None: errors.append("full-load J field is absent")
    return {"passed":not errors,"errors":errors,"final_time":final.name if final else None}

def validate(case, include_runtime=True, results=False, runtime_dir=None):
    errors=[]; warnings=[]
    if not (case/"caseManifest.json").is_file(): return {"case":str(case),"passed":False,"errors":["caseManifest.json missing"]}
    manifest=json.loads((case/"caseManifest.json").read_text()); counts=manifest["mesh_counts"]; mesh=load_mesh(case)
    expected=counts["radial"]*counts["circumferential"]*counts["axial"]
    if mesh["n_cells"]!=expected: errors.append(f"cell count {mesh['n_cells']} != {expected}")
    expected_patches={"innerWall":counts["circumferential"]*counts["axial"],"outerWall":counts["circumferential"]*counts["axial"],"bottom":counts["radial"]*counts["circumferential"],"top":counts["radial"]*counts["circumferential"]}
    actual={p:d["nFaces"] for p,d in mesh["boundary"].items()}
    if actual!=expected_patches: errors.append(f"patch sizes {actual} != {expected_patches}")
    geometry=manifest["geometry_m"]; radii=[math.hypot(p[0],p[1]) for p in mesh["points"]]; zs=[p[2] for p in mesh["points"]]
    if abs(min(radii)-geometry["inner_radius"])>1e-12 or abs(max(radii)-geometry["outer_radius"])>1e-12 or abs(min(zs))>1e-12 or abs(max(zs)-geometry["height"])>1e-12: errors.append("geometry tolerance failed")
    for rel in ("0/D","0/p","0/pointD","0/f0","0/f1","0/f0f","0/f1f","constant/solidProperties","constant/mechanicalProperties","system/fvSchemes","system/fvSolution","petscOptions.split"):
        if not (case/rel).is_file(): errors.append(f"missing {rel}")
    mechanical=(case/"constant/mechanicalProperties").read_text() if (case/"constant/mechanicalProperties").is_file() else ""
    solid=(case/"constant/solidProperties").read_text() if (case/"constant/solidProperties").is_file() else ""
    checks={"law":"type GultekinTwoFibreElastic;" in mechanical,"solid_model":"nonLinearGeometryTotalLagrangianTotalDisplacement;" in solid,"bulk":"5000000" in mechanical,"fibre_angle":abs(manifest["fibre_angle_degrees"]-40)<1e-12,"pressure":abs(manifest["final_loads"]["pressure_Pa"]-66661.184)<1e-9,"extension":abs(manifest["final_loads"]["extension_m"]-.002)<1e-15,"rotation":abs(manifest["final_loads"]["rotation_degrees"]-60)<1e-12,"schemes":"leastSquaresS4f" in (case/"system/fvSchemes").read_text() if (case/"system/fvSchemes").is_file() else False,"petsc":("-pc_type fieldsplit" in (case/"petscOptions.split").read_text() if (case/"petscOptions.split").is_file() else False)}
    errors += [f"input check failed: {name}" for name,passed in checks.items() if not passed]
    fibre=json.loads((case/"fibreMetadata.json").read_text()) if (case/"fibreMetadata.json").is_file() else {}
    if fibre.get("maximum_unit_error",1)>1e-12 or fibre.get("maximum_radial_component",1)>1e-12 or abs(fibre.get("angle_min_degrees",0)-40)>1e-10 or abs(fibre.get("angle_max_degrees",0)-40)>1e-10: errors.append("fibre normalization/angle gate failed")
    runtime=runtime_audit(runtime_dir) if include_runtime else {"passed":True,"skipped":True}
    quality=json.loads((case/"meshQuality.json").read_text()) if (case/"meshQuality.json").is_file() else {"passed":False}
    if not quality.get("passed"): errors.append("mesh quality gate failed or metadata is absent")
    result={"case":str(case),"mesh":case.name,"passed":not errors and (not include_runtime or runtime["passed"]),"errors":errors,"warnings":warnings,"checks":checks,"mesh_quality":quality,"runtime":runtime,"mesh_hash":mesh_hash(case)}
    if include_runtime and not runtime["passed"]:
        result["errors"].append("runtime provenance/model preflight failed")
    if results:
        result["results"]=result_gate(case); result["passed"]=result["passed"] and result["results"]["passed"]
    (case/"validation.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--meshes",default="mesh1,mesh2,mesh3,mesh4,mesh5,mesh6,mesh7"); ap.add_argument("--offline",action="store_true"); ap.add_argument("--results",action="store_true"); ap.add_argument("--runtime-dir")
    args=ap.parse_args(); rows=[validate(ROOT/"runs"/x.strip(),not args.offline,args.results,args.runtime_dir) for x in args.meshes.split(",") if x.strip()]
    (ROOT/"runs/validation_summary.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8"); print(json.dumps(rows,indent=2)); raise SystemExit(0 if all(x["passed"] for x in rows) else 1)
