#!/usr/bin/env python3
"""Public workflow behind Allrun."""
from __future__ import annotations
import argparse, json, os, re, shutil, subprocess, sys
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path

try:
    import fcntl
except ImportError:  # pragma: no cover - OpenFOAM is normally run on Unix.
    fcntl = None

ROOT=Path(__file__).resolve().parents[1]; SCRIPTS=ROOT/"scripts"
NAMES=("mesh1","mesh2","mesh3","mesh4","mesh5","mesh6","mesh7")


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def last_time(log):
    if not log.exists():
        return None
    matches = re.findall(r"^Time = ([^\r\n]+)", log.read_text(encoding="utf-8", errors="replace"), re.MULTILINE)
    return matches[-1].strip() if matches else None


class CaseLock:
    """Prevent generation, solving, and postprocessing from overlapping."""

    def __init__(self, name):
        self.name = name
        self.case = ROOT / "runs" / name
        self.path = self.case / ".solids4Foam.lock"
        self.stream = None

    def __enter__(self):
        if fcntl is None:
            raise SystemExit("Case locking requires a Unix-like system with fcntl")
        self.case.mkdir(parents=True, exist_ok=True)
        self.stream = self.path.open("a+", encoding="utf-8")
        try:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.stream.close()
            self.stream = None
            raise SystemExit(
                f"{self.name} is already being generated, solved, or postprocessed; "
                f"wait for that run to finish before starting another"
            )
        self.stream.seek(0)
        self.stream.truncate()
        self.stream.write(json.dumps({"pid": os.getpid(), "started_utc": utc_now()}) + "\n")
        self.stream.flush()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.stream is not None:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_UN)
            self.stream.close()

def selected(args):
    if args.mesh and args.meshes: raise SystemExit("use --mesh or --meshes, not both")
    raw=args.meshes or args.mesh or "mesh2"
    names=[x.strip() for x in raw.split(",") if x.strip()]
    unknown=sorted(set(names)-set(NAMES))
    if unknown: raise SystemExit("Unknown mesh names: "+",".join(unknown))
    return names

def call(script,*args):
    command=[sys.executable,str(SCRIPTS/script),*map(str,args)]
    print("+ "+" ".join(command),flush=True); result=subprocess.run(command,cwd=ROOT)
    if result.returncode: raise SystemExit(result.returncode)

def runtime_info(runtime_dir=None):
    command=[sys.executable,str(SCRIPTS/"validate_case.py"),"--meshes","mesh2","--offline"]
    # The real runtime gate is performed by validate_case's Python API below.
    sys.path.insert(0,str(SCRIPTS)); from validate_case import runtime_audit
    return runtime_audit(runtime_dir)

def run_case(name,runtime_dir=None,dry=False):
    info=runtime_info(runtime_dir)
    if not info["passed"]:
        print(json.dumps({"runtime_gate":"failed","details":info},indent=2),file=sys.stderr)
        raise SystemExit("runtime gate failed before solving")
    case=ROOT/"runs"/name
    exe=Path(info["executable"])
    env=os.environ.copy(); env["FOAM_CASE"]=str(case)
    loader_dirs=[]
    if info.get("library"): loader_dirs.append(str(Path(info["library"]).parent))
    for variable in ("FOAM_LD_LIBRARY_PATH","DYLD_LIBRARY_PATH","LD_LIBRARY_PATH"):
        loader_dirs.extend(item for item in env.get(variable,"").split(":") if item)
    # Preserve order while removing duplicates. FOAM_LD_LIBRARY_PATH is the
    # variable OpenFOAM's bashrc uses on macOS; DYLD/LD are used by loaders.
    loader_dirs=list(dict.fromkeys(loader_dirs))
    if loader_dirs:
        loader_path=":".join(loader_dirs)
        env["FOAM_LD_LIBRARY_PATH"]=loader_path
        env["DYLD_LIBRARY_PATH"]=loader_path
        env["LD_LIBRARY_PATH"]=loader_path
    control=case/"system/controlDict"
    selected_control=case/"system/controlDict.dryRun" if dry else control
    backup=control.read_bytes()
    log=case/("log.solids4Foam.dryRun" if dry else "log.solids4Foam")
    state_path=case/"runState.json"
    state={
        "mesh": name,
        "status": "running",
        "pid": os.getpid(),
        "started_utc": utc_now(),
        "log": str(log.relative_to(ROOT)),
        "executable": str(exe),
    }
    write_json(state_path, state)
    try:
        if dry: shutil.copy2(selected_control,control)
        with log.open("w",encoding="utf-8") as stream:
            result=subprocess.run([str(exe),"-case",str(case)],cwd=case,env=env,stdout=stream,stderr=subprocess.STDOUT)
        state.update(
            status="completed" if result.returncode == 0 else "failed",
            return_code=result.returncode,
            final_time=last_time(log),
            finished_utc=utc_now(),
        )
        write_json(state_path, state)
        if result.returncode:
            raise SystemExit(f"solids4Foam failed for {name}; see {log}")
    except BaseException as error:
        if state["status"] == "running":
            state.update(status="failed", error=str(error), finished_utc=utc_now(), final_time=last_time(log))
            write_json(state_path, state)
        raise
    finally:
        if dry: control.write_bytes(backup)

def main():
    ap=argparse.ArgumentParser(prog="./Allrun")
    ap.add_argument("--mesh"); ap.add_argument("--meshes"); ap.add_argument("--runtime-dir")
    ap.add_argument("--generate-only",action="store_true"); ap.add_argument("--validate-only",action="store_true"); ap.add_argument("--run-only",action="store_true"); ap.add_argument("--postprocess-only",action="store_true"); ap.add_argument("--plot-only",action="store_true"); ap.add_argument("--dry-run",action="store_true"); ap.add_argument("--diagnostics",action="store_true")
    args=ap.parse_args(); names=selected(args)
    stages=sum(bool(getattr(args,x)) for x in ("generate_only","validate_only","run_only","postprocess_only","plot_only"))
    if stages>1: raise SystemExit("workflow stage flags are mutually exclusive")
    with ExitStack() as locks:
        for name in names:
            locks.enter_context(CaseLock(name))
        if args.generate_only:
            call("generate_case.py","--meshes",",".join(names)); return
        if args.validate_only:
            call("validate_case.py","--meshes",",".join(names),*([] if args.runtime_dir is None else ["--runtime-dir",args.runtime_dir])); return
        if args.postprocess_only:
            call("write_reference_basis_stress.py","--meshes",",".join(names),*( ["--diagnostics"] if args.diagnostics else [])); return
        if args.plot_only:
            call("plot_results.py","--meshes",",".join(names)); return
        if not args.run_only:
            call("generate_case.py","--meshes",",".join(names))
            call("validate_case.py","--meshes",",".join(names),"--offline")
        for name in names: run_case(name,args.runtime_dir,args.dry_run)
        call("write_reference_basis_stress.py","--meshes",",".join(names),*( ["--diagnostics"] if args.diagnostics else []))
        call("validate_case.py","--meshes",",".join(names),"--results",*([] if args.runtime_dir is None else ["--runtime-dir",args.runtime_dir]))

if __name__=="__main__": main()
