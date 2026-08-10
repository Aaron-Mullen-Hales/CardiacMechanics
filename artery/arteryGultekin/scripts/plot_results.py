#!/usr/bin/env python3
"""Create a compact mesh summary plot when matplotlib is available."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--meshes",default="mesh1,mesh2,mesh3,mesh4,mesh5,mesh6,mesh7"); args=ap.parse_args()
    rows=[]
    for name in [x.strip() for x in args.meshes.split(",") if x.strip()]:
        p=ROOT/"runs"/name/"referenceBasisStress.json"
        if p.is_file(): rows.append(json.loads(p.read_text()))
    out=ROOT/"runs/plot_data.csv"; out.parent.mkdir(exist_ok=True)
    with out.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=["mesh","final_time","maximum_inplane_trace_identity_error_kPa","maximum_full_trace_identity_error_kPa"]); writer.writeheader()
        writer.writerows(rows)
    try:
        import matplotlib.pyplot as plt
        if rows:
            x=[r["mesh"] for r in rows]; y=[r["maximum_full_trace_identity_error_kPa"] for r in rows]
            plt.figure(); plt.plot(x,y,"o-"); plt.ylabel("max trace identity error (kPa)"); plt.xlabel("mesh"); plt.tight_layout(); plt.savefig(ROOT/"runs/stress_trace_identity.png",dpi=160); plt.close()
    except ImportError:
        pass
    print(json.dumps({"rows":len(rows),"csv":str(out)},indent=2))
if __name__=="__main__": main()
