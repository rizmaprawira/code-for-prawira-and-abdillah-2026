#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

RUNS = [
    "CTRL",
    "EXP_B",
    "EXP_S-123",
    "EXP_B_S-123",
    "EXP_B_S-1",
    "EXP_B_S-12",
]

def run_cmd(args, cwd):
    print(f"Running command: {' '.join(args)} in {cwd}")
    res = subprocess.run(args, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"Error executing command: {res.stderr}")
        return False
    return True

def main():
    results_dir = BASE_DIR / "results"
    t15_dir = results_dir / "t15"
    t29_dir = results_dir / "t29"

    t15_dir.mkdir(parents=True, exist_ok=True)
    t29_dir.mkdir(parents=True, exist_ok=True)

    print(f"Ensured directories exist:\n  {t15_dir}\n  {t29_dir}")

    for run in RUNS:
        run_path = BASE_DIR / run
        if not run_path.exists():
            print(f"Warning: Experiment folder {run} does not exist. Skipping.")
            continue

        ctl_path = run_path / "post" / "linear.t21l20.moist.ctl"
        if not ctl_path.exists():
            print(f"Notice: Output CTL for {run} does not exist yet. Run model first.")
            continue

        # Check/Generate psi_q plots if needed
        psi_q_dir = run_path / "plots" / "psi_q"
        existing_pngs = list(psi_q_dir.glob("*.png")) if psi_q_dir.exists() else []
        if len(existing_pngs) < 29:
            print(f"Generating psi_q plots for {run}...")
            run_cmd([
                "python3", "scripts/plot_psi_q.py",
                "--ctl", str(ctl_path),
                "--outdir", "plots/psi_q"
            ], cwd=str(run_path))

        # Copy t15
        for var_name, sub, prefix in [("psi850", "psi850", "psi_850"), ("psi500", "psi500", "psi_500"), ("z850", "z850", "z_850"), ("z500", "z500", "z_500")]:
            src_t15 = run_path / "plots" / sub / f"{prefix}_t015.png"
            if src_t15.exists():
                shutil.copy2(src_t15, t15_dir / f"{run}_{prefix}_t015.png")

            src_t29 = run_path / "plots" / sub / f"{prefix}_t029.png"
            if src_t29.exists():
                shutil.copy2(src_t29, t29_dir / f"{run}_{prefix}_t029.png")

        print(f"Gathered results for {run}")

    print("All gather operations complete.")

if __name__ == "__main__":
    main()
