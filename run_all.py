#!/usr/bin/env python3
"""
Run all meta-analysis scripts sequentially.

GLP-1 / Dual GIP/GLP-1 Agonists — Lean Body Mass Meta-Analysis
Alayci et al. 2026

Usage:
    python run_all.py
"""

import subprocess
import sys
import os

SCRIPTS = [
    ("01_primary_meta_analysis.py",         "Primary DL meta-analysis, subgroups, leave-one-out"),
    ("02_publication_bias.py",              "Doi/LFK, Peters, Begg, Egger, funnel plot"),
    ("03_advanced_sensitivity.py",          "HKSJ, prediction interval, Bayesian (overall)"),
    ("04_placebo_subgroup_advanced.py",     "HKSJ, prediction interval, Bayesian (placebo k=2)"),
    ("05_tau2_estimators_rve.py",           "REML, Paule-Mandel, profile likelihood, RVE, 3-level"),
]

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    print("=" * 70)
    print("  GLP-1 / Dual Agonist — LBM Meta-Analysis")
    print("  Running all analyses...")
    print("=" * 70)
    
    failed = []
    for script, description in SCRIPTS:
        path = os.path.join(script_dir, script)
        print(f"\n{'─' * 70}")
        print(f"  Running: {script}")
        print(f"  {description}")
        print(f"{'─' * 70}")
        
        result = subprocess.run(
            [sys.executable, path],
            capture_output=False
        )
        
        if result.returncode != 0:
            print(f"  ⚠ FAILED: {script}")
            failed.append(script)
        else:
            print(f"  ✓ Completed: {script}")
    
    print(f"\n{'=' * 70}")
    if failed:
        print(f"  ⚠ {len(failed)} script(s) failed: {', '.join(failed)}")
    else:
        print(f"  ✓ All {len(SCRIPTS)} scripts completed successfully.")
        print(f"  Figures saved to figures/ directory.")
    print(f"{'=' * 70}")

if __name__ == "__main__":
    main()
