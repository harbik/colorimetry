#!/usr/bin/env python3
"""
Validate the colorimetry Rust library against LuxPy's TM-30 implementation.

LuxPy reference:  https://github.com/ksmet1977/luxpy
Standard:         ANSI/IES TM-30-20 / CIE 224:2017

Run with:
    python3 validate_tm30.py
"""

import numpy as np
import luxpy as lx
from luxpy.color.cri.iestm30.metrics import spd_to_ies_tm30_metrics

# ── helpers ──────────────────────────────────────────────────────────────────

def tm30(spd):
    """Compute TM-30 metrics for a given SPD array (row 0 = wavelengths)."""
    return spd_to_ies_tm30_metrics(spd, cri_type='ies-tm30')

def scalar(x):
    """Squeeze a 1-element array to a plain float."""
    return float(np.squeeze(x))

def vec16(x):
    """Squeeze a (1,16) or (16,) array to a list of 16 floats."""
    return np.squeeze(x).tolist()

def print_section(title):
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)

def compare(label, rust, python, tol):
    ok = abs(rust - python) <= tol
    flag = "  OK" if ok else "  FAIL <<<<<"
    print(f"  {label:<35} Rust={rust:8.3f}  LuxPy={python:8.3f}  Δ={rust-python:+.3f}{flag}")

def compare_bins(label, rust_vals, py_vals, tol):
    print(f"  {label}")
    any_fail = False
    for j, (r, p) in enumerate(zip(rust_vals, py_vals)):
        ok = abs(r - p) <= tol
        flag = "" if ok else "  <<< FAIL"
        if not ok:
            any_fail = True
        print(f"    bin {j:2d}  Rust={r:8.3f}  LuxPy={p:8.3f}  Δ={r-p:+.3f}{flag}")
    if not any_fail:
        print(f"    All 16 bins within ±{tol}")

# ── get LuxPy's built-in CIE illuminant SPDs ─────────────────────────────────

# LuxPy stores illuminants in lx._CIE_ILLUMINANTS
# Each SPD is shape (2, N): row 0 = wavelengths, row 1 = power values.
D65  = lx._CIE_ILLUMINANTS['D65']
F1   = lx._CIE_ILLUMINANTS['F1']
F2   = lx._CIE_ILLUMINANTS['F2']
F12  = lx._CIE_ILLUMINANTS['F12']

illuminants = [
    ("D65",  D65),
    ("F1",   F1),
    ("F2",   F2),
    ("F12",  F12),
]

# ── Rust reference values (from passing tests in cfi.rs) ─────────────────────
# Format: {name: {metric: value}}
rust_ref = {
    "D65":  {"cct": 6504.0, "duv": 0.003,  "Rf": 100.0,  "Rg": 100.0},
    "F1":   {"cct": 6425.0, "duv": 0.0072, "Rf": 80.68,  "Rg": 89.83},
    "F2":   {"cct": 4225.0, "duv": 0.0019, "Rf": 70.21,  "Rg": 86.44},
    "F12":  {"cct": 3003.0, "duv": 0.0001, "Rf": 77.7,   "Rg": 102.4},
}

# ── run comparisons ───────────────────────────────────────────────────────────

for name, spd in illuminants:
    print_section(f"CIE {name}")
    d = tm30(spd)

    py_cct  = scalar(d['cct'])
    py_duv  = scalar(d['duv'])
    py_Rf   = scalar(d['Rf'])
    py_Rg   = scalar(d['Rg'])
    py_Rfhj  = vec16(d['Rfhj'])
    py_Rcshj = vec16(d['Rcshj'])
    py_Rhshj = vec16(d['Rhshj'])

    r = rust_ref[name]

    print()
    print("  Scalar metrics (tolerance in parentheses):")
    compare("CCT  (K, ±10)",  r["cct"],  py_cct,  10.0)
    compare("Duv  (±0.001)",  r["duv"],  py_duv,   0.001)
    compare("Rf   (±0.6)",    r["Rf"],   py_Rf,    0.6)
    compare("Rg   (±0.6)",    r["Rg"],   py_Rg,    0.6)

    print()
    print("  Local metrics from LuxPy (Rust not yet implemented — for reference):")
    print(f"  {'Bin':>4}  {'Rf,hj':>8}  {'Rcs,hj %':>10}  {'Rhs,hj':>8}")
    for j in range(16):
        angle = j * 22.5
        print(f"  {j:4d}  {py_Rfhj[j]:8.2f}  {py_Rcshj[j]:10.2f}  {py_Rhshj[j]:8.4f}   ({angle:.1f}°–{angle+22.5:.1f}°)")

print()
print("=" * 70)
print("  Done. Rows marked <<< FAIL exceed tolerance.")
print("  Local metrics (Rf,hj / Rcs,hj / Rhs,hj) are LuxPy reference values")
print("  to use when implementing those metrics in the Rust library.")
print("=" * 70)
