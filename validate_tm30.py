#!/usr/bin/env python3
"""
Validate the colorimetry Rust library against LuxPy's TM-30 implementation.

LuxPy reference:  https://github.com/ksmet1977/luxpy  (tested with v1.12.5)
Standard:         ANSI/IES TM-30-20 / CIE 224:2017

Note on rust_ref values
-----------------------
This script cannot invoke the Rust binary directly, so `rust_ref` contains the
expected values that the Rust library should produce — taken from the assertions
in the passing Rust unit tests in src/illuminant/cfi.rs.  The comparison is
therefore:

    LuxPy (independent Python implementation) ↔ expected Rust output

If a metric is flagged FAIL, it means LuxPy disagrees with the Rust library's
expected output by more than the stated tolerance.

For F12 local metrics (Rf,hj / Rcs,hj / Rhs,hj), the expected Rust values come
from the IES TM-30 Spectral Calculator JSON (an independent reference from the
standard authors).  For D65/F1/F2 local metrics, the Rust tests use LuxPy itself
as their baseline, so those comparisons are circular and appear here mainly as
a consistency check.

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

# LuxPy stores illuminants in lx._CIE_ILLUMINANTS (internal API, verified on
# v1.12.5 — pin your LuxPy version if this breaks across upgrades).
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

# ── Rust reference values ─────────────────────────────────────────────────────
# Global metrics (CCT/Duv/Rf/Rg): from passing Rust unit tests in cfi.rs.
# Local metrics:
#   D65  — trivially perfect (test source == reference illuminant)
#   F1   — Rust tests use LuxPy as baseline (comparison below is circular)
#   F2   — same as F1
#   F12  — Rust tests use IES TM-30 Spectral Calculator JSON as baseline
#           (independent from LuxPy; comparison below is a genuine cross-check)
#
# Rcs,hj is a dimensionless fraction (not percent). LuxPy's 'Rcshj' key also
# returns fractions, so no unit conversion is needed.
rust_ref = {
    "D65": {
        "cct": 6504.0, "duv": 0.003, "Rf": 100.0, "Rg": 100.0,
        "Rfhj":  [100.0] * 16,
        "Rcshj": [0.0]   * 16,
        "Rhshj": [0.0]   * 16,
    },
    "F1": {
        "cct": 6425.0, "duv": 0.0072, "Rf": 80.68, "Rg": 89.83,
        "Rfhj":  [64.52, 76.27, 70.52, 81.72, 86.43, 91.88, 88.93, 81.30,
                  86.87, 79.56, 82.99, 90.61, 86.35, 76.16, 69.85, 76.14],
        "Rcshj": [-0.20, -0.13, -0.07,  0.02,  0.07,  0.01, -0.05, -0.10,
                  -0.11, -0.07, -0.01,  0.04,  0.07,  0.02, -0.09, -0.10],
        "Rhshj": [-0.0248,  0.0843,  0.1564,  0.1126,  0.0579, -0.0427,
                  -0.0474, -0.0511,  0.0252,  0.0975,  0.0985,  0.0302,
                  -0.0833, -0.1430, -0.2770, -0.1015],
    },
    "F2": {
        "cct": 4225.0, "duv": 0.0019, "Rf": 70.21, "Rg": 86.44,
        "Rfhj":  [60.20, 61.28, 52.52, 68.30, 79.50, 87.51, 76.74, 72.71,
                  76.13, 62.27, 69.57, 76.48, 81.32, 71.36, 63.60, 65.27],
        "Rcshj": [-0.25, -0.18, -0.09,  0.05,  0.11,  0.04, -0.08, -0.15,
                  -0.17, -0.15, -0.04,  0.05,  0.11,  0.07, -0.06, -0.16],
        "Rhshj": [-0.0221,  0.1400,  0.2444,  0.1963,  0.0915, -0.0673,
                  -0.1226, -0.0848,  0.0061,  0.1659,  0.1906,  0.1147,
                  -0.0806, -0.1467, -0.2635, -0.1703],
    },
    "F12": {
        "cct": 3003.0, "duv": 0.0001, "Rf": 77.7, "Rg": 102.4,
        # Source: IES TM-30 Spectral Calculator JSON (independent from LuxPy)
        "Rfhj":  [78.222, 86.843, 80.612, 68.521, 72.478, 79.408, 72.551, 79.812,
                  82.356, 77.682, 74.540, 80.448, 82.391, 76.766, 79.243, 76.808],
        "Rcshj": [-0.08510, -0.03160, -0.01216,  0.09247,  0.18373,  0.13567,
                   0.10277, -0.03677, -0.04746, -0.12181, -0.12371,  0.01845,
                   0.08819,  0.04632,  0.04186, -0.03594],
        "Rhshj": [-0.04646,  0.03206,  0.09473,  0.17904,  0.12972,  0.00023,
                  -0.14516, -0.12447, -0.09830,  0.03149,  0.14218,  0.08155,
                  -0.02598, -0.06301, -0.09517, -0.10676],
    },
}

# ── run comparisons ───────────────────────────────────────────────────────────

for name, spd in illuminants:
    print_section(f"CIE {name}")
    d = tm30(spd)

    py_cct   = scalar(d['cct'])
    py_duv   = scalar(d['duv'])
    py_Rf    = scalar(d['Rf'])
    py_Rg    = scalar(d['Rg'])
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
    if name in ("F1", "F2"):
        print("  Local metrics — Rust baseline is LuxPy itself (circular check):")
    elif name == "F12":
        print("  Local metrics — Rust baseline is IES TM-30 Spectral Calculator:")
    else:
        print("  Local metrics — D65 is its own reference (trivially perfect):")
    compare_bins("Rf,hj  (±10.0)",       r["Rfhj"],  py_Rfhj,  10.0)
    compare_bins("Rcs,hj (±0.05)",       r["Rcshj"], py_Rcshj,  0.05)
    compare_bins("Rhs,hj (±0.04 rad)",   r["Rhshj"], py_Rhshj,  0.04)

print()
print("=" * 70)
print("  Done. Rows marked <<< FAIL exceed tolerance.")
print("  F12 local metrics are cross-checked against the IES TM-30 Spectral")
print("  Calculator (independent reference). F1/F2 use LuxPy as their baseline")
print("  (circular, but confirms internal LuxPy consistency).")
print("=" * 70)
