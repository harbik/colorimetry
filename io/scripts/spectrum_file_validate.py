"""
spectrum_file_validate.py
--------------------------
Validate UV-Vis spectral JSON files against the spectrum_file_schema.json schema.

Usage:
    python spectrum_file_validate.py myfile.json
    python spectrum_file_validate.py myfile.json --schema path/to/spectrum_file_schema.json
    python spectrum_file_validate.py *.json
"""

import json
import sys
import argparse
from pathlib import Path

try:
    import jsonschema
    from jsonschema import Draft202012Validator
except ImportError:
    print("ERROR: jsonschema not installed. Run: pip install jsonschema")
    sys.exit(1)


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def validate_wavelength_data_length(data: dict) -> list[str]:
    """
    Cross-field check: wavelength_axis.values_nm and spectral_data.values
    must have the same length. JSON Schema can't enforce this natively.
    Also checks uncertainty length if present.
    """
    errors = []

    spectra = []
    if data.get("file_type") == "single":
        spectra = [data.get("spectrum", {})]
    elif data.get("file_type") == "batch":
        spectra = data.get("spectra", [])

    for sp in spectra:
        sid = sp.get("id", "<unknown>")
        wl_axis = sp.get("wavelength_axis", {})
        vals = sp.get("spectral_data", {}).get("values", [])
        unc = sp.get("spectral_data", {}).get("uncertainty")

        # Resolve wavelength count from either values_nm or range_nm
        if "values_nm" in wl_axis:
            wl = wl_axis["values_nm"]
            n_wl = len(wl)
        elif "range_nm" in wl_axis:
            r = wl_axis["range_nm"]
            try:
                n_wl = round((r["end"] - r["start"]) / r["interval"]) + 1
            except (KeyError, ZeroDivisionError):
                n_wl = 0
            wl = []  # no list to check monotonicity against for range_nm
        else:
            n_wl = 0
            wl = []

        if n_wl != len(vals):
            errors.append(
                f"Spectrum '{sid}': wavelength_axis has {n_wl} points "
                f"but spectral_data.values has {len(vals)} points — must match."
            )
        if unc is not None and len(unc) != len(vals):
            errors.append(
                f"Spectrum '{sid}': spectral_data.uncertainty has {len(unc)} points "
                f"but spectral_data.values has {len(vals)} points — must match."
            )

        # Check wavelengths are strictly increasing (only for explicit values_nm)
        if wl and any(wl[i] >= wl[i+1] for i in range(len(wl)-1)):
            errors.append(
                f"Spectrum '{sid}': wavelength_axis.values_nm is not strictly increasing."
            )

        # Check scale vs value range
        scale = sp.get("spectral_data", {}).get("scale", "fractional")
        mtype = sp.get("metadata", {}).get("measurement_type", "")
        if mtype in ("reflectance", "transmittance") and scale == "fractional":
            out_of_range = [v for v in vals if not (0.0 <= v <= 1.0)]
            if out_of_range:
                errors.append(
                    f"Spectrum '{sid}': measurement_type='{mtype}', scale='fractional' "
                    f"but {len(out_of_range)} value(s) fall outside [0, 1]. "
                    f"First offender: {out_of_range[0]}"
                )

    return errors


def validate_file(file_path: Path, schema: dict) -> bool:
    print(f"\n{'='*60}")
    print(f"Validating: {file_path}")
    print(f"{'='*60}")

    # Load data
    try:
        data = load_json(file_path)
    except json.JSONDecodeError as e:
        print(f"  ✗ Invalid JSON: {e}")
        return False

    # JSON Schema validation
    validator = Draft202012Validator(schema)
    schema_errors = list(validator.iter_errors(data))

    if schema_errors:
        print(f"  ✗ Schema validation failed ({len(schema_errors)} error(s)):")
        for err in schema_errors:
            path = " → ".join(str(p) for p in err.absolute_path) or "(root)"
            print(f"    • [{path}] {err.message}")
    else:
        print("  ✓ Schema validation passed")

    # Cross-field checks
    cross_errors = validate_wavelength_data_length(data)
    if cross_errors:
        print(f"  ✗ Cross-field validation failed ({len(cross_errors)} error(s)):")
        for err in cross_errors:
            print(f"    • {err}")
    else:
        print("  ✓ Cross-field checks passed")

    # Summary
    all_ok = not schema_errors and not cross_errors
    if all_ok:
        ftype = data.get("file_type", "?")
        if ftype == "batch":
            n = len(data.get("spectra", []))
            print(f"  ✓ Valid batch file — {n} spectrum/spectra")
        else:
            sid = data.get("spectrum", {}).get("id", "?")
            print(f"  ✓ Valid single-spectrum file — id: '{sid}'")
    else:
        print("  ✗ File is INVALID")

    return all_ok


def main():
    parser = argparse.ArgumentParser(description="Validate UV-Vis spectral JSON files.")
    parser.add_argument("files", nargs="+", help="JSON file(s) to validate")
    parser.add_argument(
        "--schema",
        default=Path(__file__).parent / "spectrum_file_schema.json",
        type=Path,
        help="Path to spectrum_file_schema.json (default: same directory as this script)"
    )
    args = parser.parse_args()

    # Load schema
    try:
        schema = load_json(args.schema)
        print(f"Schema loaded: {args.schema}")
    except FileNotFoundError:
        print(f"ERROR: Schema file not found: {args.schema}")
        sys.exit(1)

    # Validate each file
    results = []
    for pattern in args.files:
        paths = list(Path(".").glob(pattern)) if "*" in str(pattern) else [Path(pattern)]
        for p in paths:
            results.append(validate_file(p, schema))

    # Final summary
    total = len(results)
    passed = sum(results)
    failed = total - passed
    print(f"\n{'='*60}")
    print(f"Summary: {passed}/{total} file(s) valid", end="")
    if failed:
        print(f", {failed} invalid")
    else:
        print(" ✓")
    print(f"{'='*60}")

    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
