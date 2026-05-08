<!-- cargo-rdme start -->

# colorimetry-io

Reader and validator for spectral data files following the
`spectrum_file_schema.json` v1.0.0 specification — a JSON format designed
for UV-Vis and visible-range spectral measurements suitable for color-science
calculations, data sharing, and long-term archiving.

## Quick start

```rust
use colorimetry_io::SpectrumFile;

let file = SpectrumFile::from_path("spectrum.json").expect("could not load file");
for sp in file.spectra() {
    println!("{}: {} points", sp.id, sp.n_points());
}
```

## File format

Files are JSON objects with a `schema_version` (semver string) and a
`file_type` of either `"single"` or `"batch"`.

### Single file

Contains exactly one spectrum under the key `"spectrum"`:

```json
{
  "schema_version": "1.0.0",
  "file_type": "single",
  "spectrum": { ... }
}
```

### Batch file

Contains one or more spectra under `"spectra"`, with an optional
`"batch_metadata"` block that holds metadata common to the whole set
(title, operator, instrument, measurement conditions):

```json
{
  "schema_version": "1.0.0",
  "file_type": "batch",
  "batch_metadata": { "title": "Munsell chips set A", "date": "2026-04-01" },
  "spectra": [ { ... }, { ... } ]
}
```

### SpectrumRecord object

Each spectrum has four required sections and two optional ones.

#### `metadata` (required)

| Field | Type | Notes |
|---|---|---|
| `measurement_type` | string enum | `reflectance`, `transmittance`, `absorbance`, `radiance`, `irradiance` |
| `date` | string | ISO 8601 date (`YYYY-MM-DD`) |
| `title` | string | optional human-readable name |
| `sample_id` | string | optional sample identifier |
| `operator` | string | optional name or ID of the operator |
| `instrument` | object | optional: `manufacturer`, `model`, `serial_number`, `detector_type`, `light_source` |
| `measurement_conditions` | object | optional: `integration_time_ms`, `averaging`, `temperature_celsius`, `geometry`, `specular_component`, `spectral_resolution_nm` |
| `tags` | string[] | optional free-form search/filter tags |
| `custom` | object | optional user-defined key/value pairs |

#### `wavelength_axis` (required)

Exactly one of `values_nm` or `range_nm` must be present — not both, not neither.

| Field | Type | Notes |
|---|---|---|
| `values_nm` | number[] | explicit wavelength list in nm, min 2 entries, strictly increasing; use for irregular grids |
| `range_nm` | object | evenly-spaced grid defined by `start`, `end`, and `interval` (all in nm); use for regular grids |

#### `spectral_data` (required)

| Field | Type | Notes |
|---|---|---|
| `values` | number[] | measured values, one per entry in `values_nm` |
| `uncertainty` | number[] | optional 1-σ per-point uncertainty, same length as `values` |
| `scale` | string enum | `"fractional"` (0–1, default) or `"percent"` (0–100) |

#### `color_science` (optional)

Metadata needed for CIE colorimetric calculations, plus optional pre-computed results:

| Field | Type | Notes |
|---|---|---|
| `illuminant` | string enum | CIE illuminant code: `D65`, `D50`, `A`, `F1`–`F12`, `LED-*`, or `"custom"` |
| `illuminant_custom_sd` | object | required when `illuminant` is `"custom"`; contains `wavelengths_nm` and `values` arrays |
| `cie_observer` | string enum | `"CIE 1931 2 degree"` (default) or `"CIE 1964 10 degree"` |
| `white_reference` | object | optional calibration tile info and spectral reflectance values |
| `results` | object | optional pre-computed colorimetric values (see below) |

All `results` fields are optional and informational — the spectral data is always the
authoritative source. Any subset of the following may be present:

| Field | Type | Notes |
|---|---|---|
| `XYZ` | `[number, number, number]` | CIE tristimulus values [X, Y, Z] |
| `xy` | `[number, number]` | CIE 1931 chromaticity coordinates [x, y] |
| `uv_prime` | `[number, number]` | CIE 1976 UCS chromaticity coordinates [u′, v′] |
| `Lab` | `[number, number, number]` | CIELAB coordinates [L\*, a\*, b\*] |
| `CCT_K` | number | Correlated color temperature in Kelvin |
| `Duv` | number | Distance from the Planckian locus (signed) in the CIE 1960 UCS |

#### `provenance` (optional)

Processing history: `software`, `software_version`, `source_file`,
`source_format`, `notes`, and an ordered `processing_steps` array
(each step has a `step` name, `description`, and optional `parameters` object).

### Full single-spectrum example

```json
{
  "schema_version": "1.0.0",
  "file_type": "single",
  "spectrum": {
    "id": "chip-5R-4-2",
    "metadata": {
      "measurement_type": "reflectance",
      "date": "2026-04-01",
      "title": "Munsell 5R 4/2",
      "instrument": { "manufacturer": "Konica Minolta", "model": "CM-700d" },
      "measurement_conditions": { "geometry": "d:8", "specular_component": "excluded" }
    },
    "wavelength_axis": {
      "range_nm": { "start": 380, "end": 780, "interval": 10 }
    },
    "spectral_data": {
      "values": [0.048, 0.051, 0.054, 0.058, 0.063],
      "scale": "fractional"
    },
    "color_science": {
      "illuminant": "D65",
      "cie_observer": "CIE 1931 2 degree",
      "results": {
        "XYZ": [17.35, 9.12, 1.18],
        "xy": [0.629, 0.330],
        "Lab": [36.1, 55.7, 37.2]
      }
    }
  }
}
```

## Validation

[`SpectrumFile::from_path`] and [`SpectrumFile::from_str`] run two validation
passes before returning:

1. **Schema validation** — checks required fields, correct types, and that
   enum fields (`measurement_type`, `illuminant`, `cie_observer`,
   `scale`) contain only allowed values.
2. **Cross-field validation** — checks that `values_nm` and `values` have
   equal length; that `uncertainty` (if present) has the same length; that
   wavelengths are strictly increasing; that reflectance/transmittance values
   lie in \[0, 1\] when `scale` is `"fractional"`; and that a custom
   illuminant is accompanied by `illuminant_custom_sd`.

Use [`SpectrumFile::from_str_unchecked`] to skip all validation when the
source is fully trusted.

## Converting to `Spectrum`

With the default `colorimetry` feature enabled, [`SpectrumRecord`] provides four
methods to convert spectral data to the fixed 380–780 nm / 1 nm grid that
[`colorimetry::spectrum::Spectrum`] requires:

| Method | When to use |
|---|---|
| [`SpectrumRecord::to_spectrum_linear`] | Any axis; simple fallback |
| [`SpectrumRecord::to_spectrum_sprague`] | Equidistant (`range_nm`) axis; accurate for coarse regular grids |
| [`SpectrumRecord::to_spectrum_smooth`] | Sub-nm pixels with a known slit width; Gaussian-weighted, uses all pixel data |
| [`SpectrumRecord::to_spectrum_binned`] | Sub-nm pixels; boxcar average to slit-width bins then Sprague |

`to_spectrum_sprague` returns an error if the axis uses `values_nm` (irregular grid).
`to_spectrum_smooth` and `to_spectrum_binned` are the right choice when the pixel
spacing is much finer than the instrument slit width — they aggregate all pixel data
rather than picking isolated samples.

<!-- cargo-rdme end -->
