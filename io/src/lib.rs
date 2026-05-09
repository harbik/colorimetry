//! # colorimetry-io
//!
//! Reader and validator for spectral data files following the
//! `spectrum_file_schema.json` v1.0.0 specification — a JSON format designed
//! for UV-Vis and visible-range spectral measurements suitable for color-science
//! calculations, data sharing, and long-term archiving.
//!
//! ## Quick start
//!
//! ```no_run
//! use colorimetry_io::SpectrumFile;
//!
//! let file = SpectrumFile::from_path("spectrum.json").expect("could not load file");
//! for sp in file.spectra() {
//!     println!("{}: {} points", sp.id, sp.n_points());
//! }
//! ```
//!
//! ## File format
//!
//! Files are JSON objects with a `schema_version` (semver string) and a
//! `file_type` of either `"single"` or `"batch"`.
//!
//! ### Single file
//!
//! Contains exactly one spectrum under the key `"spectrum"`:
//!
//! ```json
//! {
//!   "schema_version": "1.0.0",
//!   "file_type": "single",
//!   "spectrum": { ... }
//! }
//! ```
//!
//! ### Batch file
//!
//! Contains one or more spectra under `"spectra"`, with an optional
//! `"batch_metadata"` block that holds metadata common to the whole set
//! (title, operator, instrument, measurement conditions):
//!
//! ```json
//! {
//!   "schema_version": "1.0.0",
//!   "file_type": "batch",
//!   "batch_metadata": { "title": "Munsell chips set A", "date": "2026-04-01" },
//!   "spectra": [ { ... }, { ... } ]
//! }
//! ```
//!
//! ### SpectrumRecord object
//!
//! Each spectrum has four required sections and two optional ones.
//!
//! #### `metadata` (required)
//!
//! | Field | Type | Notes |
//! |---|---|---|
//! | `measurement_type` | string enum | `reflectance`, `transmittance`, `absorbance`, `radiance`, `irradiance` |
//! | `date` | string | ISO 8601 date (`YYYY-MM-DD`) |
//! | `title` | string | optional human-readable name |
//! | `sample_id` | string | optional sample identifier |
//! | `operator` | string | optional name or ID of the operator |
//! | `instrument` | object | optional: `manufacturer`, `model`, `serial_number`, `detector_type`, `light_source` |
//! | `measurement_conditions` | object | optional: `integration_time_ms`, `averaging`, `temperature_celsius`, `geometry`, `specular_component`, `spectral_resolution_nm` |
//! | `tags` | string[] | optional free-form search/filter tags |
//! | `custom` | object | optional user-defined key/value pairs |
//!
//! #### `wavelength_axis` (required)
//!
//! Exactly one of `values_nm` or `range_nm` must be present — not both, not neither.
//!
//! | Field | Type | Notes |
//! |---|---|---|
//! | `values_nm` | number[] | explicit wavelength list in nm, min 2 entries, strictly increasing; use for irregular grids |
//! | `range_nm` | object | evenly-spaced grid defined by `start`, `end`, and `interval` (all in nm); use for regular grids |
//!
//! #### `spectral_data` (required)
//!
//! | Field | Type | Notes |
//! |---|---|---|
//! | `values` | number[] | measured values, one per entry in `values_nm` |
//! | `uncertainty` | number[] | optional 1-σ per-point uncertainty, same length as `values` |
//! | `scale` | string enum | `"fractional"` (0–1, default) or `"percent"` (0–100) |
//!
//! #### `color_science` (optional)
//!
//! Metadata needed for CIE colorimetric calculations, plus optional pre-computed results:
//!
//! | Field | Type | Notes |
//! |---|---|---|
//! | `illuminant` | string enum | CIE illuminant code: `D65`, `D50`, `A`, `F1`–`F12`, `LED-*`, or `"custom"` |
//! | `illuminant_custom_sd` | object | required when `illuminant` is `"custom"`; contains `wavelengths_nm` and `values` arrays |
//! | `cie_observer` | string enum | `"CIE 1931 2 degree"` (default), `"CIE 1964 10 degree"`, `"CIE 2015 2 degree"`, or `"CIE 2015 10 degree"` |
//! | `white_reference` | object | optional calibration tile info and spectral reflectance values |
//! | `results` | object | optional pre-computed colorimetric values (see below) |
//!
//! All `results` fields are optional and informational — the spectral data is always the
//! authoritative source. Any subset of the following may be present:
//!
//! | Field | Type | Notes |
//! |---|---|---|
//! | `XYZ` | `[number, number, number]` | CIE tristimulus values [X, Y, Z] |
//! | `xy` | `[number, number]` | CIE 1931 chromaticity coordinates [x, y] |
//! | `uv_prime` | `[number, number]` | CIE 1976 UCS chromaticity coordinates [u′, v′] |
//! | `Lab` | `[number, number, number]` | CIELAB coordinates [L\*, a\*, b\*] |
//! | `CCT_K` | number | Correlated color temperature in Kelvin |
//! | `Duv` | number | Distance from the Planckian locus (signed) in the CIE 1960 UCS |
//!
//! #### `provenance` (optional)
//!
//! Processing history: `software`, `software_version`, `source_file`,
//! `source_format`, `notes`, and an ordered `processing_steps` array
//! (each step has a `step` name, `description`, and optional `parameters` object).
//!
//! ### Full single-spectrum example
//!
//! ```json
//! {
//!   "schema_version": "1.0.0",
//!   "file_type": "single",
//!   "spectrum": {
//!     "id": "chip-5R-4-2",
//!     "metadata": {
//!       "measurement_type": "reflectance",
//!       "date": "2026-04-01",
//!       "title": "Munsell 5R 4/2",
//!       "instrument": { "manufacturer": "Konica Minolta", "model": "CM-700d" },
//!       "measurement_conditions": { "geometry": "d:8", "specular_component": "excluded" }
//!     },
//!     "wavelength_axis": {
//!       "range_nm": { "start": 380, "end": 780, "interval": 10 }
//!     },
//!     "spectral_data": {
//!       "values": [0.048, 0.051, 0.054, 0.058, 0.063],
//!       "scale": "fractional"
//!     },
//!     "color_science": {
//!       "illuminant": "D65",
//!       "cie_observer": "CIE 1931 2 degree",
//!       "results": {
//!         "XYZ": [17.35, 9.12, 1.18],
//!         "xy": [0.629, 0.330],
//!         "Lab": [36.1, 55.7, 37.2]
//!       }
//!     }
//!   }
//! }
//! ```
//!
//! ## Validation
//!
//! [`SpectrumFile::from_path`] and [`SpectrumFile::from_str`] run two validation
//! passes before returning:
//!
//! 1. **Schema validation** — checks required fields, correct types, and that
//!    enum fields (`measurement_type`, `illuminant`, `cie_observer`,
//!    `scale`) contain only allowed values.
//! 2. **Cross-field validation** — checks that `values_nm` and `values` have
//!    equal length; that `uncertainty` (if present) has the same length; that
//!    wavelengths are strictly increasing; that reflectance/transmittance values
//!    lie in \[0, 1\] when `scale` is `"fractional"`; and that a custom
//!    illuminant is accompanied by `illuminant_custom_sd`.
//!
//! Use [`SpectrumFile::from_str_unchecked`] to skip all validation when the
//! source is fully trusted.
//!
//! ## Converting to `Spectrum`
//!
//! With the default `colorimetry` feature enabled, [`SpectrumRecord`] provides four
//! methods to convert spectral data to the fixed 380–780 nm / 1 nm grid that
//! [`colorimetry::spectrum::Spectrum`] requires:
//!
//! | Method | When to use |
//! |---|---|
//! | [`SpectrumRecord::to_spectrum_linear`] | Any axis; simple fallback |
//! | [`SpectrumRecord::to_spectrum_sprague`] | Equidistant (`range_nm`) axis; accurate for coarse regular grids |
//! | [`SpectrumRecord::to_spectrum_smooth`] | Sub-nm pixels with a known slit width; Gaussian-weighted, uses all pixel data |
//! | [`SpectrumRecord::to_spectrum_binned`] | Sub-nm pixels; boxcar average to slit-width bins then Sprague |
//!
//! `to_spectrum_sprague` returns an error if the axis uses `values_nm` (irregular grid).
//! `to_spectrum_smooth` and `to_spectrum_binned` are the right choice when the pixel
//! spacing is much finer than the instrument slit width — they aggregate all pixel data
//! rather than picking isolated samples.

use serde::{Deserialize, Serialize};
use std::path::Path;
use thiserror::Error;

// ─────────────────────────────────────────────────────────────────────────────
// Error type
// ─────────────────────────────────────────────────────────────────────────────

/// All errors that can occur while loading or validating a spectrum file.
#[derive(Debug, Error)]
pub enum SpectrumFileError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("JSON parse error: {0}")]
    Json(#[from] serde_json::Error),

    /// Structural schema violation (wrong type, missing required field,
    /// value not in allowed enum set, etc.)
    #[error("Schema validation failed:\n{0}")]
    SchemaValidation(String),

    /// Cross-field constraint violation (array length mismatch,
    /// non-monotonic wavelengths, value out of physical range, etc.)
    #[error("Cross-field validation failed:\n{0}")]
    CrossFieldValidation(String),
}

pub type Result<T> = std::result::Result<T, SpectrumFileError>;

// ─────────────────────────────────────────────────────────────────────────────
// Top-level file enum
// ─────────────────────────────────────────────────────────────────────────────

/// The top-level structure of a spectrum JSON file.
/// Tagged by `file_type`: either `"single"` or `"batch"`.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "file_type", rename_all = "snake_case")]
pub enum SpectrumFile {
    Single {
        schema_version: String,
        spectrum: SpectrumRecord,
    },
    Batch {
        schema_version: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        batch_metadata: Option<BatchMetadata>,
        spectra: Vec<SpectrumRecord>,
    },
}

impl SpectrumFile {
    // ── Constructors ──────────────────────────────────────────────────────────

    /// Load and fully validate a UV-Vis JSON file from a file path.
    /// Runs structural schema validation then cross-field checks.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let raw = std::fs::read_to_string(path)?;
        Self::from_str(&raw)
    }

    /// Load and fully validate a UV-Vis JSON file from a JSON string.
    pub fn from_str(json: &str) -> Result<Self> {
        // 1. Parse into untyped Value for structural checks
        let value: serde_json::Value = serde_json::from_str(json)?;

        // 2. Structural / schema-level validation
        validate_schema(&value)?;

        // 3. Deserialise into typed structs
        let file: SpectrumFile = serde_json::from_value(value)?;

        // 4. Cross-field validation
        file.validate_cross_fields()?;

        Ok(file)
    }

    /// Deserialise without any validation. Useful when you fully trust the source.
    pub fn from_str_unchecked(json: &str) -> Result<Self> {
        Ok(serde_json::from_str(json)?)
    }

    // ── Accessors ─────────────────────────────────────────────────────────────

    /// Returns all spectra in the file (works for both single and batch).
    pub fn spectra(&self) -> Vec<&SpectrumRecord> {
        match self {
            SpectrumFile::Single { spectrum, .. } => vec![spectrum],
            SpectrumFile::Batch { spectra, .. } => spectra.iter().collect(),
        }
    }

    /// The schema version declared in the file.
    pub fn schema_version(&self) -> &str {
        match self {
            SpectrumFile::Single { schema_version, .. } => schema_version,
            SpectrumFile::Batch { schema_version, .. } => schema_version,
        }
    }

    /// Batch metadata, if this is a batch file.
    pub fn batch_metadata(&self) -> Option<&BatchMetadata> {
        match self {
            SpectrumFile::Batch { batch_metadata, .. } => batch_metadata.as_ref(),
            _ => None,
        }
    }

    // ── Cross-field validation ────────────────────────────────────────────────

    fn validate_cross_fields(&self) -> Result<()> {
        let mut errors: Vec<String> = Vec::new();

        for sp in self.spectra() {
            let id = &sp.id;
            let wl = sp.wavelength_axis.wavelengths_nm();
            let vals = &sp.spectral_data.values;

            // wavelength count == value count
            if wl.len() != vals.len() {
                errors.push(format!(
                    "SpectrumRecord '{id}': wavelength_axis has {} points \
                     but spectral_data.values has {} — must match.",
                    wl.len(),
                    vals.len()
                ));
            }

            // uncertainty length == value count
            if let Some(u) = &sp.spectral_data.uncertainty {
                if u.len() != vals.len() {
                    errors.push(format!(
                        "SpectrumRecord '{id}': spectral_data.uncertainty has {} points \
                         but spectral_data.values has {} — must match.",
                        u.len(),
                        vals.len()
                    ));
                }
            }

            // wavelengths strictly increasing
            if wl.windows(2).any(|w| w[0] >= w[1]) {
                errors.push(format!(
                    "SpectrumRecord '{id}': wavelength_axis is not strictly increasing."
                ));
            }

            // reflectance / transmittance in [0,1] when scale = fractional
            let scale = sp.spectral_data.scale.as_deref().unwrap_or("fractional");
            let is_bounded = matches!(
                sp.metadata.measurement_type,
                MeasurementType::Reflectance | MeasurementType::Transmittance
            );
            if is_bounded && scale == "fractional" {
                let bad: Vec<f64> = vals
                    .iter()
                    .copied()
                    .filter(|&v| !(0.0..=1.0).contains(&v))
                    .collect();
                if !bad.is_empty() {
                    errors.push(format!(
                        "SpectrumRecord '{id}': measurement_type={:?}, scale='fractional' \
                         but {} value(s) fall outside [0,1]. First offender: {}",
                        sp.metadata.measurement_type,
                        bad.len(),
                        bad[0]
                    ));
                }
            }

            // custom illuminant requires illuminant_custom_sd
            if let Some(cs) = &sp.color_science {
                if cs.illuminant.as_deref() == Some("custom") && cs.illuminant_custom_sd.is_none() {
                    errors.push(format!(
                        "SpectrumRecord '{id}': color_science.illuminant is 'custom' \
                         but illuminant_custom_sd is missing."
                    ));
                }
                if let Some(csd) = &cs.illuminant_custom_sd {
                    if csd.wavelengths_nm.len() != csd.values.len() {
                        errors.push(format!(
                            "SpectrumRecord '{id}': illuminant_custom_sd.wavelengths_nm ({}) \
                             and .values ({}) must have equal length.",
                            csd.wavelengths_nm.len(),
                            csd.values.len()
                        ));
                    }
                }
            }
        }

        if errors.is_empty() {
            Ok(())
        } else {
            Err(SpectrumFileError::CrossFieldValidation(errors.join("\n")))
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Structs — mirror the JSON schema
// ─────────────────────────────────────────────────────────────────────────────

/// A single spectral measurement.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectrumRecord {
    pub id: String,
    pub metadata: SpectrumMetadata,
    pub wavelength_axis: WavelengthAxis,
    pub spectral_data: SpectralData,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub color_science: Option<ColorScience>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub provenance: Option<Provenance>,
}

impl SpectrumRecord {
    /// Returns all `(wavelength_nm, value)` pairs.
    pub fn points(&self) -> Vec<(f64, f64)> {
        self.wavelength_axis
            .wavelengths_nm()
            .into_iter()
            .zip(self.spectral_data.values.iter().copied())
            .collect()
    }

    /// Wavelength range as `(min_nm, max_nm)`, or `None` if the axis is empty.
    pub fn wavelength_range_nm(&self) -> Option<(f64, f64)> {
        let wl = self.wavelength_axis.wavelengths_nm();
        Some((*wl.first()?, *wl.last()?))
    }

    /// Number of spectral data points.
    pub fn n_points(&self) -> usize {
        self.spectral_data.values.len()
    }
}

#[cfg(feature = "colorimetry")]
impl SpectrumRecord {
    /// Convert to a [`colorimetry::spectrum::Spectrum`] using linear interpolation.
    ///
    /// Works with both `values_nm` and `range_nm` axes. Linear interpolation is
    /// robust but less accurate than Sprague for smooth spectra on a regular grid.
    pub fn to_spectrum_linear(
        &self,
    ) -> std::result::Result<colorimetry::spectrum::Spectrum, colorimetry::Error> {
        let wl = self.wavelength_axis.wavelengths_nm();
        let vals = self.normalized_values();
        colorimetry::spectrum::Spectrum::linear_interpolate(&wl, &vals)
    }

    /// Convert to a [`colorimetry::spectrum::Spectrum`] using Sprague (5th-order) interpolation.
    ///
    /// Sprague interpolation is more accurate than linear for smooth spectra, but
    /// requires equidistant input. Returns an error if the axis uses `values_nm`
    /// (irregular grid) instead of `range_nm`.
    pub fn to_spectrum_sprague(
        &self,
    ) -> std::result::Result<colorimetry::spectrum::Spectrum, colorimetry::Error> {
        let range = self.wavelength_axis.range_nm.as_ref().ok_or_else(|| {
            colorimetry::Error::ErrorString(
                "Sprague interpolation requires equidistant data (range_nm); \
                 use to_spectrum_linear for values_nm axes"
                    .into(),
            )
        })?;
        let vals = self.normalized_values();
        colorimetry::spectrum::Spectrum::sprague_interpolate([range.start, range.end], &vals)
    }

    /// Convert to a [`colorimetry::spectrum::Spectrum`] using Gaussian-weighted kernel regression.
    ///
    /// For each output wavelength (380–780 nm at 1 nm steps), computes a Gaussian-weighted
    /// average of all input pixels. σ is derived from `spectral_resolution_nm` (FWHM) via
    /// σ = FWHM / 2.355. Pixels further than 4σ from an output wavelength do not contribute.
    ///
    /// This is the recommended method when pixel spacing is much finer than the slit width,
    /// as it uses all available data and the Gaussian weighting matches the typical slit profile.
    /// Output wavelengths with no data within 4σ are skipped; the result is extrapolated by
    /// `linear_interpolate` for any gap at the edges of the 380–780 nm range.
    pub fn to_spectrum_smooth(
        &self,
        spectral_resolution_nm: f64,
    ) -> std::result::Result<colorimetry::spectrum::Spectrum, colorimetry::Error> {
        if spectral_resolution_nm <= 0.0 {
            return Err(colorimetry::Error::ErrorString(
                "spectral_resolution_nm must be positive".into(),
            ));
        }
        let sigma = spectral_resolution_nm / 2.355_f64;
        let cutoff = 4.0 * sigma;

        let vals = self.normalized_values();
        let wls = self.wavelength_axis.wavelengths_nm();

        let mut wl_out: Vec<f64> = Vec::with_capacity(401);
        let mut v_out: Vec<f64> = Vec::with_capacity(401);

        for w in 380_u32..=780 {
            let w_f = w as f64;
            let mut weighted_sum = 0f64;
            let mut weight_total = 0f64;

            for (&wl, &v) in wls.iter().zip(vals.iter()) {
                let d = wl - w_f;
                if d.abs() <= cutoff {
                    let weight = (-0.5 * (d / sigma).powi(2)).exp();
                    weighted_sum += weight * v;
                    weight_total += weight;
                }
            }

            if weight_total > 1e-10 {
                wl_out.push(w_f);
                v_out.push(weighted_sum / weight_total);
            }
        }

        if wl_out.len() < 2 {
            return Err(colorimetry::Error::ErrorString(
                "spectral data does not cover the 380–780 nm output range".into(),
            ));
        }

        colorimetry::spectrum::Spectrum::linear_interpolate(&wl_out, &v_out)
    }

    /// Convert to a [`colorimetry::spectrum::Spectrum`] using boxcar binning followed by
    /// Sprague interpolation.
    ///
    /// Input pixels are averaged into equidistant bins of width `bin_width_nm` spanning
    /// the data's wavelength range. The binned values are then Sprague-interpolated onto
    /// the 380–780 nm / 1 nm output grid. Set `bin_width_nm` to the instrument's optical
    /// resolution. Falls back to linear interpolation if any bins contain no pixels.
    /// Returns an error if fewer than 7 bins are populated (Sprague minimum).
    pub fn to_spectrum_binned(
        &self,
        bin_width_nm: f64,
    ) -> std::result::Result<colorimetry::spectrum::Spectrum, colorimetry::Error> {
        if bin_width_nm <= 0.0 {
            return Err(colorimetry::Error::ErrorString(
                "bin_width_nm must be positive".into(),
            ));
        }

        let vals = self.normalized_values();
        let wls = self.wavelength_axis.wavelengths_nm();

        let (wl_min, wl_max) = match (wls.first(), wls.last()) {
            (Some(&a), Some(&b)) => (a, b),
            _ => return Err(colorimetry::Error::ErrorString("no spectral data".into())),
        };

        let n_bins = ((wl_max - wl_min) / bin_width_nm).round() as usize + 1;
        let mut sums = vec![0f64; n_bins];
        let mut counts = vec![0u32; n_bins];

        for (&wl, &v) in wls.iter().zip(vals.iter()) {
            let idx = ((wl - wl_min) / bin_width_nm).round() as usize;
            if idx < n_bins {
                sums[idx] += v;
                counts[idx] += 1;
            }
        }

        let mut wl_out: Vec<f64> = Vec::with_capacity(n_bins);
        let mut v_out: Vec<f64> = Vec::with_capacity(n_bins);
        for i in 0..n_bins {
            if counts[i] > 0 {
                wl_out.push(wl_min + i as f64 * bin_width_nm);
                v_out.push(sums[i] / counts[i] as f64);
            }
        }

        if v_out.len() < 7 {
            return Err(colorimetry::Error::ProvideAtLeastNValues(7));
        }

        if wl_out.len() == n_bins {
            // All bins populated — use Sprague on the equidistant grid
            colorimetry::spectrum::Spectrum::sprague_interpolate([wl_min, wl_max], &v_out)
        } else {
            // Gaps present — fall back to linear on the populated subset
            colorimetry::spectrum::Spectrum::linear_interpolate(&wl_out, &v_out)
        }
    }

    fn normalized_values(&self) -> Vec<f64> {
        let mut vals = self.spectral_data.values.clone();
        if self.spectral_data.scale.as_deref() == Some("percent") {
            vals.iter_mut().for_each(|v| *v /= 100.0);
        }
        vals
    }
}

/// Descriptive metadata for one spectrum.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectrumMetadata {
    pub measurement_type: MeasurementType,
    /// ISO 8601 date (YYYY-MM-DD).
    pub date: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub title: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sample_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub time: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub operator: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub instrument: Option<Instrument>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub measurement_conditions: Option<MeasurementConditions>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tags: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub custom: Option<serde_json::Value>,
}

/// The physical quantity measured.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MeasurementType {
    Reflectance,
    Transmittance,
    Absorbance,
    Radiance,
    Irradiance,
}

/// Minimal instrument identification.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Instrument {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub manufacturer: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub model: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub serial_number: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub detector_type: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub light_source: Option<String>,
}

/// Physical conditions under which the measurement was made.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeasurementConditions {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub integration_time_ms: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub averaging: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature_celsius: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub geometry: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub specular_component: Option<SpecularComponent>,
    /// Optical (spectral) resolution of the instrument in nm, typically the FWHM of the slit function.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub spectral_resolution_nm: Option<f64>,
}

/// Whether the specular component is included or excluded.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpecularComponent {
    Included,
    Excluded,
    #[serde(rename = "not applicable")]
    NotApplicable,
}

/// The wavelength axis of the measurement. All values are in nm.
///
/// Exactly one of `values_nm` or `range_nm` must be present.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WavelengthAxis {
    /// Explicit wavelength list in nm. Use for irregular grids.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub values_nm: Option<Vec<f64>>,
    /// Evenly-spaced grid descriptor. Use for regular grids.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub range_nm: Option<WavelengthRange>,
}

impl WavelengthAxis {
    /// Returns the wavelength values in nm, expanding `range_nm` if that variant is used.
    pub fn wavelengths_nm(&self) -> Vec<f64> {
        if let Some(v) = &self.values_nm {
            v.clone()
        } else if let Some(r) = &self.range_nm {
            r.expand()
        } else {
            vec![]
        }
    }
}

/// Evenly-spaced wavelength grid defined by start, end, and interval (all in nm).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WavelengthRange {
    pub start: f64,
    pub end: f64,
    pub interval: f64,
}

impl WavelengthRange {
    /// Expands the range into an explicit list of wavelength values in nm.
    pub fn expand(&self) -> Vec<f64> {
        let n = ((self.end - self.start) / self.interval).round() as usize + 1;
        (0..n)
            .map(|i| self.start + i as f64 * self.interval)
            .collect()
    }
}

/// The measured spectral values.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectralData {
    pub values: Vec<f64>,
    /// Optional per-point uncertainty (1 standard deviation), same length as `values`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub uncertainty: Option<Vec<f64>>,
    /// `"fractional"` (0–1) or `"percent"` (0–100). Defaults to `"fractional"`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scale: Option<String>,
}

/// Metadata required for CIE colorimetry and color-science calculations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColorScience {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub illuminant: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub illuminant_custom_sd: Option<CustomIlluminantSd>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cie_observer: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub white_reference: Option<WhiteReference>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub results: Option<ColorScienceResults>,
}

/// Pre-computed colorimetric results derived from the spectral data.
///
/// All fields are optional and informational — the spectral data is always the authoritative
/// source. Any subset may be present.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColorScienceResults {
    /// CIE tristimulus values [X, Y, Z].
    #[serde(rename = "XYZ", skip_serializing_if = "Option::is_none")]
    pub xyz: Option<[f64; 3]>,
    /// CIE 1931 chromaticity coordinates [x, y].
    #[serde(skip_serializing_if = "Option::is_none")]
    pub xy: Option<[f64; 2]>,
    /// CIE 1976 UCS chromaticity coordinates [u′, v′].
    #[serde(skip_serializing_if = "Option::is_none")]
    pub uv_prime: Option<[f64; 2]>,
    /// CIELAB coordinates [L*, a*, b*].
    #[serde(rename = "Lab", skip_serializing_if = "Option::is_none")]
    pub lab: Option<[f64; 3]>,
    /// Correlated color temperature in Kelvin.
    #[serde(rename = "CCT_K", skip_serializing_if = "Option::is_none")]
    pub cct_k: Option<f64>,
    /// Distance from the Planckian locus (signed) in the CIE 1960 UCS.
    #[serde(rename = "Duv", skip_serializing_if = "Option::is_none")]
    pub duv: Option<f64>,
}

/// Spectral power distribution for a custom illuminant.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CustomIlluminantSd {
    pub wavelengths_nm: Vec<f64>,
    pub values: Vec<f64>,
}

/// White reference / calibration tile.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WhiteReference {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub manufacturer: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub serial_number: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub calibration_date: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reference_values: Option<Vec<f64>>,
}

/// Processing history and software trail.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Provenance {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub software: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub software_version: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_file: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_format: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub processing_steps: Option<Vec<ProcessingStep>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub notes: Option<String>,
}

/// A single processing step applied to the raw data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingStep {
    pub step: String,
    pub description: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub parameters: Option<serde_json::Value>,
}

/// Optional metadata common to all spectra in a batch file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchMetadata {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub title: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub operator: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub date: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub instrument: Option<Instrument>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub measurement_conditions: Option<MeasurementConditions>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Structural schema validator (pure serde_json, no external crate)
// ─────────────────────────────────────────────────────────────────────────────
//
// Checks that are enforced here (equivalent to JSON Schema):
//   - Required top-level fields present and correct type
//   - file_type is "single" or "batch"
//   - schema_version matches semver pattern
//   - Each spectrum has required fields (id, metadata, wavelength_axis, spectral_data)
//   - measurement_type is one of the allowed enum values
//   - scale, if present, is "fractional" or "percent"
//   - values_nm has at least 2 entries
//   - All numeric arrays contain only numbers

const ALLOWED_MEASUREMENT_TYPES: &[&str] = &[
    "reflectance",
    "transmittance",
    "absorbance",
    "radiance",
    "irradiance",
];

const ALLOWED_ILLUMINANTS: &[&str] = &[
    "D65", "D50", "D55", "D75", "A", "B", "C", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8",
    "F9", "F10", "F11", "F12", "LED-B1", "LED-B2", "LED-B3", "LED-B4", "LED-B5", "LED-BH1",
    "LED-RGB1", "LED-V1", "LED-V2", "custom",
];

const ALLOWED_OBSERVERS: &[&str] = &[
    "CIE 1931 2 degree",
    "CIE 1964 10 degree",
    "CIE 2015 2 degree",
    "CIE 2015 10 degree",
];

fn validate_schema(v: &serde_json::Value) -> Result<()> {
    let mut errors: Vec<String> = Vec::new();

    // Top-level must be an object
    let obj = match v.as_object() {
        Some(o) => o,
        None => {
            return Err(SpectrumFileError::SchemaValidation(
                "Root value must be a JSON object.".into(),
            ))
        }
    };

    // schema_version: required, string, semver-ish
    match obj.get("schema_version") {
        None => errors.push("Missing required field: schema_version".into()),
        Some(sv) => {
            if !sv.is_string() {
                errors.push("schema_version must be a string".into());
            } else {
                let s = sv.as_str().unwrap();
                let parts: Vec<&str> = s.split('.').collect();
                if parts.len() != 3 || parts.iter().any(|p| p.parse::<u32>().is_err()) {
                    errors.push(format!(
                        "schema_version '{s}' does not look like semver (e.g. 1.0.0)"
                    ));
                }
            }
        }
    }

    // file_type: required, "single" or "batch"
    let file_type = match obj.get("file_type") {
        None => {
            errors.push("Missing required field: file_type".into());
            None
        }
        Some(ft) => match ft.as_str() {
            Some(s @ "single") | Some(s @ "batch") => Some(s.to_string()),
            Some(other) => {
                errors.push(format!(
                    "file_type must be 'single' or 'batch', got '{other}'"
                ));
                None
            }
            None => {
                errors.push("file_type must be a string".into());
                None
            }
        },
    };

    match file_type.as_deref() {
        Some("single") => {
            match obj.get("spectrum") {
                None => errors.push("Single file must have a 'spectrum' field".into()),
                Some(sp) => validate_spectrum(sp, "spectrum", &mut errors),
            }
            if obj.contains_key("spectra") {
                errors.push(
                    "Single file must not have a 'spectra' array (use file_type='batch')".into(),
                );
            }
        }
        Some("batch") => {
            match obj.get("spectra") {
                None => errors.push("Batch file must have a 'spectra' array".into()),
                Some(arr) => match arr.as_array() {
                    None => errors.push("'spectra' must be an array".into()),
                    Some(items) => {
                        if items.is_empty() {
                            errors.push("'spectra' array must not be empty".into());
                        }
                        for (i, sp) in items.iter().enumerate() {
                            validate_spectrum(sp, &format!("spectra[{i}]"), &mut errors);
                        }
                    }
                },
            }
            if obj.contains_key("spectrum") {
                errors.push("Batch file must not have a 'spectrum' field".into());
            }
        }
        _ => {} // already reported
    }

    if errors.is_empty() {
        Ok(())
    } else {
        Err(SpectrumFileError::SchemaValidation(errors.join("\n")))
    }
}

fn validate_spectrum(v: &serde_json::Value, path: &str, errors: &mut Vec<String>) {
    let obj = match v.as_object() {
        Some(o) => o,
        None => {
            errors.push(format!("{path}: must be an object"));
            return;
        }
    };

    // id: required string
    require_string(obj, "id", path, errors);

    // metadata: required object
    if let Some(meta) = require_object(obj, "metadata", path, errors) {
        validate_metadata(meta, &format!("{path}.metadata"), errors);
    }

    // wavelength_axis: required object
    if let Some(wa) = require_object(obj, "wavelength_axis", path, errors) {
        validate_wavelength_axis(wa, &format!("{path}.wavelength_axis"), errors);
    }

    // spectral_data: required object
    if let Some(sd) = require_object(obj, "spectral_data", path, errors) {
        validate_spectral_data(sd, &format!("{path}.spectral_data"), errors);
    }

    // color_science: optional
    if let Some(cs) = obj.get("color_science") {
        if let Some(cso) = cs.as_object() {
            validate_color_science(cso, &format!("{path}.color_science"), errors);
        } else {
            errors.push(format!("{path}.color_science must be an object"));
        }
    }
}

fn validate_metadata(
    obj: &serde_json::Map<String, serde_json::Value>,
    path: &str,
    errors: &mut Vec<String>,
) {
    // measurement_type: required, enum
    match obj.get("measurement_type") {
        None => errors.push(format!("{path}: missing required field 'measurement_type'")),
        Some(mt) => match mt.as_str() {
            None => errors.push(format!("{path}.measurement_type must be a string")),
            Some(s) if !ALLOWED_MEASUREMENT_TYPES.contains(&s) => errors.push(format!(
                "{path}.measurement_type '{s}' is not allowed. Must be one of: {}",
                ALLOWED_MEASUREMENT_TYPES.join(", ")
            )),
            _ => {}
        },
    }
    // date: required string
    require_string(obj, "date", path, errors);
}

fn validate_wavelength_axis(
    obj: &serde_json::Map<String, serde_json::Value>,
    path: &str,
    errors: &mut Vec<String>,
) {
    let has_values = obj.contains_key("values_nm");
    let has_range = obj.contains_key("range_nm");

    match (has_values, has_range) {
        (false, false) => {
            errors.push(format!(
                "{path}: exactly one of 'values_nm' or 'range_nm' must be present (neither found)"
            ));
            return;
        }
        (true, true) => {
            errors.push(format!(
                "{path}: exactly one of 'values_nm' or 'range_nm' must be present (both found)"
            ));
            return;
        }
        _ => {}
    }

    if has_values {
        match obj.get("values_nm").and_then(|v| v.as_array()) {
            None => errors.push(format!("{path}.values_nm must be an array")),
            Some(items) => {
                if items.len() < 2 {
                    errors.push(format!("{path}.values_nm must have at least 2 elements"));
                }
                if items.iter().any(|x| !x.is_number()) {
                    errors.push(format!("{path}.values_nm must contain only numbers"));
                }
            }
        }
    } else {
        match obj.get("range_nm").and_then(|v| v.as_object()) {
            None => errors.push(format!("{path}.range_nm must be an object")),
            Some(r) => {
                for field in ["start", "end", "interval"] {
                    match r.get(field) {
                        None => errors
                            .push(format!("{path}.range_nm: missing required field '{field}'")),
                        Some(v) if !v.is_number() => {
                            errors.push(format!("{path}.range_nm.{field} must be a number"))
                        }
                        _ => {}
                    }
                }
                if let Some(iv) = r.get("interval").and_then(|v| v.as_f64()) {
                    if iv <= 0.0 {
                        errors.push(format!("{path}.range_nm.interval must be positive"));
                    }
                }
            }
        }
    }
}

fn validate_spectral_data(
    obj: &serde_json::Map<String, serde_json::Value>,
    path: &str,
    errors: &mut Vec<String>,
) {
    // values: required, array of numbers, min 2
    match obj.get("values") {
        None => errors.push(format!("{path}: missing required field 'values'")),
        Some(arr) => match arr.as_array() {
            None => errors.push(format!("{path}.values must be an array")),
            Some(items) => {
                if items.len() < 2 {
                    errors.push(format!("{path}.values must have at least 2 elements"));
                }
                if items.iter().any(|x| !x.is_number()) {
                    errors.push(format!("{path}.values must contain only numbers"));
                }
            }
        },
    }

    // uncertainty: optional array of non-negative numbers
    if let Some(unc) = obj.get("uncertainty") {
        match unc.as_array() {
            None => errors.push(format!("{path}.uncertainty must be an array")),
            Some(items) => {
                if items.iter().any(|x| !x.is_number()) {
                    errors.push(format!("{path}.uncertainty must contain only numbers"));
                } else if items.iter().any(|x| x.as_f64().unwrap_or(0.0) < 0.0) {
                    errors.push(format!("{path}.uncertainty values must be non-negative"));
                }
            }
        }
    }

    // scale: optional, enum
    if let Some(sc) = obj.get("scale") {
        match sc.as_str() {
            None => errors.push(format!("{path}.scale must be a string")),
            Some(s) if s != "fractional" && s != "percent" => errors.push(format!(
                "{path}.scale must be 'fractional' or 'percent', got '{s}'"
            )),
            _ => {}
        }
    }
}

fn validate_color_science(
    obj: &serde_json::Map<String, serde_json::Value>,
    path: &str,
    errors: &mut Vec<String>,
) {
    // illuminant: optional, enum
    if let Some(il) = obj.get("illuminant") {
        match il.as_str() {
            None => errors.push(format!("{path}.illuminant must be a string")),
            Some(s) if !ALLOWED_ILLUMINANTS.contains(&s) => errors.push(format!(
                "{path}.illuminant '{s}' is not a recognised CIE illuminant"
            )),
            _ => {}
        }
    }

    // cie_observer: optional, enum
    if let Some(obs) = obj.get("cie_observer") {
        match obs.as_str() {
            None => errors.push(format!("{path}.cie_observer must be a string")),
            Some(s) if !ALLOWED_OBSERVERS.contains(&s) => errors.push(format!(
                "{path}.cie_observer '{s}' not recognised. Must be one of: {}",
                ALLOWED_OBSERVERS.join(", ")
            )),
            _ => {}
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn require_string(
    obj: &serde_json::Map<String, serde_json::Value>,
    key: &str,
    path: &str,
    errors: &mut Vec<String>,
) {
    match obj.get(key) {
        None => errors.push(format!("{path}: missing required field '{key}'")),
        Some(v) if !v.is_string() => errors.push(format!("{path}.{key} must be a string")),
        _ => {}
    }
}

fn require_object<'a>(
    obj: &'a serde_json::Map<String, serde_json::Value>,
    key: &str,
    path: &str,
    errors: &mut Vec<String>,
) -> Option<&'a serde_json::Map<String, serde_json::Value>> {
    match obj.get(key) {
        None => {
            errors.push(format!("{path}: missing required field '{key}'"));
            None
        }
        Some(v) => match v.as_object() {
            None => {
                errors.push(format!("{path}.{key} must be an object"));
                None
            }
            Some(o) => Some(o),
        },
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_single(mtype: &str, wls: &[f64], vals: &[f64]) -> String {
        let wl_s: Vec<String> = wls.iter().map(|w| w.to_string()).collect();
        let v_s: Vec<String> = vals.iter().map(|v| v.to_string()).collect();
        format!(
            r#"{{"schema_version":"1.0.0","file_type":"single","spectrum":{{"id":"t1",
            "metadata":{{"measurement_type":"{mtype}","date":"2026-04-29"}},
            "wavelength_axis":{{"values_nm":[{wl}]}},
            "spectral_data":{{"values":[{v}]}}}}}}"#,
            mtype = mtype,
            wl = wl_s.join(","),
            v = v_s.join(","),
        )
    }

    fn wls_41() -> Vec<f64> {
        (0..41).map(|i| 380.0 + i as f64 * 10.0).collect()
    }
    fn vals_41() -> Vec<f64> {
        (0..41).map(|i| i as f64 / 100.0).collect()
    }

    #[test]
    fn valid_single_spectrum() {
        let file =
            SpectrumFile::from_str(&make_single("reflectance", &wls_41(), &vals_41())).unwrap();
        let spectra = file.spectra();
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].n_points(), 41);
        assert_eq!(file.schema_version(), "1.0.0");
    }

    #[test]
    fn valid_batch_file() {
        let json = r#"{"schema_version":"1.0.0","file_type":"batch","spectra":[
            {"id":"a","metadata":{"measurement_type":"reflectance","date":"2026-04-29"},
             "wavelength_axis":{"values_nm":[380,390,400]},
             "spectral_data":{"values":[0.1,0.2,0.3]}},
            {"id":"b","metadata":{"measurement_type":"transmittance","date":"2026-04-29"},
             "wavelength_axis":{"values_nm":[380,390,400]},
             "spectral_data":{"values":[0.5,0.6,0.7]}}
        ]}"#;
        let file = SpectrumFile::from_str(json).unwrap();
        assert_eq!(file.spectra().len(), 2);
    }

    #[test]
    fn missing_measurement_type_is_schema_error() {
        let json = r#"{"schema_version":"1.0.0","file_type":"single","spectrum":{"id":"x",
            "metadata":{"date":"2026-04-29"},
            "wavelength_axis":{"values_nm":[380,390,400]},
            "spectral_data":{"values":[0.1,0.2,0.3]}}}"#;
        assert!(matches!(
            SpectrumFile::from_str(json),
            Err(SpectrumFileError::SchemaValidation(_))
        ));
    }

    #[test]
    fn invalid_measurement_type_is_schema_error() {
        let json = make_single("fluorescence", &[380.0, 390.0], &[0.1, 0.2]);
        assert!(matches!(
            SpectrumFile::from_str(&json),
            Err(SpectrumFileError::SchemaValidation(_))
        ));
    }

    #[test]
    fn wavelength_value_length_mismatch() {
        let wls = vec![380.0, 390.0, 400.0];
        let vals = vec![0.1, 0.2]; // too short
        assert!(matches!(
            SpectrumFile::from_str(&make_single("reflectance", &wls, &vals)),
            Err(SpectrumFileError::CrossFieldValidation(_))
        ));
    }

    #[test]
    fn non_monotonic_wavelengths() {
        let wls = vec![380.0, 370.0, 400.0];
        let vals = vec![0.1, 0.2, 0.3];
        assert!(matches!(
            SpectrumFile::from_str(&make_single("reflectance", &wls, &vals)),
            Err(SpectrumFileError::CrossFieldValidation(_))
        ));
    }

    #[test]
    fn reflectance_out_of_range() {
        let wls = vec![380.0, 390.0, 400.0];
        let vals = vec![0.1, 1.5, 0.3];
        assert!(matches!(
            SpectrumFile::from_str(&make_single("reflectance", &wls, &vals)),
            Err(SpectrumFileError::CrossFieldValidation(_))
        ));
    }

    #[test]
    fn absorbance_above_one_is_ok() {
        // Absorbance is not bounded by [0,1]
        let wls = vec![380.0, 390.0, 400.0];
        let vals = vec![0.1, 1.8, 2.5];
        assert!(SpectrumFile::from_str(&make_single("absorbance", &wls, &vals)).is_ok());
    }

    #[test]
    fn custom_illuminant_missing_sd() {
        let json = r#"{"schema_version":"1.0.0","file_type":"single","spectrum":{"id":"x",
            "metadata":{"measurement_type":"reflectance","date":"2026-04-29"},
            "wavelength_axis":{"values_nm":[380,390,400]},
            "spectral_data":{"values":[0.1,0.2,0.3]},
            "color_science":{"illuminant":"custom"}}}"#;
        assert!(matches!(
            SpectrumFile::from_str(json),
            Err(SpectrumFileError::CrossFieldValidation(_))
        ));
    }

    #[test]
    fn points_iterator_correct() {
        let wls = vec![380.0, 390.0, 400.0];
        let vals = vec![0.1, 0.2, 0.3];
        let file = SpectrumFile::from_str(&make_single("reflectance", &wls, &vals)).unwrap();
        let pts = file.spectra()[0].points();
        assert_eq!(pts, vec![(380.0, 0.1), (390.0, 0.2), (400.0, 0.3)]);
    }

    #[test]
    fn wavelength_range_accessor() {
        let file =
            SpectrumFile::from_str(&make_single("reflectance", &wls_41(), &vals_41())).unwrap();
        assert_eq!(
            file.spectra()[0].wavelength_range_nm(),
            Some((380.0, 780.0))
        );
    }

    #[test]
    fn invalid_scale_value() {
        let json = r#"{"schema_version":"1.0.0","file_type":"single","spectrum":{"id":"x",
            "metadata":{"measurement_type":"reflectance","date":"2026-04-29"},
            "wavelength_axis":{"values_nm":[380,390,400]},
            "spectral_data":{"values":[0.1,0.2,0.3],"scale":"ratio"}}}"#;
        assert!(matches!(
            SpectrumFile::from_str(json),
            Err(SpectrumFileError::SchemaValidation(_))
        ));
    }

    // ── WavelengthAxis and WavelengthRange unit tests ─────────────────────────

    #[test]
    fn wavelength_axis_values_nm_variant() {
        let axis = WavelengthAxis {
            values_nm: Some(vec![380.0, 450.0, 550.0, 700.0]),
            range_nm: None,
        };
        assert_eq!(axis.wavelengths_nm(), vec![380.0, 450.0, 550.0, 700.0]);
    }

    #[test]
    fn wavelength_axis_range_nm_variant() {
        let axis = WavelengthAxis {
            values_nm: None,
            range_nm: Some(WavelengthRange {
                start: 380.0,
                end: 400.0,
                interval: 10.0,
            }),
        };
        let wls = axis.wavelengths_nm();
        assert_eq!(wls.len(), 3);
        assert!((wls[0] - 380.0).abs() < 1e-10);
        assert!((wls[1] - 390.0).abs() < 1e-10);
        assert!((wls[2] - 400.0).abs() < 1e-10);
    }

    #[test]
    fn wavelength_range_expand_direct() {
        let r = WavelengthRange {
            start: 380.0,
            end: 780.0,
            interval: 10.0,
        };
        let wls = r.expand();
        assert_eq!(wls.len(), 41);
        assert!((wls[0] - 380.0).abs() < 1e-10);
        assert!((wls[40] - 780.0).abs() < 1e-10);
    }

    // ── Cross-field validation edge cases ─────────────────────────────────────

    #[test]
    fn uncertainty_length_mismatch_is_error() {
        let json = r#"{
            "schema_version": "1.0.0",
            "file_type": "single",
            "spectrum": {
                "id": "x",
                "metadata": {"measurement_type": "reflectance", "date": "2026-04-29"},
                "wavelength_axis": {"values_nm": [380, 390, 400]},
                "spectral_data": {"values": [0.1, 0.2, 0.3], "uncertainty": [0.01, 0.01]}
            }
        }"#;
        assert!(matches!(
            SpectrumFile::from_str(json),
            Err(SpectrumFileError::CrossFieldValidation(_))
        ));
    }

    #[test]
    fn illuminant_custom_sd_length_mismatch_is_error() {
        let json = r#"{
            "schema_version": "1.0.0",
            "file_type": "single",
            "spectrum": {
                "id": "x",
                "metadata": {"measurement_type": "reflectance", "date": "2026-04-29"},
                "wavelength_axis": {"values_nm": [380, 390, 400]},
                "spectral_data": {"values": [0.1, 0.2, 0.3]},
                "color_science": {
                    "illuminant": "custom",
                    "illuminant_custom_sd": {
                        "wavelengths_nm": [380, 390, 400],
                        "values": [1.0, 1.1]
                    }
                }
            }
        }"#;
        assert!(matches!(
            SpectrumFile::from_str(json),
            Err(SpectrumFileError::CrossFieldValidation(_))
        ));
    }

    // ── from_path and from_str_unchecked ──────────────────────────────────────

    #[test]
    fn from_path_loads_single_example() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/scripts/example_single.json");
        let file = SpectrumFile::from_path(path).unwrap();
        assert_eq!(file.spectra().len(), 1);
        assert_eq!(file.spectra()[0].id, "sample-001");
    }

    #[test]
    fn from_path_loads_batch_example() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/scripts/example_batch.json");
        let file = SpectrumFile::from_path(path).unwrap();
        assert_eq!(file.spectra().len(), 2);
    }

    #[test]
    fn from_str_unchecked_skips_cross_field_validation() {
        // 3 wavelengths but only 2 values — cross-field check rejects this, unchecked accepts it
        let json = r#"{
            "schema_version": "1.0.0",
            "file_type": "single",
            "spectrum": {
                "id": "x",
                "metadata": {"measurement_type": "reflectance", "date": "2026-04-29"},
                "wavelength_axis": {"values_nm": [380, 390, 400]},
                "spectral_data": {"values": [0.1, 0.2]}
            }
        }"#;
        assert!(SpectrumFile::from_str_unchecked(json).is_ok());
        assert!(matches!(
            SpectrumFile::from_str(json),
            Err(SpectrumFileError::CrossFieldValidation(_))
        ));
    }

    // ── batch_metadata accessor ───────────────────────────────────────────────

    #[test]
    fn batch_metadata_fields_accessible() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/scripts/example_batch.json");
        let file = SpectrumFile::from_path(path).unwrap();
        let meta = file
            .batch_metadata()
            .expect("batch file must have metadata");
        assert_eq!(
            meta.title.as_deref(),
            Some("Ceramic tile color survey - April 2026")
        );
        assert_eq!(meta.operator.as_deref(), Some("J. Smith"));
    }

    #[test]
    fn batch_metadata_returns_none_for_single_file() {
        let file =
            SpectrumFile::from_str(&make_single("reflectance", &wls_41(), &vals_41())).unwrap();
        assert!(file.batch_metadata().is_none());
    }

    // ── Conversion methods ────────────────────────────────────────────────────

    #[cfg(feature = "colorimetry")]
    mod conversion_tests {
        use super::*;

        fn vals_41() -> Vec<f64> {
            (0..41).map(|i| i as f64 / 40.0).collect()
        }

        fn make_range_record(vals: Vec<f64>) -> SpectrumRecord {
            SpectrumRecord {
                id: "test".into(),
                metadata: SpectrumMetadata {
                    measurement_type: MeasurementType::Reflectance,
                    date: "2026-04-29".into(),
                    title: None,
                    description: None,
                    sample_id: None,
                    time: None,
                    operator: None,
                    instrument: None,
                    measurement_conditions: None,
                    tags: None,
                    custom: None,
                },
                wavelength_axis: WavelengthAxis {
                    values_nm: None,
                    range_nm: Some(WavelengthRange {
                        start: 380.0,
                        end: 780.0,
                        interval: 10.0,
                    }),
                },
                spectral_data: SpectralData {
                    values: vals,
                    uncertainty: None,
                    scale: None,
                },
                color_science: None,
                provenance: None,
            }
        }

        fn make_values_record(wls: Vec<f64>, vals: Vec<f64>) -> SpectrumRecord {
            SpectrumRecord {
                id: "test".into(),
                metadata: SpectrumMetadata {
                    measurement_type: MeasurementType::Reflectance,
                    date: "2026-04-29".into(),
                    title: None,
                    description: None,
                    sample_id: None,
                    time: None,
                    operator: None,
                    instrument: None,
                    measurement_conditions: None,
                    tags: None,
                    custom: None,
                },
                wavelength_axis: WavelengthAxis {
                    values_nm: Some(wls),
                    range_nm: None,
                },
                spectral_data: SpectralData {
                    values: vals,
                    uncertainty: None,
                    scale: None,
                },
                color_science: None,
                provenance: None,
            }
        }

        #[test]
        fn to_spectrum_linear_range_nm() {
            let rec = make_range_record(vals_41());
            assert!(rec.to_spectrum_linear().is_ok());
        }

        #[test]
        fn to_spectrum_linear_values_nm() {
            let wls: Vec<f64> = (0..41).map(|i| 380.0 + i as f64 * 10.0).collect();
            let rec = make_values_record(wls, vals_41());
            assert!(rec.to_spectrum_linear().is_ok());
        }

        #[test]
        fn to_spectrum_sprague_happy_path() {
            let rec = make_range_record(vals_41());
            assert!(rec.to_spectrum_sprague().is_ok());
        }

        #[test]
        fn to_spectrum_sprague_error_on_values_nm() {
            let wls: Vec<f64> = (0..41).map(|i| 380.0 + i as f64 * 10.0).collect();
            let rec = make_values_record(wls, vals_41());
            assert!(rec.to_spectrum_sprague().is_err());
        }

        #[test]
        fn to_spectrum_smooth_happy_path() {
            let rec = make_range_record(vals_41());
            // FWHM=15 nm → cutoff≈25.5 nm, so every 380–780 nm output point has coverage
            // from the 10 nm spaced input pixels
            assert!(rec.to_spectrum_smooth(15.0).is_ok());
        }

        #[test]
        fn to_spectrum_smooth_negative_resolution_is_error() {
            let rec = make_range_record(vals_41());
            assert!(rec.to_spectrum_smooth(-1.0).is_err());
        }

        #[test]
        fn to_spectrum_smooth_zero_resolution_is_error() {
            let rec = make_range_record(vals_41());
            assert!(rec.to_spectrum_smooth(0.0).is_err());
        }

        #[test]
        fn to_spectrum_binned_happy_path() {
            let rec = make_range_record(vals_41());
            // bin_width matches pixel spacing → all 41 bins populated, Sprague path taken
            assert!(rec.to_spectrum_binned(10.0).is_ok());
        }

        #[test]
        fn to_spectrum_binned_negative_bin_width_is_error() {
            let rec = make_range_record(vals_41());
            assert!(rec.to_spectrum_binned(-1.0).is_err());
        }

        #[test]
        fn to_spectrum_binned_too_few_bins_is_error() {
            // 5 bins populated — fewer than the 7-point Sprague minimum
            let wls: Vec<f64> = (0..5).map(|i| 380.0 + i as f64 * 50.0).collect();
            let vals = vec![0.1, 0.2, 0.3, 0.4, 0.5];
            let rec = make_values_record(wls, vals);
            assert!(rec.to_spectrum_binned(50.0).is_err());
        }

        #[test]
        fn normalized_values_percent_scale() {
            // Values 0..100 (percent) should be divided by 100 before interpolation
            let percent_vals: Vec<f64> = (0..41).map(|i| i as f64 * 2.5).collect();
            let mut rec = make_range_record(percent_vals);
            rec.spectral_data.scale = Some("percent".into());
            assert!(rec.to_spectrum_linear().is_ok());
        }
    }
}
