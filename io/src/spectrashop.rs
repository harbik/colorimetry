//! Parser for the SpectraShop tab-separated text export format.
//!
//! ## Format licensing
//!
//! The SpectraShop format is proprietary to Robin Myers Imaging. The format
//! specification is published at:
//! <https://www.chromaxion.com/spectral_library/SpectraShop_Import-Export_Format.pdf>
//! specifically to allow third-party readers and writers.
//!
//! ## Data-file licensing (Chromaxion Spectral Library)
//!
//! Spectral data files distributed by Robin Myers Imaging / Chromaxion are free
//! for personal, scientific, and teaching use. Redistribution requires attribution
//! to Chromaxion.com or Robin Myers. Commercial sale requires express written
//! permission. All data © Robin D. Myers, all rights reserved worldwide.
//! See <https://www.chromaxion.com/spectral-library.php> for the full terms.

use super::ALLOWED_ILLUMINANTS;
use super::*;

// Accumulated file-level metadata parsed from a SpectraShop header block.
#[derive(Default, Clone)]
struct SsMeta {
    file_descriptor: Option<String>,
    spectrum_type: Option<String>,
    created: Option<String>,
    originator: Option<String>,
    acquire_note: Option<String>,
    note: Option<String>,
    illuminant: Option<String>,
    observer: Option<String>,
    instrumentation: Option<String>,
    instrument_serial: Option<String>,
    measurement_geometry: Option<String>,
    measurement_source: Option<String>,
    measurement_aperture: Option<String>,
    measurement_filter: Option<String>,
    sample_backing: Option<String>,
    surface: Option<String>,
    manufacturer: Option<String>,
    material: Option<String>,
    model_num: Option<String>,
    prod_date: Option<String>,
    serial: Option<String>,
    nmeasure: Option<u32>,
    sw_version: Option<String>,
}

pub(super) fn ss_parse(input: &str, source_file: Option<&str>) -> Result<SpectrumFile> {
    let mut meta = SsMeta::default();
    let mut in_data_format = false;
    let mut in_data = false;

    // A SpectraShop file can contain multiple BEGIN_DATA_FORMAT/BEGIN_DATA block
    // pairs (e.g., different wavelength ranges in the same file).  Each pair is
    // processed independently when END_DATA is encountered.
    let mut current_format: Vec<String> = Vec::new();
    let mut current_data: Vec<String> = Vec::new();
    let mut all_records: Vec<SpectrumRecord> = Vec::new();
    let mut global_idx: usize = 0;
    let mut saw_data_block = false;
    let mut any_records_parsed = false;

    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // Section markers (case-insensitive exact match on the trimmed line)
        let upper = trimmed.to_uppercase();
        if upper == "BEGIN_DATA_FORMAT" {
            current_format.clear();
            in_data_format = true;
            continue;
        }
        if upper == "END_DATA_FORMAT" {
            in_data_format = false;
            continue;
        }
        if upper == "BEGIN_DATA" {
            current_data.clear();
            in_data = true;
            saw_data_block = true;
            continue;
        }
        if upper == "END_DATA" {
            in_data = false;
            // Process the completed format+data section.
            if !current_format.is_empty() && !current_data.is_empty() {
                any_records_parsed = true;
                let n_fields = current_format.len();
                if current_data.len() % n_fields != 0 {
                    return Err(SpectrumFileError::SchemaValidation(format!(
                        "SpectraShop: data token count ({}) is not a multiple of \
                         field count ({}) in section ending near record {}",
                        current_data.len(),
                        n_fields,
                        global_idx
                    )));
                }
                let n_section = current_data.len() / n_fields;
                for r in 0..n_section {
                    let row = &current_data[r * n_fields..(r + 1) * n_fields];
                    all_records.push(ss_build_record(
                        &current_format,
                        row,
                        &meta,
                        global_idx,
                        source_file,
                    )?);
                    global_idx += 1;
                }
            }
            continue;
        }

        if in_data_format {
            for field in trimmed.split_whitespace() {
                current_format.push(field.to_uppercase());
            }
            continue;
        }
        if in_data {
            ss_tokenize_into(trimmed, &mut current_data);
            continue;
        }

        // Header / inter-section metadata line: KEYWORD<tab>value.
        // These can appear before the first section or between sections.
        if let Some((key, value)) = ss_split_kv(trimmed) {
            ss_apply_kv(&mut meta, &key, &value);
        }
    }

    if !saw_data_block {
        return Err(SpectrumFileError::SchemaValidation(
            "SpectraShop: no BEGIN_DATA/END_DATA block found".into(),
        ));
    }
    if !any_records_parsed {
        return Err(SpectrumFileError::SchemaValidation(
            "SpectraShop: BEGIN_DATA/END_DATA block present but contained no parseable records \
             (check that BEGIN_DATA_FORMAT is also present and non-empty)"
                .into(),
        ));
    }
    if all_records.is_empty() {
        return Err(SpectrumFileError::SchemaValidation(
            "SpectraShop: file contains no spectral records".into(),
        ));
    }

    let schema_version = "1.0.0".to_string();
    if all_records.len() == 1 {
        Ok(SpectrumFile::Single {
            schema_version,
            spectrum: all_records.into_iter().next().unwrap(),
        })
    } else {
        Ok(SpectrumFile::Batch {
            schema_version,
            batch_metadata: None,
            spectra: all_records,
        })
    }
}

// Apply a parsed KEYWORD+VALUE pair to the SsMeta accumulator.
fn ss_apply_kv(meta: &mut SsMeta, key: &str, value: &str) {
    match key {
        // "SpectraShop\t5.0" first line — treat version as sw_version
        "SPECTRASHOP" => {
            if meta.sw_version.is_none() {
                meta.sw_version = ss_str(value);
            }
        }
        "FILE_DESCRIPTOR" => meta.file_descriptor = ss_str(value),
        "SPECTRUM_TYPE" => meta.spectrum_type = ss_str(value),
        "CREATED" => meta.created = ss_str(value),
        "ORIGINATOR" => meta.originator = ss_str(value),
        "ACQUIRE_NOTE" => meta.acquire_note = ss_str(value),
        "NOTE" => meta.note = ss_str(value),
        "ILLUMINANT" => meta.illuminant = ss_str(value),
        "OBSERVER" => meta.observer = ss_str(value),
        "INSTRUMENTATION" => meta.instrumentation = ss_str(value),
        "INSTRUMENT_SERIAL" | "INSTRUMENT_SERIAL_NUMBER" => meta.instrument_serial = ss_str(value),
        "MEASUREMENT_GEOMETRY" => meta.measurement_geometry = ss_str(value),
        "MEASUREMENT_SOURCE" => meta.measurement_source = ss_str(value),
        "MEASUREMENT_APERTURE" => meta.measurement_aperture = ss_str(value),
        "MEASUREMENT_FILTER" => meta.measurement_filter = ss_str(value),
        "SAMPLE_BACKING" => meta.sample_backing = ss_str(value),
        "SURFACE" => meta.surface = ss_str(value),
        "MANUFACTURER" => meta.manufacturer = ss_str(value),
        "MATERIAL" => meta.material = ss_str(value),
        "MODEL" | "MODEL_NUM" => meta.model_num = ss_str(value),
        "PROD_DATE" => meta.prod_date = ss_str(value),
        "SERIAL" => meta.serial = ss_str(value),
        "NMEASUREMENTS" | "NMEASURE" => {
            meta.nmeasure = ss_str(value)
                .as_deref()
                .unwrap_or(value)
                .trim()
                .parse::<u32>()
                .ok()
        }
        "SW_VERSION" => meta.sw_version = ss_str(value),
        // Informational / ignored
        "NUMBER_OF_SETS" | "NUMBER_OF_FIELDS" | "RGB_SPACE" => {}
        _ => {}
    }
}

// Split "KEYWORD<whitespace>rest" → (UPPERCASE_KEY, rest).
fn ss_split_kv(line: &str) -> Option<(String, String)> {
    let sep = line
        .find('\t')
        .or_else(|| line.find(|c: char| c.is_whitespace()))?;
    let key = line[..sep].trim().to_uppercase();
    let value = line[sep..].trim().to_string();
    if key.is_empty() {
        return None;
    }
    Some((key, value))
}

// Strip surrounding double-quotes and return None for empty strings.
fn ss_str(s: &str) -> Option<String> {
    let s = s.trim();
    let s = s
        .strip_prefix('"')
        .and_then(|s| s.strip_suffix('"'))
        .unwrap_or(s);
    if s.is_empty() {
        None
    } else {
        Some(s.to_string())
    }
}

// Parse a decimal number; accepts comma as decimal separator.
fn ss_float(s: &str) -> Option<f64> {
    s.trim().replace(',', ".").parse::<f64>().ok()
}

// Parse "4 mm" → Some(4.0); "LAV" → None.
fn ss_aperture_mm(s: &str) -> Option<f64> {
    ss_float(s.trim().split_whitespace().next().unwrap_or(s))
}

// Map SpectraShop SPECTRUM_TYPE to MeasurementType.
fn ss_measurement_type(s: &str) -> MeasurementType {
    let lower = s.trim().to_lowercase();
    if lower.contains("transmis") {
        MeasurementType::Transmittance
    } else if lower.contains("emissive") || lower.contains("irradiance") {
        MeasurementType::Irradiance
    } else {
        MeasurementType::Reflectance
    }
}

// Map SpectraShop OBSERVER string to a CIE observer enum value.
fn ss_observer(s: &str) -> Option<String> {
    let lower = s.trim().to_lowercase();
    if lower.contains("10") {
        Some("CIE 1964 10 degree".into())
    } else if lower.contains("2") {
        Some("CIE 1931 2 degree".into())
    } else {
        None
    }
}

// Map a SpectraShop ILLUMINANT string to an allowed illuminant code, if recognised.
fn ss_illuminant(s: &str) -> Option<String> {
    if ALLOWED_ILLUMINANTS.contains(&s) && s != "custom" {
        Some(s.to_string())
    } else {
        None
    }
}

// Build a WavelengthAxis from a list of wavelengths, using range_nm for regular grids.
fn ss_wavelength_axis(wls: &[f64]) -> WavelengthAxis {
    if wls.len() >= 2 {
        let interval = wls[1] - wls[0];
        let regular = interval > 0.0
            && wls
                .windows(2)
                .all(|w| (w[1] - w[0] - interval).abs() < 0.001);
        if regular {
            return WavelengthAxis {
                values_nm: None,
                range_nm: Some(WavelengthRange {
                    start: wls[0],
                    end: *wls.last().unwrap(),
                    interval,
                }),
            };
        }
    }
    WavelengthAxis {
        values_nm: Some(wls.to_vec()),
        range_nm: None,
    }
}

// Parse a date from ISO 8601 (YYYY-MM-DD...) or MM/DD/YYYY formats.
fn ss_parse_date(s: &str) -> Option<String> {
    let s = s.trim();
    if s.len() >= 10 && s.as_bytes().get(4) == Some(&b'-') && s.as_bytes().get(7) == Some(&b'-') {
        return Some(s[..10].to_string());
    }
    let parts: Vec<&str> = s.splitn(3, '/').collect();
    if parts.len() == 3 {
        let month = parts[0].trim();
        let day = parts[1].trim();
        let year = parts[2].trim().split_whitespace().next().unwrap_or("");
        if month.parse::<u32>().is_ok()
            && day.parse::<u32>().is_ok()
            && year.len() == 4
            && year.parse::<u32>().is_ok()
        {
            return Some(format!("{year}-{month:0>2}-{day:0>2}"));
        }
    }
    None
}

// Tokenize a SpectraShop data line into tab-delimited fields, stripping the
// outermost pair of quotes from each field. This correctly handles sample names
// that contain embedded quotes (e.g. `"United Drug Co. "Graph Blue""`), because
// only the first and last `"` are stripped rather than the first matching pair.
fn ss_tokenize_into(line: &str, tokens: &mut Vec<String>) {
    for field in line.split('\t') {
        if field.is_empty() {
            continue; // trailing or adjacent tabs — skip
        }
        let s = field.trim();
        if s.starts_with('"') && s.ends_with('"') && s.len() >= 2 {
            tokens.push(s[1..s.len() - 1].to_string());
        } else if !s.is_empty() {
            tokens.push(s.to_string());
        }
    }
}

fn ss_build_record(
    fields: &[String],
    values: &[String],
    meta: &SsMeta,
    idx: usize,
    source_file: Option<&str>,
) -> Result<SpectrumRecord> {
    // SAMPLE_ID1 = machine/file id (becomes SpectrumRecord.id)
    // SAMPLE_ID2 / SAMPLE_NAME = human label (becomes metadata.title; falls back to SAMPLE_ID1)
    // SAMPLE_ID3 = rarely used; goes to custom if non-empty
    let mut ss_id1: Option<String> = None;
    let mut ss_id2: Option<String> = None;
    let mut wls: Vec<f64> = Vec::new();
    let mut spectral_vals: Vec<f64> = Vec::new();
    let mut custom_map: serde_json::Map<String, serde_json::Value> = serde_json::Map::new();

    let mut i = 0;
    while i < fields.len() {
        let field = fields[i].as_str();
        let val = &values[i];

        match field {
            "SAMPLE_ID" | "SAMPLE_ID1" => {
                ss_id1 = ss_str(val);
                i += 1;
            }
            "SAMPLE_NAME" | "SAMPLE_ID2" => {
                ss_id2 = ss_str(val);
                i += 1;
            }
            "SAMPLE_ID3" => {
                if let Some(s) = ss_str(val) {
                    custom_map.insert("sample_id3".into(), serde_json::Value::String(s));
                }
                i += 1;
            }
            "SPECTRAL_NM" => {
                // Consume the paired SPECTRAL_NM + SPECTRAL_VAL together.
                if i + 1 < fields.len() && fields[i + 1] == "SPECTRAL_VAL" {
                    if let (Some(wl), Some(sv)) = (ss_float(val), ss_float(&values[i + 1])) {
                        wls.push(wl);
                        spectral_vals.push(sv);
                    }
                    i += 2;
                } else {
                    i += 1;
                }
            }
            "SPECTRAL_VAL" => i += 1,
            _ => {
                if let Some(s) = ss_str(val) {
                    custom_map.insert(field.to_lowercase(), serde_json::Value::String(s));
                }
                i += 1;
            }
        }
    }

    if wls.is_empty() {
        return Err(SpectrumFileError::SchemaValidation(format!(
            "SpectraShop: record {idx} has no SPECTRAL_NM/SPECTRAL_VAL pairs"
        )));
    }

    // id comes from SAMPLE_ID / SAMPLE_ID1.
    // title prefers SAMPLE_NAME / SAMPLE_ID2 (the human label) and falls back to
    // SAMPLE_ID / SAMPLE_ID1 when no separate name field is present (e.g. Smarties).
    let id = ss_id1.clone().unwrap_or_else(|| format!("{}", idx + 1));
    let title = ss_id2.clone().or_else(|| ss_id1.clone());
    let date = meta
        .created
        .as_deref()
        .and_then(ss_parse_date)
        .unwrap_or_else(|| "1970-01-01".to_string());
    let measurement_type = meta
        .spectrum_type
        .as_deref()
        .map(ss_measurement_type)
        .unwrap_or(MeasurementType::Reflectance);

    let instrument =
        (meta.instrumentation.is_some() || meta.instrument_serial.is_some()).then(|| Instrument {
            manufacturer: None,
            model: meta.instrumentation.clone(),
            serial_number: meta.instrument_serial.clone(),
            detector_type: None,
            light_source: None,
        });

    let aperture_mm = meta
        .measurement_aperture
        .as_deref()
        .and_then(ss_aperture_mm);
    let filter = meta.measurement_filter.as_deref().and_then(|s| {
        if s.eq_ignore_ascii_case("none") {
            None
        } else {
            Some(s.to_string())
        }
    });
    let mc_any = meta.measurement_geometry.is_some()
        || aperture_mm.is_some()
        || filter.is_some()
        || meta.nmeasure.is_some();
    let measurement_conditions = mc_any.then(|| MeasurementConditions {
        integration_time_ms: None,
        averaging: meta.nmeasure,
        temperature_celsius: None,
        geometry: meta.measurement_geometry.clone(),
        specular_component: None,
        spectral_resolution_nm: None,
        measurement_aperture_mm: aperture_mm,
        measurement_filter: filter,
    });

    let illuminant = meta.illuminant.as_deref().and_then(ss_illuminant);
    let observer = meta.observer.as_deref().and_then(ss_observer);
    let color_science = (illuminant.is_some() || observer.is_some()).then(|| ColorScience {
        illuminant,
        illuminant_custom_sd: None,
        cie_observer: observer,
        white_reference: None,
        results: None,
    });

    // File-level fields that don't map to explicit SpectrumRecord fields go in custom.
    for (k, v) in [
        ("manufacturer", &meta.manufacturer),
        ("material", &meta.material),
        ("model_num", &meta.model_num),
        ("serial", &meta.serial),
        ("prod_date", &meta.prod_date),
        ("measurement_source", &meta.measurement_source),
    ] {
        if let Some(s) = v {
            custom_map
                .entry(k)
                .or_insert_with(|| serde_json::Value::String(s.clone()));
        }
    }
    let custom = (!custom_map.is_empty()).then(|| serde_json::Value::Object(custom_map));

    let notes: Vec<&str> = [&meta.acquire_note, &meta.note]
        .iter()
        .filter_map(|o| o.as_deref())
        .collect();
    let notes_str = (!notes.is_empty()).then(|| notes.join("; "));

    Ok(SpectrumRecord {
        id,
        metadata: SpectrumMetadata {
            measurement_type,
            date,
            title,
            description: meta.file_descriptor.clone(),
            sample_id: None,
            time: None,
            operator: meta.originator.clone(),
            instrument,
            measurement_conditions,
            surface: meta.surface.clone(),
            sample_backing: meta.sample_backing.clone(),
            tags: None,
            copyright: None,
            custom,
        },
        wavelength_axis: ss_wavelength_axis(&wls),
        spectral_data: SpectralData {
            values: spectral_vals,
            uncertainty: None,
            scale: None,
        },
        color_science,
        provenance: Some(Provenance {
            software: Some("SpectraShop".into()),
            software_version: meta.sw_version.clone(),
            source_file: source_file.map(str::to_string),
            source_format: Some("SpectraShop Text".into()),
            processing_steps: None,
            notes: notes_str,
        }),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ss_input_single() -> String {
        // Tab-separated SpectraShop header + data block with 3 wavelength points.
        [
            "FILE_DESCRIPTOR\tTest file",
            "SPECTRUM_TYPE\tReflective",
            "CREATED\t2014-05-15",
            "ORIGINATOR\tTest User",
            "ILLUMINANT\tD65",
            "OBSERVER\t2 degree",
            "INSTRUMENTATION\tCM-700d",
            "INSTRUMENT_SERIAL\tABC123",
            "MEASUREMENT_GEOMETRY\td:8",
            "MEASUREMENT_FILTER\tNone",
            "SAMPLE_BACKING\tBlack",
            "NMEASURE\t3",
            "SW_VERSION\t5.0.0",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSAMPLE_NAME\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "s1\tRedTile\t380\t0.048\t390\t0.051\t400\t0.054",
            "END_DATA",
        ]
        .join("\n")
    }

    #[test]
    fn spectrashop_single_record() {
        let file = SpectrumFile::from_spectrashop_str(&ss_input_single()).unwrap();
        assert!(matches!(file, SpectrumFile::Single { .. }));
        let spectra = file.spectra();
        assert_eq!(spectra.len(), 1);
        let sp = spectra[0];
        assert_eq!(sp.id, "s1");
        assert_eq!(sp.metadata.title.as_deref(), Some("RedTile"));
        assert_eq!(sp.metadata.description.as_deref(), Some("Test file"));
        assert_eq!(sp.metadata.operator.as_deref(), Some("Test User"));
        assert_eq!(sp.metadata.date, "2014-05-15");
        assert_eq!(sp.metadata.sample_backing.as_deref(), Some("Black"));
        assert!(matches!(
            sp.metadata.measurement_type,
            MeasurementType::Reflectance
        ));
        assert_eq!(sp.n_points(), 3);
        let cs = sp.color_science.as_ref().unwrap();
        assert_eq!(cs.illuminant.as_deref(), Some("D65"));
        assert_eq!(cs.cie_observer.as_deref(), Some("CIE 1931 2 degree"));
        let mc = sp.metadata.measurement_conditions.as_ref().unwrap();
        assert_eq!(mc.geometry.as_deref(), Some("d:8"));
        assert!(mc.measurement_filter.is_none()); // "None" filtered out
        assert_eq!(mc.averaging, Some(3));
        let prov = sp.provenance.as_ref().unwrap();
        assert_eq!(prov.software.as_deref(), Some("SpectraShop"));
        assert_eq!(prov.software_version.as_deref(), Some("5.0.0"));
        assert_eq!(prov.source_format.as_deref(), Some("SpectraShop Text"));
        let instr = sp.metadata.instrument.as_ref().unwrap();
        assert_eq!(instr.model.as_deref(), Some("CM-700d"));
        assert_eq!(instr.serial_number.as_deref(), Some("ABC123"));
    }

    #[test]
    fn spectrashop_batch_two_records() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "1\t380\t0.1\t390\t0.2",
            "2\t380\t0.3\t390\t0.4",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        assert!(matches!(file, SpectrumFile::Batch { .. }));
        assert_eq!(file.spectra().len(), 2);
        assert_eq!(file.spectra()[0].id, "1");
        assert_eq!(file.spectra()[1].id, "2");
    }

    #[test]
    fn spectrashop_regular_grid_uses_range_nm() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "s1\t380\t0.1\t390\t0.2\t400\t0.3",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        let sp = &file.spectra()[0];
        assert!(sp.wavelength_axis.range_nm.is_some());
        assert!(sp.wavelength_axis.values_nm.is_none());
        let r = sp.wavelength_axis.range_nm.as_ref().unwrap();
        assert!((r.start - 380.0).abs() < 1e-9);
        assert!((r.end - 400.0).abs() < 1e-9);
        assert!((r.interval - 10.0).abs() < 1e-9);
    }

    #[test]
    fn spectrashop_irregular_grid_uses_values_nm() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "s1\t380\t0.1\t400\t0.2\t450\t0.3",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        let sp = &file.spectra()[0];
        assert!(sp.wavelength_axis.values_nm.is_some());
        assert!(sp.wavelength_axis.range_nm.is_none());
    }

    #[test]
    fn spectrashop_date_mm_dd_yyyy() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "CREATED\t05/15/2014 13:40:29",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "s1\t380\t0.1\t390\t0.2",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        assert_eq!(file.spectra()[0].metadata.date, "2014-05-15");
    }

    #[test]
    fn spectrashop_transmissive_type() {
        let input = [
            "SPECTRUM_TYPE\tTransmissive",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "s1\t380\t0.1\t390\t0.2",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        assert!(matches!(
            file.spectra()[0].metadata.measurement_type,
            MeasurementType::Transmittance
        ));
    }

    #[test]
    fn spectrashop_10_degree_observer() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "OBSERVER\t10 degree",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "s1\t380\t0.1\t390\t0.2",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        let cs = file.spectra()[0].color_science.as_ref().unwrap();
        assert_eq!(cs.cie_observer.as_deref(), Some("CIE 1964 10 degree"));
    }

    #[test]
    fn spectrashop_aperture_mm_parsed() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "MEASUREMENT_GEOMETRY\td:8",
            "MEASUREMENT_APERTURE\t4 mm",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "s1\t380\t0.1\t390\t0.2",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        let mc = file.spectra()[0]
            .metadata
            .measurement_conditions
            .as_ref()
            .unwrap();
        assert!((mc.measurement_aperture_mm.unwrap() - 4.0).abs() < 1e-9);
    }

    #[test]
    fn spectrashop_lav_aperture_is_none() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "MEASUREMENT_GEOMETRY\td:8",
            "MEASUREMENT_APERTURE\tLAV",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "s1\t380\t0.1\t390\t0.2",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        let mc = file.spectra()[0]
            .metadata
            .measurement_conditions
            .as_ref()
            .unwrap();
        assert!(mc.measurement_aperture_mm.is_none());
    }

    #[test]
    fn spectrashop_quoted_sample_name_with_spaces() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSAMPLE_NAME\tSPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "s1\t\"Red Munsell\"\t380\t0.1\t390\t0.2",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        assert_eq!(
            file.spectra()[0].metadata.title.as_deref(),
            Some("Red Munsell")
        );
    }

    #[test]
    fn spectrashop_missing_data_format_is_error() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "BEGIN_DATA",
            "s1\t380\t0.1\t390\t0.2",
            "END_DATA",
        ]
        .join("\n");

        assert!(matches!(
            SpectrumFile::from_spectrashop_str(&input),
            Err(SpectrumFileError::SchemaValidation(_))
        ));
    }

    #[test]
    fn spectrashop_empty_data_block_is_error() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "BEGIN_DATA_FORMAT",
            "SAMPLE_ID\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "END_DATA",
        ]
        .join("\n");

        assert!(matches!(
            SpectrumFile::from_spectrashop_str(&input),
            Err(SpectrumFileError::SchemaValidation(_))
        ));
    }

    #[test]
    fn spectrashop_id_fallback_when_no_sample_id_field() {
        let input = [
            "SPECTRUM_TYPE\tReflective",
            "BEGIN_DATA_FORMAT",
            "SPECTRAL_NM\tSPECTRAL_VAL\tSPECTRAL_NM\tSPECTRAL_VAL",
            "END_DATA_FORMAT",
            "BEGIN_DATA",
            "380\t0.1\t390\t0.2",
            "END_DATA",
        ]
        .join("\n");

        let file = SpectrumFile::from_spectrashop_str(&input).unwrap();
        // No SAMPLE_ID field → fallback to "1"
        assert_eq!(file.spectra()[0].id, "1");
    }
}
