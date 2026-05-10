#![cfg(feature = "spectrashop")]
/// Integration test: parse every SpectraShop .txt file in io/data/spectrashop/.
///
/// For each file we assert:
///   - Parsing succeeds (no error)
///   - At least one spectrum was produced
///   - Every spectrum has at least 2 wavelength points
///   - Every spectrum has a non-empty id
///   - Wavelengths are strictly increasing (cross-field invariant)
use colorimetry_io::SpectrumFile;
use std::path::PathBuf;

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/spectrashop")
}

fn collect_txt_files(dir: &std::path::Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    if let Ok(entries) = std::fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                files.extend(collect_txt_files(&path));
            } else if path.extension().and_then(|e| e.to_str()) == Some("txt") {
                files.push(path);
            }
        }
    }
    files.sort();
    files
}

#[test]
fn parse_all_spectrashop_files() {
    let dir = data_dir();
    let files = collect_txt_files(&dir);
    assert!(
        !files.is_empty(),
        "No .txt files found under {}",
        dir.display()
    );

    let mut failures: Vec<String> = Vec::new();

    for path in &files {
        let rel = path.strip_prefix(&dir).unwrap_or(path);

        let file = match SpectrumFile::from_spectrashop_path(path) {
            Ok(f) => f,
            Err(e) => {
                failures.push(format!("{}: PARSE ERROR — {e}", rel.display()));
                continue;
            }
        };

        let spectra = file.spectra();
        if spectra.is_empty() {
            failures.push(format!("{}: no spectra produced", rel.display()));
            continue;
        }

        for sp in &spectra {
            if sp.id.is_empty() {
                failures.push(format!("{}: spectrum with empty id", rel.display()));
            }
            if sp.n_points() < 2 {
                failures.push(format!(
                    "{}: spectrum '{}' has fewer than 2 points",
                    rel.display(),
                    sp.id
                ));
            }
            let wls = sp.wavelength_axis.wavelengths_nm();
            if wls.windows(2).any(|w| w[0] >= w[1]) {
                failures.push(format!(
                    "{}: spectrum '{}' wavelengths not strictly increasing",
                    rel.display(),
                    sp.id
                ));
            }
        }
    }

    if !failures.is_empty() {
        let total = files.len();
        let n_fail = failures.len();
        panic!(
            "{n_fail} issue(s) across {total} file(s):\n{}",
            failures.join("\n")
        );
    }

    println!("All {} SpectraShop files parsed successfully.", files.len());
}

#[test]
fn spectrashop_smarties_spot_check() {
    let path = data_dir().join("candies/Smarties.txt");
    let file = SpectrumFile::from_spectrashop_path(&path).unwrap();
    let spectra = file.spectra();

    // File says NUMBER_OF_SETS 8
    assert_eq!(spectra.len(), 8, "expected 8 Smarties records");

    // First record is "Red"
    assert_eq!(spectra[0].id, "Red");
    assert_eq!(spectra[0].metadata.title.as_deref(), Some("Red"));

    // Wavelength grid should be regular (range_nm)
    assert!(
        spectra[0].wavelength_axis.range_nm.is_some(),
        "expected range_nm for regular grid"
    );

    // Observer and illuminant should be mapped
    let cs = spectra[0]
        .color_science
        .as_ref()
        .expect("color_science missing");
    assert_eq!(cs.illuminant.as_deref(), Some("D65"));
    assert_eq!(cs.cie_observer.as_deref(), Some("CIE 1931 2 degree"));

    // Software provenance
    let prov = spectra[0].provenance.as_ref().expect("provenance missing");
    assert_eq!(prov.software.as_deref(), Some("SpectraShop"));
}

#[test]
fn spectrashop_colorchecker_spot_check() {
    let path = data_dir().join("charts/ColorChecker 1977 #1.txt");
    let file = SpectrumFile::from_spectrashop_path(&path).unwrap();
    assert_eq!(file.spectra().len(), 24, "ColorChecker has 24 patches");
}

#[test]
fn spectrashop_filters_transmissive() {
    let path = data_dir().join("filters/Wratten Filters.txt");
    let file = SpectrumFile::from_spectrashop_path(&path).unwrap();
    assert!(!file.spectra().is_empty());
    // SPECTRUM_TYPE = Transmissive
    assert!(matches!(
        file.spectra()[0].metadata.measurement_type,
        colorimetry_io::MeasurementType::Transmittance
    ));
}

#[test]
fn spectrashop_monitor_irradiance() {
    let path = data_dir().join("monitors/Apple 13 inch.txt");
    let file = SpectrumFile::from_spectrashop_path(&path).unwrap();
    // SPECTRUM_TYPE = Emissive-monitor → Irradiance
    assert!(matches!(
        file.spectra()[0].metadata.measurement_type,
        colorimetry_io::MeasurementType::Irradiance
    ));
    assert_eq!(file.spectra().len(), 4); // R, G, B, W
}
