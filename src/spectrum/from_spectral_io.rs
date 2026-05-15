use spectral_io::SpectrumRecord;

use crate::spectrum::into_spectrum::{IntoSpectrum, SpectralSample};

impl IntoSpectrum for SpectrumRecord {
    fn spectral_sample(&self) -> SpectralSample {
        let mut values = self.spectral_data.values.clone();
        if self.spectral_data.scale.as_deref() == Some("percent") {
            values.iter_mut().for_each(|v| *v /= 100.0);
        }
        SpectralSample {
            kind: self.metadata.measurement_type,
            wavelengths_nm: self.wavelength_axis.wavelengths_nm(),
            values,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spectral_io::{
        MeasurementType, SpectralData, SpectrumMetadata, SpectrumRecord, WavelengthAxis,
        WavelengthRange,
    };

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
                surface: None,
                sample_backing: None,
                tags: None,
                copyright: None,
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
                surface: None,
                sample_backing: None,
                tags: None,
                copyright: None,
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
        assert!(make_range_record(vals_41()).to_spectrum_linear().is_ok());
    }

    #[test]
    fn to_spectrum_linear_values_nm() {
        let wls: Vec<f64> = (0..41).map(|i| 380.0 + i as f64 * 10.0).collect();
        assert!(make_values_record(wls, vals_41())
            .to_spectrum_linear()
            .is_ok());
    }

    #[test]
    fn to_spectrum_sprague_happy_path() {
        assert!(make_range_record(vals_41()).to_spectrum_sprague().is_ok());
    }

    #[test]
    fn to_spectrum_sprague_accepts_equidistant_values_nm() {
        // Evenly-spaced values_nm is now accepted: equidistance is detected
        // from the wavelength array itself, not from the axis variant.
        let wls: Vec<f64> = (0..41).map(|i| 380.0 + i as f64 * 10.0).collect();
        assert!(make_values_record(wls, vals_41())
            .to_spectrum_sprague()
            .is_ok());
    }

    #[test]
    fn to_spectrum_sprague_error_on_irregular_wavelengths() {
        let wls = vec![
            380.0, 390.0, 405.0, 420.0, 440.0, 460.0, 490.0, 520.0, 560.0, 600.0, 640.0, 680.0,
            720.0, 760.0, 780.0,
        ];
        let vals: Vec<f64> = (0..15).map(|i| i as f64 / 14.0).collect();
        assert!(make_values_record(wls, vals).to_spectrum_sprague().is_err());
    }

    #[test]
    fn to_spectrum_smooth_happy_path() {
        assert!(make_range_record(vals_41())
            .to_spectrum_smooth(15.0)
            .is_ok());
    }

    #[test]
    fn to_spectrum_smooth_negative_resolution_is_error() {
        assert!(make_range_record(vals_41())
            .to_spectrum_smooth(-1.0)
            .is_err());
    }

    #[test]
    fn to_spectrum_smooth_zero_resolution_is_error() {
        assert!(make_range_record(vals_41())
            .to_spectrum_smooth(0.0)
            .is_err());
    }

    #[test]
    fn to_spectrum_binned_happy_path() {
        assert!(make_range_record(vals_41())
            .to_spectrum_binned(10.0)
            .is_ok());
    }

    #[test]
    fn to_spectrum_binned_negative_bin_width_is_error() {
        assert!(make_range_record(vals_41())
            .to_spectrum_binned(-1.0)
            .is_err());
    }

    #[test]
    fn to_spectrum_binned_too_few_bins_is_error() {
        let wls: Vec<f64> = (0..5).map(|i| 380.0 + i as f64 * 50.0).collect();
        assert!(make_values_record(wls, vec![0.1, 0.2, 0.3, 0.4, 0.5])
            .to_spectrum_binned(50.0)
            .is_err());
    }

    #[test]
    fn to_spectrum_binned_non_divisible_range_does_not_drop_last_point() {
        let wls: Vec<f64> = (0..81).map(|i| 380.0 + i as f64 * 5.0).collect();
        let vals: Vec<f64> = (0..81).map(|i| i as f64 / 80.0).collect();
        assert!(make_values_record(wls, vals)
            .to_spectrum_binned(15.0)
            .is_ok());
    }

    #[test]
    fn normalized_values_percent_scale() {
        let percent_vals: Vec<f64> = (0..41).map(|i| i as f64 * 2.5).collect();
        let mut rec = make_range_record(percent_vals);
        rec.spectral_data.scale = Some("percent".into());
        assert!(rec.to_spectrum_linear().is_ok());
    }

    #[test]
    fn spectral_sample_preserves_measurement_type() {
        use spectral_io::MeasurementType;
        for kind in [
            MeasurementType::Reflectance,
            MeasurementType::Transmittance,
            MeasurementType::Absorbance,
            MeasurementType::Radiance,
            MeasurementType::Irradiance,
            MeasurementType::Emission,
        ] {
            let mut rec = make_range_record(vals_41());
            rec.metadata.measurement_type = kind;
            assert_eq!(rec.spectral_sample().kind, kind);
        }
    }
}
