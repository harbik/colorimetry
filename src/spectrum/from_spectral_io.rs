use spectral_io::SpectrumRecord;

use crate::{spectrum::Spectrum, Error};

fn normalized_values(record: &SpectrumRecord) -> Vec<f64> {
    let mut vals = record.spectral_data.values.clone();
    if record.spectral_data.scale.as_deref() == Some("percent") {
        vals.iter_mut().for_each(|v| *v /= 100.0);
    }
    vals
}

/// Extension trait that converts a `spectral_io::SpectrumRecord` into a
/// [`Spectrum`] using one of four interpolation / aggregation strategies.
pub trait IntoSpectrum {
    /// Linear interpolation onto the 380–780 nm / 1 nm grid.
    ///
    /// Works with both `values_nm` and `range_nm` axes. Simple and robust;
    /// less accurate than Sprague for smooth spectra on a regular grid.
    fn to_spectrum_linear(&self) -> Result<Spectrum, Error>;

    /// Sprague (5th-order) interpolation onto the 380–780 nm / 1 nm grid.
    ///
    /// More accurate than linear for smooth spectra but requires equidistant
    /// input. Returns an error if the axis uses `values_nm` (irregular grid).
    fn to_spectrum_sprague(&self) -> Result<Spectrum, Error>;

    /// Gaussian-weighted kernel regression onto the 380–780 nm / 1 nm grid.
    ///
    /// `spectral_resolution_nm` is the instrument FWHM; σ = FWHM / 2.355.
    /// Pixels further than 4σ from an output wavelength do not contribute.
    /// Best choice when pixel spacing is much finer than the slit width.
    fn to_spectrum_smooth(&self, spectral_resolution_nm: f64) -> Result<Spectrum, Error>;

    /// Boxcar binning (width = `bin_width_nm`) then Sprague interpolation.
    ///
    /// Set `bin_width_nm` to the instrument's optical resolution. Falls back
    /// to linear interpolation if any bins are empty. Returns an error if
    /// fewer than 7 bins are populated (Sprague minimum).
    fn to_spectrum_binned(&self, bin_width_nm: f64) -> Result<Spectrum, Error>;
}

impl IntoSpectrum for SpectrumRecord {
    fn to_spectrum_linear(&self) -> Result<Spectrum, Error> {
        let wl = self.wavelength_axis.wavelengths_nm();
        let vals = normalized_values(self);
        Spectrum::linear_interpolate(&wl, &vals)
    }

    fn to_spectrum_sprague(&self) -> Result<Spectrum, Error> {
        let range = self.wavelength_axis.range_nm.as_ref().ok_or_else(|| {
            Error::ErrorString(
                "Sprague interpolation requires equidistant data (range_nm); \
                 use to_spectrum_linear for values_nm axes"
                    .into(),
            )
        })?;
        let vals = normalized_values(self);
        Spectrum::sprague_interpolate([range.start, range.end], &vals)
    }

    fn to_spectrum_smooth(&self, spectral_resolution_nm: f64) -> Result<Spectrum, Error> {
        if spectral_resolution_nm <= 0.0 {
            return Err(Error::ErrorString(
                "spectral_resolution_nm must be positive".into(),
            ));
        }
        let sigma = spectral_resolution_nm / 2.355_f64;
        let cutoff = 4.0 * sigma;

        let vals = normalized_values(self);
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
            return Err(Error::ErrorString(
                "spectral data does not cover the 380–780 nm output range".into(),
            ));
        }

        Spectrum::linear_interpolate(&wl_out, &v_out)
    }

    fn to_spectrum_binned(&self, bin_width_nm: f64) -> Result<Spectrum, Error> {
        if bin_width_nm <= 0.0 {
            return Err(Error::ErrorString("bin_width_nm must be positive".into()));
        }

        let vals = normalized_values(self);
        let wls = self.wavelength_axis.wavelengths_nm();

        let (wl_min, wl_max) = match (wls.first(), wls.last()) {
            (Some(&a), Some(&b)) => (a, b),
            _ => return Err(Error::ErrorString("no spectral data".into())),
        };

        let n_bins = ((wl_max - wl_min) / bin_width_nm + 1e-9).floor() as usize + 1;
        let mut sums = vec![0f64; n_bins];
        let mut counts = vec![0u32; n_bins];

        for (&wl, &v) in wls.iter().zip(vals.iter()) {
            let idx = ((wl - wl_min) / bin_width_nm + 1e-9).floor() as usize;
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
            return Err(Error::ProvideAtLeastNValues(7));
        }

        if wl_out.len() == n_bins {
            let last_bin_nm = wl_min + (n_bins - 1) as f64 * bin_width_nm;
            Spectrum::sprague_interpolate([wl_min, last_bin_nm], &v_out)
        } else {
            Spectrum::linear_interpolate(&wl_out, &v_out)
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
    fn to_spectrum_sprague_error_on_values_nm() {
        let wls: Vec<f64> = (0..41).map(|i| 380.0 + i as f64 * 10.0).collect();
        assert!(make_values_record(wls, vals_41())
            .to_spectrum_sprague()
            .is_err());
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
        // 380–780 nm at 5 nm steps (81 points), binned with 15 nm.
        // 400 / 15 = 26.67, so the last bin does not align exactly with 780 nm.
        // Before the fix, round() assigned idx=27 for 780 nm which equalled n_bins
        // (27) and was silently dropped.
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
}
