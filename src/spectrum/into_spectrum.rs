use crate::{spectrum::Spectrum, Error};

/// The physical quantity represented by a spectral measurement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeasurementKind {
    Reflectance,
    Transmittance,
    Absorbance,
    Radiance,
    Irradiance,
}

/// Raw spectral samples returned by [`IntoSpectrum::spectral_data`].
///
/// Values must be in the natural scale for the measurement kind: reflectance
/// and transmittance in `[0, 1]` (fractional), radiance and irradiance in
/// their physical units. Implementors are responsible for any normalisation
/// (e.g. converting percent to fractional) before returning.
pub struct SpectralSample {
    pub kind: MeasurementKind,
    pub wavelengths_nm: Vec<f64>,
    pub values: Vec<f64>,
}

/// Convert spectral data into a [`Spectrum`] on the 380–780 nm / 1 nm grid.
///
/// Implementors provide only [`spectral_data`](IntoSpectrum::spectral_data);
/// the four interpolation / aggregation methods are provided as defaults.
pub trait IntoSpectrum {
    /// Return the raw spectral samples: measurement kind, wavelengths (nm),
    /// and corresponding values normalised to the natural scale.
    fn spectral_data(&self) -> SpectralSample;

    /// Linear interpolation onto the 380–780 nm / 1 nm grid.
    ///
    /// Works with both evenly-spaced and irregular wavelength axes.
    fn to_spectrum_linear(&self) -> Result<Spectrum, Error> {
        let SpectralSample {
            wavelengths_nm,
            values,
            ..
        } = self.spectral_data();
        Spectrum::linear_interpolate(&wavelengths_nm, &values)
    }

    /// Sprague (5th-order) interpolation onto the 380–780 nm / 1 nm grid.
    ///
    /// Requires equidistant input; returns an error if the wavelength axis is
    /// irregular. Use `to_spectrum_linear` for non-uniform grids.
    fn to_spectrum_sprague(&self) -> Result<Spectrum, Error> {
        let SpectralSample {
            wavelengths_nm,
            values,
            ..
        } = self.spectral_data();
        let range = equidistant_range(&wavelengths_nm).ok_or_else(|| {
            Error::ErrorString(
                "Sprague interpolation requires equidistant data; \
                 use to_spectrum_linear for irregular grids"
                    .into(),
            )
        })?;
        Spectrum::sprague_interpolate(range, &values)
    }

    /// Gaussian-weighted kernel regression onto the 380–780 nm / 1 nm grid.
    ///
    /// `spectral_resolution_nm` is the instrument FWHM; σ = FWHM / 2.355.
    /// Pixels further than 4σ from an output wavelength do not contribute.
    fn to_spectrum_smooth(&self, spectral_resolution_nm: f64) -> Result<Spectrum, Error> {
        if spectral_resolution_nm <= 0.0 {
            return Err(Error::ErrorString(
                "spectral_resolution_nm must be positive".into(),
            ));
        }
        let SpectralSample {
            wavelengths_nm,
            values,
            ..
        } = self.spectral_data();
        let sigma = spectral_resolution_nm / 2.355_f64;
        let cutoff = 4.0 * sigma;

        let mut wl_out: Vec<f64> = Vec::with_capacity(401);
        let mut v_out: Vec<f64> = Vec::with_capacity(401);

        for w in 380_u32..=780 {
            let w_f = w as f64;
            let mut weighted_sum = 0f64;
            let mut weight_total = 0f64;

            for (&wl, &v) in wavelengths_nm.iter().zip(values.iter()) {
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

    /// Boxcar binning (width = `bin_width_nm`) then Sprague interpolation.
    ///
    /// Falls back to linear interpolation if any bins are empty. Returns an
    /// error if fewer than 7 bins are populated (Sprague minimum).
    fn to_spectrum_binned(&self, bin_width_nm: f64) -> Result<Spectrum, Error> {
        if bin_width_nm <= 0.0 {
            return Err(Error::ErrorString("bin_width_nm must be positive".into()));
        }

        let SpectralSample {
            wavelengths_nm,
            values,
            ..
        } = self.spectral_data();

        let (wl_min, wl_max) = match (wavelengths_nm.first(), wavelengths_nm.last()) {
            (Some(&a), Some(&b)) => (a, b),
            _ => return Err(Error::ErrorString("no spectral data".into())),
        };

        let n_bins = ((wl_max - wl_min) / bin_width_nm + 1e-9).floor() as usize + 1;
        let mut sums = vec![0f64; n_bins];
        let mut counts = vec![0u32; n_bins];

        for (&wl, &v) in wavelengths_nm.iter().zip(values.iter()) {
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

fn equidistant_range(wls: &[f64]) -> Option<[f64; 2]> {
    let (&first, rest) = wls.split_first()?;
    let &last = wls.last()?;
    if rest.is_empty() {
        return None;
    }
    let step = (last - first) / (wls.len() - 1) as f64;
    let ok = wls
        .windows(2)
        .all(|w| (w[1] - w[0] - step).abs() < 1e-6 * step.abs().max(1.0));
    ok.then_some([first, last])
}
