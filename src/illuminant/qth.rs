// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

//! Quartz Tungsten Halogen (QTH) illuminant spectral power models.
//!
//! This module provides functions to compute the spectral power distribution (SPD) of Quartz Tungsten Halogen (QTH) lamps.
//! QTH lamps are commonly used in photography, microscopy, and other applications requiring a stable and continuous light source.
//! The SPD is modeled using Planck's law of black-body radiation, adjusted for the specific characteristics of QTH lamps.

use std::borrow::Cow;

use crate::{
    illuminant::planck::Planck,
    spectrum::{wavelength, Spectrum, NS, SPECTRUM_WAVELENGTH_RANGE},
    traits::Light,
};

/// Quartz Tungsten Halogen Illuminant, modeled according to the work
/// of M. Ojanen et all., Metrologia 49 (2012) S53–S58.
///
/// # Notes
/// - The spectral emissivity Tungsten data this model uses, is based on the measurements
///   of De Vos (1954), covering a temperature range from 1600K to 2800K. For higher temperatures
///   the model extrapolates the emissivity values.
/// - The model includes a recycling factor to account for the effect of filament reflections
///   inside the lamp envelope, which increases the effective emissivity of the filament.
///   A typical value for the recycling factor is 0.5 for a double coiled FEL lamp, frequently
///   used a calibration lamp for spectrometers, as estimated by Ojanen et all
///
/// # References
/// * M. Ojanen et all., Metrologia 49 (2012) S53–S58
/// * J.C De Vos, Physica 20, 690 (1954)
/// * M. Pon and Jan P. Hessler, Applied Optics 23, 7, 975-978 (1984).
pub struct QTH {
    filament_temperature: f64, // actual filament temperature in Kelvin, not the correlated color temperature
    scale: f64,                //
    recycling_factor: f64, // fraction of light reflected back to the filament, typically 0.5 for FEL lamps
    include_envelope: bool, // whether to include the quartz envelope transmission
}

impl QTH {
    /// Coefficients for a third degree polynomial converting correlated color temperature (CCT) to
    /// filament temperature (T).  This polynomial was fitted by calculating the CCT of a QTH lamp
    /// for a range of filament temperatures, and then fitting a polynomial to the results, with the
    /// filament temperature as dependent parameter.
    ///
    /// # Note:
    /// - The coefficients were derived from a lamp with the default recycling factor of 0.5 and including
    ///   the quartz envelope transmission. Different coefficients may be needed for other configurations,
    ///   but the differences are expected to be very small. Please report if you find significant deviations.
    /// - Correlated color temperatures are about 10 to 40K higher than the actual filament temperature,
    ///   with the larger differences at higher temperatures, as their color points are a bit above the Planckian locus,
    ///   with tint values in the range frpom 0.06 for a temperature of 1600K to 0.38 for a temperature of 3200K,
    const CCT2T: [f64; 4] = [-189.0, 82.8, -13.5, 0.768];

    /// Create a new QTH illuminant model with the given correlated color temperature
    /// The correlated color should be in the range from 1600 to 3200K
    pub fn new(cct: f64) -> Self {
        assert!(
            (1600.0..=3250.0).contains(&cct),
            "Temperature should be between 2500K and 3200K for QTH lamps."
        );
        let rt = 10000.0 / cct;
        let delta_temp = Self::CCT2T
            .into_iter()
            .rev()
            .fold(0.0, |acc, coeff| acc * rt + coeff);
        Self::from_filament_temperature(cct + delta_temp)
    }

    pub fn from_filament_temperature(temperature: f64) -> Self {
        assert!(
            (1559.0..=3200.0).contains(&temperature),
            "Temperature should be between 1600K and 3200K for QTH lamps."
        );
        let scale = 1.0; // Scale can be adjusted based on lamp power and efficiency
        let recycling_factor = 0.5;
        QTH {
            filament_temperature: temperature,
            scale,
            recycling_factor,
            include_envelope: true,
        }
    }

    pub fn with_recyclingr(mut self, recycling_factor: f64) -> Self {
        assert!(
            (0.0..1.0).contains(&recycling_factor),
            "Recycling factor should be between 0.0 and 1.0"
        );
        self.recycling_factor = recycling_factor;
        self
    }

    pub fn without_recycling(mut self) -> Self {
        self.recycling_factor = 0.0;
        self
    }

    pub fn with_scale(mut self, scale: f64) -> Self {
        assert!(scale > 0.0, "Scale should be positive");
        self.scale = scale;
        self
    }

    pub fn with_envelope(mut self) -> Self {
        self.include_envelope = true;
        self
    }

    pub fn without_envelope(mut self) -> Self {
        self.include_envelope = false;
        self
    }

    /// Empirical model for tungsten emissivity as a function of wavelength.
    /// Source: Robert D. Larrabee, *Journal of the Optical Society of America*,
    /// Vol. 49, No. 6, pp. 619–623 (1959).
    ///
    /// Notes:
    /// - The wavelength parameter `wl` may be given in nanometers or meters.
    /// - Equation (6) from the paper is used. Although the text refers to units of
    ///   “millimicrons” (an older name for nanometers, 1 mµ = 1 nm), the calculation
    ///   requires the wavelength to be expressed in micrometers for correct results.
    /// - Table 6 provides emissivity values for tungsten across wavelengths
    ///   from 310 to 800 nm and temperatures from 1600 K to 2400 K.
    ///   The author notes, however, that the actual measured wavelength ranges
    ///   were narrower than the tabulated values.
    pub fn tungsten_emissivity_larrabee(&self, wl: f64) -> f64 {
        let l = wavelength(wl); // wavelenght in m
        let wavelength_um = wavelength(l) * 1E6; // convert to micrometers
        if wavelength_um < 450.0 {
            0.6075 - 0.3 * wavelength_um - 0.3265E-4 * self.filament_temperature
                + 0.59E-4 * wavelength_um * self.filament_temperature
        } else if wavelength_um < 680.0 {
            0.4655 + 0.01558 * wavelength_um + 0.2675E-4 * self.filament_temperature
                - 0.7305E-4 * wavelength_um * self.filament_temperature
        } else {
            0.6552 - 0.2633 * wavelength_um - 0.7333E-4 * self.filament_temperature
                + 0.7417E-4 * wavelength_um * self.filament_temperature
        }
    }

    /// Empirical model for tungsten emissivity as a function of wavelength and temperature.
    /// This model is based on the work of J.C. De Vos (1954) and further refined by M. Pon and Jan P. Hessler (1984).
    /// # Inputs
    /// * `wl`: Wavelength in nanometers or meters.
    ///
    /// # Notes
    /// * Use a recycling factor of 0.0 to get the original De Vos emissivity.
    ///   The default value of 0.5 is the value as estimated by Ojanen et all.
    pub fn tungsten_emissivity_devos(&self, wl: f64) -> f64 {
        const COEFFS: [[f64; 8]; 8] = [
            [
                380.0, 0.47245, -0.0155, -0.0086, -0.0229, 0.0000, -2.860, 0.000,
            ],
            [
                450.0, 0.46361, -0.0172, -0.1304, 0.0000, 0.0000, 0.520, 0.000,
            ],
            [
                530.0, 0.45549, -0.0173, -0.1150, 0.0000, 0.0000, -0.500, 0.000,
            ],
            [
                610.0, 0.44297, -0.0177, -0.1482, 0.0000, 0.0000, 0.723, 0.000,
            ],
            [
                700.0, 0.43151, -0.0207, -0.1441, -0.0551, 0.0000, -0.278, -0.190,
            ],
            [
                850.0, 0.40610, -0.0259, -0.1889, 0.0087, 0.0290, -0.126, 0.246,
            ],
            [
                1270.0, 0.32835, 0.0000, -0.1686, 0.0737, 0.0000, 0.046, 0.016,
            ],
            [
                2100.0, 0.22631, 0.0431, -0.0829, 0.0241, 0.0000, 0.040, -0.026,
            ],
        ];

        let ew_lt = |coef: [f64; 8], wl_um: f64| -> f64 {
            const T0: f64 = 2200.0;
            let dt = (self.filament_temperature - T0) / 1000.0;
            let [l0, a0, a1, b0, b1, b2, c0, c1] = coef;
            let dwl = wl_um - l0 / 1000.0; // in micrometers
            let ew =
                a0 + a1 * dt + (b0 + b1 * dt + b2 * dt * dt) * dwl + (c0 + c1 * dt) * dwl * dwl;

            // account for filament reflections
            ew / (1.0 - self.recycling_factor * (1.0 - ew)) // account for recycling factor
        };

        let wl = wavelength(wl); // wavelenght in m
        let wl_um = (wl) * 1E6; // convert to nanometers
        let wl_nm = (wl) * 1E9; // convert to nanometers

        if wl_nm < 340.0 {
            panic!("Wavelength too low for tungsten_emissivity_devos")
        };
        if wl_nm < 420.0 {
            return ew_lt(COEFFS[0], wl_um);
        };
        if wl_nm < 480.0 {
            return ew_lt(COEFFS[1], wl_um);
        };
        if wl_nm < 580.0 {
            return ew_lt(COEFFS[2], wl_um);
        };
        if wl_nm < 640.0 {
            return ew_lt(COEFFS[3], wl_um);
        };
        if wl_nm < 760.0 {
            return ew_lt(COEFFS[4], wl_um);
        };
        if wl_nm < 940.0 {
            return ew_lt(COEFFS[5], wl_um);
        };
        if wl_nm < 1600.0 {
            return ew_lt(COEFFS[6], wl_um);
        };
        if wl_nm <= 2600.0 {
            return ew_lt(COEFFS[7], wl_um);
        };
        panic!("Wavelength too high for tungsten_emissivity_devos");
    }

    /// Empirical model for the transmission of the quartz envelope as a function of wavelength.
    /// # Inputs
    /// * `wl`: Wavelength in nanometers or meters.
    ///
    /// # Notes
    /// * From Fig 1, from Ojanen (2012).
    ///   The model is valid for wavelengths from 400 nm to 800 nm.
    ///   Temperature dependence is ignored.
    /// * The transmission data in Fig 1 was digitized using WebPlotDigitizer, and
    ///   a 4th order polynomial was fitted to the data using a least squares fit.
    ///
    /// # References
    /// * M. Ojanen et all., Metrologia 49 (2012) S53–S58
    /// * [WebPlotDigitizer](https://automeris.io/WebPlotDigitizer)
    ///
    pub fn envelope_transmission(&self, wl: f64) -> f64 {
        const CF: [f64; 5] = [0.6, 2.17, -4.49, 4.25, -1.53];
        let l = wavelength(wl); // wavelength in m
        let wl_um = wavelength(l) * 1E6; // convert to micrometers

        // Horner's method
        CF.into_iter()
            .rev()
            .fold(0.0, |acc, coeff| acc * wl_um + coeff)
    }
}

impl Light for QTH {
    fn spectrum(&self) -> Cow<'_, Spectrum> {
        let planck = Planck::new(self.filament_temperature);
        let mut values = [0.0; NS];

        for (i, l) in SPECTRUM_WAVELENGTH_RANGE.enumerate() {
            let wl = wavelength(l);
            let radiance = planck.at_wavelength(wl); // convert nm to m
            let emissivity = self.tungsten_emissivity_devos(wl);
            if self.include_envelope {
                let transmission = self.envelope_transmission(wl);
                values[i] = radiance * emissivity * transmission * self.scale;
            } else {
                values[i] = radiance * emissivity * self.scale;
            }
        }
        Cow::Owned(Spectrum::new(values))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fel_emissivity() {
        // Table II Larrabee 1959
        let wants = [
            [1600.0, 310.0, 0.479],
            [1800.0, 390.0, 0.473],
            [2000.0, 500.0, 0.453],
            [2200.0, 600.0, 0.440],
            [2400.0, 700.0, 0.419],
            [1600.0, 800.0, 0.422],
        ];
        for [t, wl, want] in wants {
            let fel = QTH::from_filament_temperature(t);
            let emissivity = fel.tungsten_emissivity_larrabee(wl);
            approx::assert_abs_diff_eq!(emissivity, want, epsilon = 4E-2);
        }
    }

    #[test]
    fn test_fel_emissivity_devos() {
        // Table II DeVos 1980
        let wants = [
            [1600.0, 1275.5, 0.32725],
            [1800.0, 1275.5, 0.32725],
            [2000.0, 1275.5, 0.32725],
            [1600.0, 2600.0, 0.16567],
            [1800.0, 2600.0, 0.1754],
            [2000.0, 2600.0, 0.18513],
            [2200.0, 2600.0, 0.19486],
            [2400.0, 2600.0, 0.20459],
            [2600.0, 2600.0, 0.21432],
            [2800.0, 2600.0, 0.22405],
        ];
        for [temp, wl, want] in wants {
            let fel = QTH::from_filament_temperature(temp)
                .without_envelope()
                .without_recycling();
            let emissivity = fel.tungsten_emissivity_devos(wl);
            approx::assert_abs_diff_eq!(emissivity, want, epsilon = 1E-4);
        }
    }

    #[test]
    fn test_transmission() {
        let fel = QTH::from_filament_temperature(2800.0);
        let wants = [
            (400.0, 0.9806),
            (500.0, 0.9951),
            (600.0, 1.0026),
            (700.0, 1.0063),
            (800.0, 1.0084),
        ];
        for (wl, want) in wants {
            let t = fel.envelope_transmission(wl);
            approx::assert_abs_diff_eq!(t, want, epsilon = 4E-3);
        }
    }

    #[test]
    #[cfg(feature = "cct")]
    // Verify that the CCT calculated from the QTH model matches the input filament temperature
    // and that the polynomial approximation is accurate enough
    fn test_qth_cct() {
        use crate::traits::Light;
        for temp in (1600..=3200).step_by(100) {
            let fel = QTH::from_filament_temperature(temp as f64);
            let cct = fel.cct().unwrap();
            println!("{temp}, {:.1}, {:.2}", cct.t(), cct.tint());
            let fel_cct = QTH::new(cct.t());
            approx::assert_abs_diff_eq!(fel_cct.filament_temperature, temp as f64, epsilon = 2.0);
        }
    }
}
