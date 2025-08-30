// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

use nalgebra::Vector3;

use super::{Cam, M16, MCAT02, MCAT02INV, MHPE};

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[derive(Clone, Copy, Debug)]
/// CIECAM viewing conditions.
///
/// The ViewConditions as recommended by CIE248:2022 are provided for various scenarios as constants, and are included as:
/// - [`CIE248_CABINET`] Viewing a surface in a cabinet
/// - [`CIE248_HOME_SCREEN`] Viewing a self-luminous display at home
/// - [`CIE248_PROJECTED_DARK`] Viewing projected images in a darkened room
/// - [`CIE248_OFFICE_SCREEN`] Viewing a self-luminous display under office illumination
///  
/// The TM30 and Color Fidelity ViewConditions are provided as [`TM30VC`].
pub struct ViewConditions {
    /// Adaptation Luminance, in cd/m2.
    /// La = Lw/5, with Lw: luminance of a perfect white object
    la: f64,
    // The Y value of the “grey surround” immediately around the stimulus, expressed in the same scale as Yw.
    yb: f64,
    c: f64,
    nc: f64,
    f: f64,
    dopt: Option<f64>,
}

impl ViewConditions {
    /// Creates a new instance of `ViewConditions`.
    ///
    /// # Arguments
    /// * `la` - Adaptation luminance in cd/m².
    /// * `yb` - Relative background luminance, as a fraction of the white point’s Yn (cd/m²) value. Typically yb = 20, if Yn = 100.
    /// * `c` - Impact of surround on chroma/contrast.
    /// * `nc` - Chromatic induction factor.
    /// * `f` - Surround factor, representing the degree of adaptation.
    /// * `dopt` - Degree of adaptation, if omitted, formula 4.3 of CIE248:2022 is used.
    ///
    /// # Returns
    /// A new instance of `ViewConditions`.
    pub const fn new(
        la: f64,
        yb: f64,
        c: f64,
        nc: f64,
        f: f64,
        dopt: Option<f64>,
    ) -> ViewConditions {
        ViewConditions {
            yb,
            f,
            nc,
            c,
            la,
            dopt,
        }
    }
    /// Creates a new `ViewConditions` instance with average surround settings,
    /// as defined in CIE 248:2022.
    ///
    /// # Arguments
    /// * `la` – Adaptation luminance in cd/m².
    ///
    /// # Returns
    /// A `ViewConditions` configured for an average surround environment.
    ///
    /// # Notes
    /// The average surround condition corresponds to a surround ratio greater than, as specified in
    /// Annex B of CIE 248:2022. It is typically used in controlled environments such as color viewing
    /// cabinets, where lighting is consistent and neutral.  This can also be used for general lighting
    /// applications, looking at a color patch not lit by a spot light.
    ///
    /// This method sets `c = 0.69`, `nc = 1.0`, and `f = 1.0`, in line with CIE 248:2022
    /// recommendations for average surround viewing conditions.
    pub const fn average_surround(la: f64) -> ViewConditions {
        ViewConditions {
            la,
            yb: 20.0,
            c: 0.69,
            nc: 1.0,
            f: 1.0,
            dopt: None,
        }
    }

    /// Creates a new `ViewConditions` instance with dim surround settings,
    /// as defined in CIE 248:2022.
    ///
    /// # Arguments
    /// * `la` – Adaptation luminance in cd/m².
    ///
    /// # Returns
    /// A `ViewConditions` configured for a dim surround environment.
    ///
    /// # Notes
    /// The dim surround condition corresponds to a surround ratio between 0.0 and 0.2,
    /// as described in Annex B of CIE 248:2022. It is typically used for viewing
    /// self-luminous displays in home settings.
    ///
    /// For example, if a display emits 80 cd/m² (white luminance), and the background wall
    /// is illuminated to about 30 lux—roughly 10 cd/m²—the surround ratio is 10 / 80 = 0.125,
    /// which falls within the dim range.
    ///
    /// This method sets `c = 0.59`, `nc = 0.9`, and `f = 0.9`, following the dim surround
    /// parameters specified by CIE 248:2022 for color appearance modeling.
    pub const fn dim_surround(la: f64) -> ViewConditions {
        ViewConditions {
            la,
            yb: 20.0,
            c: 0.59,
            nc: 0.9,
            f: 0.9,
            dopt: None,
        }
    }

    /// Creates a new `ViewConditions` instance with dark surround settings,
    /// as defined in CIE 248:2022.
    ///
    /// # Arguments
    /// * `la` – Adaptation luminance in cd/m².
    ///
    /// # Returns
    /// A `ViewConditions` configured for a dark surround environment.
    ///
    /// # Notes
    /// The dark surround condition corresponds to a surround ratio approaching zero, as described in
    /// Annex B of CIE 248:2022.  This setting is used for viewing projected or cinematic images in
    /// darkened rooms where ambient luminance is minimal.
    ///
    /// For example, if the display luminance is 80 cd/m² and the room walls appear black, with a
    /// luminance less than 1 cd/m² of luminance, the surround ratio would be around 0.0125.
    ///
    /// This method sets `c = 0.525`, `nc = 0.8`, and `f = 0.8`,
    /// matching the parameters recommended by CIE 248:2022
    /// for color appearance modeling in dark environments.
    pub const fn dark_surround(la: f64) -> ViewConditions {
        ViewConditions {
            la,
            yb: 20.0,
            c: 0.525,
            nc: 0.8,
            f: 0.8,
            dopt: None,
        }
    }

    /// Sets the degree of adaptation.
    /// This method allows you to specify the degree of adaptation directly.
    /// The value should be between 0.0 and 1.0, where 0.0 represents no adaptation
    /// and 1.0 represents full adaptation.
    pub fn set_degree_of_adaptation(&self, d: f64) -> ViewConditions {
        ViewConditions {
            dopt: Some(d.clamp(0.0, 1.0)),
            ..*self
        }
    }

    ///  Impact of surround on chroma/contrast  
    pub fn impact_of_surround(&self) -> f64 {
        self.c
    }

    ///  Chromatic induction factor
    pub fn chromatic_induction_factor(&self) -> f64 {
        self.nc
    }

    ///  Adapting luminance in cd/m²  
    pub fn adapting_luminance(&self) -> f64 {
        self.la
    }

    ///  relative background luminance, as a fraction of the white point’s Yn (cd/m2) value  
    pub fn relative_background_luminance(&self) -> f64 {
        self.yb
    }

    ///  surround factor, representing the degree of adaptation
    pub fn surround_factor(&self) -> f64 {
        self.f
    }

    pub fn reference_values(&self, xyzn: Vector3<f64>, cam: Cam) -> super::ReferenceValues {
        let mut rgb_w = match cam {
            Cam::CieCam16 => M16 * xyzn,
            Cam::CieCam02 => MCAT02 * xyzn,
        };

        let vcd = self.degree_of_adaptation();
        let yw = xyzn[1];
        let d_rgb = rgb_w.map(|v| vcd * yw / v + 1.0 - vcd);
        let n = self.yb / yw;
        let z = n.sqrt() + 1.48;
        let nbb = 0.725 * n.powf(-0.2);
        let ncb = nbb;
        rgb_w.component_mul_assign(&Vector3::from(d_rgb));
        let qu = 150f64.max(rgb_w[0].max(rgb_w[1]).max(rgb_w[2]));
        match cam {
            Cam::CieCam02 => {
                rgb_w = MCAT02INV * rgb_w;
                rgb_w = MHPE * rgb_w;
                rgb_w.apply(|v| self.lum_adapt02(v));
            }
            Cam::CieCam16 => {
                rgb_w.apply(|q| self.lum_adapt16(q, 0.26, qu));
            }
        }
        let aw = super::achromatic_rsp(rgb_w, nbb);

        super::ReferenceValues {
            n,
            z,
            nbb,
            ncb,
            d_rgb: d_rgb.into(),
            aw,
            qu,
        }
    }

    #[inline]
    pub fn k(&self) -> f64 {
        1. / (5. * self.la + 1.)
    }

    #[inline]
    pub fn luminance_level_adaptation_factor(&self) -> f64 {
        let k = self.k();
        k.powi(4) * self.la + (1. - k.powi(4)).powi(2) / 10. * (5.0 * self.la).powf(1. / 3.)
    }

    /// Degree of Adaptation, if omitted, formula 4.3 of CIE248:2022 is used.``
    pub fn degree_of_adaptation(&self) -> f64 {
        if let Some(d) = self.dopt {
            d.clamp(0.0, 1.0)
        } else {
            (self.f * (1.0 - (1.0 / 3.6) * ((-self.la - 42.0) / 92.0).exp())).clamp(0.0, 1.0)
        }
    }

    pub fn lum_adapt02(&self, v: &mut f64) {
        let t = (self.luminance_level_adaptation_factor() * *v / 100.).powf(0.42);
        *v = v.signum() * 400. * t / (27.13 + t) + 0.1;
    }

    /// Hyperbolic post-adaptation response compression function
    ///
    /// As used in CIECAM02 and CAM16
    ///
    /// Modified hyperbolic post-adaptation response compression function
    ///
    /// Updated to fix some issues with the function as used in CIECAM02 and CAM16.
    /// See Section 3.2 in CIE248:2022.
    /// CIE recommends to use:
    /// l = 0.26;
    /// u = max(150.0, Rwc, Gwc, Bwc);
    pub fn lum_adapt16(&self, q: &mut f64, ql: f64, qu: f64) {
        let fl = self.luminance_level_adaptation_factor();

        // CIECAM16 eq 3.4
        let f = |q: f64| -> f64 {
            let t = (fl * q / 100.).powf(0.42);
            400.0 * t / (27.13 + t)
        };

        // CIECAM16 eq 3.5
        let fp = |qu: f64| -> f64 {
            let t = fl * qu / 100.;
            let den = 1.68 * 27.13 * fl * t.powf(-0.58);
            let nom = (27.13 + t.powf(0.42)).powi(2);
            den / nom
        };

        // CIECAM16 eq 3.3
        if *q < ql {
            *q = f(ql) * *q / ql;
        } else if *q > qu {
            *q = f(qu) * fp(qu) * (*q - qu);
        } else {
            *q = f(*q);
        };

        // CIECAM16 eq 3.2
        *q += 0.1;
    }
}

impl Default for ViewConditions {
    fn default() -> Self {
        Self::average_surround(100.0)
    }
}

/// ViewConditions as used for TM30 and Color FIdelity Calculations
/// (ANSI/TM-30-20-2020)
pub const TM30VC: ViewConditions = ViewConditions::new(100.0, 20.0, 0.69, 1.0, 1.0, Some(1.0));

/// Values according to Table 1, CIE 248:2022
///    
pub const CIE248_CABINET: ViewConditions = ViewConditions::average_surround(
    63.7, // adaptation luminance (cd·m⁻²) = 0.2 * Lw (perfect white object)
);

/// 2) Viewing a self-luminous display at home
pub const CIE248_HOME_SCREEN: ViewConditions = ViewConditions::dim_surround(
    16.0, // adaptation luminance (cd·m⁻²)
);

/// 3) Viewing projected images in a darkened room
pub const CIE248_PROJECTED_DARK: ViewConditions = ViewConditions::dark_surround(
    30.0, // adaptation luminance (cd·m⁻²)
);

/// 4) Viewing a self-luminous display under office illumination
pub const CIE248_OFFICE_SCREEN: ViewConditions = ViewConditions::average_surround(
    16.0, // adaptation luminance (cd·m⁻²)
);
