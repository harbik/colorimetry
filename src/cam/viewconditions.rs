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
    dopt: Option<f64>,

    f: f64,

    /// Adaptation Luminance, in cd/m2.
    /// La = Lw/5, with Lw: luminance of a perfect white object
    la: f64,
    nc: f64,
    yb: f64,
    c: f64,
}

impl ViewConditions {
    pub fn new(yb: f64, f: f64, nc: f64, c: f64, la: f64, dopt: Option<f64>) -> ViewConditions {
        ViewConditions {
            yb,
            f,
            nc,
            c,
            la,
            dopt,
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

    ///  relative background luminance, as a fraction of the white point’s Y value  
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
            (self.f * (1.0 - (1.0 / 3.6) * ((-1.0 * self.la - 42.0) / 92.0).exp())).clamp(0.0, 1.0)
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
        Self {
            yb: 20.0,
            c: 0.69,
            nc: 1.0,
            f: 1.0,
            la: 100.0,
            dopt: None,
        }
    }
}

/// ViewConditions as used for TM30 and Color FIdelity Calculations
/// (ANSI/TM-30-20-2020)
pub const TM30VC: ViewConditions = ViewConditions {
    yb: 20.0,
    c: 0.69,
    nc: 1.0,
    f: 1.0,
    la: 100.0,
    dopt: Some(1.0),
};

/// Values according to Table 1, CIE 248:2022
///    
/// 1) Surface colour evaluation in a color viewing cabinet
///
///   – Adopted white point: Viewing-cabinet white ⇒ Yₙ = 318.3 cd·m⁻²  
///   – Adaptation luminance: L_A = 63.7 cd·m⁻²  
///   – Surround: Average ⇒ (f = 1.0, c = 0.69, nc = 1.0)  
///   – Surround ratio: S_R = 1 ⇒ background luminance Y_b = 1 × 318.3 = 318.3 cd·m⁻²
///
pub const CIE248_CABINET: ViewConditions = ViewConditions {
    yb: 318.3,  // = 1 × 318.3
    f: 1.0,     // “average” surround
    nc: 1.0,    // “average” surround
    c: 0.69,    // “average” surround
    la: 63.7,   // adaptation luminance (cd·m⁻²)
    dopt: None, // let the library compute D via CIE248:2022
};

/// 2) Viewing a self-luminous display at home
///
///   – Adopted white point: between display white and ambient white ⇒ Yₙ = 80 cd·m⁻²  
///   – Adaptation luminance: L_A = 16 cd·m⁻²  
///   – Surround: Dim ⇒ (f = 0.9, c = 0.59, nc = 0.9)  
///   – Surround ratio: S_R = 0.15 ⇒ Y_b = 0.15 × 80 = 12 cd·m⁻²
///
pub const CIE248_HOME_SCREEN: ViewConditions = ViewConditions {
    yb: 12.0, // = 0.15 × 80
    f: 0.9,   // “dim” surround
    nc: 0.9,  // “dim” surround
    c: 0.59,  // “dim” surround
    la: 16.0, // adaptation luminance (cd·m⁻²)
    dopt: None,
};

/// 3) Viewing projected images in a darkened room
///
///   – Adopted white point: between projector white and ambient white ⇒ Yₙ = 150 cd·m⁻²  
///   – Adaptation luminance: L_A = 30 cd·m⁻²  
///   – Surround: Dark ⇒ (f = 0.8, c = 0.525, nc = 0.8)  
///   – Surround ratio: S_R = 0 ⇒ Y_b = 0 × 150 = 0 cd·m⁻²
///
pub const CIE248_PROJECTED_DARK: ViewConditions = ViewConditions {
    yb: 0.0,  // = 0 × 150
    f: 0.8,   // “dark” surround
    nc: 0.8,  // “dark” surround
    c: 0.525, // “dark” surround
    la: 30.0, // adaptation luminance (cd·m⁻²)
    dopt: None,
};

/// 4) Viewing a self-luminous display under office illumination
///
///   – Adopted white point: between display white and office white ⇒ Yₙ = 80 cd·m⁻²  
///   – Adaptation luminance: L_A = 16 cd·m⁻²  
///   – Surround: Average ⇒ (f = 1.0, c = 0.69, nc = 1.0)  
///   – Surround ratio: S_R = 2 ⇒ Y_b = 2 × 80 = 160 cd·m⁻²
///
pub const CIE248_OFFICE_SCREEN: ViewConditions = ViewConditions {
    yb: 160.0, // = 2 × 80
    f: 1.0,    // “average” surround
    nc: 1.0,   // “average” surround
    c: 0.69,   // “average” surround
    la: 16.0,  // adaptation luminance (cd·m⁻²)
    dopt: None,
};
