use nalgebra::Vector3;
use wasm_bindgen::prelude::wasm_bindgen;

use super::{achromatic_rsp, Cam, M16, MCAT02, MCAT02INV, MHPE};

#[wasm_bindgen]
#[derive(Clone, Copy, Debug)]
pub struct ViewConditions {
    /// Degree of Adaptation, if omitted, formula 4.3 of CIE248:2022 is used.``
    pub dopt: Option<f64>,

    pub f: f64,

    /// Adaptation Luminance, in cd/m2.
    /// La = Lw/5, with Lw: luminance of a perfect white object
    pub la: f64,
    pub nc: f64,
    pub yb: f64,
    pub c: f64,
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

    pub fn reference_values(&self, xyzn: Vector3<f64>, cam: Cam) -> super::ReferenceValues {
        let mut rgb_w = match cam {
            Cam::CieCam16 => M16 * xyzn,
            Cam::CieCam02 => MCAT02 * xyzn,
        };

        let vcd = self.dd();
        let yw = xyzn[1];
        let d_rgb = rgb_w.map(|v| vcd * yw / v + 1.0 - vcd);
        let n = self.yb / yw;
        let z = n.sqrt() + 1.48;
        let nbb = 0.725 * n.powf(-0.2);
        let ncb = nbb;
        rgb_w.component_mul_assign(&Vector3::from(d_rgb));
        if cam == Cam::CieCam02 {
            rgb_w = MCAT02INV * rgb_w;
            rgb_w = MHPE * rgb_w;
        }
        let qu = 150f64.max(rgb_w[0].max(rgb_w[1]).max(rgb_w[2]));
        rgb_w.apply(|q| self.lum_adapt(q, 0.26, qu));
        let aw = achromatic_rsp(rgb_w, nbb);

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
    pub fn f_l(&self) -> f64 {
        let k = self.k();
        k.powi(4) * self.la + (1. - k.powi(4)).powi(2) / 10. * (5.0 * self.la).powf(1. / 3.)
    }

    pub fn dd(&self) -> f64 {
        if let Some(d) = self.dopt {
            d.clamp(0.0, 1.0)
        } else {
            (self.f * (1.0 - (1.0 / 3.6) * ((-1.0 * self.la - 42.0) / 92.0).exp())).clamp(0.0, 1.0)
        }
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
    pub fn lum_adapt(&self, q: &mut f64, ql: f64, qu: f64) {
        let fl = self.f_l();

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

pub const TM30VC: ViewConditions = ViewConditions {
    yb: 20.0,
    c: 0.69,
    nc: 1.0,
    f: 1.0,
    la: 100.0,
    dopt: Some(1.0),
};

/// Values according to Table 1, record 2, CIE 248:2022
pub const CIE_HOME_DISPLAY: ViewConditions = ViewConditions {
    yb: 20.0,
    c: 0.59,
    nc: 0.9,
    f: 0.8,
    la: 16.0,
    dopt: None,
};
