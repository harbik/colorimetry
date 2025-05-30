use wasm_bindgen::prelude::wasm_bindgen;


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
