#[derive(Clone)]
/// General representation of the RGB encoding and decoding functions used by color spaces.
/// The function type is determined by the number of parameters supplied.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, PartialEq)]
pub struct GammaCurve {
    pub(super) category: usize,
    pub(super) data: [f64; 7],
}

impl GammaCurve {
    pub const fn new(category: usize, data: [f64; 7]) -> Self {
        Self { category, data }
    }

    // from rgb coordinates to xyz, gamma > 1.0
    pub fn decode(&self, x: f64) -> f64 {
        if !(0.0..=1.0).contains(&x) {
            f64::NAN
        } else {
            match self.category {
                0 => x,
                1 => {
                    let g = self.data[0];
                    x.powf(g)
                }
                3 => {
                    let [g, a, b, ..] = self.data;
                    if x >= -b / a {
                        (a * x + b).powf(g)
                    } else {
                        0.0
                    }
                }
                4 => {
                    let [g, a, b, c, ..] = self.data;
                    if x >= -b / a {
                        (a * x + b).powf(g) + c
                    } else {
                        c
                    }
                }
                5 => {
                    let [g, a, b, c, d, ..] = self.data;
                    if x >= d {
                        (a * x + b).powf(g)
                    } else {
                        c * x
                    }
                }
                7 => {
                    let [g, a, b, c, d, e, f] = self.data;
                    if x >= d {
                        (a * x + b).powf(g) + e
                    } else {
                        c * x + f
                    }
                }
                _ => f64::NAN,
            }
        }
    }

    // from xyz coordinates to rgb, gamma < 1.0
    pub fn encode(&self, x: f64) -> f64 {
        let x = x.clamp(0.0, 1.0);
        match self.category {
            0 => x, // no gamma correction
            1 => {
                let g = self.data[0];
                x.powf(1.0 / g)
            }
            3 => {
                let [g, a, b, ..] = self.data;
                if x >= 0.0 {
                    (x.powf(1.0 / g) - b) / a
                } else {
                    f64::NAN
                }
            }
            5 => {
                let [g, a, b, c, d, ..] = self.data;
                if x >= d * c {
                    (x.powf(1.0 / g) - b) / a
                } else {
                    x / c
                }
            }
            _ => f64::NAN,
        }
    }

    pub fn encode_u8(&self, v: f64, vmin: Option<f64>, vmax: Option<f64>) -> u8 {
        match (v, vmin, vmax) {
            (v, None, None) => (self.encode(v) * 255.0).round() as u8, // clamped
            (v, Some(vmin), None) => (self.encode((v - vmin) / (1.0 - vmin)) * 255.0).round() as u8,
            (v, Some(vmin), Some(vmax)) => {
                (self.encode((v - vmin) / (vmax - vmin)) * 255.0).round() as u8
            }
            (v, None, Some(vmax)) => (self.encode(v / vmax) * 255.0).round() as u8,
        }
    }
}

#[test]
fn test_gamma_srgb() {
    let gc = GammaCurve::new(
        5,
        [
            2.4,
            1.0 / 1.055,
            0.055 / 1.055,
            1.0 / 12.92,
            0.04045,
            0.0,
            0.0,
        ],
    );
    println!("{} {}", gc.encode(0.0), gc.encode_u8(0.0, None, None));
    println!("{} {}", gc.encode(1.0), gc.encode_u8(1.0, None, None));
    println!(
        "{} {}",
        gc.encode(1.0000001),
        gc.encode_u8(1.000001, None, None)
    );
    println!(
        "{} {}",
        gc.encode(-0.0000001),
        gc.encode_u8(-0.000001, None, None)
    );
    println!(
        "{} {}",
        gc.encode(1.0000001),
        gc.encode_u8(1.000001, Some(0.0), Some(1.0))
    );
}
