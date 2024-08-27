

#[derive(Clone)]
/// General representation of the RGB encoding and decoding functions used by color spaces.
/// The function type is determined by the number of parameters supplied.
pub struct GammaCurve {
    p: Vec<f64>
}

impl GammaCurve{
    pub fn new(p: Vec<f64>) -> Self { Self { p } 
}


    // from rgb coordinates to xyz, gamma > 1.0
    pub fn decode(&self, x: f64) -> f64 { 
        if x<0.0 || x>1.0 { 
            f64::NAN
        } else {
            match self.p.len() {
                1 => {
                    let g = self.p[0];
                    x.powf(g)
                },
                3 => {
                    let [g, a, b] = self.p[..] else {panic!()}; // never reached
                    if x>= -b/a {
                        (a*x + b).powf(g)
                    } else {
                        0.0
                    }
                },
                4 => {
                    let [g, a, b, c] = self.p[..] else {panic!()};
                    if x>= -b/a {
                        (a*x + b).powf(g) + c
                    } else {
                       c 
                    }
                },
                5 => {
                    let [g, a, b, c, d] = self.p[..] else {panic!()};
                    if x>= d {
                        (a*x + b).powf(g)
                    } else {
                        c*x
                    }
                },
                7 => {
                    let [g,a,b, c, d, e, f] = self.p[..] else {panic!()};
                    if x>= d {
                        (a*x + b).powf(g) + e
                    } else {
                        c*x + f
                    }
                },
                _ => f64::NAN
            }
        }
    }

    // from xyz coordinates to rgb, gamma < 1.0
    pub fn encode(&self, x:f64) -> f64 {
        let x = x.clamp(0.0, 1.0);
        match self.p.len() {
            1 => {
                let g = self.p[0];
                x.powf(1.0/g)
            },
            3 => {
                let [g, a, b] = self.p[..] else {panic!()}; // never reached
                todo!()
            },
            4 => {
                let [g, a, b, c] = self.p[..] else {panic!()};
                todo!()
            },
            5 => {
                let [g, a, b, c, d] = self.p[..] else {panic!()};
                if x>= d * c {
                    (x.powf(1.0/g) - b) / a
                } else {
                    x / c
                }
            },
            7 => {
                let [g,a,b, c, d, e, f] = self.p[..] else {panic!()};
                todo!()
            },
            _ => f64::NAN

        }

    }

    pub fn encode_u8(&self, v: f64, vmin: Option<f64>, vmax: Option<f64>) -> u8 {
        match (v, vmin, vmax) {
            (v, None, None) =>  (self.encode(v) * 255.0).round() as u8, // clamped
            (v, Some(vmin), None) => (self.encode((v-vmin)/(1.0-vmin)) * 255.0).round() as u8,
            (v, Some(vmin), Some(vmax)) =>  (self.encode((v-vmin)/(vmax-vmin)) * 255.0).round() as u8,
            (v, None, Some(vmax)) =>  (self.encode(v/vmax) * 255.0).round() as u8,
        }
    }

}

#[test]
fn test_gamma_srgb(){
    let gc = GammaCurve::new(vec![2.4, 1.0/1.055, 0.055/1.055, 1.0/12.92, 0.04045]);
    println!("{} {}", gc.encode(0.0), gc.encode_u8(0.0, None, None));
    println!("{} {}", gc.encode(1.0), gc.encode_u8(1.0, None, None));
    println!("{} {}", gc.encode(1.0000001), gc.encode_u8(1.000001, None, None));
    println!("{} {}", gc.encode(-0.0000001), gc.encode_u8(-0.000001, None, None));
    println!("{} {}", gc.encode(1.0000001), gc.encode_u8(1.000001, Some(0.0), Some(1.0)));

}
