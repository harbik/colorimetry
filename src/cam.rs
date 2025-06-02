mod viewconditions;
pub use viewconditions::ViewConditions;

mod cam16;
pub use crate::cam::cam16::CieCam16;

mod cam02;
pub use crate::cam::cam02::CieCam02;

use crate::{observer::Observer, xyz::XYZ};
use nalgebra::{matrix, vector, SMatrix, Vector3};
use std::f64::consts::PI;

#[derive(Debug)]
pub struct CamJCh {
    /// Colorimetric Observer used.
    /// The standard requires use of the CIE1931 standard observer.
    observer: Observer,

    /// Correlates of Lightness, Chroma, and hue-angle
    jch: Vector3<f64>,

    /// Tristimulus values of the reference white beng use
    xyzn: Vector3<f64>,

    /// Viewing Conditions
    vc: ViewConditions,
}

/// Values used in CieCam forward and backward transformations, only dependent on the view
/// conditions and the reference white.
/// For a set of colors, viewed under the same conditions, these have to be calculated only once.
#[derive(Clone, Copy, Debug)]
pub struct ReferenceValues {
    n: f64,
    z: f64,
    nbb: f64,
    ncb: f64,
    d_rgb: [f64; 3],
    aw: f64,
    qu: f64, // see lum_adapt
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cam {
    CieCam02,
    CieCam16,
}

impl CamJCh {
    const P1C: f64 = 50_000.0 / 13.0;
    const P3: f64 = 21.0 / 20.0;
    const NOM: f64 = 1.0; // is listed in the standard as (2.0+Self::P3)*460.0/C16_3; // this is 1!! But his is in Luo step 3 as part of nominators
    const RCPR_9: f64 = 1.0 / 0.9;
    const DEN1: f64 = ((2.0 + Self::P3) * 220.0) / 1403.0;
    const DEN2: f64 = (Self::P3 * 6300.0 - 27.0) / 1403.0;

    pub fn new(
        observer: Observer,
        jch: Vector3<f64>,
        xyzn: Vector3<f64>,
        vc: ViewConditions,
    ) -> Self {
        Self {
            observer,
            jch,
            xyzn,
            vc,
        }
    }

    pub fn from_xyz(
        xyz: XYZ,
        xyzn: XYZ,
        vc: ViewConditions,
        cam: Cam,
    ) -> Result<Self, crate::Error> {
        if xyz.observer != xyzn.observer {
            return Err(crate::Error::RequireSameObserver);
        }
        let xyz_vec = xyz.xyz;
        let xyzn_vec = xyzn.xyz;

        let ReferenceValues {
            n,
            z,
            nbb,
            ncb,
            d_rgb,
            aw,
            qu,
        } = vc.reference_values(xyzn.xyz, cam);
        //  let vcdd = vc.dd();
        //  let vcfl = vc.f_l();

        let mut rgb = match cam {
            Cam::CieCam16 => M16 * xyz_vec,
            Cam::CieCam02 => MCAT02 * xyz_vec,
        };

        // rgb_C
        rgb.component_mul_assign(&Vector3::from(d_rgb));

        // rgb_pw & rgb_pa
        match cam {
            Cam::CieCam16 => {
                rgb.apply(|v| vc.lum_adapt16(v, 0.26, qu));
            }
            Cam::CieCam02 => {
                rgb = MCAT02INV * rgb;
                rgb = MHPE * rgb;
                rgb.apply(|v| vc.lum_adapt02(v));
            }
        }

        let ca = rgb[0] - 12.0 / 11.0 * rgb[1] + rgb[2] / 11.0;
        let cb = (rgb[0] + rgb[1] - 2. * rgb[2]) / 9.0;

        // calculate h in radians
        let mut h = cb.atan2(ca);
        if h < 0.0 {
            h += 2.0 * PI;
        }

        // calculate J = jj
        let jj = 100.0 * (achromatic_rsp(rgb, nbb) / aw).powf(vc.c * z);

        // calculate C = cc
        let et = 0.25f64 * ((h + 2.0).cos() + 3.8);
        let t = (50000.0 / 13.0 * ncb * vc.nc * et * (ca * ca + cb * cb).sqrt())
            / (rgb[0] + rgb[1] + 21.0 / 20.0 * rgb[2]);
        let cc = t.powf(0.9) * (jj / 100.).sqrt() * (1.64 - (0.29f64).powf(n)).powf(0.73);

        Ok(Self {
            vc,
            jch: Vector3::new(jj, cc, h * 180.0 / PI),
            observer: xyz.observer,
            xyzn: xyzn_vec,
        })
    }

    pub fn xyz(
        &self,
        opt_xyzn: Option<XYZ>,
        opt_viewconditions: Option<ViewConditions>,
        cam: Cam,
    ) -> Result<XYZ, crate::Error> {
        let vc = opt_viewconditions.unwrap_or_default();
        let xyzn = if let Some(white) = opt_xyzn {
            if white.observer == self.observer {
                white.xyz
            } else {
                return Err(crate::Error::RequireSameObserver);
            }
        } else {
            self.xyzn
        };
        let ReferenceValues {
            n,
            z,
            nbb,
            ncb,
            d_rgb,
            aw,
            qu: _,
        } = vc.reference_values(xyzn, cam);
        let d_rgb_vec = Vector3::from(d_rgb);
        let &[lightness, chroma, hue_angle] = self.jch.as_ref();
        let t = (chroma / ((lightness / 100.0).sqrt() * (1.64 - 0.29f64.powf(n)).powf(0.73)))
            .powf(Self::RCPR_9);
        let p1 = (Self::P1C * vc.nc * ncb * eccentricity(hue_angle)) / t; // NaN if t=0, but OK, as check on t==0.0 if used
        let p2 = achromatic_response_from_lightness(aw, vc.c, z, lightness) / nbb + 0.305;
        let (a, b) = match hue_angle.to_radians().sin_cos() {
            (_, _) if t.is_nan() || t == 0.0 => (0.0, 0.0),
            (hs, hc) if hs.abs() >= hc.abs() => {
                let b = p2 * Self::NOM / (p1 / hs + Self::DEN1 * hc / hs + Self::DEN2);
                (b * hc / hs, b)
            }
            (hs, hc) => {
                let a = p2 * Self::NOM / (p1 / hc + Self::DEN1 + Self::DEN2 * hs / hc);
                (a, a * hs / hc)
            }
        };

        // rgb_a
        let m = matrix![ 460.0, 451.0, 288.0; 460.0, -891.0, -261.0; 460.0, -220.0, -6_300.0; ]
            / 1_403.0;
        let rgb_p = (m * vector![p2, a, b]).map(|x| inv_cone_adaptation(vc.f_l(), x)); // Step 4 & 5
        let rgb = rgb_p.component_div(&d_rgb_vec);

        let xyz = match cam {
            Cam::CieCam16 => M16INV * rgb,
            Cam::CieCam02 => {
                let rgb_c = (MCAT02 * MHPEINV * rgb_p).component_div(&d_rgb_vec); // Step 6 & 7
                MCAT02INV * rgb_c
            }
        };
        Ok(XYZ::from_vecs(xyz, self.observer))
    }
}

pub trait CamTransforms {
    const UCS_KL: f64 = 1.0;
    const UCS_C1: f64 = 0.007;
    const UCS_C2: f64 = 0.0228;

    fn jch_vec(&self) -> &Vector3<f64>;
    fn view_conditions(&self) -> &ViewConditions;
    fn observer(&self) -> Observer;
    fn xyzn(&self) -> &Vector3<f64>;

    /// Returns the JCh values of the color as an array.
    ///
    /// The JCh values are a lightness, chroma, and hue angle representation of the color.
    fn jch(&self) -> [f64; 3] {
        (*self.jch_vec()).into()
    }

    /// CIECAM16 “raw Jab” coordinates, wich are analogous to CIELAB’s *a*/*b* values.
    ///
    /// **Note:**  
    /// - These raw Jab coordinates are _not_ perceptually uniform.  
    /// - Euclidean distances in this space do **not** reliably correspond to visual color differences.  
    ///
    /// For a perceptually uniform version, use `jab_prime()`, which applies the CAM16-UCS non-linear  
    /// stretching (with constants `C1 = 0.007` and `C2 = 0.0228`) to produce `(J′, a′, b′)`.
    fn jab(&self) -> [f64; 3] {
        let &[jj, cc, h] = self.jch_vec().as_ref();
        let vc = self.view_conditions();
        let m = cc * vc.f_l().powf(0.25);
        let mprime = 1.0 / Self::UCS_C2 * (1.0 + Self::UCS_C2 * m).ln();
        [
            (1.0 + 100.0 * Self::UCS_C1) * jj / (1.0 + Self::UCS_C1 * jj),
            mprime * (h * PI / 180.0).cos(),
            mprime * (h * PI / 180.0).sin(),
        ]
    }

    /// This function returns the CIECAM16 UCS "Ja'b'" values, which are a perceptually uniform representation  
    /// of color.
    ///  
    /// **Note:**  
    /// - **J′**: Lightness  
    /// - **a′**, **b′**: Chromatic components, derived from the non-linear CAM16-UCS transformation applied  
    ///   to the raw Jab values.
    ///
    /// Euclidean distance between colors in this space more reliably matches human-perceived color differences.  
    /// The constants `C1 = 0.007` and `C2 = 0.0228` are used for the transformation to UCS type.  
    ///  
    /// For the raw, non-uniform Jab values, see the `jab()` function.
    fn jab_prime(&self) -> [f64; 3] {
        let &[jj, cc, h] = self.jch_vec().as_ref();
        let vc = self.view_conditions();
        let m = cc * vc.f_l().powf(0.25);
        let mprime = 1.0 / Self::UCS_C2 * (1.0 + Self::UCS_C2 * m).ln();
        [
            (1.0 + 100.0 * Self::UCS_C1) * jj / (1.0 + Self::UCS_C1 * jj),
            mprime * (h * PI / 180.0).cos(),
            mprime * (h * PI / 180.0).sin(),
        ]
    }

    /// Calculates the CIECAM-UCS ΔE′ (prime) color difference between two colors.
    ///
    /// This method converts each color to its CAM16-UCS Ja′b′ coordinates and then
    /// computes the Euclidean distance, which closely matches perceived color differences.
    ///
    /// # Formula
    /// ```text
    /// ΔE′ = sqrt((J1′ - J2′)² + (a1′ - a2′)² + (b1′ - b2′)²)
    /// ```
    ///
    /// # Parameters
    /// - `other`: The `CieCam16` instance to compare against.
    ///
    /// # Returns
    /// The ΔE′ value as an `f64`, representing the perceptual color difference.
    ///
    /// # Errors
    /// Returns an error if the observers of the two colors do not match.
    fn de_cam(&self, other: &Self) -> Result<f64, crate::Error> {
        if self.observer() != other.observer() {
            return Err(crate::Error::RequireSameObserver);
        }

        let jabp1 = self.jab_prime();
        let jabp2 = other.jab_prime();
        Ok(Self::delta_e_prime_from_jabp(
            jabp1.as_ref(),
            jabp2.as_ref(),
        ))
    }

    /// Calculates the CIECAM16 distance between two Ja'b' values.
    /// This is the CIECAM16 distance formula, which is used to calculate the
    /// perceptual difference between two colors in the CIECAM16 color space.
    fn delta_e_prime_from_jabp(jabp1: &[f64], jabp2: &[f64]) -> f64 {
        if (Self::UCS_KL - 1.0).abs() < f64::EPSILON {
            super::math::distance(jabp1, jabp2)
        } else {
            todo!()
        }
    }
}

#[inline]
fn achromatic_rsp(rgb: Vector3<f64>, nbb: f64) -> f64 {
    (2.0 * rgb[0] + rgb[1] + rgb[2] / 20.0 - 0.305) * nbb
}

pub fn eccentricity(hue_angle: f64) -> f64 {
    0.25 * ((hue_angle.to_radians() + 2.0).cos() + 3.8)
}

pub fn achromatic_response_from_lightness(aw: f64, c: f64, z: f64, lightness: f64) -> f64 {
    aw * (lightness / 100.0).powf(1.0f64 / (c * z))
}

fn inv_cone_adaptation(f_l: f64, x: f64) -> f64 {
    let x = x - 0.1;
    let t = 27.13 * x.abs() / (400.0 - x.abs());
    x.signum() * ((100.0 * t.powf(1.0 / 0.42)) / f_l)
}

const MCAT02: SMatrix<f64, 3, 3> = matrix![
     0.7328,  0.4296,  -0.1624;
    -0.7036,  1.6975,   0.0061;
     0.0030,  0.0136,   0.9834;
];

/// Inverse CIECAT02 Chromatic Adaptation as a Matrix
const MCAT02INV: SMatrix<f64, 3, 3> = matrix![
    1.096123820835514, 		-0.2788690002182872, 	0.18274517938277304;
    0.45436904197535916,	 0.4735331543074117,	0.0720978037172291;
    -0.009627608738429353, 	-0.005698031216113419,	1.0153256399545427;
];

const MHPE: SMatrix<f64, 3, 3> = matrix![
     0.38971, 0.68898, -0.07868;
    -0.22981, 1.18340,  0.04641;
     0.00000, 0.00000,  1.00000;
];


const MHPEINV: SMatrix<f64, 3, 3> = matrix![
    1.9101968340520348, -1.1121238927878747,  0.20190795676749937;
    0.3709500882486886,  0.6290542573926132, -0.000008055142184359149;
    0.0,  				 0.0,  				  1.0;
];

const M16: SMatrix<f64, 3, 3> = matrix![
    0.401288, 0.650173, -0.051461;
    -0.250268, 1.204414, 0.045854;
    -0.002079, 0.048952, 0.953127
];

const M16INV: SMatrix<f64, 3, 3> = matrix![
    1.86206786, -1.01125463, 0.14918677;
    0.38752654, 0.62144744, -0.00897398;
    -0.01584150, -0.03412294, 1.04996444
];
