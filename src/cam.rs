use std::f64::consts::PI;

const P1C: f64 = 50_000.0 / 13.0;
const P3: f64 = 21.0 / 20.0;
const C16_3: f64 = 1_403.0;
const NOM: f64 = 1.0; // is listed in the standard as (2.0+Self::P3)*460.0/C16_3; // this is 1!! But his is in Luo step 3 as part of nominators
const DEN1: f64 = ((2.0 + P3) * 220.0) / C16_3;
const DEN2: f64 = (P3 * 6300.0 - 27.0) / C16_3;
const RCPR_9: f64 = 1.0 / 0.9;
const UCS_KL: f64 = 1.0;
const UCS_C1: f64 = 0.007;
const UCS_C2: f64 = 0.0228;
use nalgebra::{matrix, vector, Matrix3, SMatrix, Vector3};

use crate::{
    error::CmtError,
    geometry::distance,
    prelude::Observer,
    traits::{Filter, Light},
    xyz::XYZ,
};

use super::viewconditions::{ReferenceValues, ViewConditions};

#[derive(Debug)]
pub struct CieCam16 {
    /// Colorimetric Observer used.
    /// The standard requires use of the CIE1931 standard observer.
    pub(crate) observer: Observer,
    /// Correlates of Lightness, Chroma, and hue-angle
    pub(crate) jch: Vector3<f64>,
    /// Tristimulus values of the reference white beng use
    pub(crate) xyzn: Vector3<f64>,
    /// Viewing Conditions
    pub(crate) vc: ViewConditions, // 6 f64
}

impl CieCam16 {
    /// CIECAM16 coordinates for a particular set of viewing conditions.
    fn new(xyz: XYZ, xyzn: XYZ, vc: ViewConditions) -> Result<Self, CmtError> {
        let xyz_vec = xyz.xyz;
        let xyzn_vec = xyzn.xyz;
        if xyz.observer != xyzn.observer {
            return Err(CmtError::RequireSameObserver);
        }
        let ReferenceValues {
            n,
            z,
            nbb,
            ncb,
            d_rgb,
            aw,
            qu,
        } = ReferenceValues::new(xyzn_vec, vc);
        let vcdd = vc.dd();
        let vcfl = vc.f_l();
        let mut rgb = M16 * xyz_vec;
        rgb.component_mul_assign(&Vector3::from(d_rgb));
        rgb.apply(|v| vc.lum_adapt(v, 0.26, qu));

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

    fn jabp(&self) -> Vector3<f64> {
        let &[jj, cc, h] = self.jch.as_ref();
        let m = cc * self.vc.f_l().powf(0.25);
        let mprime = 1.0 / UCS_C2 * (1.0 + UCS_C2 * m).ln();
        Vector3::new(
            (1.0 + 100.0 * UCS_C1) * jj / (1.0 + UCS_C1 * jj),
            mprime * (h * PI / 180.0).cos(),
            mprime * (h * PI / 180.0).sin(),
        )
    }

    fn jchp(&self) -> Vector3<f64> {
        let &[jj, cc, h] = self.jch.as_ref();
        let m = cc * self.vc.f_l().powf(0.25);
        let mprime = 1.0 / UCS_C2 * (1.0 + UCS_C2 * m).ln();
        Vector3::new(
            (1.0 + 100.0 * UCS_C1) * jj / (1.0 + UCS_C1 * jj),
            mprime * (h * PI / 180.0).cos(),
            mprime * (h * PI / 180.0).sin(),
        )
    }

    // Static equation Observer::Cam::delta_e_prime, same for CAM16 and CAM02 values
    fn delta_e_prime(jabp1: &[f64], jabp2: &[f64]) -> f64 {
        if (UCS_KL - 1.0).abs() < f64::EPSILON {
            distance(jabp1, jabp2)
        } else {
            todo!()
        }
    }

    /// Inverse Transform, using optional different view conditions or a different reference white.
    ///
    /// This can be used to match colors with different viewing conditions, for example to match the
    /// colors as viewied in a viewing cabinet to the colors on a display.
    /// Also a different white adaptation state can be given, as a chromatic adaptation transform is
    /// included in the CIECAM model.
    /// If no different white adaptation or viewing conditionns are given, this is a straight back
    /// transform from the input parameters, which can for example be used to test the backward
    /// transform.
    fn xyz(&self, white_opt: Option<XYZ>, vc_opt: Option<ViewConditions>) -> Result<XYZ, CmtError> {
        let vc = vc_opt.unwrap_or(self.vc);
        let xyzn = if let Some(white) = white_opt {
            if white.observer == self.observer {
                white.xyz
            } else {
                return Err(CmtError::RequireSameObserver);
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
            qu,
        } = ReferenceValues::new(xyzn, vc);
        let d_rgb_vec = Vector3::from(d_rgb);
        let &[lightness, chroma, hue_angle] = self.jch.as_ref();
        let t = (chroma / ((lightness / 100.0).sqrt() * (1.64 - 0.29f64.powf(n)).powf(0.73)))
            .powf(RCPR_9);
        let p1 = (P1C * vc.nc * ncb * eccentricity(hue_angle)) / t; // NaN if t=0, but OK, as check on t==0.0 if used
        let p2 = achromatic_response_from_lightness(aw, vc.c, z, lightness) / nbb + 0.305;
        let (a, b) = match hue_angle.to_radians().sin_cos() {
            (_, _) if t.is_nan() || t == 0.0 => (0.0, 0.0),
            (hs, hc) if hs.abs() >= hc.abs() => {
                let b = p2 * NOM / (p1 / hs + DEN1 * hc / hs + DEN2);
                (b * hc / hs, b)
            }
            (hs, hc) => {
                let a = p2 * NOM / (p1 / hc + DEN1 + DEN2 * hs / hc);
                (a, a * hs / hc)
            }
        };

        // rgb_a
        let m = matrix![ 460.0, 451.0, 288.0; 460.0, -891.0, -261.0; 460.0, -220.0, -6_300.0; ]
            / 1_403.0;
        let rgb_p = (m * vector![p2, a, b]).map(|x| inv_cone_adaptation(vc.f_l(), x)); // Step 4 & 5
        let rgb = rgb_p.component_div(&d_rgb_vec);

        let xyz = M16INV * rgb;
        Ok(XYZ::from_vecs(xyz, self.observer))
    }
}

#[inline]
pub fn achromatic_rsp(rgb: Vector3<f64>, nbb: f64) -> f64 {
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

pub const MCAT02: SMatrix<f64, 3, 3> = matrix![
     0.7328,  0.4296,  -0.1624;
    -0.7036,  1.6975,   0.0061;
     0.0030,  0.0136,   0.9834;
];

/**
   Inverse CIECAT02 Chromatic Adaptation as a Matrix
*/
pub const MCAT02INV: SMatrix<f64, 3, 3> = matrix![
    1.096123820835514, 		-0.2788690002182872, 	0.18274517938277304;
    0.45436904197535916,	 0.4735331543074117,	0.0720978037172291;
    -0.009627608738429353, 	-0.005698031216113419,	1.0153256399545427;
];

pub const MHPE: SMatrix<f64, 3, 3> = matrix![
     0.38971, 0.68898, -0.07868;
    -0.22981, 1.18340,  0.04641;
     0.00000, 0.00000,  1.00000;
];

pub const MHPEINVLUO: SMatrix<f64, 3, 3> = matrix![
    1.910197, -1.112124,  0.201908;
    0.370950,  0.629054, -0.000008;
    0.000000,  0.000000,  1.000000;
];

pub const MHPEINV: SMatrix<f64, 3, 3> = matrix![
    1.9101968340520348, -1.1121238927878747,  0.20190795676749937;
    0.3709500882486886,  0.6290542573926132, -0.000008055142184359149;
    0.0,  				 0.0,  				  1.0;
];

pub const MCAT02INVLUO: SMatrix<f64, 3, 3> = matrix![
     1.096124, -0.278869, 0.182745;
     0.454369,  0.473533, 0.072098;
    -0.009628, -0.005698, 1.015326;
];

pub static M16: Matrix3<f64> = matrix![
    0.401288, 0.650173, -0.051461;
    -0.250268, 1.204414, 0.045854;
    -0.002079, 0.048952, 0.953127
];

pub static M16INV: Matrix3<f64> = matrix![
    1.86206786, -1.01125463, 0.14918677;
    0.38752654, 0.62144744, -0.00897398;
    -0.01584150, -0.03412294, 1.04996444
];

pub static MRGBAINV: Matrix3<f64> = matrix![
460.0/C16_3, 451.0/C16_3, 288.0/C16_3;
460.0/C16_3, -891.0/C16_3, -261.0/C16_3;
460.0/C16_3, -220.0/C16_3, -6_300.0/C16_3;
 ];

#[cfg(test)]
mod cam_test {
    use super::*;
    use crate::prelude::*;
    use approx::{assert_abs_diff_eq, assert_ulps_eq};

    #[test]
    fn test_m16() {
        approx::assert_abs_diff_eq!(M16INV * M16, Matrix3::identity(), epsilon = 1E-8);
    }

    #[test]
    fn test_worked_example() {
        // see section 7 CIE 248:2022
        let xyz = XYZ::new([60.70, 49.60, 10.29], Observer::Std1931);
        let xyzn = XYZ::new([96.46, 100.0, 108.62], Observer::Std1931);
        let vc = ViewConditions::new(16.0, 1.0, 1.0, 0.69, 40.0, None);
        let cam = CieCam16::new(xyz, xyzn, vc).unwrap();
        let &[j, c, h] = cam.jch.as_ref();
        // println!("J:\t{j:?}\nC:\t{c:?}\nh:\t{h:?}");
        approx::assert_abs_diff_eq!(j, 70.4406, epsilon = 1E-4);
        approx::assert_abs_diff_eq!(c, 58.6035, epsilon = 1E-4);
        approx::assert_abs_diff_eq!(h, 57.9145, epsilon = 1E-4);

        // inverse transformation, with no change in white adaptation of viewing conditions.
        let xyz_rev = cam.xyz(None, None).unwrap();
        assert_abs_diff_eq!(xyz, xyz_rev, epsilon = 1E-4);
    }
}
