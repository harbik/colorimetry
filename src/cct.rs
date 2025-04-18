/*!
### Calculation of Correlated Color Temperature and Tint

**Correlated Color Temperature (CCT)** and **Tint** are two parameters used to describe the color
*appearance of lamps. They are based on how incandescent lamps emit light through **thermal
*emission**.

Any object above 0 Kelvin emits electromagnetic radiation, mostly in the infrared range—beyond the
range of human vision. When the object gets hot enough, like a metal heated by a blacksmith, the
emitted light becomes visible. It starts as dark red, shifts to orange as it gets hotter, and
eventually appears yellow-white.

The spectral distribution of thermal radiation from an ideal black, non-reflective object is
described by **Planck’s law**. A perfect blackbody has a **spectral reflectivity** of zero and a
**spectral emissivity** of one across the entire spectrum. Real materials (often referred to as gray
bodies) have more complex behaviors: their reflectivity varies with temperature and depends on the
material’s properties, so they can’t be described by a simple model.

The **chromaticity** of a blackbody emitter depends only on its temperature and traces a curve in a
chromaticity diagram called the **Planckian Locus**. Light sources that are not perfect blackbodies
(like modern lamps or gray emitters) typically fall near, but not exactly on, this curve. Their
color can still be compared to that of a blackbody at a certain temperature—this comparison defines
their **Correlated Color Temperature (CCT)**. The deviation from the Planckian Locus is known as
**Tint**.

Originally, the CIE defined CCT as the temperature of the blackbody with the shortest distance to
the measured color point in the **CIE 1960 (u,v)** color space. Although the newer **CIE 1976
(u',v')** space is now preferred, the original metric is still used—by applying a scaling factor of
2/3 to the v' coordinate. This avoids breaking compatibility with existing industry standards and
ensures continuity, even though it shifts the definition from perception-based to a purely
mathematical one.

CCT calculations are based on the **CIE 1931 Standard Observer**.

# References

*/

use core::f64;
use std::{cmp::max, sync::OnceLock};

use approx::{assert_ulps_eq, relative_eq, ulps_eq, AbsDiffEq, RelativeEq, UlpsEq};

use crate::{geometry::distance_to_line, physics::planck, error::CmtError, observer::{Observer, ObserverData}, data::observers::CIE1931, spectrum::NS, xyz::XYZ};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CCT(f64, f64);

impl AbsDiffEq for CCT {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        // Mired temperature?
        self.0.abs_diff_eq(&other.0, epsilon) &&
        self.1.abs_diff_eq(&other.1, epsilon)
    }
}

impl RelativeEq for CCT {
    fn default_max_relative() -> Self::Epsilon {
        f64::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: Self::Epsilon, max_relative: Self::Epsilon)
        -> bool {
        self.0.relative_eq(&other.0, epsilon, max_relative) &&
        self.1.relative_eq(&other.1, epsilon, max_relative)
    }
}

impl UlpsEq for CCT {
    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.0.ulps_eq(&other.0, epsilon, max_ulps ) &&
        self.1.ulps_eq(&other.1, epsilon, max_ulps )
    }
}

impl CCT {
    /// Create from a Correlated Color Temperature and a Planckian Locus Distance.
    /// Both parameters are range restricted: 1000<cct<1_000_000, -0.05<duv<0.05
    pub fn try_new(cct: f64, duv: f64) -> Result<Self, CmtError> {
        match (cct,duv) {
            (_, d) if d>0.05 => Err(CmtError::CCTDuvHighError),
            (_, d) if d< -0.05 => Err(CmtError::CCTDuvLowError),
            (t, d) if ulps_eq!(t,im2t(0)) || ulps_eq!(t, im2t(N_STEPS-1)) => Ok(Self(t,d)),
            (t, _) if t>im2t(0) => Err(CmtError::CCTTemperatureTooHigh),
            (t, _) if t<im2t(N_STEPS-1) => Err(CmtError::CCTTemperatureTooLow),
            (t, d) => Ok(Self(t,d)),
        }
    }

    /// Create from a Correlated Color Temperature and a Tint value.
    pub fn try_new_with_tint(cct: f64, tint: f64) -> Result<Self, CmtError> {
        Self::try_new(cct, tint/1000.0)
    }

    pub fn t(&self) -> f64 {
        self.0
    }
    pub fn d(&self) -> f64 {
        self.1
    }

    pub fn tint(&self) -> f64 {
        1000.0 * self.1
    }

}

/// Get cct and duv values as an array.
impl From<CCT> for [f64;2] {
    fn from(cct: CCT) -> Self {
        [cct.0, cct.1]
    }
}

/// Number of iterations in binary search.
pub const N_DEPTH: usize = 12; // 2^12 iso temperature line

 
/// Number of iso-temperature lines used, over a range from 1 to N_MIRED_MAX, in this case library
/// corresponding to a temperature range from 1000 to 1_000_000 Kelvin.
/// 
/// An iso-temperature line is here defined by a point on the Planckian locus in the CIE 1960 UCS (u,v) space,
/// and a slope.
pub const N_STEPS:usize = 2 << (N_DEPTH -1);

/// Mired is a unit to express color temperature: M = 1_000_000K / T, with T an absolute correlated
/// color temperature in units of Kelvin.
/// A maximum mired temperature value of 1000 corresponds to a minimum temperature of 1000K.
/// Its minimum value here is 1, corresponding to a maximum correlated color temperature of 1_000_000K.
const MIRED_MAX: usize = 1000;

/// Calculate Correlated Color temperature and Planckian Deviation distance, using the CIE1960 color
/// space, using Robertson and Binary Search algorithms.
impl TryFrom<XYZ> for CCT {
    type Error = CmtError;

    fn try_from(xyz: XYZ) -> Result<Self, Self::Error> {
        if xyz.observer != Observer::Std1931 { return Err(CmtError::RequiresCIE1931XYZ); }
        let [u, v] = xyz.uv60();
        // index bounderies N_STEPS-1 length lookup table e.g. 0-4095
        let [mut imlow, mut imhigh]  = [0usize, N_STEPS-1];
        let [mut dlow, mut dhigh] = [0.0, 0.0];
       
        for _ in 0..N_DEPTH {
            let im = (imhigh+imlow)/2;
            let &[ub, vb, m] = robertson_table(im);
            let d = distance_to_line(u, v, ub, vb, m);
            if d<0.0 { // line is located left of (xyz)
                imlow = im;
                dlow = d;
            } else { // line is located right of (xyz)
                imhigh = im;
                dhigh = d;
            }
        }  
        if imlow == 0 { // xyz is in the first interval, or above high temperature limit
            let &[ub, vb, m] = robertson_table(imlow);
            let d = distance_to_line(u, v, ub, vb, m);
            match distance_to_line(u, v, ub, vb, m){
                d if ulps_eq!(d,0.0,epsilon=1E-10) => { // at low temp limit
                    dlow = 0.0;
                    dhigh = 1.0;
                    imhigh = 1;
                }
                d if d>0.0 => return Err(CmtError::CCTTemperatureTooHigh),
                d => dlow = d 
            };
        };
        if imhigh == N_STEPS - 1 { // xyz is in the last interval, or below low temperature limit
            let &[ub, vb, m] = robertson_table(imhigh);
            match distance_to_line(u, v, ub, vb, m){
                d if ulps_eq!(d,0.0,epsilon=1E-10) => { // at high temp limit
                    dhigh = -0.0;
                    imlow = N_STEPS - 2;
                }
                d if d<0.0 => return Err(CmtError::CCTTemperatureTooLow),
                d => dlow = d 
            };
        };

        let t = robertson_interpolate(im2t(imlow), dlow, im2t(imhigh), dhigh);
        let d = duv_interpolate(u, v, imlow, imhigh, dlow, dhigh);
        CCT::try_new(t, d)
    }

}
    
/// Calculate tristimulus values from a Correlated Color Temperature and a Planckian Locus Distance.
/// This can fail for lower temperatures, and positive distances, with a chromaticity outsiode the CIE 1931 gamut.
impl TryFrom<CCT> for XYZ {
    type Error = CmtError;
    fn try_from(cct: CCT) -> Result<Self, Self::Error> {
        let CCT(t, d) = cct;
        let [u0,v0,m] = iso_temp_line(t);
        let du = m.signum() * d/(m*m+1.0).sqrt();
        let dv = m * du;
        XYZ::try_from_luv60(u0+du, v0+dv, None, None)
    }
    
}

/// Calculates Robertson's Table values for a temperature value of t, in units of Kelvin.
/// These are the coordinates of the blackbody locus at temperature T, and it's line normal,
/// which is the slope of the curve at that point rotated by 90º.
fn iso_temp_line(t: f64) -> [f64;3] {
    let xyz = CIE1931.xyz_planckian_locus(t);
    let [x, y, z] = xyz.values();
    let [u, v ]= xyz.uv60();
    let [dx, dy, dz] = CIE1931.xyz_planckian_locus_slope(t).values();
    let sigma = x + 15.0 * y + 3.0 * z;
    let dsigma = dx + 15.0 * dy + 3.0 * dz;
    let den = 6.0 * y * dsigma - 6.0 * dy * sigma;
    let m = if ulps_eq!(den, 0.0) {
        f64::MAX
    } else {
        (4.0 * dx * sigma - 4.0 * x * dsigma)/den
    };
    [u,v,m]
}

/// Table row of Robertson's Iso Correlated Color Temperature lines, with 4096
/// `(u,v)`` (CIE1960) chromaticity coordinates, and Plankian locus slopes `m`.
///
/// These are used for calculating correlated color temperatures from
/// chromaticity coordinates, as implemente in `XYZ`'s cct method.
/// Index 0 corresponds to a color temperature of 1000K, and index N_ROBERTSON-1 (4095 in this case) to a
/// temperature of 1_000_00K (see function `im2t``).
/// This table is empty on start-up, and rows get filled each time a table
/// entry is requested.
///
/// For more information, see `The Improved Robertson Method for Calculating
/// Correlated Color Temperature` by Gerard Harbers.
pub fn robertson_table(im: usize) -> &'static [f64;3] {
    static ROBERTSON_TABLE: OnceLock<[OnceLock<[f64;3]>;N_STEPS]> = OnceLock::new();

    // Get reference to table, or initialize it when not done yet.
    const UVM_EMPTY: OnceLock<[f64;3]> = OnceLock::new();
    let robertson_table = ROBERTSON_TABLE.get_or_init(|| {
        [UVM_EMPTY; N_STEPS]
    });
    
    // Get table row, or calculate when not done yet.
    robertson_table[im].get_or_init(||{
            let cct = im2t(im);
            iso_temp_line(cct)
    })
}

fn im2t(im: usize) -> f64 {
    1E6/( 1.0 + ((((MIRED_MAX-1)*im)) as f64) / ((N_STEPS-1) as f64)) 
}  

#[test]
fn im2t_test(){
    let i0 = im2t(0);
    assert_ulps_eq!(i0, 1E6);
    let ins = im2t(N_STEPS-1);
    assert_ulps_eq!(ins, 1000.0);
}


fn robertson_interpolate(tp: f64, dp: f64, tn: f64, dn: f64) -> f64 {
    (tp.recip() + dp/(dp-dn) * (tn.recip() - tp.recip())).recip()
}

/// Calculate distance to the Blackbody locus
/// 
/// Finds intersection between iso-temperature lines at (ux,vx), and calculates the distance to this point for the
/// table points. Linear interpolates the value between thse two, and substracts this value of the distance from 
/// the test point to the intersection.
fn duv_interpolate(u: f64, v: f64, imp: usize, imn: usize, dp: f64, dh: f64) -> f64 {
    let [ul, vl, ml] = robertson_table(imp);
    let [uh, vh, mh] = robertson_table(imn);

    // (ux, vx) intersection point of the two iso temperature lines
    let ux =(ml * ul - vl - mh * uh + vh)/(ml - mh);
    let vx = ml * (ux - ul) + vl;

    // distances blackbody points to intersection
    let ddl = (ul - ux).hypot(vl - vx);
    let ddh = (uh - ux).hypot(vh - vx);
    let dd = (u - ux).hypot(v - vx);

    // interpolated distance at blackbody locus
    let d = ddl + dp/(dp-dh) * (ddh - ddl);
    dd - d
}
    
#[test]
fn test_ends(){
    // Temperature  at the low end, should pass...
    let cct0 = CCT::try_new(1000.0, 0.0).unwrap();
    let xyz: XYZ = cct0.try_into().unwrap();
    let cct: CCT = xyz.try_into().unwrap();
    assert_ulps_eq!(cct, cct0);

    // Temperature  at the low end, should pass...
    let cct0 = CCT::try_new(1E6, 0.0).unwrap();
    let xyz: XYZ = cct0.try_into().unwrap();
    let cct: CCT = xyz.try_into().unwrap();
    assert_ulps_eq!(cct, cct0);

    // Temperature  at the low end, should pass...
    let cct0 = CCT::try_new(1005.0, 0.0).unwrap();
    let xyz: XYZ = cct0.try_into().unwrap();
    let cct: CCT = xyz.try_into().unwrap();
    assert_ulps_eq!(cct, cct0, epsilon = 1E-5);

    // Temperature  at the low end, should pass...
    let cct0 = CCT::try_new(1E6 - 100.0, 0.0).unwrap();
    let xyz: XYZ = cct0.try_into().unwrap();
    let cct: CCT = xyz.try_into().unwrap();
    approx::assert_abs_diff_eq!(cct, cct0, epsilon=0.2);
}

#[test]
fn test_cct(){
    // Temperature too low here...
    let xyz: XYZ = CCT(999.0, 0.0).try_into().unwrap();
    let cct: Result<CCT,_> = xyz.try_into();
    assert_eq!(cct, Err(CmtError::CCTTemperatureTooLow));

    // ... and too high here.
    let xyz: XYZ = CCT(1E6+1.0, 0.0).try_into().unwrap();
    let cct: Result<CCT,_> = xyz.try_into();
    assert_eq!(cct, Err(CmtError::CCTTemperatureTooHigh));

    // DUV too high...
    let xyz: XYZ = CCT(6000.0, 0.051).try_into().unwrap();
    let cct: Result<CCT,_> = xyz.try_into();
    assert_eq!(cct, Err(CmtError::CCTDuvHighError));

    // and too low here.
    let xyz: XYZ = CCT(6000.0, -0.051).try_into().unwrap();
    let cct: Result<CCT,_> = xyz.try_into();
    assert_eq!(cct, Err(CmtError::CCTDuvLowError));


    // Test round trip random values, and check the difference by distance in uv-prime space.  For a
    // 4096 size table, distances are found to be less than 6E-5 over the full range of temperatures
    // (1000K, 1_000_000K) and duv values (-0.05, 0.05). This is a relatively slow test, as it tends
    // to fill the Robertson lookup table fully, with each entry requiring to calculate tristimulus
    // values from a Planckian spectrum. It will speed up when more than ~ 5_000 values are tested, 
    // as the table will be completely calculated.

    let seed: u64 = 27; // using a fixed seed, generating the same random numbers at every run.
    let mut rng = <rand::rngs::StdRng as rand::SeedableRng>::seed_from_u64(seed);

    for i in 0..100 {
      //  let mut rng = rand::thread_rng();
        let mired = rand::Rng::gen_range(&mut rng, 1.0..1000.0); // mired temp
        let t = 1E6/mired;
        let d = rand::Rng::gen_range(&mut rng, -0.05..0.05);
        
        // Get a CCT to test. Unwrap OK as range restricted above.
        let cct0  = CCT::try_new(t,d).unwrap();

        // calculate XYZ0, skip value if not a valid point
        let Ok(xyz0):Result<XYZ,_> = CCT::try_new(t,d).unwrap().try_into() else { continue;};
        
        // calculate the cct to test.
        let cct: CCT = xyz0.try_into().unwrap();

        // calculate XYZ, should not fail, a CCT and Duv very close to already fetted values.
        let xyz: XYZ = cct.try_into().unwrap();
        let d = xyz.uv_prime_distance(&xyz0);
        assert_ulps_eq!(d, 0.0, epsilon=6E-5);
    }

}

#[test]
#[cfg(feature = "cie-illuminants")]
fn f1_test() {
    let xyz_f1 = CIE1931.xyz_cie_table(&crate::std_illuminants::StdIlluminant::F1, None);
    // value from CIE Standard CIE15:2004 Table T8.1
    approx::assert_ulps_eq!(xyz_f1.cct().unwrap().t(), 6430.0, epsilon = 0.5);
}

#[test]
#[cfg(feature = "cie-illuminants")]
fn f3_1_test() {
    let xyz_f3_1 = CIE1931.xyz_cie_table(&crate::std_illuminants::StdIlluminant::F3_1, None);
    // value from CIE Standard CIE15:2004 Table T8.1    
    approx::assert_ulps_eq!(xyz_f3_1.cct().unwrap().t(), 2932.0, epsilon = 0.5);
}

