use core::f64;
use std::sync::OnceLock;
use wasm_bindgen::prelude::wasm_bindgen;
use nalgebra::{Matrix3, SMatrix};
use crate::{lab::Lab, physics::planck, RgbSpaceTag, spc::{Category, Spectrum, NS}, xyz::XYZ, CmError, LineAB, StdIlluminant};



/**
    Light-weight identifier added to the `XYZ` and `RGB` datasets,
    representing the colorimetric standard observer used.
 */
#[wasm_bindgen]
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub enum ObserverTag { 
    #[default] Std1931, 
    Std1976, 
    Std2015, 
    Std2015_10
}

impl ObserverTag {
    /**
        Get the `Observer` object, associated with the tag.
     */
    pub fn observer(&self) -> &'static Observer {
        match self {
            ObserverTag::Std1931 =>  &crate::data::CIE1931,
            _ => todo!()
        }
    }
}


/**
    A data structure to define Standard Observers, such as the CIE 1931 2º and
    the CIE 2015 standard observers.
    
    These are defined in the form of the three color matching functions,
    typically denoted by $\hat{x}(\lamda)$,$\hat{y}{\lambda}$, and $\hat{z}(\lambda)$.
    Traditionally, the CIE1931 Colorimetric Standard Observer is used almost exclusively,
    but is known to be not a good representation of human vision in the blue region of the
    spectrum. We also know now that the way you see color varies with age, and your healty,
    and that not everyone sees to same color.

    In this library colors are represented by spectral distributions, to allow color modelling
    with newer, and better standard observers, such as the CIE2015 Observer, derived from
    the sensitivities of the cones in the retina of your eye, the biological color receptors
    of light.

    It's main purpose is to calculate `XYZ` tristimulus values for a general stimulus,
    in from of a `Spectrum`.
*/
#[wasm_bindgen]
pub struct Observer {
    pub(crate) data: SMatrix<f64, 3, NS>,
    pub(crate) lumconst: f64,
    pub(crate) id: ObserverTag,
}

impl Observer {
    /**
        Calculates Tristimulus valus of a single standard spectrum.
        For example used with illuminants, or stimuli.
        Tristimulus values are the basis for all calculations in CIE Colorimetry.
    */
    pub fn xyz(&self, s: &Spectrum) -> XYZ {
        let t = self.data * s.data;
        XYZ {
            data:  t * self.lumconst,
            obs_id: self.id,
            yw: None
        }
    }

    /**
        Tristimulus Values for the Standard Illuminants in this library,
        normalized to a luminous value of y = yw = 100.0.

        Also provided is the absolute luminous value of the illuminant, which
        typically has little value from a colorimetric perspective, but is used
        for scaling.

        Values are calculated on first use.
    */
    pub fn xyz_std_illuminant(&self, std_illuminant: &StdIlluminant) -> &(XYZ, f64) {
        const EMPTY:OnceLock<(XYZ, f64)> = OnceLock::new();
        const XYZ_STD_ILLUMINANTS_LEN: usize = 32;
        static XYZ_STD_ILLUMINANTS : OnceLock<[OnceLock<(XYZ, f64)>;XYZ_STD_ILLUMINANTS_LEN]> = OnceLock::new();
        let xyz_std_illuminants = XYZ_STD_ILLUMINANTS.get_or_init(||[EMPTY; XYZ_STD_ILLUMINANTS_LEN]);
        xyz_std_illuminants[*std_illuminant as usize].get_or_init(||{
            let mut xyz = self.xyz(std_illuminant.spectrum());
            xyz.yw = Some(100.0);
            let yabs = xyz.data.y;
            xyz.data.iter_mut().for_each(|v|*v = *v * 100.0/yabs);
            (xyz , yabs)
        })
    }

    /// Calculates XYZ tristimulus values for a Planckian emitter
    /// for an observer. 
    pub fn xyz_planckian_locus(&self, cct: f64) -> XYZ {
        let mut l = 380.0E-9; // wavelength in unit of meter
        let step = 1E-9; // wavelength interval
        let [mut x, mut y, mut z]: [f64;3]  = Default::default();
        let mut pow = 0.0;
        for i in 0..NS {
            let p = planck(l, cct);
            pow += p;
            x +=  p * self.data[(i,0)];
            y +=  p * self.data[(i,1)];
            z +=  p * self.data[(i,2)];
            l += step;
        }
        let scale = self.lumconst/pow; // Scale for 1 W/m2 irradiance.
        XYZ::new(x * scale, y * scale, z * scale, None, self.id)
    }

    pub fn planckian_locus_slope(&self, _cct: f64) -> f64 {
        todo!()
    }
    
    /// XYZ tristimulus values for the CIE standard daylight illuminant D65.
    /// The values are calculated on first use.
    pub fn xyz_d65(&self) -> XYZ {
        self.xyz_std_illuminant(&StdIlluminant::D65).0
    }

    /// XYZ tristimulus values for the CIE standard daylight illuminant D50.
    /// The values are calculated on first use.
    pub fn xyz_d50(&self) -> XYZ {
        self.xyz_std_illuminant(&StdIlluminant::D50).0
    }

    /// XYZ tristimulus values for the CIE standard A, or incandescent lamp.
    /// The values are calculated on first use.
    pub fn xyz_a(&self) -> XYZ {
        self.xyz_std_illuminant(&StdIlluminant::A).0
    }

    /// Calulates Tristimulus values for a sample illuminated with a `StandardIlluminant`.
    /// The values are normalized for a white illuminance of 100 cd/m2.
    pub fn xyz_with_std_illuminant(&self, illuminant: &StdIlluminant, sample: &Spectrum) -> XYZ {
        let std_illuminant = illuminant.spectrum();
        let t = self.data * std_illuminant.data.component_mul(&sample.data);
        let (_, yabs) = self.xyz_std_illuminant(illuminant);
        XYZ {
            data:  t * self.lumconst * 100.0 / *yabs,
            obs_id: self.id,
            yw: Some(100.0),
        }
    }


    /// Calculates the L*a*b* CIELAB D65 values of a ColorPatch or Filter, using D65 as an illuminant.
    /// Accepts a Filter or ColorPatch Spectrum only.
    /// Returns f64::NAN's otherwise.
    pub fn lab_d65(&self, sample: &Spectrum) -> Lab {
        if sample.cat != Category::Filter && sample.cat != Category::ColorPatch { // invalid
            Lab::new(f64::NAN, f64::NAN, f64::NAN, self.xyz_d65())
        } else {
            let &[x, y, z] = self.xyz_with_std_illuminant(&StdIlluminant::D65, sample).data.as_ref();
            Lab::new(x, y, z, self.xyz_d65())
        }
    }

    /// Calculate the Spectral Locus, or (x,y) coordinates of the _horse shoe_,
    /// the boundary of area of all physical colors in a chromiticity diagram,
    /// as XYZ tristimulus values.  Tristimulus values are returned here,
    /// instead of chromaticity xy coordinates, to make use of the available XYZ
    /// transforms.
    /// 
    /// This function limits the input values to produce unique chromaticity
    /// values only.  Spectral locus points tend to freeze, or even fold back to
    /// lower wavelength values at the blue and red perimeter ends.  This can be
    /// quite anoying, for example when trying to calculate dominant wavelength,
    /// or when creating plots.  To get the allowed range, use the
    /// `spectral_wavelength_min` and `spectral_wavelength_max` methods.
    /// 
    /// See Wikipedia's [CIE 1931 Color Space](https://en.wikipedia.org/wiki/CIE_1931_color_space).
    pub fn spectral_locus_by_nm(&self, l: usize) -> Result<XYZ, CmError> {
        let min = self.spectral_locus_nm_min();
        let max = self.spectral_locus_nm_max();
        if l<380 || l>780 {
            return Err(CmError::WavelengthOutOfRange);
        };
        if l<min || l>max {
            Err(CmError::NoUniqueSpectralLocus(min, max))
        } else {
            let &[x, y, z] = self.data.column(l-380).as_ref();
            Ok(XYZ::new(x, y, z, None, self.id))

        }
    }

    /// Unrestricted, direct, access to the spectal locus data.
    /// To get unique values only please use the `spectral_locus_by_nm` function.
    pub fn spectral_locus_by_index(&self, i:usize) -> [f64;2] {
        let &[x, y, z] = self.data.column(i).as_ref();
        let s = x + y + z;
        [x/s, y/s]
    }

    /// The index value of the blue spectral locus edge.
    /// Any further spectral locus points will hover around this edge, and will not have a unique wavelength.
    pub fn spectral_locus_index_min(&self) -> usize {
        static MIN: OnceLock<usize> = OnceLock::new();
        *MIN.get_or_init(||{
            const START: usize = 100;
            let mut lp = LineAB::try_new(self.spectral_locus_by_index(START), [0.33333, 0.33333]).unwrap();
            let mut m = START - 1;
            loop {
                let l = LineAB::try_new(self.spectral_locus_by_index(m), [0.33333, 0.33333]).unwrap();
                match (m, l.angle_diff(lp)) {
                    (0, d) if d> -f64::EPSILON => break m + 1,
                    (0, _)  => break 0,
                    (1.., d) if d> -f64::EPSILON => break m,
                    _ => {
                        m = m - 1;
                        lp = l;
                    }
                }
            }
        })
    }

    pub fn spectral_locus_nm_min(&self) -> usize {
        self.spectral_locus_index_min()+380
    }

    /// The index value of the red spectral locus edge.
    /// Any further spectral locus points will hover around this edge.
    pub fn spectral_locus_index_max(&self) -> usize {
        static MAX: OnceLock<usize> = OnceLock::new();
        *MAX.get_or_init(||{
            const START: usize = 300;
            let mut lp = LineAB::try_new(self.spectral_locus_by_index(START), [0.33333, 0.33333]).unwrap();
            let mut m = START + 1;
            loop {
                let l = LineAB::try_new(self.spectral_locus_by_index(m), [0.33333, 0.33333]).unwrap();
                match (m, l.angle_diff(lp)) {
                    (400, d) if d< f64::EPSILON => break m-1,
                    (400, _)  => break 400,
                    (..400, d) if d< f64::EPSILON => break m-1,
                    _ => {
                        m = m + 1;
                        lp = l;
                    }
                }
            }
        })
    }

    pub fn spectral_locus_nm_max(&self) -> usize {
        self.spectral_locus_index_max()+380
    }


 
    /// Calculates the RGB to XYZ matrix, for a particular color space.
    /// The matrices are buffered.
    pub fn rgb2xyz(&self, rgbspace: &RgbSpaceTag) -> &'static Matrix3<f64> {
        const EMPTY:OnceLock<Matrix3<f64>> = OnceLock::new();
        const RGB2XYZ_AR_LEN: usize = 16;
        static RGB2XYZ_AR : OnceLock<[OnceLock<Matrix3<f64>>;RGB2XYZ_AR_LEN]> = OnceLock::new();
        let rgb2xyz_ar =RGB2XYZ_AR.get_or_init(||[EMPTY;RGB2XYZ_AR_LEN]);
        rgb2xyz_ar[*rgbspace as usize].get_or_init(||{
            let (space,_) = rgbspace.rgb_space();
            let mut rgb2xyz = Matrix3::from_iterator(space.primaries.iter().flat_map(|s|self.xyz(s).set_illuminance(1.0).values()));
            let xyzw = self.xyz(&space.white).set_illuminance(1.0);
            let decomp = rgb2xyz.lu();
            // unwrap: only used with library color spaces
            let rgbw = decomp.solve(&xyzw.data).unwrap();
            for (i, mut col) in rgb2xyz.column_iter_mut().enumerate() {
                col *= rgbw[i];
            }
            rgb2xyz

        })
    }

    /// Calculates the RGB to XYZ matrix, for a particular color space.
    /// The matrices are buffered.
    pub fn xyz2rgb(&self, rgbspace: RgbSpaceTag) -> &'static Matrix3<f64> {
        const EMPTY:OnceLock<Matrix3<f64>> = OnceLock::new();
        const XYZ2RGB_AR_LEN: usize = 16;
        static XYZ2RGB_AR : OnceLock<[OnceLock<Matrix3<f64>>;XYZ2RGB_AR_LEN]> = OnceLock::new();
        let xyz2rgb =XYZ2RGB_AR.get_or_init(||[EMPTY;XYZ2RGB_AR_LEN]);
        xyz2rgb[rgbspace as usize].get_or_init(||{
            // unwrap: only used with library color spaces
            self.rgb2xyz(&rgbspace).try_inverse().unwrap()
        })
    }


}




#[cfg(test)]
mod obs_test {

    use crate::CIE1931;
    use approx::assert_ulps_eq;

    #[test]
    fn test_spectral_locus(){
        let [x,y] = CIE1931.spectral_locus_by_nm(CIE1931.spectral_locus_nm_min()).unwrap().chromaticity();
        assert_ulps_eq!(x, 0.17411, epsilon=1E-5);
        assert_ulps_eq!(y, 0.00496, epsilon=1E-5);

        let [x,y] = CIE1931.spectral_locus_by_nm(CIE1931.spectral_locus_nm_max()).unwrap().chromaticity();
        assert_ulps_eq!(x, 0.73469, epsilon=1E-5);
        assert_ulps_eq!(y, 0.26531, epsilon=1E-5);
    }

    #[test]
    fn test_spectral_locus_min_max(){
        let min = CIE1931.spectral_locus_index_min();
        println!("{min}");
        let max = CIE1931.spectral_locus_index_max();
        println!("{max}");
    }

    #[test]
    // Data from http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
    // Lindbloom's values are reproduced with an accuracy of 3E-4, which is
    // a small, but significant difference.  This difference is caused by a difference in the display's white point,
    // due to wavelength domain differences.  Here we use a domain
    // from 380 to 780 with a step of 1 nanometer for the spectra, and in
    // specific for the color matching functions, as recommended by
    // CIE15:2004, and the color matching functions provided by the CIE in
    // their online dataset section. The spectra for the primaries are chosen to match the RGB primary values as given by the Color Space specifications.
    // The white point uses the standard's spectral distribution, as provided by the CIE, clipped to a domain from 380 to 780 nanometer.
    // See `colorimetry::data::D65`.
    fn test_rgb2xyz_cie1931(){
        let want =  nalgebra::Matrix3::new(0.4124564,  0.3575761,  0.1804375, 0.2126729,  0.7151522,  0.0721750, 0.0193339,  0.1191920,  0.9503041);
        let got = CIE1931.rgb2xyz(&crate::RgbSpaceTag::SRGB);
        approx::assert_ulps_eq!(want, got, epsilon = 3E-4);
    }

    #[test]
    // Check the inverse transformation.
    // See comments at `test_rgb2xyz_cie1931`.
    fn test_xyz2rgb_cie1931(){
        let want =  nalgebra::Matrix3::new(3.2404542, -1.5371385, -0.4985314, -0.9692660,  1.8760108,  0.0415560, 0.0556434, -0.2040259,  1.0572252);
        let got = CIE1931.xyz2rgb(crate::RgbSpaceTag::SRGB);
        approx::assert_ulps_eq!(want, got, epsilon = 3E-4);
    }

    #[test]
    fn test_xyz_std_illuminants(){
        let (xyz, yabs) = CIE1931.xyz_std_illuminant(&crate::StdIlluminant::D65);
        println!("{xyz:?} {yabs:.4}");
    }
    
}