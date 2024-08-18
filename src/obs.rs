use std::sync::OnceLock;
use wasm_bindgen::prelude::wasm_bindgen;
use nalgebra::SMatrix;
use crate::{lab::Lab, physics::planck, rgb::RGB, spc::{Spectrum, Category, NS}, xyz::XYZ};



#[wasm_bindgen]
#[derive(Clone, Copy, Default)]
pub enum ObsId { 
    #[default] Std1931, 
    Std1976, 
    Std2015, 
    Std2015_10
}



/// A data structure to define Standard Observers, such as the CIE 1931 2º and the CIE 2015 standard observers.
/// These are defined in the form of three discrete representations of color matching functions.
#[wasm_bindgen]
pub struct Observer {
    pub(crate) data: SMatrix<f64, 3, NS>,
    pub(crate) lumconst: f64,
    pub(crate) id: ObsId,
}


impl Observer {
    /// Calculates Tristimulus valus of a single standard spectrum.
    /// For example used with illuminants, or stimuli.
    /// Tristimulus values are the basis for all calculations in CIE Colorimetry.
    pub fn xyz(&self, s: &Spectrum) -> XYZ {
        let t = self.data * s.data;
        XYZ {
           data:  t * self.lumconst,
           obs_id: self.id
        }
    }

    /// Calulates Tristimulus values from a multiplicative combination of two standard spectra.
    /// There are no checks on spectral categories.
    /// Most commonly, one will be an illuminant, the other a Filter or ColorPatch.
    pub fn xyz2(&self, s1: &Spectrum, s2: &Spectrum) -> XYZ {

        let t = self.data * s1.data.component_mul(&s2.data);
        XYZ {
           data:  t * self.lumconst,
           obs_id: self.id
        }
    }


    /// Calculates the L*a*b* CIELAB D65 values of a ColorPatch or Filter, using D65 as an illuminant.
    /// Accepts a Filter or ColorPatch Spectrum only.
    /// Returns f64::NAN's otherwise.
    pub fn lab_d65(&self, s: &Spectrum) -> Lab {
        if s.cat != Category::Filter && s.cat != Category::ColorPatch { // invalid
            Lab::new(f64::NAN, f64::NAN, f64::NAN, self.d65())
        } else {
            let &[x, y, z] = self.xyz2(&crate::data::D65,s).data.as_ref();
            Lab::new(x, y, z, self.d65())
        }
    }

    /// Calculates XYZ tristimulus values for a Planckian emitter.
    pub fn xyz_planck(&self, cct: f64) -> XYZ {
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
        XYZ::new(x * scale, y * scale, z * scale, self.id)
    }

    pub fn planckian_slope(&self, _cct: f64) -> f64 {
        todo!()
    }

    pub fn d65(&self) -> XYZ {
        static D65XYZ: OnceLock<XYZ> = OnceLock::new();
        D65XYZ.get_or_init(||{
            self.xyz(&crate::data::D65)
        }).clone()
    }

    pub fn d50(&self) -> XYZ {
        static D50XYZ: OnceLock<XYZ> = OnceLock::new();
        D50XYZ.get_or_init(||{
            self.xyz(&crate::data::D50)
        }).clone()
    }

    pub fn s_rgb(&self) -> RGB {
        todo!()
    }

}

impl Observer {
    /// Table row of Robertson's Iso Correlated Color Temperature lines, with 4096
    /// `(u,v)`` (CIE1960) chromaticity coordinates, and Plankian locus slopes `m`.
    ///
    /// These are used for calculating correlated color temperatures from
    /// chromaticity coordinates, as implemente in `XYZ`'s cct method.
    /// Index 0 corresponds to a color temperature of 1000K, and index 4096 to a
    /// temperature of 1_000_00K (see function `im2t``).
    /// This table is empty on start-up, and rows get filled each time a table
    /// entry is requested.
    ///
    /// For more information, see `The Improved Robertson Method for Calculating
    /// Correlated Color Temperature` by Gerard Harbers.
    pub fn robertson_4k_table(&self, im: usize) -> [f64;3] {
        static ROBERTSON_TABLE: OnceLock<[OnceLock<[f64;3]>;4096]> = OnceLock::new();

        // Get reference to table, or initialize it when not done yet.
        const EMPTY_ROW: OnceLock<[f64;3]> = OnceLock::new();
        let robertson_table = ROBERTSON_TABLE.get_or_init(|| {
            [EMPTY_ROW; 4096]
        });
        
        // Get table row, or calculate when not done yet.
        let _uvm = robertson_table[im].get_or_init(||{
                let cct = im2t(im);
                let [u,v] = self.xyz_planck(cct).uv60();
                let m = self.planckian_slope(cct);
                [u,v,m]
            }
        );
        todo!()
    }

}

const MIRED_MAX: usize = 1000;
const NP:usize = 4096;
fn im2t(im: usize) -> f64 {
    1E6/( 1.0 + ((((MIRED_MAX-1)*im)) as f64) / ((NP-1) as f64)) 
}
