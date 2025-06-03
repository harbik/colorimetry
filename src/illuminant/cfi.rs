use crate::math::distance;
use crate::observer::CIE1964;
use crate::{
    cam::{CamTransforms, CieCam02, TM30VC},
    colorant::{CFI_DATA, N_CFI},
    illuminant::Illuminant,
};

use super::CCT;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
#[allow(dead_code)]
pub struct CFI {
    jabp_ts: [[f64; 3]; N_CFI], // CIE1964
    jabp_rs: [[f64; 3]; N_CFI],
    cct: CCT, // using CIE1931
}

impl CFI {
    pub fn new(illuminant: &Illuminant) -> Result<Self, crate::Error> {
        let cct = illuminant.cct()?;
        let ref_illuminant = Illuminant::cfi_reference(cct.t())?;
        let vc = TM30VC;
        let xyzn_r = CIE1964.xyz(&ref_illuminant, None).set_illuminance(100.0);
        let xyzn_t = CIE1964.xyz(illuminant, None).set_illuminance(100.0);
        let mut jabp_ts = [[0f64; 3]; N_CFI];
        let mut jabp_rs = [[0f64; 3]; N_CFI];
        for (i, cfi_ces) in CFI_DATA.iter().enumerate() {
            let xyz_t = CIE1964.xyz(illuminant, Some(cfi_ces));
            let jabp_t = CieCam02::from_xyz(xyz_t, xyzn_t, vc)?.jab_prime();
            jabp_ts[i] = jabp_t;

            let xyz_r = CIE1964.xyz(&ref_illuminant, Some(cfi_ces));
            let jabp_r = CieCam02::from_xyz(xyz_r, xyzn_r, vc)?.jab_prime();
            jabp_rs[i] = jabp_r;
        }

        Ok(CFI {
            jabp_ts,
            jabp_rs,
            cct,
        })
    }

    pub fn special_color_fidelity_indices(&self) -> [f64; N_CFI] {
        let mut cfis = [0f64; N_CFI];
        for (i, (jabd, jabr)) in self.jabp_ts.iter().zip(self.jabp_rs.iter()).enumerate() {
            let de = distance(jabd, jabr);
            cfis[i] = rf_from_de(de);
        }
        cfis
    }

    pub fn general_color_fidelity_index(&self) -> f64 {
        let mut sum = 0.0;
        for (jabp_t, jabp_r) in self.jabp_ts.iter().zip(self.jabp_rs.iter()) {
            let de = distance(jabp_t, jabp_r);
            sum += de;
        }
        rf_from_de(sum / N_CFI as f64)
    }

    pub fn cct(&self) -> CCT {
        self.cct
    }
}

const CF: f64 = 6.73; // was 7.54 in TM30-15
pub fn rf_from_de(de: f64) -> f64 {
    10.0 * (((100.0 - CF * de) / 10.0).exp() + 1.0).ln()
}

#[cfg(test)]
mod tests {
    use super::CFI;
    use crate::{
        colorant::N_CFI,
        illuminant::{D65, F1, F2, F12},
    };
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_cfi() {
        let cfi = CFI::new(&D65).unwrap();
        let cfis = cfi.special_color_fidelity_indices();
        assert_eq!(cfis.len(), N_CFI);
        for &cfi in cfis.iter() {
            assert_abs_diff_eq!(cfi, 100.0, epsilon = 0.01); // Rf for CIE D65
        }
        assert!(cfi.general_color_fidelity_index() >= 0.0);
    }

    #[test]
    fn test_fl1() {
        let cfi = CFI::new(&F1).unwrap();
        assert_abs_diff_eq!(cfi.cct().t(), 6428.0, epsilon = 3.0); // Rf for DE=1
        assert_abs_diff_eq!(cfi.general_color_fidelity_index(), 81.0, epsilon = 0.5);
        // Rf for DE=1
    }

    #[test]
    fn test_fl2() {
        let cfi = CFI::new(&F2).unwrap();
        assert_abs_diff_eq!(cfi.cct().t(), 4225.0, epsilon = 1.0); // Rf for DE=1
        assert_abs_diff_eq!(cfi.general_color_fidelity_index(), 70.0, epsilon = 0.5);
        // Rf for DE=1
    }

    #[test]
    fn test_fl12() {
        let cfi = CFI::new(&F12).unwrap();
        assert_abs_diff_eq!(cfi.cct().t(), 3000.0, epsilon = 1.0); // Rf for DE=1
        assert_abs_diff_eq!(cfi.general_color_fidelity_index(), 78.0, epsilon = 0.6);
        // Rf for DE=1
    }
}
