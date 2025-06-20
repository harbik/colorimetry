mod cielch_gamut_cie1931_d50;
mod cielch_gamut_cie1931_d65;
mod cielch_gamut_cie1964_d50;
mod cielch_gamut_cie1964_d65;
mod cielch_gamut_cie2015_10_d50;
mod cielch_gamut_cie2015_10_d65;
mod cielch_gamut_cie2015_d50;
mod cielch_gamut_cie2015_d65;

pub enum CieLChGamut {
    CielchGamutCie1931D50,
    CielchGamutCie1931D65,
    CielchGamutCie1964D50,
    CielchGamutCie1964D65,
    CielchGamutCie2015_10D50,
    CielchGamutCie2015_10D65,
    CielchGamutCie2015D50,
    CielchGamutCie2015D65,
}

use cielch_gamut_cie1931_d50::CIELCH_GAMUT_CIE1931_D50;
use cielch_gamut_cie1931_d65::CIELCH_GAMUT_CIE1931_D65;
use cielch_gamut_cie1964_d50::CIELCH_GAMUT_CIE1964_D50;
use cielch_gamut_cie1964_d65::CIELCH_GAMUT_CIE1964_D65;
use cielch_gamut_cie2015_10_d50::CIELCH_GAMUT_CIE2015_10_D50;
use cielch_gamut_cie2015_10_d65::CIELCH_GAMUT_CIE2015_10_D65;
use cielch_gamut_cie2015_d50::CIELCH_GAMUT_CIE2015_D50;
use cielch_gamut_cie2015_d65::CIELCH_GAMUT_CIE2015_D65;

impl CieLChGamut {
    pub fn data(&self) -> &[[u8; 72]; 99] {
        match self {
            CieLChGamut::CielchGamutCie1931D50 => &CIELCH_GAMUT_CIE1931_D50,
            CieLChGamut::CielchGamutCie1931D65 => &CIELCH_GAMUT_CIE1931_D65,
            CieLChGamut::CielchGamutCie1964D50 => &CIELCH_GAMUT_CIE1964_D50,
            CieLChGamut::CielchGamutCie1964D65 => &CIELCH_GAMUT_CIE1964_D65,
            CieLChGamut::CielchGamutCie2015_10D50 => &CIELCH_GAMUT_CIE2015_10_D50,
            CieLChGamut::CielchGamutCie2015_10D65 => &CIELCH_GAMUT_CIE2015_10_D65,
            CieLChGamut::CielchGamutCie2015D50 => &CIELCH_GAMUT_CIE2015_D50,
            CieLChGamut::CielchGamutCie2015D65 => &CIELCH_GAMUT_CIE2015_D65,
        }
    }

    pub fn in_gamut(&self, lab: super::CieLab) -> bool {
        let [l, c, h] = lab.lch();
        todo!()
    }
}
