use std::collections::HashMap;

use crate::illuminant::CieIlluminant;
use crate::observer::Observer;
use std::sync::LazyLock;

macro_rules! generate_relxyz_gamuts {
    ($(($identifier:ident, $observer:ident, $illuminant:ident)),* $(,)?) => {
        $(
            pub(super) static $identifier: LazyLock<RelXYZGamutData> = LazyLock::new(|| {
                let opt_colors = Observer::$observer.optimal_colors(CieIlluminant::$illuminant);
                let hashmap = opt_colors.max_luminance_per_chromaticity_bin();
                RelXYZGamutData::new(
                    Observer::$observer,
                    hashmap,
                )
            });
        )*
    };
}

generate_relxyz_gamuts!(
    (RELXYZ_GAMUT_CIE1931_D50_DATA, Cie1931, D50),
    (RELXYZ_GAMUT_CIE1931_D65_DATA, Cie1931, D65),
);

#[cfg(feature = "supplemental-observers")]
generate_relxyz_gamuts!(
    (CIEXYZ_GAMUT_CIE1964_D50_DATA, Cie1964, D50),
    (CIEXYZ_GAMUT_CIE1964_D65_DATA, Cie1964, D65),
    (CIEXYZ_GAMUT_CIE2015_D50_DATA, Cie2015, D50),
    (CIEXYZ_GAMUT_CIE2015_D65_DATA, Cie2015, D65),
    (CIEXYZ_GAMUT_CIE2015_10_D50_DATA, Cie2015_10, D50),
    (CIEXYZ_GAMUT_CIE2015_10_D65_DATA, Cie2015_10, D65),
);

pub struct RelXYZGamutData {
    observer: Observer,
    max_luminances: HashMap<[u16; 2], u16>,
}

impl RelXYZGamutData {
    pub const fn new(observer: Observer, max_luminances: HashMap<[u16; 2], u16>) -> Self {
        RelXYZGamutData {
            observer,
            max_luminances,
        }
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, strum::Display)]
/// Represents the relative XYZ gamut for different observers and reference white points.
/// This enum provides a way to access the gamut data for various CIE observers and illuminants.
pub enum RelXYZGamut {
    #[default]
    Cie1931D65,
    Cie1931D50,
    #[cfg(feature = "supplemental-observers")]
    Cie1964D50,
    #[cfg(feature = "supplemental-observers")]
    Cie1964D65,
    #[cfg(feature = "supplemental-observers")]
    Cie2015D50,
    #[cfg(feature = "supplemental-observers")]
    Cie2015D65,
    #[cfg(feature = "supplemental-observers")]
    Cie2015_10D50,
    #[cfg(feature = "supplemental-observers")]
    Cie2015_10D65,
}

impl RelXYZGamut {
    pub fn observer(&self) -> Observer {
        self.data().observer
    }

    pub fn data(&self) -> &RelXYZGamutData {
        match self {
            RelXYZGamut::Cie1931D50 => &*RELXYZ_GAMUT_CIE1931_D50_DATA,
            RelXYZGamut::Cie1931D65 => &*RELXYZ_GAMUT_CIE1931_D65_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie1964D50 => &*CIEXYZ_GAMUT_CIE1964_D50_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie1964D65 => &*CIEXYZ_GAMUT_CIE1964_D65_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie2015D50 => &*CIEXYZ_GAMUT_CIE2015_D50_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie2015D65 => &*CIEXYZ_GAMUT_CIE2015_D65_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie2015_10D50 => &*CIEXYZ_GAMUT_CIE2015_10_D50_DATA,
            #[cfg(feature = "supplemental-observers")]
            RelXYZGamut::Cie2015_10D65 => &*CIEXYZ_GAMUT_CIE2015_10_D65_DATA,
        }
    }

    pub fn max_luminance_for_bin(&self, x: u16, y: u16) -> Option<u16> {
        self.data().max_luminances.get(&[x, y]).copied()
    }
}
