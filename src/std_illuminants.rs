/*!
# Standard Illuminants

Many of the CIE Standard Illuminants are made available in this module through
`StdIlluminant` for data-defined illuminants and the
`Spectrum::cie_d_illuminant(cct: f64)` function for generic D-illuminant.

The `StdIlluminant` object gives access to all the CIE illuminants defined in this library, obtained from the datasets published in the
CIE15::2018 standard, downloaded from the [CIE Website](https://web.archive.org/web/20240314231650/https://cie.co.at/data-tables) August 2024.

As the data is compiled into the library, you can choose only to include the two basic illuminants `D65` and `D50` to limit the size of your
recent `F3_X` series included here,
executable by using the `--no-default-features` in the compiler, or the `default-features = false` in its dependence declaration of this crate in
the `cargo.toml` file of your application. The library uses the "cie-illuminants" feature flag to select the inclusion of these illuminants.

If you use this library in a JavaScript application, the default `colorimetry` package excludes all the features that allow fast load times for
lightweight applications.  Use the `colorimetry-all` package to use the library with all its features enabled.

For more detailed information on the CIE Standard Illuminant Datasets, see
[Standard illuminant](https://en.wikipedia.org/wiki/Standard_illuminant#White_points_of_standard_illuminants)
on Wikipedia.  Instead of a dash, use the `_` character to access these
illuminants by their name here, so use `StdIlluminant::LED_BH1` to use the
phosphor-converted Blue LED and Red LED standard illuminant.
The Fluorescent `F3_X` series is included here, with X ranging from 1 to 15.
*/

use std::{ops::Deref, vec};
use nalgebra::{ArrayStorage, SMatrix};
use wasm_bindgen::prelude::*;
use crate::{CmtError, Spectrum, RefWhite, Illuminant, NS};

/**
The CIE Standard Illuminants, available in the library, defined as enums.

The illuminants D65 and D50 are always included in this library, all the others will be only
included with the cie-illumiants feature flag.
A static reference to the spectra can be obtained using the "spectrum" method.

```
// print all the StdIlluminants
    use colorimetry::StdIlluminant;
    use strum::IntoEnumIterator;

    for spc in StdIlluminant::iter() {
        println!{"{spc}"};
    }
```

 */



// This macro generates the `StdIlluminant` enumerator, representing the standard illuminants
// available in the library.  It adds the illuminants defined as static data by their name as an
// identifier to an `enum`, and add a `spectrum` method to access its data.  This is somewhat
// contrived, as it is a work around for a wasm-bindgen bug, which does not obey the feature-cfg on
// a enum item.  see [issue](https://github.com/rustwasm/wasm-bindgen/issues/3297).
macro_rules! std_illuminants {
    ($($val:ident)* [$($cieval:ident)* ]) => {
        // only basic selection
        #[cfg(not(feature="cie-illuminants"))]
        #[allow(non_camel_case_types)]
        #[wasm_bindgen]
        #[derive(Clone, Copy, strum_macros::Display, strum_macros::EnumIter)]
        pub enum StdIlluminant  {
                $($val,)*
        }

        // extended selection
        #[cfg(feature="cie-illuminants")]
        #[allow(non_camel_case_types)]
        #[wasm_bindgen]
        #[derive(Clone, Copy, strum_macros::Display, strum_macros::EnumIter)]
        pub enum StdIlluminant  {
                $($val,)*
                $($cieval,)*
        }

        impl StdIlluminant {
            pub fn illuminant(&self) -> &crate::Illuminant {
                match self {
                    $(Self::$val => &crate::$val,)*
                    $(
                        #[cfg(feature="cie-illuminants")]
                        Self::$cieval => &crate::$cieval,
                    )*
                }
            }
        }
    };
}

impl From<StdIlluminant> for Illuminant {
    fn from(std_illuminant: StdIlluminant) -> Self {
        std_illuminant.illuminant().clone()
    }
}

impl RefWhite for StdIlluminant {
    fn xyzn(&self, observer: crate::Observer, illuminance: Option<f64>) -> crate::XYZ {
        observer.data().xyz_cie_table(self, illuminance)
    }
    
    fn spectrum(&self) -> &Illuminant {
        self.illuminant()
    }
}

std_illuminants!(D65 D50 [A F1 F2 F3 F4 F5 F6 F7 F8 F9 F10 F11 F12
   F3_1 F3_2 F3_3 F3_4 F3_5 F3_6 F3_7 F3_8 F3_9 F3_10 F3_11 F3_12 F3_13 F3_14 F3_15
   LED_B1 LED_B2 LED_B3 LED_B4 LED_B5 LED_BH1 LED_RGB1 LED_V1 LED_V2
   ]);




