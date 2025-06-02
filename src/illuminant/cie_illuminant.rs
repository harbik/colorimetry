//! # CIE Standard Illuminants Enumerator

use crate::{illuminant::Illuminant, spectrum::Spectrum, traits::Light};
use std::borrow::Cow;

macro_rules! std_illuminants {
    ($($val:ident)* [$($cieval:ident)* ]) => {
        // only basic selection
        #[cfg(not(feature="cie-illuminants"))]
        #[allow(non_camel_case_types)]
        #[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
        #[derive(Clone, Copy, Debug, strum_macros::Display, strum_macros::EnumIter)]
        pub enum CieIlluminant  {
                $($val,)*
        }

        // extended selection
        #[cfg(feature="cie-illuminants")]
        /// A **lightweight enum** representing the CIE standard illuminants from the CIE 15:2018 datasets
        /// (downloaded August 2024). Each variant holds a zero-cost reference to its precompiled spectrum,
        /// making it easy to include as a field in your own types without pulling in heavy data structures.
        ///
        /// This enum implements `IntoEnumIterator`, so you can **iterate through every standard illuminant**
        /// (useful for testing, batch conversions, or validation).
        ///
        /// - Use `CieIlluminant::iter()` or `CieIlluminant::spectrum()` to list or retrieve any built-in illuminant.
        /// - For a generic D-series illuminant at any correlated color temperature, use
        ///   `Spectrum::cie_d_illuminant(cct: f64)`.
        ///
        /// By default, only **D65** and **D50** are included. To pull in the full set of fluorescent “F3_X”
        /// series and other CIE illuminants, enable the `"cie-illuminants"` feature in `Cargo.toml`
        /// (or build with `--features cie-illuminants`). Omit that feature (or use `--no-default-features`)
        /// to keep your binary lean.
        ///
        /// In JavaScript/WebAssembly builds, the `colorimetry` package excludes these extra spectra by default
        /// for faster load times. To include them, use the `colorimetry-all` bundle instead.
        ///
        /// For more background, see the Wikipedia article on
        /// [Standard illuminant white points](https://en.wikipedia.org/wiki/Standard_illuminant#White_points_of_standard_illuminants).
        ///
        /// # Examples
        /// ```rust
        /// use colorimetry::prelude::*;
        /// use strum::IntoEnumIterator;
        ///
        /// // Iterate through and print all available CIE illuminants:
        /// for illum in CieIlluminant::iter() {
        ///     println!("{illum}");
        /// }
        /// ```
        #[allow(non_camel_case_types)]
        #[cfg_attr(target_arch = "wasm32", wasm_bindgen::prelude::wasm_bindgen)]
        #[derive(Clone, Debug, Copy, strum_macros::Display, strum_macros::EnumIter)]
        pub enum CieIlluminant  {
                $($val,)*
                $($cieval,)*
        }

        impl CieIlluminant {
            pub fn illuminant(&self) -> &crate::illuminant::Illuminant {
                match self {
                    $(Self::$val => &crate::illuminant::cie_data::$val,)*
                    $(
                        #[cfg(feature="cie-illuminants")]
                        Self::$cieval => &crate::illuminant::cie_data::$cieval,
                    )*
                }
            }
        }
    };
}

impl From<CieIlluminant> for Illuminant {
    fn from(std_illuminant: CieIlluminant) -> Self {
        std_illuminant.illuminant().clone()
    }
}

impl AsRef<Illuminant> for CieIlluminant {
    fn as_ref(&self) -> &Illuminant {
        self.illuminant()
    }
}

impl Light for CieIlluminant {
    fn xyzn(
        &self,
        observer: crate::observer::Observer,
        illuminance: Option<f64>,
    ) -> crate::xyz::XYZ {
        observer.data().xyz_cie_table(self, illuminance)
    }

    fn spectrum(&self) -> Cow<'_, Spectrum> {
        // Cow::Borrowed(self.illuminant())
        (&self.illuminant().0).into()
    }
}

std_illuminants!(D65 D50 [A F1 F2 F3 F4 F5 F6 F7 F8 F9 F10 F11 F12
F3_1 F3_2 F3_3 F3_4 F3_5 F3_6 F3_7 F3_8 F3_9 F3_10 F3_11 F3_12 F3_13 F3_14 F3_15
LED_B1 LED_B2 LED_B3 LED_B4 LED_B5 LED_BH1 LED_RGB1 LED_V1 LED_V2
]);
