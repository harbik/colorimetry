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


use std::vec;
use nalgebra::{ArrayStorage, SMatrix};
use wasm_bindgen::prelude::*;

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


use crate::{Category::Illuminant, CmtError, Spectrum, NS};

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
            pub fn spectrum(&self) -> &Spectrum {
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

std_illuminants!(D65 D50 [A F1 F2 F3 F4 F5 F6 F7 F8 F9 F10 F11 F12
   F3_1 F3_2 F3_3 F3_4 F3_5 F3_6 F3_7 F3_8 F3_9 F3_10 F3_11 F3_12 F3_13 F3_14 F3_15
   LED_B1 LED_B2 LED_B3 LED_B4 LED_B5 LED_BH1 LED_RGB1 LED_V1 LED_V2
   ]);


impl Spectrum {
    pub fn d_illuminant(cct: f64) -> Result<Spectrum, CmtError> {
        if cct<4000.0 || cct>25000.0 {
            Err(CmtError::OutOfRange{name:"CIE D Illuminant Temperature".to_string(), low: 4000.0, high: 25000.0})
        } else { 
            let xd = match cct {
                t if t < 7000.0 => {
                    0.244063 + 0.09911E3 / t + 2.9678E6 / t.powi(2) - 4.607E9 / t.powi(3)
                }
                t => 0.23704 + 0.24748E3 / t + 1.9018E6 / t.powi(2) - 2.0064E9 / t.powi(3),
            };
            let yd = -3. * xd.powi(2) + 2.87 * xd - 0.275;
            let m = 0.0241 + 0.2562 * xd - 0.7341 * yd;
            let m1 = (-1.3515 - 1.7703 * xd + 5.9114 * yd) / m;
            let m2 = (0.03 - 31.4424 * xd + 30.0717 * yd) / m;
            let mut v = [0.0; CIE_D_S_LEN];
            v.iter_mut().enumerate().for_each(|(i,x)|*x = CIE_D_S[(i,0)] + m1 * CIE_D_S[(i,1)] + m2 * CIE_D_S[(i,2)]);
            Spectrum::linear_interpolate(Illuminant, &[380.0, 780.0], &v)
        }
    }
}


#[test]
fn test_d_illuminant(){

    let s = Spectrum::d_illuminant(6504.0).unwrap();
    let xyz = crate::CIE1931.xyz(&s).set_illuminance(100.0);
    approx::assert_ulps_eq!(xyz, crate::CIE1931.xyz_d65(), epsilon = 2E-2);
}

const CIE_D_S_LEN: usize = 81;

static CIE_D_S: SMatrix::<f64, CIE_D_S_LEN, 3> = SMatrix::from_array_storage(ArrayStorage([
    [63.40, 64.60, 65.80, 80.30, 94.80, 99.80, 104.80, 105.35, 105.90, 101.35, 96.80, 105.35, 113.90, 119.75, 125.60, 125.55, 125.50, 123.40, 121.30, 121.30,
    121.30, 117.40, 113.50, 113.30, 113.10, 111.95, 110.80, 108.65, 106.50, 107.65, 108.80, 107.05, 105.30, 104.85, 104.40, 102.20, 100.00, 98.00, 96.00, 95.55,
    95.10, 92.10, 89.10, 89.80, 90.50, 90.40, 90.30, 89.35, 88.40, 86.20, 84.00, 84.55, 85.10, 83.50, 81.90, 82.25, 82.60, 83.75, 84.90, 83.10, 81.30, 76.60,
    71.90, 73.10, 74.30, 75.35, 76.40, 69.85, 63.30, 67.50, 71.70, 74.35, 77.00, 71.10, 65.20, 56.45, 47.70, 58.15, 68.60, 66.80, 65.00],
    [38.50, 36.75, 35.00, 39.20, 43.40, 44.85, 46.30, 45.10, 43.90, 40.50, 37.10, 36.90, 36.70, 36.30, 35.90, 34.25, 32.60, 30.25, 27.90, 26.10, 24.30, 22.20,
    20.10, 18.15, 16.20, 14.70, 13.20, 10.90, 8.60, 7.35, 6.10, 5.15, 4.20, 3.05, 1.90, 0.95, 0.00, -0.80, -1.60, -2.55, -3.50, -3.50, -3.50, -4.65, -5.80,
    -6.50, -7.20, -7.90, -8.60, -9.05, -9.50, -10.20, -10.90, -10.80, -10.70, -11.35, -12.00, -13.00, -14.00, -13.80, -13.60, -12.80, -12.00, -12.65, -13.30,
    -13.10, -12.90, -11.75, -10.60, -11.10, -11.60, -11.90, -12.20, -11.20, -10.20, -9.00, -7.80, -9.50, -11.20, -10.80, -10.40],
    [3.00, 2.10, 1.20, 0.05, -1.10, -0.80, -0.50, -0.60, -0.70, -0.95, -1.20, -1.90, -2.60, -2.75, -2.90, -2.85, -2.80, -2.70, -2.60, -2.60, -2.60, -2.20,
    -1.80, -1.65, -1.50, -1.40, -1.30, -1.25, -1.20, -1.10, -1.00, -0.75, -0.50, -0.40, -0.30, -0.15, 0.00, 0.10, 0.20, 0.35, 0.50, 1.30, 2.10, 2.65, 3.20,
    3.65, 4.10, 4.40, 4.70, 4.90, 5.10, 5.90, 6.70, 7.00, 7.30, 7.95, 8.60, 9.20, 9.80, 10.00, 10.20, 9.25, 8.30, 8.95, 9.60, 9.05, 8.50, 7.75, 7.00, 7.30,
    7.60, 7.80, 8.00, 7.35, 6.70, 5.95, 5.20, 6.30, 7.40, 7.10, 6.80]
]));

