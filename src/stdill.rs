/**!
The CIE Standard Illuminants, available in the library, defined as enums.

The illuminants D65 and D50 are always included in this library, all the others will be only
included with the cie-illumiants feature flag.
A static reference to the spectra can be obtained using the "spectrum" method.

 */


use crate::Spectrum;

// This macro 
macro_rules! std_illuminants {
    ($($val:ident)* [$($cieval:ident)* ]) => {
        #[allow(non_camel_case_types)]
        #[derive(Clone, Copy)]
        pub enum StdIlluminant  {
                $($val,)*
                $(
                    #[cfg(feature="cie-illuminants")]
                    $cieval,
                )*
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
   F3_1 F3_2 F3_3 F3_4 F3_5 F3_6 F3_7 F3_8 F3_9 F3_10 F3_11 F3_12 F3_13  F3_14  F3_15
   LED_B1 LED_B2 LED_B3 LED_B4 LED_B5 LED_BH1 LED_RGB1 LED_V1 LED_V2
   ]);