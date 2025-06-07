/*!
This is a Rust library for working with color and light — great for projects in lighting, imaging, or anything that needs accurate color handling.
It includes tools to simulate how colors look under different lights, convert between color spaces, and follow well-known standards from groups like the CIE, ICC, and IES.
You can use it to build lighting tools, visualize spectra, or get the right transforms for your color profiles.

It also has early support for JavaScript and WebAssembly, so you can run it in the browser, and use it with JavaScript Runtimes such as Deno.

# Usage
To use this library in a Rust application, run the command:
 ```bash
    cargo add colorimetry
```
or add this line to the dependencies in your Cargo.toml file:
```text
    colorimetry = "0.0.6"
```

# Examples
<details>
<summary><i>Calculate Tristimulus Values for Illuminants</i></summary>

This example calculates the XYZ tristimulus values of the D65 illuminant for both the CIE 1931 2º standard observer and the CIE 2015 10º observer.

```
    use colorimetry::prelude::*;
    use colorimetry::illuminant::D65;
    use approx::assert_abs_diff_eq as check;

    // D65 Tristimulus values, using the CIE1931 standard observer by default
    let xyz_d65 = D65.xyz(None).set_illuminance(100.0);

    let [x, y, z] = xyz_d65.values();
    check!([x, y, z].as_ref(), [95.04, 100.0, 108.86].as_ref(),  epsilon = 5E-3);

# #[cfg(feature = "supplemental-observers")]
# {
    // D65 Tristimulus values using the CIE2015 10º observer
    // This requires the `supplemental-observers` feature (enabled by default)
    let xyz_d65_10 = D65
        .xyz(Some(Observer::Cie2015_10))
        .set_illuminance(100.0);

    let [x_10, y_10, z_10] = xyz_d65_10.values();
   check!([x_10, y_10, z_10].as_ref(), [94.72, 100.0, 107.143].as_ref(), epsilon = 5E-3);
# }
```
</details>

<details>
<summary><i>Calculate Correlated Color Temperature and Tint</i></summary>

The correlated color temperature (CCT) of an illuminant, typically expressed in kelvin (K),
describes whether a light source appears warm (low CCT) or cool (high CCT). It is a key parameter
for characterizing the visual appearance of white light .
This example calculates both the correlated color temperature and the deviation from the Planckian
locus, often referred to as the tint.

```
    use colorimetry::prelude::*;
# #[cfg(feature="cie-illuminants")]
    use colorimetry::illuminant::A;
    use approx::assert_abs_diff_eq as check;

    // Calculate CCT and Duv for the A illuminant
    // Requires `cct`, and `cie-illuminants` features
# #[cfg(all(feature="cct", feature="cie-illuminants"))]
    let [cct, duv] = A.cct().unwrap().values();
# #[cfg(all(feature="cct", feature="cie-illuminants"))]
    check!([cct, duv].as_ref(), [2855.4977, 0.0].as_ref(),  epsilon = 5E-4);
```
</details>

<details>
<summary><i>Calculate Color Fidelity Index for Illuminants</i></summary>

The CIE has announced that the Color Fidelity Index (CFI) will replace the Color Rendering Index
(CRI) as the standard metric for evaluating color rendering. Both indices aim to quantify how
accurately a light source reproduces the colors of illuminated objects. However, the CFI offers a
significant improvement in accuracy by using 99 reference color samples and more advanced color
difference metrics, compared to the CRI’s use of only 8 samples.
Below is an example calculation of the general Color Fidelity Index for the CIE F2 illuminant:

```
    use colorimetry::prelude::*;
# #[cfg(feature = "cie-illuminants")]
    use colorimetry::illuminant::F2;
    use approx::assert_abs_diff_eq as check;

# #[cfg(all(feature = "cfi", feature = "cie-illuminants"))]
# {
    // Calculate the Color Fidelity Index of the CIE F2 standard illuminant
    // Requires `cfi`, and `cie-illuminants` features
    let cf_f2 = F2.cfi().unwrap();
    let cf = cf_f2.general_color_fidelity_index();
    check!(cf, 70.3,  epsilon = 1E-1);
# }
```
</details>

<details>

<summary><i>Calculate Spectral Locus to Plot in Diagrams</i></summary>

The spectral locus is the boundary in a chromaticity diagram that encloses all perceivable,
physically realizable colors. Due to its shape, it is sometimes informally referred to as the
"horseshoe."
Below, we compute the chromaticity coordinates that define the spectral locus.

```
    use colorimetry::observer::Observer::Cie1931;
    let mut locus = Vec::new();
    let wavelength_range = Cie1931.spectral_locus_wavelength_range();
    for wavelength in wavelength_range {
        // unwrap OK because nm is in range
        let xyz = Cie1931.xyz_at_wavelength(wavelength).unwrap();
        let chromaticity = xyz.chromaticity();
        locus.push([wavelength as f64, chromaticity.x(), chromaticity.y()]);
    }
    println!("{locus:?}");
```
</details>

<details>
<summary><i>Calculate XYZ/RGB Transformation Matrices, for any Observer, for use in Color Profiles</i></summary>
This is usually done with the CIE 1931 Standard Observer, but this library supports any observer—as long as both the color space and the data use the same one.  
Instead of fixed XYZ values, it computes conversions from the spectral definitions of the primaries to be able to do so.
Here, we compute transformation matrices for the `DisplayP3` color space using both the `Cie1931` and `Cie2015` observers.  

```
    use approx::assert_abs_diff_eq as check;
    use colorimetry::observer::Observer;
    use colorimetry::rgb::RgbSpace::DisplayP3;

    let xyz2rgb_31 = Observer::Cie1931.xyz2rgb(DisplayP3);
    let want31 = nalgebra::Matrix3::new(
         2.4933, -0.9313,	-0.4027,
        -0.8298,  1.7629,	 0.0236,
         0.0355, -0.076,	 0.9574
    );
    check!(xyz2rgb_31, want31, epsilon=5E-4);

    let rgb2xyz_31 = Observer::Cie1931.rgb2xyz(DisplayP3);
    let want31inv = nalgebra::Matrix3::new(
        0.4866, 0.2656, 0.1981,
        0.2291, 0.6917, 0.0792,
        0.0001, 0.0451, 1.0433,

    );
    check!(rgb2xyz_31, want31inv, epsilon=5E-4);
    # #[cfg(feature = "supplemental-observers")]
    # {
    // requires `supplemental-observers` 
    use colorimetry::observer::Observer::Cie2015;

    let xyz2rgb_15 = Cie2015.xyz2rgb(DisplayP3);
    let want15 = nalgebra::Matrix3::new(
        2.5258,  -1.0009,	-0.3649,
        -0.9006,   1.8546,	-0.0011,
        0.0279,  -0.0574,	 0.95874
    );
    check!(xyz2rgb_15, want15, epsilon=5E-4);
    # }
```
</details>

<details>
<summary><i>Find the closest spectral match in the Munsell Color Book</i></summary>

This example finds the best matching color in the Munsell Color Book for a given sample—in this case, the R9 test color used in the CRI color rendering method.
It uses the `CieCam16::de_ucs` color difference metric and the `Cie2015_10` standard observer to calculate perceptual similarity.

The closest match identified is Munsell "5R 5/14", a vivid red hue, with a color difference of just 3 ΔE.  
In practical terms, a ΔE of 3 is considered a close match—just at the threshold where most observers might start to notice a difference under controlled viewing conditions.

```
# #[cfg(all(feature= "cri", feature = "supplemental-observers", feature = "munsell"))]
# {
    use approx::assert_abs_diff_eq as check;
    use colorimetry::observer::Observer::Cie2015_10;
    use colorimetry::colorant::MunsellCollection;
    use colorimetry::colorant::TCS;

    let cri_r9 = &TCS[8];
    let (key, delta_e) = MunsellCollection::match_ciecam16(
        cri_r9,
        None,
        None,
        Some(Cie2015_10),
    )
    .unwrap();
    assert_eq!(key, "5R4/14");
    check!(delta_e, 3.0, epsilon = 5e-2);
# }
```
</details>

# Features

- [`Spectrum`] Standard fixed-grid spectral representation over a wavelength range from 380 to 780 nanometers, with 1-nanometer intervals.  
  - Uses [`nalgebra`] vector and matrix types for fast integration and transformations, with access to a wide range of numerical algorithms.  
  - Supports interpolation from irregular spectral data using [`Spectrum::linear_interpolate`] and [`Spectrum::sprague_interpolate`].  
  - Optional smoothing using a Gaussian filter via [`Spectrum::smooth`].

- Generate spectral distributions from analytical models:  
  - [`Illuminant::planckian`] Planck’s law for blackbody radiators  
  - [`Illuminant::led`] Spectral power distribution of a LED  
  - [`Colorant::gaussian`] Gaussian color filters  
  - [`Stimulus::from_rgb`] Spectral distribution of an RGB color pixel using Gaussian spectral primaries

- CIE Standard Illuminants:  
  - Daylight: [`D65`], [`D50`]  
  - Incandescent: [`A`]  
  - Fluorescent: [`F1`], [`F2`], [`F3`], [`F4`], [`F5`], [`F6`], [`F7`], [`F8`], [`F9`], [`F10`], [`F11`], [`F12`]  
  - Extended fluorescent set: [`F3_1`], [`F3_2`], [`F3_3`], [`F3_4`], [`F3_5`], [`F3_6`], [`F3_7`], [`F3_8`], [`F3_9`], [`F3_10`], [`F3_11`], [`F3_12`], [`F3_13`], [`F3_14`], [`F3_15`]  
  - LED: [`LED_B1`], [`LED_B2`], [`LED_B3`], [`LED_B4`], [`LED_B5`], [`LED_BH1`], [`LED_RGB1`], [`LED_V1`]

- Illuminant metrics:  
  - [`CCT`] Correlated color temperature, including distance to the blackbody locus for tint indication[^1]
  - [`CRI`] Color rendering index[^2]  
  - [`CFI`] Color fidelity index[^3]

- Advanced color (appearance) models:  
  - [`CieLab`], [`CieCam02`], [`CieCam16`]  
  - Color difference methods: [`CieLab::ciede`], [`CieLab::ciede2000`], [`CieCam02::de_ucs`], [`CieCam16::de_ucs`]

- Spectrally based RGB color spaces, with support for non-CIE 1931 observers and generic transformations between [`Rgb`] and [`XYZ`]:  
  - RGB Color Spaces [`RgbSpace::SRGB`],  [`RgbSpace::Adobe`], [`RgbSpace::DisplayP3`]
  - Define RGB colors using spectral primaries, allowing switching observers

- Multiple CIE Standard Observers (see [`Observer`]), with transformations to [`XYZ`]:  
  - [`Observer::Cie1931`] the CIE 1931 2º standard observer  
  - [`Observer::Cie1964`] the CIE 1964 10º standard observer  
  - [`Observer::Cie2015`] the CIE 2015 2º cone fundamentals-based observer  
  - [`Observer::Cie2015_10`] the CIE 2015 10º cone fundamentals-based observer

# Feature Flags

- `supplemental-observers`  
  Adds addiational observers such as `Cie1964`, `Cie2015`, and `Cie2015_10`. Enabled by default.


- `cie-illuminants`  
  The `D65` and `D50` illuminants are always included - if you want to use one of the other CIE illuminants, set this feature flag.
- `munsell`  
  Include reflection spectra for Munsell colors.

- `cct`  
  Calculates correlated color temperatures (CCT) for illuminants.
  Generates a 4096-entry lookup table (each entry containing three `f64` values).
  Memory is reserved at compile time but computed on demand.
  (Included automatically if the `cri` feature is enabled).

- `cri`  
  Enables Color Rendering Index (CRI) calculations, providing Ra and R1–R14 values for illuminants.  
  Loads an additional 14 test color sample spectra.

- `cfi`  
  Enables Color Fidelity Index (CRI) calculations, providing R<sub>f</sub> and R<sub>f,1</sub>–R<sub>f,99</sub> values for illuminants.  
  Loads the 99 CES test color sample spectra.

<details>
<summary>How to enable a feature?</summary>

To enable a feature, such as `cri` and `munsell`, use 

```bash
cargo add colorimetry -F cri,munsell
```
or
```bash
cargo add colorimetry --features cri,munsell
```

Alternatively, configure features manually in your `Cargo.toml`:

```toml
colorimetry = { version = "0.0.6", features = ["cri", "munsell"] }
```

</details>

## License
All content &copy;2025 Harbers Bik LLC, and licensed under either of the

 * Apache License, Version 2.0,
   ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>), or the
 * MIT license
   [LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>,

at your option.

## Contribution

Unless you explicitly state otherwise, any Contribution intentionally submitted
for inclusion in the Work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

[`nalgebra`]:https://docs.rs/nalgebra/latest/nalgebra/ 
[`Spectrum`]: https://docs.rs/colorimetry/latest/colorimetry/spectrum/struct.Spectrum.html
[`Spectrum::linear_interpolate`]: https://docs.rs/colorimetry/latest/colorimetry/spectrum/struct.Spectrum.html#method.linear_interpolate
[`Spectrum::sprague_interpolate`]: https://docs.rs/colorimetry/latest/colorimetry/spectrum/struct.Spectrum.html#method.sprague_interpolate
[`Spectrum::smooth`]: https://docs.rs/colorimetry/latest/colorimetry/spectrum/struct.Spectrum.html#method.smooth
[`Illuminant`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.Illuminant.html
[`Illuminant::planckian`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.Illuminant.html#method.planckian
[`Illuminant::led`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.Illuminant.html#method.led
[`CieLab`]: https://docs.rs/colorimetry/latest/colorimetry/lab/struct.CieLab.html
[`CieLab::ciede`]: https://docs.rs/colorimetry/latest/colorimetry/lab/struct.CieLab.html#method.ciede
[`CieLab::ciede2000`]: https://docs.rs/colorimetry/latest/colorimetry/lab/struct.CieLab.html#method.ciede2000
[`CieCam02`]: https://docs.rs/colorimetry/latest/colorimetry/cam/struct.CieCam02.html
[`CieCam02::de_ucs`]: https://docs.rs/colorimetry/latest/colorimetry/cam/struct.CieCam02.html#method.de_ucs
[`CieCam16`]: https://docs.rs/colorimetry/latest/colorimetry/cam/struct.CieCam16.html
[`CieCam16::de_ucs`]: https://docs.rs/colorimetry/latest/colorimetry/cam/struct.CieCam16.html#method.de_ucs
[`CCT`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CCT.html
[`CRI`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CRI.html
[`CFI`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html
[`Colorant`]: https://docs.rs/colorimetry/latest/colorimetry/colorant/struct.Colorant.html
[`Colorant::gaussian`]: https://docs.rs/colorimetry/latest/colorimetry/colorant/struct.Colorant.html#method.gaussian
[`Stimulus`]: https://docs.rs/colorimetry/latest/colorimetry/stimulus/struct.Stimulus.html
[`Stimulus::from_srgb`]: https://docs.rs/colorimetry/latest/colorimetry/stimulus/struct.Stimulus.html#method.from_srgb
[`Stimulus::from_rgb`]: https://docs.rs/colorimetry/latest/colorimetry/stimulus/struct.Stimulus.html#method.from_rgb
[Colorimetric Observers]: https://docs.rs/colorimetry/latest/colorimetry/observer/index.html
[`Observer`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html
[`Observer::Cie1931`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html#variant.Cie1931
[`Observer::Cie1964`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html#variant.Cie1964
[`Observer::Cie2015`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html#variant.Ci2015e
[`Observer::Cie2015_10`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html#variant.Cie2015_10
[`ObserverData`]:https://docs.rs/colorimetry/latest/colorimetry/observer/enum.ObserverData.html 
[`Observer.xyz`]: https://docs.rs/colorimetry/latest/colorimetry/observer/struct.ObserverData.html#method.xyz
[`CIE1931`]: https://docs.rs/colorimetry/latest/colorimetry/data/observers/static.CIE1931.html
[`CIE1964`]: https://docs.rs/colorimetry/latest/colorimetry/data/observers/static.CIE1964.html
[`CIE2015`]: https://docs.rs/colorimetry/latest/colorimetry/data/observers/static.CIE2015.html
[`CIE2015_10`]: https://docs.rs/colorimetry/latest/colorimetry/data/observers/static.CIE2015_10.html
[`XYZ`]: https://docs.rs/colorimetry/latest/colorimetry/xyz/struct.XYZ.html
[`Rgb`]: https://docs.rs/colorimetry/latest/colorimetry/rgb/struct.RGB.html
[`WideRgb`]: https://docs.rs/colorimetry/latest/colorimetry/widergb/struct.WideRgb.html
[`RgbSpace`]: https://docs.rs/colorimetry/latest/colorimetry/rgbspace/enum.RgbSpace.html
[`RgbSpace::SRGB`]: https://docs.rs/colorimetry/latest/colorimetry/rgbspace/enum.RgbSpace.html#variant.SRGB
[`RgbSpace::Adobe`]: https://docs.rs/colorimetry/latest/colorimetry/rgbspace/enum.RgbSpace.html#variant.Adobe
[`RgbSpace::DisplayP3`]: https://docs.rs/colorimetry/latest/colorimetry/rgbspace/enum.RgbSpace.html#variant.DisplayP3
[`RgbSpaceData`]: https://docs.rs/colorimetry/latest/colorimetry/rgbspace/struct.RgbSpaceData.html

[`D50`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.D50.html 
[`D65`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.D65.html 
[`A`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.A.html 
[`F1`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F1.html 
[`F2`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F2.html 
[`F3`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3.html 
[`F4`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F4.html 
[`F5`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F5.html 
[`F6`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F6.html 
[`F7`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F7.html 
[`F8`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F8.html 
[`F9`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F9.html 
[`F10`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F10.html 
[`F11`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F11.html 
[`F12`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F12.html 
[`F3_1`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_1.html 
[`F3_2`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_2.html 
[`F3_3`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_3.html 
[`F3_4`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_4.html 
[`F3_5`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_5.html 
[`F3_6`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_6.html 
[`F3_7`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_7.html 
[`F3_8`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_8.html 
[`F3_9`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_9.html 
[`F3_10`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_10.html 
[`F3_11`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_11.html 
[`F3_12`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_12.html 
[`F3_13`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_13.html 
[`F3_14`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_14.html 
[`F3_15`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.F3_15.html 
[`LED_B1`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.LED_B1.html 
[`LED_B2`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.LED_B2.html 
[`LED_B3`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.LED_B3.html 
[`LED_B4`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.LED_B4.html 
[`LED_B5`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.LED_B5.html 
[`LED_BH1`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.LED_BH1.html 
[`LED_RGB1`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.LED_RGB1.html 
[`LED_V1`]: https://docs.rs/colorimetry/latest/colorimetry/data/illuminants/static.LED_V1.html 

[^1]: Commission Internationale de l'Éclairage. (2004). *CIE 015:2004: Colorimetry* (3rd ed.). Vienna.  
[^2]: Commission Internationale de l'Éclairage, (1995). *CIE 13.3-1995: Method of Measuring and Specifying Colour Rendering Properties of Light Sources*, Vienna.  
[^3]: Commission Internationale de l'Éclairage. (2017). *CIE 224:2017: Colour Fidelity Index for accurate scientific use*. Vienna.  

*/

// This library defines many floating-point constants for colorimetry.
// Clippy's `approx_constant` lint would otherwise generate numerous false positives
// by flagging constants close to standard values.
// To suppress these, `#![allow(clippy::approx_constant)]` is applied.
#![allow(clippy::approx_constant)]


pub mod cam;
pub mod colorant;
mod error;
pub mod illuminant;
pub mod lab;
pub mod math;
pub mod observer;
pub mod prelude;
pub mod rgb;
pub mod spectrum;
pub mod stimulus;
pub mod traits;
pub mod xyz;

pub use error::Error;
