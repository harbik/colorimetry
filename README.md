# Colorimetry

![Build Status](https://github.com/harbik/colorimetry/actions/workflows/build-and-test.yml/badge.svg)

<!-- cargo-rdme start -->

This is a Rust library for spectral colorimetry — the science of measuring and predicting color based on the spectral composition of light.
It implements the core calculation methods defined in **CIE 15:2018 Colorimetry**[^cie15] and the color rendering metrics defined in **CIE 224:2017 Colour Fidelity Index**[^cie224] / **ANSI/IES TM-30**.
Use it to build lighting-quality tools, simulate how colors shift under different light sources, convert between color spaces, or compute accurate observer-weighted tristimulus values.

It also has early support for JavaScript and WebAssembly, so you can run it in the browser and use it with JavaScript runtimes such as Deno.

## Usage

To use this library in a Rust application, run the command:

 ```bash
    cargo add colorimetry
```

or add this line to the dependencies in your Cargo.toml file:

```text
    colorimetry = "0.0.8"
```

## Examples

<details>
<summary><strong>Calculate Tristimulus Values for Illuminants</strong></summary>

This example calculates the XYZ tristimulus values of the D65 illuminant for both the CIE 1931 2º standard observer and the CIE 2015 10º observer.

```rust
  use colorimetry::illuminant::D65;

  // D65 Tristimulus values, using the CIE1931 standard observer by default
  let xyz_d65 = D65.xyz(None).set_illuminance(100.0);

  let [x, y, z] = xyz_d65.to_array();
  // [95.04, 100.0, 108.86]

  // D65 Tristimulus values using the CIE2015 10º observer
  use colorimetry::observer::Observer::Cie2015_10;
  let xyz_d65_10 = D65
    .xyz(Some(Cie2015_10)).set_illuminance(100.0);

  let [x_10, y_10, z_10] = xyz_d65_10.to_array();
  //[94.72, 100.0, 107.143]
```

</details>

<details>
<summary><strong>Calculate Correlated Color Temperature and Tint</strong></summary>

The correlated color temperature (CCT) of an illuminant, typically expressed in kelvin (K),
describes whether a light source appears warm (low CCT) or cool (high CCT). It is a key parameter
for characterizing the visual appearance of white light .
This example calculates both the correlated color temperature and the deviation from the Planckian
locus, often referred to as the tint.

```rust
  use colorimetry::illuminant::A;

  // Calculate CCT and Duv for the A illuminant
  // Requires `cct`, and `cie-illuminants` features
  let [cct, duv] = A.cct().unwrap().to_array();
  // [2855.4977, 0.0]
```

</details>

<details>
<summary><strong>Calculate Color Fidelity Index for Illuminants</strong></summary>

The **CIE 2017 Colour Fidelity Index** (**R<sub>f</sub>**, CIE 224:2017[^cie224]) is the modern replacement for the
older Color Rendering Index (**R<sub>a</sub>**, CRI). Both indices quantify how accurately a light source
reproduces the colors of illuminated objects relative to a reference illuminant. The CFI is significantly
more accurate: it uses 99 Color Evaluation Samples (CES) covering a broad range of real-world colors and
applies the CIECAM02-UCS perceptual color space, compared to CRI’s 8 pastel samples and the outdated
CIE 1964 (U\*V\*W\*) space. This library implements the version harmonized with **ANSI/IES TM-30-20/24**
(scaling constant CF = 6.73).

Below is an example calculation of the general Colour Fidelity Index for the CIE F2 illuminant:

```rust
  use colorimetry::illuminant::F2;

  // Calculate the Color Fidelity Index of the CIE F2 standard illuminant
  // Requires `cfi`, and `cie-illuminants` features
  let cf_f2 = F2.cfi().unwrap();
  let cf = cf_f2.general_color_fidelity_index();
  // 70.3
```

</details>

<details>

<summary><strong>Calculate Spectral Locus to Plot in Diagrams</strong></summary>

The spectral locus is the boundary in a chromaticity diagram that encloses all perceivable,
physically realizable colors. Due to its shape, it is sometimes informally referred to as the
"horseshoe."
Below, we compute the chromaticity coordinates that define the spectral locus.

```rust
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
<summary><strong>Calculate XYZ/RGB Transformation Matrices, for any Observer, for use in Color Profiles</strong></summary>

This is usually done with the CIE 1931 Standard Observer, but this library supports any observer—as long as both the color space and the data use the same one.
Instead of fixed XYZ values, it computes conversions from the spectral definitions of the primaries to be able to do so.
Here, we compute transformation matrices for the `DisplayP3` color space using both the `Cie1931` and `Cie2015` observers.

```rust
  use colorimetry::observer::Observer;
  use colorimetry::rgb::RgbSpace::DisplayP3;

  let xyz2rgb_31 = Observer::Cie1931.xyz2rgb_matrix(DisplayP3);
  //  2.4933, -0.9313, -0.4027,
  // -0.8298,  1.7629,  0.0236,
  //  0.0355, -0.076,   0.9574

  let rgb2xyz_31 = Observer::Cie1931.rgb2xyz_matrix(DisplayP3);
  // 0.4866, 0.2656, 0.1981,
  // 0.2291, 0.6917, 0.0792,
  // 0.0001, 0.0451, 1.0433,

  use colorimetry::observer::Observer::Cie2015;

  let xyz2rgb_15 = Cie2015.xyz2rgb_matrix(DisplayP3).clone();
  //  2.5258,  -1.0009, -0.3649,
  // -0.9006,   1.8546, -0.0011,
  //  0.0279,  -0.0574,  0.95874
```

</details>

<details>
<summary><strong>Find the closest spectral match in the Munsell Color Book</strong></summary>

This example finds the best matching color in the Munsell Color Book for a given sample—in this case, the R9 test color used in the CRI color rendering method.
It uses the `CieCam16::de_ucs` color difference metric and the `Cie2015_10` standard observer to calculate perceptual similarity.

The closest match identified is Munsell "5R 5/14", a vivid red hue, with a color difference of just 3 ΔE.
In practical terms, a ΔE of 3 is considered a close match—just at the threshold where most observers might start to notice a difference under controlled viewing conditions.

```rust
  // requires `cri` and `munsell` features
  use colorimetry::observer::Observer::Cie2015_10;
  use colorimetry::colorant::{MunsellCollection, TCS};

  let cri_r9 = &TCS[8];
  let (key, delta_e) = MunsellCollection::match_ciecam16(
    cri_r9,
    None,
    None,
    Some(Cie2015_10),
  ).unwrap();
  // ("5R4/14", 2.85)
```

</details>

<details>
<summary><strong>Match a paint color to a <i>sRGB</i> display value under realistic viewing conditions</strong></summary>

This example matches the Munsell paint chip <i>5 BG 5/8</i>—a teal/blue-green color—
to its nearest <i>sRGB</i>
<span style="display: inline-block; width: 1em; height: 1em; background-color: rgb(0, 113, 138);
border-radius: 50%; vertical-align: middle; border: 1px solid #000;"></span>
<span>rgb(0, 113, 138)</span>
equivalent, mimicking real-world viewing conditions.

Instead of the traditional <i>CIE 1931 2°</i> observer, this match uses the <i>CIE 2015 10° observer</i>,
which more accurately reflects how paint colors appear on walls. The illumination is based on the warm-white
<i>LED_B2</i> standard illuminant (≈ 3000 K). Together, these adjustments help the display color reflect
what you'd actually see on a freshly painted surface.

```rust
  // requires `cie-illuminants`, and `munsell` features
  use colorimetry::{
    cam::{ViewConditions, CIE248_HOME_SCREEN},
    colorant::Munsell,
    illuminant::LED_B2,
    observer::Observer::{Cie1931, Cie2015_10},
    rgb::RgbSpace::SRGB,
  };

  let paint = Munsell::try_new("5BG5/8").unwrap();
  let vc = ViewConditions::average_surround(6.0);
  let cam_paint = Cie2015_10.ciecam16(&LED_B2, &paint, vc);
  let rgb_2015 = cam_paint
    .rgb(SRGB, Some(CIE248_HOME_SCREEN))
    .unwrap()
    .compress();

  // Use a spectral representation of the Cie2015_10 RGB pixel, using the `Rgb`'s Light trait,
  // and calculate its XYZ tristimulus and RGB values for the CIE 1931 standard observer, the
  // observer
  // required for the sRGB color space.
  let xyz_1931 = Cie1931.xyz(&rgb_2015, None);
  let rgb_1931 = xyz_1931.rgb(SRGB).compress();
  let [r, g, b]: [u8; 3] = rgb_1931.into();
  //  (0, 113, 138)
```

</details>

## Capabilities

- Uses a Standard fixed-grid spectral representation [`Spectrum`], with a wavelength domain ranging from 380 to 780 nanometers, with 1-nanometer intervals, ensuring high precision in capturing fine spectral details critical for accurate colorimetry calculations.

  - Uses [`nalgebra`] vector and matrix types for fast integration and transformations, chosen for its high performance, numerical stability, and compatibility with other mathematical libraries.
  - Supports interpolation from irregular spectral data using [`Spectrum::linear_interpolate`] and [`Spectrum::sprague_interpolate`].
  - Optional smoothing using a Gaussian filter via [`Spectrum::smooth`], which is typically used to reduce noise in spectral data or to smooth out irregularities in measured spectra for better analysis.

- Generate spectral distributions from analytical models:

  - [`Illuminant::planckian`] Planck’s law for blackbody radiators
  - [`Illuminant::led`] Spectral power distribution of a LED
  - [`Colorant::gaussian`] Gaussian color filters
  - [`Stimulus::from_rgb`] Spectral distribution of an RGB color pixel using Gaussian spectral primaries

- Includes Spectral Representations CIE Standard Illuminants (optional):

  - Daylight: [`D65`], [`D50`]
  - Incandescent: [`A`]
  - Fluorescent: [`F1`], [`F2`], [`F3`], [`F4`], [`F5`], [`F6`], [`F7`], [`F8`], [`F9`], [`F10`], [`F11`], [`F12`]
  - Extended fluorescent set: [`F3_1`], [`F3_2`], [`F3_3`], [`F3_4`], [`F3_5`], [`F3_6`], [`F3_7`], [`F3_8`], [`F3_9`], [`F3_10`], [`F3_11`], [`F3_12`], [`F3_13`], [`F3_14`], [`F3_15`]
  - LED: [`LED_B1`], [`LED_B2`], [`LED_B3`], [`LED_B4`], [`LED_B5`], [`LED_BH1`], [`LED_RGB1`], [`LED_V1`]

- Includes Various Colorant Collections (optional):

  - Munsell Color System [`MunsellCollection`], with over 1,000 colors
  - Color Evaluation Samples [`CES`], a set of 99 test colors used in the Color Fidelity Index (CFI) calculations

- Calculate Illuminant metrics:

  - [`CCT`] Correlated color temperature (CIE 15:2018[^cie15] §9), including distance to the blackbody locus (Duv) for tint indication[^1]
  - [`CRI`] Color Rendering Index R<sub>a</sub> / R<sub>1</sub>–R<sub>14</sub>[^2]
  - [`CFI`] CIE 2017 Colour Fidelity Index R<sub>f</sub> / R<sub>f,1</sub>–R<sub>f,99</sub> (CIE 224:2017[^cie224])[^3]:
    - [`CFI::general_color_fidelity_index`] — overall R<sub>f</sub> score (0–100)
    - [`CFI::general_color_gamut_index`] — gamut index R<sub>g</sub>, area of the 16-bin a′b′ polygon relative to the reference (ANSI/IES TM-30)
    - [`CFI::rf_hj`] — per-hue-bin fidelity R<sub>f,hj</sub> for each of the 16 hue sectors
    - [`CFI::rcs_hj`] — per-hue-bin chroma shift R<sub>cs,hj</sub> (fraction; positive = saturation boost)
    - [`CFI::rhs_hj`] — per-hue-bin hue shift R<sub>hs,hj</sub> (radians, wrapped to (−π, π])

- Use Advanced color (appearance) models:

  - [`CieLab`], [`CieCam02`], [`CieCam16`]
  - Color difference methods: [`CieLab::ciede`], [`CieLab::ciede2000`], [`CieCam02::de_ucs`], [`CieCam16::de_ucs`]

- Work with Spectrally based RGB color spaces, with support for non-CIE 1931 observers and generic transformations between [`Rgb`] and [`XYZ`]:

  - RGB Color Spaces [`RgbSpace::SRGB`],  [`RgbSpace::Adobe`], [`RgbSpace::DisplayP3`]
  - Define RGB colors using spectral primaries, allowing switching observers

- Includes Multiple CIE Standard Observers (see [`Observer`]), with transformations to [`XYZ`]:

  - [`Observer::Cie1931`] the CIE 1931 2º standard observer
  - [`Observer::Cie1964`] the CIE 1964 10º standard observer (optional, enabled by default)1
  - [`Observer::Cie2015`] the CIE 2015 2º cone fundamentals-based observer (optional, enabled by default)
  - [`Observer::Cie2015_10`] the CIE 2015 10º cone fundamentals-based observer (optional, enabled by default)

- Compute gamut boundaries for CIE XYZ ([`RelXYZGamut`]) and CIELAB ([`CieLChGamut`]) using optimal colors ([`OptimalColors`]).

  - Estimate the total number of perceptually distinct colors.
  - Plot CIE LCh iso-hue lines and iso-lightness contours on chromaticity diagrams.
  - Plot maximum luminance value contours on chromaticity diagrams.
  - Use in hue-neutral gamut-compression algorithms.

[^cie15]: CIE 15:2018 *Colorimetry*, 4th edition. Commission Internationale de l'Éclairage. ISBN 978-3-902842-13-8. Available at <https://cie.co.at/publications/colorimetry-4th-edition>.
[^cie224]: CIE 224:2017 *Colour Fidelity Index for Accurate Scientific Use*. Commission Internationale de l'Éclairage. ISBN 978-3-902842-55-8. Available at <https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use>.
[^1]: McCamy, C. S. (1992). Correlated color temperature as an explicit function of chromaticity coordinates. *Color Research & Application*, 17(2), 142–144. The Ohno (2014) method is used internally for higher accuracy.
[^2]: CIE 13.3-1995 *Method of Measuring and Specifying Colour Rendering Properties of Light Sources*. The R<sub>a</sub> (CRI) metric is based on 8 test color samples in CIE 1964 (U\*V\*W\*) space; see Wyszecki & Stiles (2000) for background.
[^3]: CIE 224:2017[^cie224] defines R<sub>f</sub> using 99 CES in CIECAM02-UCS J′a′b′ space with a softplus formula (CF = 6.73). This library also implements the per-bin metrics R<sub>f,hj</sub>, R<sub>cs,hj</sub>, R<sub>hs,hj</sub> and the gamut index R<sub>g</sub> from ANSI/IES TM-30-20/24.

## Standards

This library implements calculations defined in the following CIE and IES standards.
The documents are not bundled with the library — each developer must obtain their own copy from the issuing body.

| Standard | Topic | Purchase |
|---|---|---|
| **CIE 15:2018** | Colorimetry, 4th edition — tristimulus values, standard observers, chromaticity, CCT, CRI | <https://cie.co.at/publications/colorimetry-4th-edition> |
| **CIE 224:2017** | Colour Fidelity Index for accurate scientific use — R<sub>f</sub>, R<sub>g</sub>, CES, CIECAM02-UCS | <https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use> |
| **ANSI/IES TM-30-20/24** | IES method for evaluating light source color rendition — per-bin R<sub>f,hj</sub>, R<sub>cs,hj</sub>, R<sub>hs,hj</sub>, CVG | <https://www.ies.org/store/> |

## Features

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
  Enables Color Fidelity Index (CFI / R<sub>f</sub>) calculations following CIE 224:2017 / ANSI/IES TM-30.
  Provides the overall R<sub>f</sub>, per-sample R<sub>f,1</sub>–R<sub>f,99</sub>, gamut index R<sub>g</sub>,
  and the 16-bin metrics R<sub>f,hj</sub> / R<sub>cs,hj</sub> / R<sub>hs,hj</sub>.
  Loads the 99 CES test color sample spectra at compile time.

<details>
<summary>How to enable a feature?</summary>

To enable a feature, such as `cri` and `munsell`, use

```bash
cargo add colorimetry -F cri,munsell
```

or, if you prefer to use the `cargo add` command with the `--features` flag, you can run:

```bash
cargo add colorimetry --features cri,munsell
```

Alternatively, configure features manually in your `Cargo.toml`:

```toml
colorimetry = { version = "0.0.8", features = ["cri", "munsell"] }
```

</details>

## Command Line Tool

This library has an associated command-line tool, named `color`, contained in the "colorimetry-cli" crate.
It provides a convenient way to perform colorimetric calculations, convert spectral data, and make color plots directly from the terminal.

To use it you have to install it using `cargo`, which in turn requires that you have Rust and Cargo installed.
Install them first by following the instructions at <https://www.rust-lang.org/tools/install>.
After having done this, run the following command:

```bash
cargo install colorimetry-cli
```

That should be it, you can now use the `color` command in your terminal.
Run the following command to see the available options and features for the `colorimetry` tool:
```bash
color --help
```

## Color Plots

The `colorimetry` library includes a plotting module, in an associated `colorimetry-plot` crate
that can be used to generate a number of color diagrams, spectral plots, and color rendering visulizations.
It's ouput is in SVG format, which can be viewed in any modern web browser or vector graphics editor, such as Inkscape.
The `colorimetry-plot` crate is not included by default, so you need to add it to your project using the following command:

```bash
cargo add colorimetry-plot
```
You can then use the `colorimetry-plot` crate to generate color plots.

## Developer Tasks with `xtask`

This project uses a Rust-based `xtask` utility for common development tasks:

- `cargo xtask check` – to run clippy, fmt, check, and readme verification,
- `cargo xtask test` – test various feature configurations,
- `cargo xtask doc` – build and open project documentation (fails on warnings), and,
- `cargo xtask wasm` – to generate the web-assembly files in `pkg` (requires `wasm-pack` and `wasm-opt`)

More commands will be added as the project evolves.

## License

All content &copy;2025 Harbers Bik LLC, and licensed under either of the

- Apache License, Version 2.0,
  ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>), or the
- MIT license
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
[`CFI::general_color_fidelity_index`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.general_color_fidelity_index
[`CFI::general_color_gamut_index`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.general_color_gamut_index
[`CFI::rf_hj`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.rf_hj
[`CFI::rcs_hj`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.rcs_hj
[`CFI::rhs_hj`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.rhs_hj
[`Colorant::gaussian`]: https://docs.rs/colorimetry/latest/colorimetry/colorant/struct.Colorant.html#method.gaussian
[`Stimulus::from_rgb`]: https://docs.rs/colorimetry/latest/colorimetry/stimulus/struct.Stimulus.html#method.from_rgb
[`Observer`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html
[`Observer::Cie1931`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html#variant.Cie1931
[`Observer::Cie1964`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html#variant.Cie1964
[`Observer::Cie2015`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html#variant.Cie2015
[`Observer::Cie2015_10`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html#variant.Cie2015_10
[`XYZ`]: https://docs.rs/colorimetry/latest/colorimetry/xyz/struct.XYZ.html
[`Rgb`]: https://docs.rs/colorimetry/latest/colorimetry/rgb/struct.RGB.html
[`RgbSpace::SRGB`]: https://docs.rs/colorimetry/latest/colorimetry/rgb/enum.RgbSpace.html#variant.SRGB
[`RgbSpace::Adobe`]: https://docs.rs/colorimetry/latest/colorimetry/rgb/enum.RgbSpace.html#variant.Adobe
[`RgbSpace::DisplayP3`]: https://docs.rs/colorimetry/latest/colorimetry/rgb/enum.RgbSpace.html#variant.DisplayP3
[`CieLChGamut`]: https://docs.rs/colorimetry/latest/colorimetry/lab/struct.CieLchGamut.html
[`RelXYZGamut`]: https://docs.rs/colorimetry/latest/colorimetry/xyz/struct.RelXYZGamut.html
[`OptimalColors`]: https://docs.rs/colorimetry/latest/colorimetry/observer/struct.OptimalColors.html

[`CES`]: https://docs.rs/colorimetry/latest/colorimetry/colorant/static.CES.html
[`MunsellCollection`]: https://docs.rs/colorimetry/latest/colorimetry/colorant/struct.MunsellCollection.html

[`D65`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.D65.html
[`D50`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.D50.html
[`A`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.A.html
[`F1`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F1.html
[`F2`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F2.html
[`F3`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3.html
[`F4`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F4.html
[`F5`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F5.html
[`F6`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F6.html
[`F7`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F7.html
[`F8`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F8.html
[`F9`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F9.html
[`F10`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F10.html
[`F11`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F11.html
[`F12`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F12.html
[`F3_1`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_1.html
[`F3_2`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_2.html
[`F3_3`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_3.html
[`F3_4`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_4.html
[`F3_5`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_5.html
[`F3_6`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_6.html
[`F3_7`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_7.html
[`F3_8`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_8.html
[`F3_9`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_9.html
[`F3_10`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_10.html
[`F3_11`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_11.html
[`F3_12`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_12.html
[`F3_13`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_13.html
[`F3_14`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_14.html
[`F3_15`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.F3_15.html
[`LED_B1`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.LED_B1.html
[`LED_B2`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.LED_B2.html
[`LED_B3`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.LED_B3.html
[`LED_B4`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.LED_B4.html
[`LED_B5`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.LED_B5.html
[`LED_BH1`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.LED_BH1.html
[`LED_RGB1`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.LED_RGB1.html
[`LED_V1`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/static.LED_V1.html

<!-- cargo-rdme end -->
