
Colorimetry
===========

A Rust library for color modeling in illumination and engineering projects, with early JavaScript/WebAssembly support.  
Algorithms follow standards from the CIE, ICC, and IES.

This **spectral** library represents colors, illuminants, filters, and patches across 380–780 nm at 1 nm intervals, with support for importing data defined over different or irregular domains.  
It supports advanced colorimetric observers beyond the outdated CIE 1931 standard.

The library includes a broad range of spectral data by default or via feature flags.  
Custom datasets are supported, with linear, Sprague, and spline interpolation available for domain remapping.

## Installation

To use this library in a Rust application, run the command:
 ```bash
    cargo add colorimetry
```
or add this line to the dependencies in your Cargo.toml file:
```toml
    colorimetry = "0.0.4"
```
The easiest way to use the objects and functions library is through its prelude.


## Features

This library includes a range of spectral data collections, with only a minimal set enabled by default.  
Additional functionality can be activated using Cargo feature flags.

<details>
<summary><strong>Enable Features</strong></summary>

To enable a feature when adding the library:

```bash
cargo add colorimetry -F cri
```
or
```bash
cargo add colorimetry --features cri
```

Alternatively, configure features manually in your `Cargo.toml`:

```toml
// add `cri` and `color-fidelity` features
colorimetry = { version = "0.0.4", features = ["cri", "color-fidelity"] }
```

</details>

<details>
<summary><strong>Default Features</strong></summary>

- **`cie-illuminants`**  
  Includes a large collection of standard illuminants (e.g., Fluorescent and LED series).  
  Enabled by default. To disable, use:

  ```bash
  cargo add colorimetry --no-default-features
  ```

- **`supplemental-observers`**  
  Adds several standard and experimental colorimetric observers beyond the CIE 1931 Standard Observer, which is always included.  
  Enabled by default.

</details>

<details>
<summary><strong>Optional Features</strong></summary>

- **`cie-illuminants`**
  Only **D50** and **D65** are included by default.
  Use this feature to include the **A**, **F_x**, **F3_x** and **LED_x** illuminants.

- **`munsell`**  
  Includes reflection spectra for Munsell colors.  
  _Note: significantly increases executable size._

- **`charts`**  
  Adds reflection spectra for various standard test charts.

- **`cri`**  
  Enables Color Rendering Index (CRI) calculations, providing Ra and R1–R14 values for illuminants.  
  Loads an additional 14 test color sample spectra.

- **`cct`**  
  Calculates correlated color temperatures (CCT) for illuminants.  
  Generates a 4096-entry lookup table (each entry containing three `f64` values).  
  Memory is reserved at compile time but computed on demand.  
  _Included automatically when the `cri` feature is enabled._

- **`color-fidelity`**  
  Calculates the CIE 224:2017 Color Fidelity Index and related metrics.  
  Includes 99 test color sample spectra.

</details>

## Concepts and Examples

- While the core data structure of the library is `Spectrum`, all colorimetric models operate on higher-level abstractions: `Illuminant`, `Colorant`, or `Stimulus`, each encapsulating spectral data in a way that aligns with standard colorimetric principles.
- Depending on the spectral type, the library provides specialized methods—for example:
  - Correlated color temperature **CCT** and color rendering index **CRI** calculations for `Illuminant`s.
  - **CIELAB** and **CIECAM** color appearance model computations for `Colorant`s.
  - **RGB** conversions and colorimetric projections for `Stimulus` data.
- The library includes multiple CIE standard observers as `Observer` instances, including the CIE 1931 2º, CIE 1964 10º, and the cone fundamentals–based CIE 2015 2º and 10º observers.

As a first step in almost any colorimetric model, an `Observer` is used to compute the tristimulus values `XYZ`, which indirectly represent the responses of the eye's cone receptors.

- Advanced color models such as **CIELAB** and **CIECAM** build on these tristimulus values to describe our color perceptions—how we see color—taking into account the state of visual adaptation and viewing conditions.
- These models also enable the estimation of perceived color differences between stimuli, making it possible to compare how colors appear across displays, printed materials, and real-world objects.

## Spectral Distributions

The [`Spectrum`](crate::spectrum::Spectrum) struct underpins all spectral calculations in this library. It stores data in a `nalgebra::SVector<f64, NS>` with length `NS = 401`, representing wavelengths from **380 nm to 780 nm** in **1 nm** steps.

To perform colorimetric calculations, use either [`Illuminant`](crate::illuminant::Illuminant) for light source spectra, or [`Colorant`](crate::colorant::Colorant) for modeling surface colors, such as paints and printed inks.
A [`Stimulus`](crate::stimulus::Stimulus) is used to model pixels in displays, where a combination of red, green, and blue sub-pixels are controlled to create sensations of color, directly viewed by looking at them.

<details>
<summary><strong>Using data with other wavelength domains</strong></summary>
If you have spectral data using another wavelength domain, two mapping functions are available to create a `Spectrum` from your data:

- **Linear interpolation** The [`linear_interpolate`](crate::spectrum::Spectrum::linear_interpolate)
  constructor takes a slice of wavelengths and a slice of spectral data as arguments, and produces a `Spectrum` if both slices have the same length.

- **Sprague interpolation**

</details>

<details>
<summary><strong>References</strong></summary>
This spectral domain aligns with standards such as:

- [CIE 15:2004 – Colorimetry](https://archive.org/details/gov.law.cie.15.2004)
- [IES LM-79-08 – Electrical and Photometric Measurements of Solid-State Lighting Products](https://webstore.ansi.org/preview-pages/IESNA/preview_IESNA%2BLM-79-08.pdf)

This 380–780 nm range is also the default domain used by the [IES TM-30 Spectral Calculator](https://www.ies.org/standards/standards-toolbox/tm-30-spectral-calculator/).
</details>


## Illuminants
An [`Illuminant`](crate::illuminant::Illuminant) is a spectral representation of a light which hits an object, which, upon scattering on its way to our eyes, creates the sensations of color we experience.
The spectral composition of an illuminant influences the colors we see.

Illuminants can be created using spectral data, in form of a `Spectrum` instance, or can be generated from various spectral models.
Alternatively, the library includes the CIE standard illuminants.

<details>
<summary><strong>From Spectral Data</strong></summary>
To get an `Illuminant` from your spectral data, first create a `Spectrum`, for example by using one of the interpolation methods, or directly using an array.

```rust
    use colorimetry::prelude::*;

    // create equal energy spectrum from an array, with values of 1.0.
    let spectrum = Spectrum::new([1.0; 401]);
    let illuminant = Illuminant::new(spectrum);
    // Use None for default CIE 1931 2º standard observer
    let [x, y] = illuminant.xyz(None).chromaticity();
    approx::assert_abs_diff_eq!(x, 0.3333, epsilon=1E-4);
    approx::assert_abs_diff_eq!(y, 0.3333, epsilon=1E-4);
```

</details>

<details>
<summary><strong>Illuminant Models</strong></summary>

- **Planckian illuminant**, a pure thermal emission based spectrum.
  Uses Plank's law, and takes an absolute temperature, in Kelvin, as argument.
  ```rust
      use crate::colorimetry::prelude::*;

      // Plankian illuminant with a temperature of 3000 Kelvin
      let p3000 = Illuminant::planckian(3000.0);
      let [x, y] = CIE1931.xyz(&p3000, None).chromaticity();

      approx::assert_abs_diff_eq!( x, 0.436_935, epsilon = 1E-6);
      approx::assert_abs_diff_eq!( y, 0.404_083, epsilon = 1E-6);
  ```

- Generic Daylight **CIE D-illuminant,** generating a daylight spectrum with a characteristic
  correlated color temperature in the range from 4000 to 25_000 Kelvin.

- **LED illuminant**, with a spectral distribution described by an analytical function,
  as proposed by Yoshi Ohno, as published in _Optical Engineering 44(11)_, 2005.

- **Equal Energy Illuminant**, with a uniform spectral distribution with an irradiance of 1 watt per square meter.
</details>

<details>
<summary><strong>CIE Standard Illuminants </strong></summary>
<i>Daylight Illuminants</i>

[`D65`](crate::std_illuminants::StdIlluminant::D65), 
[`D50`](crate::std_illuminants::StdIlluminant::D65)
</details>

<details>
<summary><strong>Additional CIE Standard Illuminants (use `cie-illuminants` feature) </strong></summary>
<i>Standard Incandescent Lamp</i>

[`A`](crate::std_illuminants::StdIlluminant::A),

<i>Fluorescent Lamps, Standard Series</i>

[`F1`](crate::std_illuminants::StdIlluminant::F1),
[`F2`](crate::std_illuminants::StdIlluminant::F2),
[`F3`](crate::std_illuminants::StdIlluminant::F3),
[`F4`](crate::std_illuminants::StdIlluminant::F4),
[`F5`](crate::std_illuminants::StdIlluminant::F5),
[`F6`](crate::std_illuminants::StdIlluminant::F6),
[`F7`](crate::std_illuminants::StdIlluminant::F7),
[`F8`](crate::std_illuminants::StdIlluminant::F8),
[`F9`](crate::std_illuminants::StdIlluminant::F9),
[`F10`](crate::std_illuminants::StdIlluminant::F10),
[`F11`](crate::std_illuminants::StdIlluminant::F11),
[`F12`](crate::std_illuminants::StdIlluminant::F12)

<i>Fluorescent Lamps, F3 Series</i>

[`F3_1`](crate::std_illuminants::StdIlluminant::F3_1),
[`F3_2`](crate::std_illuminants::StdIlluminant::F3_2),
[`F3_3`](crate::std_illuminants::StdIlluminant::F3_3),
[`F3_4`](crate::std_illuminants::StdIlluminant::F3_4),
[`F3_5`](crate::std_illuminants::StdIlluminant::F3_5),
[`F3_6`](crate::std_illuminants::StdIlluminant::F3_6),
[`F3_7`](crate::std_illuminants::StdIlluminant::F3_7),
[`F3_8`](crate::std_illuminants::StdIlluminant::F3_8),
[`F3_9`](crate::std_illuminants::StdIlluminant::F3_9),
[`F3_10`](crate::std_illuminants::StdIlluminant::F3_10),
[`F3_11`](crate::std_illuminants::StdIlluminant::F3_11),
[`F3_12`](crate::std_illuminants::StdIlluminant::F3_12),
[`F3_13`](crate::std_illuminants::StdIlluminant::F3_13),
[`F3_14`](crate::std_illuminants::StdIlluminant::F3_14),
[`F3_15`](crate::std_illuminants::StdIlluminant::F3_15),

<i>LED Illuminants</i>

[`LED_B1`](crate::std_illuminants::StdIlluminant::LED_B1),
[`LED_B2`](crate::std_illuminants::StdIlluminant::LED_B2),
[`LED_B3`](crate::std_illuminants::StdIlluminant::LED_B3),
[`LED_B4`](crate::std_illuminants::StdIlluminant::LED_B4),
[`LED_B5`](crate::std_illuminants::StdIlluminant::LED_B5),
[`LED_BH1`](crate::std_illuminants::StdIlluminant::LED_BH1),
[`LED_RGB1`](crate::std_illuminants::StdIlluminant::LED_RGB1),
[`LED_V1`](crate::std_illuminants::StdIlluminant::LED_V1),

</details>

<details>
<summary><strong>Correlated Color Temperature (CCT)</strong></summary>

Illuminants are typically characterized by their **correlated color temperature (CCT)**, expressed in Kelvin (K), and by their **tint**, which describes the chromaticity deviation from the Planckian (blackbody) locus.

The **CCT** is defined as the temperature of the Planckian (ideal blackbody) radiator whose perceived color most closely matches that of the test light source, when viewed under identical conditions. Because many real-world light sources (e.g., fluorescent or LED lamps) do not emit light that exactly matches any blackbody radiator, their color temperature is termed *correlated* rather than exact.

CCT is not derived directly from spectral data, but is calculated using the chromaticity coordinates (typically in the CIE 1931 (x, y) color space) by finding the closest point on the Planckian locus—usually by minimizing the Euclidean or perceptual distance in color space.

In this library, an advanced, high accuracy, iterative Robertson's method is used to calculate both values.

**Reference**

Commission Internationale de l'Éclairage. (2004). *CIE 015:2004: Colorimetry* (3rd ed.). Vienna: CIE.

Here we us Plank's law, to create an illuminant spectrum, and check its temperature and tint.
  ```rust
      # #[cfg(feature = "cct")]{
      // this example requires `cct` feature enabled
      use crate::colorimetry::prelude::*;

      // Plankian illuminant with a temperature of 3000 Kelvin
      let p3000 = Illuminant::planckian(3000.0);

      // calculate CCT and Duv for this illuminant
      // unwrap OK as we know values should be approximately 3000.0, and 0.0
      let [cct, duv] = p3000.cct().unwrap().values();

      approx::assert_abs_diff_eq!( cct, 3000.0, epsilon = 1E-4);
      approx::assert_abs_diff_eq!( duv, 0.0, epsilon = 1E-6);
    # }
  ```


</details>

<details>
<summary><strong>Correlated Color Rendering Index (CRI)</strong></summary>

The CIE Color Rendering Index (CRI), including the general color rendering index (Rₐ) and the individual special color rendering indices (R₁ through R₁₅), can be calculated using the `cri` method, which follows the procedure specified in *CIE 13.3-1995: Method of Measuring and Specifying Colour Rendering Properties of Light Sources* (Commission Internationale de l'Éclairage, 1995).


  ```rust
    # #[cfg(all(feature = "cri", feature = "cie-illuminants"))]{
    // this example requires `cri` and `cie-illuminants` features enabled

    use crate::colorimetry::prelude::*;

    let f3_11 = StdIlluminant::F3_11.illuminant();
    let cri = f3_11.cri().unwrap();

    let expected_ra = 78.0;
    approx::assert_abs_diff_eq!(cri.ra(), expected_ra, epsilon = 1.0);

    let expected_values = [
        90.0, 86.0, 49.0, 82.0, 81.0, 70.0, 85.0, 79.0, 24.0, 34.0, 64.0, 50.0, 90.0, 67.0,
    ];

    approx::assert_abs_diff_eq!(
        cri.values().as_ref(),
        expected_values.as_ref(),
        epsilon = 1.0
    );

    # }
  ```

</details>

## Colorants

A `Colorant` represents a color filter or surface (such as a color patch) defined by spectral values  
ranging from `0.0` to `1.0`. Each value corresponds to the proportion of light at a given wavelength  
that is **not absorbed**:

- `0.0` → Full absorption (no transmission or reflection)
- `1.0` → Full transmission or reflection (no absorption)

This model is commonly used in color science to describe the spectral behavior of materials and follows  
conventions used in CIE colorimetry.

In color models, a `Colorant` spectrum in not used directly, but is always associated with an illuminant;
without an illuminant, objects appear black.

<details>
<summary><strong>Create from spectral data</strong></summary>
A `Colorant` can be created from a [`Spectrum`](crate::spectrum::Spectrum)` which, besides using a
direct array, include various interpolation constructors and smoothing methods.

</details>

<details>
<summary><strong>Colorant Models</strong></summary>
The library defines different model based constructors.
Here are a couple of examples.

```rust

// Create a perfect white color patch.
let white = Colorant::white()); 

// Create a gray neutral colorant with 30% reflectance at all wavelengths
let gray = Colorant::gray(0.3); 

// A perfect absorber (black)
let black = Colorant::black();



```

</details>


## Stimuli

A `Stimulus` represents the spectral power distribution (SPD) of light reaching the eye —  
the physical input that gives rise to color perception.

It encapsulates a [`Spectrum`](crate::spectrum::Spectrum), which contains the  
spectral data (e.g., from a light source, reflected surface, or transmitted medium)  
as a function of wavelength.

In colorimetric terms, the `Stimulus` models the energy that interacts with the  
human visual system. When evaluated using a standard [`Observer`](crate::observer::Observer),  
it yields CIE XYZ tristimulus values that quantify the perceived color.

For example, a `Stimulus` might represent:
- Emitted light from a display pixel  
- Reflected light from a surface under an illuminant  
- Transmitted light through a colored filter

<details>
<summary><strong>Constructors</strong></summary>
</details>


## Standard Observers

In colorimetry, color perception is modeled as the response of the human visual system to spectral stimuli. Human vision is trichromatic, based on the relative excitations of three cone types (L, M, and S) in the retina. These physiological responses are abstracted in the CIE XYZ color space as X, Y, and Z tristimulus values.

Tristimulus values are computed by integrating a spectral power distribution with a set of three **color matching functions** (CMFs), which represent the average spectral sensitivity of the human eye. A set of CMFs defines a **standard observer**.

The primary observer used in most applications is the **CIE 1931 2º Standard Observer**, derived from color matching experiments with foveal (central) vision. This observer is represented in this library by a static [`Observer`](crate::observer::Observer) instance: [`CIE1931`](crate::data::observers::CIE1931).

With the **`supplemental-observers`** feature enabled, the library also includes:
- `CIE1964_10`: the CIE 1964 10º standard observer for larger visual fields,
- `CIE2015_2` and `CIE2015_10`: cone fundamentals–based observers defined over 390–830 nm (CIE 170-2:2015).

The core method [`Observer::xyz`](crate::observer::Observer::xyz) maps a spectral distribution to a [`XYZ`](crate::xyz::XYZ) tristimulus value. These serve as the basis for advanced color appearance models such as **CIELAB** and **CIECAM**, which incorporate adaptation state, luminance level, and surround context.




## Advanced Colorimetry
Here are some examples of advanced colorimetry task, facilitated by this library (_Under Development_).

<details>
<summary><strong>Color Perception Differences Between Different Observers</strong></summary>

[`Stimulus::srgb`](crate::stimulus::Stimulus::srgb), and [`Stimulus::rgb`](crate::stimulus::Stimulus::rgb), which create a `Stimulus` of a set of RGB pixel values.
The first takes three `u8` arguments, while the second uses a [`RGB`](crate::rgb::RGB) object as argument.
This function allows calculating the perceived color difference between different observers, from the perspective of a single observer.

```rust
    use colorimetry::prelude::*;
    let red = Stimulus::srgb(255, 0, 0);
    approx::assert_abs_diff_eq!(
        CIE1931.xyz(&red, None).chromaticity().as_ref(),
        &[0.64, 0.33].as_ref(),
        epsilon = 1E-5
    );
```


</details>



## License
All content &copy;2025 Harbers Bik LLC, and licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
 * MIT license
   [LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>

at your option.

## Contribution

Unless you explicitly state otherwise, any Contribution intentionally submitted
for inclusion in the Work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.