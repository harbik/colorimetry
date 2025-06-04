
![Build Status](https://github.com/harbik/colorimetry/actions/workflows/build-and-test.yml/badge.svg)

A Rust library for color modeling in illumination and engineering projects, with early JavaScript/WebAssembly support.
Algorithms follow standards from the CIE, ICC, and IES.
It intends to provide a comprehensive framework for spectral colorimetry.

# Usage
To use this library in a Rust application, run the command:
 ```bash
    cargo add colorimetry
```
or add this line to the dependencies in your Cargo.toml file:
```toml
    colorimetry = "0.0.6"
```

# Examples
<details>
<summary>Tristimulus Values</summary>

```rust
    use colorimetry::prelude::*;
    use colorimetry::illuminant::D65;
    use approx::assert_abs_diff_eq as check;

    // D65 Tristimulus values, using the CIE1931 standard observer by default
    let xyz_d65 = D65.xyz(None).set_illuminance(100.0);

    let [x, y, z] = xyz_d65.values();
    check!([x, y, z].as_ref(), [95.04, 100.0, 108.86].as_ref(),  epsilon = 5E-3);

    # #[cfg(feature = "supplemental-observers")]
    # {
    // bring in the 2015‐10º observer only if that feature is enabled
    use colorimetry::observer::Observer::Std2015_10;

    // D65 Tristimulus values using the CIE2015 10º observer
    let xyz_d65_10 = D65
        .xyz(Some(Std2015_10))
        .set_illuminance(100.0);

    let [x_10, y_10, z_10] = xyz_d65_10.values();
    check!([x_10, y_10, z_10].as_ref(), [94.72, 100.0, 107.143].as_ref(), epsilon = 5E-3);
    # }
```
</details>

<details>
<summary>Correlated Color Temperature</summary>

```rust
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
<summary>Color Fidelity Index</summary>

Here is an example calculating the general color fidelity index of the CIE F2 illuminant:
```rust
    use colorimetry::prelude::*;
    # #[cfg(feature="cie-illuminants")]
    use colorimetry::illuminant::F2;
    use approx::assert_abs_diff_eq as check;

    // Calculate the Color Fidelity Index of the CIE F2 standard illuminant
    // Requires `cfi`, and `cie-illuminants` features
    # #[cfg(all(feature="cfi", feature="cie-illuminants"))]
    let cf_f2 = F2.cfi().unwrap();
    # #[cfg(all(feature="cfi", feature="cie-illuminants"))]
    let cf = cf_f2.general_color_fidelity_index();
    # #[cfg(all(feature="cfi", feature="cie-illuminants"))]
    check!(cf, 70.3,  epsilon = 1E-1);
```
</details>


# Features

- [`Spectrum`] Standard fixed grid spectral representation, over a wavelength domain from 380 to 780 nanometers with 1 nanometer intervals.
- Transformations from irregular spectral data using [`Spectrum::linear_interpolate`] and [`Spectrum::sprague_interpolate`] interpolation, and optional smoothing using a Gaussian filter [`Spectrum::smooth`].
- Uses [`nalgebra`] vector/matrix definitions, for fast integration and transformations with access to a large set of algorithms.
- Generate spectral distributions from analytical models
  - [`Illuminant::planckian`] Planck’s law for blackbody radiators  
  - [`Illuminant::led`] Spectral power distribution of a LED
  - [`Colorant::gaussian`] Gaussian color filters
  - [`Stimulus::from_rgb`] Specral distribution of an RGB color pixel using Gaussian spectral primaries.
- CIE Standard Illuminants
  - Daylight Illuminants: [`D65`], [`D50`]
  - Incandescent lamp: [`A`]
  - Fluorescent standard lamps: [`F1`], [`F2`], [`F3`], [`F4`], [`F5`], [`F6`], [`F7`], [`F8`], [`F9`], [`F10`], [`F11`], [`F12`] fluorescent lamps
  - Extended set of Fluorescent lamps: [`F3_1`], [`F3_2`], [`F3_3`], [`F3_4`], [`F3_5`], [`F3_6`], [`F3_7`], [`F3_8`], [`F3_9`], [`F3_10`], [`F3_11`], [`F3_12`], [`F3_13`], [`F3_14`], [`F3_15`]
  - Standard LED Illuminants: [`LED_B1`], [`LED_B2`], [`LED_B3`], [`LED_B4`], [`LED_B5`], [`LED_BH1`], [`LED_RGB1`], [`LED_V1`]
- Illuminant metrics
  - Correlated color temperature [`CCT`], including distance to blackbody locus for indication of its tint
  - Color rendering index [`CRI`],
  - Color fidelity index [`CFI`].
- Advanced color (appearance) models
  - [`CieLab`],
  - [`CieCam02`]
  - [`CieCam16`]
- <u>Spectral based red, green, and blue color spaces</u>, with support for non-CIE1931 observers, and generic transformations between  [`Rgb`] and [`XYZ`] color models:
  - sRGB
  - Adobe RGB
  - DisplayP3
- <u>CIE Standard Observers</u> support (see [`Observer`] and [`ObserverData`]) with transformations to [`XYZ`] tristimulus values;
  - CIE 1931 2º,
  - CIE 1964 10º,
  - CIE 2015 2º,
  - CIE 2015 10º,
- Color Difference Metrics

# Usage
To use this library in a Rust application, run the command:
 ```bash
    cargo add colorimetry
```
or add this line to the dependencies in your Cargo.toml file:
```toml
    colorimetry = "0.0.6"
```
The easiest way to use the objects and functions in this library is through its prelude.

```rust
    use colorimetry::prelude::*;
    use colorimetry::illuminant::D65;

    // D65 Tristimulus values, using the CIE1931 standard observer by default
    let [x, y, z] = D65.xyz(None).set_illuminance(100.0).values();

    approx::assert_ulps_eq!(x, 95.04, epsilon = 5E-3);
    approx::assert_ulps_eq!(y, 100.0, epsilon = 5E-3);
    approx::assert_ulps_eq!(z, 108.86, epsilon = 5E-3);
```

# Flags

- `supplemental-observers` (enabled by default)
- `cie-illuminants`
- `munsell` Include reflection spectra for Munsell colors.
- `cct`
  Calculates correlated color temperatures (CCT) for illuminants.
  Generates a 4096-entry lookup table (each entry containing three `f64` values).
  Memory is reserved at compile time but computed on demand.
  (Included automatically if the `cri` feature is enabled).
- `cri`
  Enables Color Rendering Index (CRI) calculations, providing Ra and R1–R14 values for illuminants.  
  Loads an additional 14 test color sample spectra.

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

# Spectral Data
All spectra—illuminants, filters, surface reflectance, or stimuli—are defined on a fixed internal grid of 401 samples spanning 380 nm to 780 nm in 1 nm increments, using the [`Spectrum`] object, and the [`nalgebra`] linear-algebra library.
This allows for efficient vector and matrix operations on spectral data, and provides access to advanced mathematical routines.  

For importing measurements from **arbitrary or irregular** wavelength grids, data can be resampling using interpolation methods (linear or Sprague–Karup), with optional **smoothing** for oversampled datasets.  

Being a spectral colorimetry library, it contains predefined spectral datasets, and facilities to generate spectral data from parameters.

Data sets which are contained within the library are implemented are:
- measurements of the **Munsell** color book patches (using the `munsell` featurej),
- the - **CRI TCS** (Test Color Samples, with `cri` feature), 
- and the 99 **TM30/CFI CES** test color samples (`cfi` feature).

Furthermore the library contains all the **CIE Illuminants** such as:
- **D65** and **D50** daylight illuminants,
- the fluorescent **Fx**, **F3x**-series,
- and the ***LEDxx** series.

Besides predefined datasets, spectral data can be generated from models, such as **Planckian** thermal emitters and **LED




- **Advanced color models**
  -  **[CIELAB]**
  -  **[CIECAM16]**


## Installation



## Features

This library includes a range of spectral data collections, with only a minimal set enabled by default.  
Additional functionality can be activated using Cargo feature flags.



<details>
<summary><strong>Default Features</strong></summary>
These features are enabled by default. To disable, use:

```bash
cargo add colorimetry --no-default-features
```

- **`cie-illuminants`**  
  Adds a large collection of standard illuminants (e.g., Fluorescent and LED series) beyond the **D50** and **D65**, which are always included.
    

- **`supplemental-observers`**  
  Adds the following CIE standard colorimetric observers beyond the **CIE 1931** 2º Standard Observer:
  - **CIE 1964** 10° Standard Observer  
  - **CIE 2015** Cone Fundamental based Standard Observers (2° & 10°)

</details>


## Spectral Distributions

The [`Spectrum`] struct underpins all spectral calculations in this library. It stores spectral data in a `nalgebra::SVector<f64, NS>` vector, with length `NS = 401`, over a wavelength domain ranging from **380 nm to 780 nm** in **1 nm** steps.

To perform colorimetric calculations, use either [`Illuminant`] for light source spectra, or [`Colorant`] for modeling surface colors, such as paints and printed inks.
A [`Stimulus`] is used to model pixels in displays, where a combination of red, green, and blue sub-pixels are controlled to create sensations of color, directly viewed by looking at them.

<details>
<summary><strong>Intialize from Array</strong></summary>

```rust
    use colorimetry::prelude::*;
    use colorimetry::spectrum::NS;

    // a black hole stimulus spectrum, using NS = 401 zero values
    let black_hole_spectrum = Spectrum::new([0.0; NS]);
    
    // the stimulus, reaching our eyes, when looking at a black hole:
    let black_hole_stimulus = Stimulus::new(black_hole_spectrum);
```

</details>

<details>
<summary><strong>Using data with other wavelength domains</strong></summary>

If you have spectral data defined over a wavelength domain different from the _380-780-1 nanometers_ as used in this library, you can use two interpolation methods converting your data into a `Spectrum`:

- **Linear interpolation**  
  The [`Spectrum::linear_interpolate`] constructor takes a slice of wavelengths and a corresponding slice of spectral values. It returns a `Spectrum` if both slices are of equal length and the wavelengths are ordered.

- **Sprague interpolation**  
  For smoother interpolation, [`Spectrum::sprague_interpolate`] implements the [Sprague-Karup interpolation method](https://www.sciencedirect.com/science/article/pii/0771050X75900273), commonly used in color science. This method requires that the input wavelengths be evenly spaced. It takes the domain bounds and a slice of spectral values as input and produces a high-resolution `Spectrum` aligned with the internal wavelength grid.

</details>

<details>
<summary><strong>References</strong></summary>
This spectral domain aligns with standards such as:

- [CIE 15:2004 – Colorimetry](https://archive.org/details/gov.law.cie.15.2004)
- [IES LM-79-08 – Electrical and Photometric Measurements of Solid-State Lighting Products](https://webstore.ansi.org/preview-pages/IESNA/preview_IESNA%2BLM-79-08.pdf)

This 380–780 nm range is also the default domain used by the [IES TM-30 Spectral Calculator](https://www.ies.org/standards/standards-toolbox/tm-30-spectral-calculator/).
</details>


## Illuminants
An [`Illuminant`] is a spectral representation of a light which hits an object, which, upon scattering on its way to our eyes, creates the sensations of color we experience.
The spectral composition of an illuminant influences the colors we see.

Illuminants can be created using spectral data, in form of a `Spectrum` instance, or can be generated from various spectral models.
Alternatively, the library includes the CIE standard illuminants.

<details>
<summary><strong>Initialize from Spectrum</strong></summary>
To get an `Illuminant` from your spectral data, first create a `Spectrum`, for example by using one of the interpolation methods, or directly using an array.

```rust
    use colorimetry::prelude::*;

    // create equal energy spectrum from an array, with values of 1.0.
    let spectrum = Spectrum::new([1.0; 401]);
    let illuminant = Illuminant::new(spectrum);
    
    // calculate chromaticity coordinates as used in the CIE 1931 chromaticity diagram
    // use `None` as argument to used the default CIE 1931 2º standard observer
    let chromaticity = illuminant.xyz(None).chromaticity();
    
    // check the values
    approx::assert_abs_diff_eq!(chromaticity.x(), 0.3333, epsilon=1E-4);
    approx::assert_abs_diff_eq!(chromaticity.y(), 0.3333, epsilon=1E-4);
```

</details>

<details>
<summary><strong>Factory Functions</strong></summary>

- **Planckian illuminant**, a pure thermal emission based spectrum.
  Uses Plank's law, and takes an absolute temperature, in Kelvin, as argument.
  ```rust
      use crate::colorimetry::prelude::*;

      // Plankian illuminant with a temperature of 3000 Kelvin
      let p3000 = Illuminant::planckian(3000.0);
      let chromaticity = CIE1931.xyz(&p3000, None).chromaticity();

      approx::assert_abs_diff_eq!( chromaticity.x(), 0.436_935, epsilon = 1E-6);
      approx::assert_abs_diff_eq!( chromaticity.y(), 0.404_083, epsilon = 1E-6);
  ```

- Generic Daylight **CIE D-illuminant,** generating a daylight spectrum with a characteristic
  correlated color temperature in the range from 4000 to 25_000 Kelvin.

- **LED illuminant**, with a spectral distribution described by an analytical function,
  as proposed by Yoshi Ohno, as published in _Optical Engineering 44(11)_, 2005.

- **Equal Energy Illuminant**, with a uniform spectral distribution with an irradiance of 1 watt per square meter.
</details>

<details><summary><strong>CIE Standard Illuminants </strong></summary>

The following standard illuminants are available in the library using the _cie-illuminants_ feature, which is enabled by default.


</details>

<details>
<summary><strong>Correlated Color Temperature (CCT)</strong></summary>

Illuminants are typically characterized by their **correlated color temperature (CCT)**, expressed in Kelvin (K), and by their **tint**, which describes the chromaticity deviation from the Planckian (blackbody) locus.

The **CCT** is defined as the temperature of the Planckian (ideal blackbody) radiator whose perceived color most closely matches that of the test light source, when viewed under identical conditions. Because many real-world light sources (e.g., fluorescent or LED lamps) do not emit light that exactly matches any blackbody radiator, their color temperature is termed *correlated* rather than exact.

CCT is not derived directly from spectral data, but is calculated using the chromaticity coordinates by finding the closest point on the Planckian locus—usually by minimizing the Euclidean or perceptual distance in color space[^3].

In this library, an advanced, high accuracy, iterative Robertson's method is used to calculate both values.

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

The CIE Color Rendering Index (CRI), including the general color rendering index (Rₐ) and the individual special color rendering indices (R₁ through R₁₅), can be calculated using the `cri` method, which follows the procedure specified by the CIE[^2].


  ```rust
    # #[cfg(all(feature = "cri", feature = "cie-illuminants"))]{
    // this example requires `cri` and `cie-illuminants` features enabled

    use crate::colorimetry::prelude::*;

    let f3_11 = CieIlluminant::F3_11.illuminant();
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

A [`Colorant`] represents a color filter or surface (such as a color patch), defined by dimensionless spectral values ranging from 0.0 to 1.0.

Each value corresponds to the proportion of light at a given wavelength that is **not absorbed**:

- `0.0` → Full absorption (no transmission or reflection)
- `1.0` → Full transmission or reflection (no absorption)

This model is commonly used in color science to describe the spectral behavior of materials and follows conventions used in CIE colorimetry.

In color models, a [`Colorant`] spectrum in not used directly, but is always associated with an illuminant;
without an illuminant, objects appear black.

<details>
<summary><strong>Initialize from Spectrum</strong></summary>

A [`Colorant`] is a wrapper around Spectrum and can be created using its new method, which accepts a Spectrum and ensures that all values lie within the range 0.0 to 1.0.
If any value falls outside this range, the constructor returns an error.

</details>

<details>
<summary><strong>Factory Functions</strong></summary>

The library defines different model based factory functions.
Here are a couple of examples.

```rust
use crate::colorimetry::prelude::*;

// Create a perfect white `Colorant` or color patch, with no absorption.
let white = Colorant::white(); 

// Create a gray neutral colorant with 30% reflectance at all wavelengths
let gray = Colorant::gray(0.3); 

// A perfect absorber absorbing all light.
let black = Colorant::black();

// A `top_hat` colorant or rectangular bandfilter, defined by a center wavelength,
// and a width, both expressed in units of nanometer or meters.
let green_mono = Colorant::top_hat(550.0, 1.0);

// A `gaussian` shaped colorant, defined by a center values, and a standard deviation 
// `sigma` value, with a peak value of 1.0. 
let red = Colorant::gaussian(610.0, 5.0);

```
</details>

<details>
<summary><strong>Mixing and Adding</strong></summary>

`Colorant`s support several operations useful for simulating physical interactions with light:

- **Multiplication** models **subtractive mixing**, such as placing multiple filters in sequence or layering ink and pigment. For example, multiplying two `Colorant`s simulates how their combined spectral transmissions reduce the overall light passing through.
  
- **Addition** is helpful in constructing synthetic spectral functions—such as building up a custom filter shape by combining multiple Gaussians. Any resulting values above `1.0` are clipped to `1.0`.

- **Scalar multiplication** adjusts the transmission intensity of a `Colorant`. This can be used to simulate partial transparency or adjust the concentration of a dye. Resulting values are clamped to the valid `[0.0, 1.0]` range: values above `1.0` become `1.0`, and those below `0.0` become `0.0`.

</details>

<details>
<summary><strong>CIELAB</strong></summary>
The [`Colorant::cielab`] method calculates a colorant's CIELAB values.
Here is an example calculating the CIELAB coordinates for a perfect white colorant:

```rust
  use crate::colorimetry::prelude::*;

  let colorant = Colorant::white();

  // use default (None) D65 illuminant  and default CIE 1931 standard observer (second None)
  let [l, a, b] = colorant.cielab(None, None).values();

  approx::assert_abs_diff_eq!(l, 100.0, epsilon = 1E-4); // L* should be 100 for white
  approx::assert_abs_diff_eq!(a, 0.0, epsilon = 1E-4); // a* should be 0 for white
  approx::assert_abs_diff_eq!(b, 0.0, epsilon = 1E-4); // b* should be 0 for white
```

</details>


## Stimuli

A [`Stimulus`] wraps a [`Spectrum`] representing the spectral power distribution of light as it arrives at an observer or sensor. This could be a pixel in a camera sensor, or a set of photoreceptors in the human eye. The spectral data is expressed in physical radiometric terms, with units corresponding to **luminance**: _candelas per square meter_ (cd/m²), integrated over the visible range.

A [`Stimulus`] may describe:
- Emitted light from a self-luminous source, such as a display pixel  
- Reflected light from an object surface element illuminated by a known source  
- Transmitted light after passing through a wavelength-selective medium, such as a colored filter element

In all cases, the stimulus encapsulates the final spectral signal available for visual or digital perception, after any combination of emission, reflection, or transmission events.

<details>
<summary><strong>Initialize from Spectrum</strong></summary>

A [`Stimulus`] is a wrapper around Spectrum and can be created using its new method, which accepts a Spectrum and ensures that all values lie within the range 0.0 to 1.0.
If any value falls outside this range, the constructor returns an error.

</details>

<details>
<summary><strong>Factory functions</strong></summary>

- [`Stimulus::from_srgb`], and [`Stimulus::from_rgb`], create a `Stimulus` of a set of RGB pixel values.
  The first takes three `u8` arguments, while the second uses a [`Rgb`] object as argument.
  This function allows calculating the perceived color difference between different observers, from the perspective of a single observer.

  ```rust
  use colorimetry::prelude::*;
  let red = Stimulus::from_srgb(255, 0, 0);
  let red_chromaticity = CIE1931.xyz(&red, None).chromaticity();
  approx::assert_abs_diff_eq!(
      red_chromaticity.to_array().as_ref(),
      &[0.64, 0.33].as_ref(),
      epsilon = 1E-5
  );
    ```
</details>


## Standard Observers

In colorimetry, color perception is modeled as the response of the human visual system to spectral stimuli. Human vision is trichromatic, based on the relative excitations of three cone types (L, M, and S) in the retina. These physiological responses are abstracted in the CIE XYZ color space as X, Y, and Z tristimulus values.

Tristimulus values are computed by integrating a spectral power distribution with a set of three **color matching functions** (CMFs), which represent the average spectral sensitivity of the human eye. A set of CMFs defines a **standard observer**.

The primary observer used in most applications is the **CIE 1931 2º Standard Observer**, derived from color matching experiments with foveal (central) vision. This observer is represented in this library by a static [`Observer`] instance: [`CIE1931`].

With the **`supplemental-observers`** feature enabled, the library also includes:
- [`CIE1964`] and the [`CIE1964`] standard observer for larger visual fields,
- [`CIE2015`] and [`CIE2015_10`] cone fundamentals–based observers (CIE 170-2:2015).




[`nalgebra`]:https://docs.rs/nalgebra/latest/nalgebra/ 
[`Spectrum`]: https://docs.rs/colorimetry/latest/colorimetry/spectrum/struct.Spectrum.html
[`Spectrum::linear_interpolate`]: https://docs.rs/colorimetry/latest/colorimetry/spectrum/struct.Spectrum.html#method.linear_interpolate
[`Spectrum::sprague_interpolate`]: https://docs.rs/colorimetry/latest/colorimetry/spectrum/struct.Spectrum.html#method.sprague_interpolate
[`Spectrum::smooth`]: https://docs.rs/colorimetry/latest/colorimetry/spectrum/struct.Spectrum.html#method.smooth
[`Illuminant`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.Illuminant.html
[`Illuminant::planckian`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.Illuminant.html#method.planckian
[`Illuminant::led`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.Illuminant.html#method.led
[`CCT`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CCT.html
[`CRI`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CRI.html
[`Colorant`]: https://docs.rs/colorimetry/latest/colorimetry/colorant/struct.Colorant.html
[`Colorant::gaussian`]: https://docs.rs/colorimetry/latest/colorimetry/colorant/struct.Colorant.html#method.gaussian
[`Stimulus`]: https://docs.rs/colorimetry/latest/colorimetry/stimulus/struct.Stimulus.html
[`Stimulus::from_srgb`]: https://docs.rs/colorimetry/latest/colorimetry/stimulus/struct.Stimulus.html#method.from_srgb
[`Stimulus::from_rgb`]: https://docs.rs/colorimetry/latest/colorimetry/stimulus/struct.Stimulus.html#method.from_rgb
[Colorimetric Observers]: https://docs.rs/colorimetry/latest/colorimetry/observer/index.html
[`Observer`]: https://docs.rs/colorimetry/latest/colorimetry/observer/enum.Observer.html
[`ObserverData`]:https://docs.rs/colorimetry/latest/colorimetry/observer/enum.ObserverData.html 
[`Observer.xyz`]: https://docs.rs/colorimetry/latest/colorimetry/observer/struct.ObserverData.html#method.xyz
[`CIE1931`]: https://docs.rs/colorimetry/latest/colorimetry/data/observers/static.CIE1931.html
[`CIE1964`]: https://docs.rs/colorimetry/latest/colorimetry/data/observers/static.CIE1964.html
[`CIE2015`]: https://docs.rs/colorimetry/latest/colorimetry/data/observers/static.CIE2015.html
[`CIE2015_10`]: https://docs.rs/colorimetry/latest/colorimetry/data/observers/static.CIE2015_10.html
[`XYZ`]: https://docs.rs/colorimetry/latest/colorimetry/xyz/struct.XYZ.html
[`Rgb`]: https://docs.rs/colorimetry/latest/colorimetry/rgb/struct.RGB.html
[`WideRgb`]: https://docs.rs/colorimetry/latest/colorimetry/widergb/struct.WideRgb.html
[CIECAM16]: https://docs.rs/colorimetry/latest/colorimetry/cam/struct.CieCam16.html
[CIELAB]: https://docs.rs/colorimetry/latest/colorimetry/lab/struct.CieLab.html
[`RgbSpace`]: https://docs.rs/colorimetry/latest/colorimetry/rgbspace/enum.RgbSpace.html
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

[^1]: “ColorChecker” is a registered trademark of X-Rite, Incorporated
[^2]: CIE 13.3-1995: Method of Measuring and Specifying Colour Rendering Properties of Light Sources* (Commission Internationale de l'Éclairage, 1995).
[^3]: Commission Internationale de l'Éclairage. (2004). *CIE 015:2004: Colorimetry* (3rd ed.). Vienna: CIE.
