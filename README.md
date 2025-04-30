
# Overview

A Rust library for color modeling in illumination and engineering projects, with early JavaScript/WebAssembly support.  
Algorithms follow standards from the CIE, ICC, and IES.

This **spectral** library represents colors, illuminants, filters, and patches across 380–780 nm at 1 nm intervals, with support for importing data defined over different or irregular domains.  
It supports advanced colorimetric observers beyond the outdated CIE 1931 standard.

The library includes a broad range of spectral data by default or via feature flags.  
Custom datasets are supported, with linear, Sprague, and spline interpolation available for domain remapping.

# Installation

To use this library in a Rust application, run the command:
 ```bash
    cargo add colorimetry
```
or add this line to the dependencies in your Cargo.toml file:
```toml
    colorimetry = "0.0.4"
```
The easiest way to use the objects and functions library is through its prelude.

This example calculates the chromaticity values of the CIE D65 illuminant.
```rust
    use colorimetry::prelude::*;
    // use CIE 1931 standard observer as defaul
    let [x, y] = D65.xyz(None).chromaticity();
    approx::assert_abs_diff_eq!(x, 0.3127, epsilon=5E-5);
    approx::assert_abs_diff_eq!(y, 0.3291, epsilon=5E-5);
```

And here we calculate the CIELAB values of a Gaussian filter, with center wavelength of 550 nanometers, and with a _FWHM_ (Full Width Half Maximum) of 25 nanometers:
```rust
    use colorimetry::prelude::*;
    // using D65 lluminant and CIE 1931 standard observer as defaults:
    let [l, a, b] = Colorant::gaussian(550.0, 25.0).cielab(None, None).values();
    approx::assert_abs_diff_eq!(l, 77.26, epsilon=5E-3);
    approx::assert_abs_diff_eq!(a, -72.93, epsilon=5E-3);
    approx::assert_abs_diff_eq!(b, 104.58, epsilon=5E-3);
```


## Features

This library includes a range of spectral data collections, with only a minimal set enabled by default.  
Additional functionality can be activated using Cargo feature flags.

<details>
<summary><strong>How to Enable Features</strong></summary>

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
colorimetry = { version = "0.0.4", features = ["cri"] }
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

## Spectral Distributions

The [`Spectrum`](crate::spectrum::Spectrum) struct underpins all spectral calculations in this library. It stores data in a `nalgebra::Vector<f64>` of length 401, representing wavelengths from **380 nm to 780 nm** in **1 nm** steps.

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
    // create spectrum from array, consisting of spectral values of 1.0.
    use colorimetry::prelude::*;
    let spectrum = Spectrum::new([1.0; 401]);
    let illuminant = Illuminant::new(spectrum);
    let xy = illuminant.xyz(None).chromaticity();
    approx::assert_abs_diff_eq!(xy.as_ref(), [0.3333, 0.3333].as_ref(), epsilon=1E-4)
```
</details>

<details>
<summary><strong>Illuminant Models</strong></summary>

- **Planckian illuminant**, a pure thermal emission based spectrum.
  Uses Plank's law, and takes an absolute temperature, in Kelvin, as argument.
  ```rust
      use crate::colorimetry::prelude::*;

      let p3000 = Illuminant::planckian(3000.0);
      let xy = CIE1931.xyz(&p3000, None).chromaticity();

      approx::assert_abs_diff_eq!( xy.as_ref(), [0.436_935,0.404_083].as_ref(), epsilon = 1E-6);
  ```

- Generic Daylight **CIE D-illuminant,** generating a daylight spectrum with a characteristic
  correlated color temperature in the range from 4000 to 25_000 Kelvin.

- **LED illuminant**, with a spectral distribution described by an analytical function,
  as proposed by Yoshi Ohno, as published in _Optical Engineering 44(11)_, 2005.

- **Equal Energy Illuminant**, with a uniform spectral distribution with an irradiance of 1 watt per square meter.
</details>

<details>
<summary><strong>CIE Standard Illuminants</strong></summary>
The CIE has defined a set of standard illuminants, in form of spectral data values.
They represent various light sources, and are included in this library using 
[`StdIlluminant](crate::std_illuminants::StdIlluminant), which is a `enum` type.

They are all included by default, and can be 
To use the others, use the `

<details>
<summary><i>Daylight Illuminants</i></summary>

- [`D65`](crate::std_illuminants::StdIlluminant::D65)
- [`D50`](crate::std_illuminants::StdIlluminant::D65)

</details>

<details>
<summary><i>Standard Fluorescent Lamps</i></summary>

- **F1**
- **F2**
- **F3**
- **F4**
- **F5**
- **F6**
- **F7**
- **F8**
- **F9**
- **F10**
- **F11**
- **F12**

</details>





</details>


`Illuminant` and `StdIlluminant` both implement the [`Light`](crate::traits::Light) trait, which is used as generic input for color models.

## Stimuli
Other interesting constructors are the [`Stimulus::srgb`](crate::stimulus::Stimulus::srgb), and [`Stimulus::rgb`](crate::stimulus::Stimulus::rgb), which create a spectrum of a set of RGB pixel values.
The first takes three `u8` arguments, while the second uses a [`RGB`](crate::rgb::RGB) object as argument.

```rust
    use colorimetry::prelude::*;
    let red = Stimulus::srgb(255, 0, 0);
    approx::assert_abs_diff_eq!(
        CIE1931.xyz(&red, None).chromaticity().as_ref(),
        &[0.64, 0.33].as_ref(),
        epsilon = 1E-5
    );
```


## The CIE Standard Colorimetric Observer
What we perceive as color are sensations in the virtual cortex located in the back of our brain.
These sensations are triggered by stimuli from the photosensitive layer in the back of our eyes, called retina.
Human vision is trichromatic, which means that light, when entering our eyes, is classified by three stimuli.
In colorimetry, at the physiological level, these stimuli are represented by the X, Y, and Z tristimulus values, using the CIE XYZ color model, and using three types of color sensitivity functions, called color matching functions, for an average or standard observer.
Color matching functions are indirect representations of the spectral sensitivities of the _L_, _M_, and _S_ cones in our retinas, as function of spectral stimuli entering our eyes.

The first and currently still almost exclusively used standard observer is the CIE 1931 Colorimetric Standard Observer.
The definition of this observer by the CIE launched the field of colorimetry.
In this library it is represented by a static instance of the [`Observer`](crate::observer::Observer) class called `CIE1931` and is always available.
With the default **supplemental-observers** feature also other observers are included, such as the CIE 1976 10º, s, the CIE 2015 2º and CIE 2015 10º observers.

The primary function of a [`Observer`](crate::observer::Observer), such as the [`CIE1931`](crate::data::observers::CIE1931) colorimetric standard observer, is the [`CIE1931.xyz`] method, which takes a spectral distribution as a single argument, and produces a [`XYZ`](crate::xyz::XYZ) object, encapsulating the CIE 1931 X, Y, and Z tristimulus values.
These tristimulus values are used by more advanced color models, such as the CIELAB and CIECAM, to describe the sensations of color in our minds.

## Paints, Dyes, and Inks
These are represented by a [`Colorant`](crate::colorant::Colorant) object, which encapsulates a reflection or transmission spectrum, defined over a domain from 380 to 780 nanometers with 1 nanometer steps, and values within a range from 0.0 to 1.0.
They implement the [`Filter`](crate::traits::Filter) trait, which is used as input for many of the color models.

## XYZ Tristimulus Values
These tristimulus values are a representation of the response of each of the three cones, and an inproduct of the spectrum and the color matching functions.
All color models are using the tristimulus values of a stimulus, essentially a light ray being detected by a set of cones, as a basis.

Although they can be initiated in this library directly using the [`XYZ::new`](crate::xyz::XYZ::new) constructor, they are typical produced by using the [`Observer::xyz`](crate::observer::Observer.xyz) function, which takes a generic `Light` and an optional `Filter` argument.
Examples of object which implement the `Light` trait are `StdIlluminant` and `Illuminant`.


## CieLab Color Model
Likewise, the `lab_d65` and `lab_d50` methods can be used to get CIELAB coordinates for a spectrum measured from a color sample, as an instance of the [`CieLab`](crate::lab::CieLab) class.

## Color

## [`RGB`](crate::rgb::RGB) Color Values, and [`RgbSpace`](crate::rgbspace::RgbSpace) Color Spaces.


## Correlated Color Temperature

## Color Rendering Metrics



# Use with Deno/TypeScript



# Use in Web Applications

# License
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