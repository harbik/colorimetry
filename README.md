
# Overview
A Colorimetry Library for color modelling in illumination and engineering projects.
Besides being a Rust library, it also provides JavaScript/WebAssembly interfaces, although still rudimentary at this stage.
The algorithms implemented try to follow the recommendations of the International Commission on Illumination (CIE),
the International Color Consortium (ICC), and the Illumination Engineering Society (IES).

It is a spectral library, using spectral representations of color, allowing to use more advanced
colorimetric observers besides the outdated and flawed CIE 1931 standard observer.
It uses spectral representations of observers, illuminants, filters, and color patches defined over
a wavelength domain from 380 to 780 nanometers, with 1 nanometer steps.
Linear, Sprague, and Spline interpolation are available if your data uses other domains.

# Installation

To use this library in a Rust application, run the command:
 ```bash
    cargo add colorimetry
```
or add this line to the dependencies in your Cargo.toml file:
```toml
    colorimetry = "0.0.3"
```
The easiest way to use the objects and functions library is through its prelude.
This example calculates the chromaticity values of the CIE D65 illuminant.
```rust
    use colorimetry::prelude::*;
    let xy = CIE1931.xyz(&D65, None).chromaticity();
    approx::assert_abs_diff_eq!(xy.as_ref(), [0.3127, 0.3291].as_ref(), epsilon=5E-5);
```

## Features 

The library includes many data collections, of which only a minimal set is included by default.
A feature is enabled by including using the flag `-F`.
For example, to include "cri" Color Rendering Index module illuminants, use:
```bash
    cargo add colorimetry -F cri
```
This can also be set manually, in your cargo.toml file, using
```toml
    colorimetry = {version = "0.0.3", features = ["cri"]}
```

The current features in this library are:

- **cie-illuminants** _default_
    Include a large collection of standard illuminants such as the Fluorescent and LED series.
    Included by default. 
- **supplemental-observers** _default_
    The CIE 1931 Standard Observer is always included, but with feature several other standard and experimental
    colorimetric observers are included as well.
    Included by default.
- **munsell**
    Include reflection spectra for the Munsell colors.
    This will increase the size of your executable quite a bit, and is not included by default.
- **charts**
    Include reflection spectra for various test charts.
- **cri** 
    Include the color rendering index module, which calculates the Ra and R1 to R14 values for illuminants.
    This loads an additional 14 test color sample spectra.
- **cct**
    Calculate correlated color temperature for illuminants.
    Builds a 4096 length lookup table, with each row consisting of 3*f64 values.
    The table rows are only calculated when required, but table space is reserved in the executable.
    This module is also included with the "cri" feature.
- **color-fidelity**
    Calculates CIE 224:2017 Color Fidelity Index, and associated values.
    Contains 99 test color samples.

## Spectral Composition
All spectral calculations in this library use the [Spectrum] class as a base, which contains the spectral data.

For practical considerations, it uses a wavelength domain from 380 to 780 nanometers, with 1 nanometer intervals, as recommended in the [CIE15:2004](https://archive.org/details/gov.law.cie.15.2004) standard.
[Spectrum] uses a [`nalgebra::Vector<f64>`] type, with a length of 401 elements, to capture these spectral data.
Historically, different wavelength domains have been recommended and used by the CIE, such as ranges from 300 to 830 nanometers, and an interval size of 5 nanometers.
The choice of domain has a small impact on calculated colorimetric values, and the reference values calculated here can differ a bit from the ones published by the CIE in the past.


## Illuminants
We need light to see.
Objects 'get color' only when they are illuminated.
In this library an [Illuminant] is a spectral representation of the light which hits an object.

The most common illuminant is daylight.
The CIE has defined the D65 standard illuminant, and recommends using this as default daylight illuminant.
This library uses [StdIlluminant] for the CIE recommended standard illuminants, in particular [StdIlluminant::D65] for default daylight.
Daylight is not constant - it varies with time of day, season, and weather.


Another source of light are electric lamps, such as incandescent light bulbs.
They generate light by thermal emission from a very hot tungsten filament in a glass envelope.
In physics, the spectral properties of thermal emission is described by Planck's law.
For incandescent light bulbs, the CIE recommends using the A-illuminant, in this library available as [StdIlluminant::A].


This example calculates the Illuminance and CIE 1931 (x, y) chromaticity
coordinates for a Planckian (thermal emission-based) illuminator with a
Correlated Color Temperature of 3000 Kelvin using the CIE 1931 standard observer.

```rust
    use crate::colorimetry::prelude::*;

    let p3000 = Illuminant::planckian(3000.0);
    let xy = CIE1931.xyz(&p3000, None).chromaticity();

    approx::assert_abs_diff_eq!(xy.as_ref(), [0.436_935,0.404_083].as_ref(), epsilon = 1E-6);
```

Besides the [Illuminant::planckian] constructor, [Illuminant] has many other constructors.
For example, [Illuminant::d65] and [Illuminant::d50] provide spectral distributions of the CIE D65, and D50 standard illuminants, defined by the CIE in tabular form.
Many other Standard Illuminants can be used, such as the A, Fluorescent, and LED Illuminants defined by the CIE, when the library is compiled with the "cie-illuminants" feature.
This feature is a default feature, but can be disabled when not used and compact binaries are required.
For example, to get the A illuminant's chromaticity:
```rust
    use colorimetry::prelude::*;

    let xy_a = CIE1931.xyz(&StdIlluminant::A, None).chromaticity();
    // see <https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_A>
    approx::assert_abs_diff_eq!(xy_a.as_ref(), [0.44758, 0.40745].as_ref(), epsilon=1E-5)
```

Other interesting constructors are the [`Stimulus::srgb`], and [`Stimulus::rgb`], which create a spectrum of a set of [RGB] pixel values.
The first takes three `u8` arguments, while the second uses an [RGB] object as argument.

```rust
    use colorimetry::prelude::*;
    let red = Stimulus::srgb(255, 0, 0);
    approx::assert_abs_diff_eq!(
        CIE1931.xyz(&red, None).chromaticity().as_ref(),
        &[0.64, 0.33].as_ref(),
        epsilon = 1E-5
    );
```

When dealing with spectral data defined over a domain not matching with the one used in this library you can use [Spectrum::linear_interpolate].
This function can also be used to define arbitrary spectral shapes by providing a set of wavelength and spectral value points.

## The CIE Standard Colorimetric [Observer]
[CIE1931] is a static instance of the [Observer] class representing colorimetric standard observers, and is always included.
With the default **supplemental-observers** feature also other observers are included, such as the CIE 1976 10º, s, the CIE 2015 2º and CIE 2015 10º observers.
An observer is represented by three functions, called color matching functions, which are supposed to be an indirect representation of the spectral sensitivities of the _L_, _M_, and _S_ cones in the back of eyes.

The primary function of an [Observer], such as the [CIE1931] colorimetric standard observer, is the [`CIE1931.xyz`] method, which takes a spectral distribution as a single argument, and produces an [XYZ] object, encapsulating the CIE 1931 X, Y, and Z tristimulus values.

## [XYZ] Tristimulus Values
These tristimulus values are a representation of the response of each of the three cones, and an inproduct of the spectrum and the color matching functions.
All color models are using the tristimulus values of a stimulus, essentially a light ray being detected by a set of cones, as a basis.

## [`CieLab`] Color Model
Likewise, the [CIE1931.lab_d65] and [CIE1931.lab_d50] methods can be used to get CIELAB coordinates for a spectrum measured from a color sample, as an instance of the [`CieLab`](crate::lab::CieLab) class.

## Color

## [RGB] Color Values, and [RgbSpace] Color Spaces.


## Correlated Color Temperature

## Color Rendering Metrics



# Use with Deno/TypeScript



# Use in Web Applications

# License
All content &copy;2024 Harbers Bik LLC, and licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
 * MIT license
   [LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>

at your option.

## Contribution

Unless you explicitly state otherwise, any Contribution intentionally submitted
for inclusion in the Work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.