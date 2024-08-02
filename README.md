
# Overview
The Colorimetry Library is a library for color calculations in illumination and engineering projects.
It can be used for Rust projects and provides JavaScript/WebAssembly interfaces.
The algorithms implemented try to follow the recommendations by the International Commission on Illumination,
the International Color Consortium, the Illumination Engineering Society, and many others.

# The `Spectrum` class is used for all spectral calculations
All spectral calculations in this library use the `Spectrum` class, which contains the spectral data and spectrum type.
For practical considerations, it uses a wavelength domain from 380 to 780 nanometers, with 1 nanometer intervals, as recommended in [CIE15:2004](https://archive.org/details/gov.law.cie.15.2004).
This results in some inconsistencies with older data, as in the past, other wavelength domains were often being used for integration.
In particular, chromaticity coordinates for the standard illuminants D65, D50, and A differ slightly from published values.

This example calculates the Illuminance and CIE 1931 (x, y) chromaticity
coordinates for a Planckian (thermal emission-based) illuminator with a
Correlated Color Temperature of 3000 Kelvin using the CIE 1931 standard observer.

```rust
use crate::colorimetry::{Spectrum, CIE1931};
use approx::assert_ulps_eq;

let p3000 = Spectrum::planckian(3000.0);
let [l, x, y] = CIE1931.xyz(&p3000).lxy();

assert_ulps_eq!(l, 20.668_927, epsilon = 1E-6);
assert_ulps_eq!(x, 0.436_935, epsilon = 1E-6);
assert_ulps_eq!(y, 0.404_083, epsilon = 1E-6);
```

Besides the `Spectrum::planck` constructor, `Spectrum` has many other constructors.
For example, `Spectrum::d65`, `Spectrum::d50`, and `Spectrm::a` provide spectral distributions of the CIE D65, D50, and A standard illuminants.

# The CIE Standard Colorimetric `Observer`
`CIE1931` is an instance of the `Observer` class representing colorimetric standard observers and also includes the CIE 1976 10º standard observers and the CIE 2015 2º and 10º cone fundamental derived observers.
Its primary function `xyz`(s: &Spectrum) -> XYZ` takes a spectral distribution as a single argument with an instance of the `XYZ` class, encapsulating the X, Y, and Z tristimulus values.

# `XYZ` Tristimulus Values

# `Lab` Color Model

# `RGB` Color Spaces

# License
All content &copy;2024 Harbers Bik LLC, and licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>?)

at your option.

## Contribution

Unless you explicitly state otherwise, any Contribution intentionally submitted
for inclusion in the Work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.