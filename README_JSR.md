
# Overview
The Colorimetry Library is used for color calculations in illumination and color engineering projects.
It is written in Rust but provides WebAssembly JavaScript/WebAssembly interface as well, as published here.
The algorithms implemented try to follow the recommendations of the International Commission on Illumination,
the International Color Consortium, the Illumination Engineering Society, and many others.

JavaScript is mainly known as the programming language for web applications, but it can also be used to write scripts to run on your computer, similar to Python.
To do that, you need a _runtime_, a program that reads and executes your script.
[Node](node.js) is such a runtime which has been around for a while.
Recently, alternative runtimes, such as [Deno](deno.com), have been developed, which are much easier to use and do not need any configuration.
Deno even allows the compilation of a script into an executable program, allowing it to be run on computers without installing the runtime itself.

Here is a brief overview of the main objects in this library, with some introductory examples for use in Deno and Web Applications.

# Use in Deno
[Deno](deno.com) is, per the Deno site,

> Deno is the open-source JavaScript runtime for the modern web. Built on web
> standards with zero-config TypeScript, unmatched security, and a complete
> built-in toolchain, Deno is the easiest, most productive way to JavaScript.


For Deno installation instructions, see the following: [install Deno](https://docs.deno.com/runtime/manual/getting_started/installation/).

# `XYZ` Tristimulus Values
Tristimulus values, typically represented by the symbols X, Y, and Z, form the basis for all colorimetric models.

Here is an example of initializing a tri-stimulus value set from two chromaticity coordinates.
```ts
// Copy paste this code in a new file, and save as: `xyz.ts`
// Open a new command line terminal, and change to the directory where you saved the file.
// Type `deno run --allow-read xyz.ts` to see this script executed.

import init, * as cmt from "jsr:@harbik/colorimetry";
await init();

// Create a new XYZ object using D65 CIE 1931 chromaticity coordinates
const xyz = new cmt.XYZ(0.31272, 0.32903);

// Get and check the corresponding tristimulus values, with a luminous value
// of 100.0
const [x, y, z] = xyz.values();
console.log("D65 Tristimulus Values are", x, y, z);

// and get back the orgiinal chromaticity coordinates:
const [xc, yc] = xyz.chromaticity();
console.log("D65 Chromaticitiy values are", xc, yc);


// to get the luminous value:
const l = xyz.luminousValue();
console.log("Confirm that the Luminous Value is 100");

```
!include license.txt

To run this example, the `--allow-read` permission has to be used, as the script fetches the required WebAssembly library files from this repository.

## Use `Spectrum` for spectral data
Tristimulus values are calculated from spectral distributions, either directly measured with a _spectrometer_, calculated from physics models, or obtained from spectral libraries.

All spectral calculations in this library use the `Spectrum` class, which contains the spectral data and spectrum type.
It uses a wavelength domain from 380 to 780 nanometers, with 1 nanometer intervals, as recommended in [CIE15:2004](https://archive.org/details/gov.law.cie.15.2004).
This results in inconsistencies with older CIE data, as other wavelength domains and intervals were often used for recommended calculations.
In particular, chromaticity coordinates for the standard illuminants D65, D50, and A differ slightly from historic published values.

This example calculates the Illuminance and CIE 1931 (x, y) chromaticity
coordinates for a Planckian (thermal emission-based) illuminator with a
Correlated Color Temperature of 3000 Kelvin using the CIE 1931 standard observer.

```rust
use crate::colorimetry::{Spectrum, CIE1931};
use approx::assert_ulps_eq;

let p3000 = Spectrum::planckian_illuminant(3000.0);
let [l, x, y] = CIE1931.xyz(&p3000).lxy();

assert_ulps_eq!(l, 20.668_927, epsilon = 1E-6);
assert_ulps_eq!(x, 0.436_935, epsilon = 1E-6);
assert_ulps_eq!(y, 0.404_083, epsilon = 1E-6);
```

Besides the `Spectrum::planck` constructor, `Spectrum` has many other constructors.
For example, `Spectrum::d65`, `Spectrum::d50`, and `Spectrm::a` provide spectral distributions of the CIE D65, D50, and A standard illuminants.

## The CIE Standard Colorimetric `Observer`
`CIE1931` is an instance of the `Observer` class representing colorimetric standard observers and also includes the CIE 1976 10º standard observers and the CIE 2015 2º and 10º cone fundamental derived observers.
Other instances are `CIE1976`, for the CIE 10º standard observer and `CIE2015``, and `CIE2015_10` for the 2º and 10º observers derived from the Cone Fundamentals.

Its primary function `CIE1931.xyz` takes a spectral distribution as a single argument with an instance of the `XYZ` class, encapsulating the X, Y, and Z tristimulus values.
Likewise, the `CIE1931.lab_d65` and `CIE1931.lab_d50` methods can be used to get CIELAB coordinates for a spectrum measured from a color sample.
These result in an instance of the `Lab` class.


## `XYZ` Tristimulus Values

## `Lab` Color Model

## `RGB` Color Spaces


# Use with Deno/TypeScript



# Use in Web Applications

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