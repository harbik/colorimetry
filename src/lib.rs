// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

/*!
This is a Rust library for **spectral colorimetry** — the science of measuring and predicting color
perception based on the full spectral power distribution of light sources, spectral reflectivity of surfaces,
and the spectral responsivity of human vision.
It implements the core methods defined in **CIE 15:2018 Colorimetry**[^cie15] and the color-quality
metrics defined in **CIE 224:2017 Colour Fidelity Index**[^cie224] / **ANSI/IES TM-30**.

## Why spectral?

Every light source emits a unique mix of wavelengths — its *spectral power distribution* (SPD).
Two light sources can look equally white yet render object colors very differently,
because their SPDs differ where the reflectance of an illuminated surface is sensitive.

We don't see the same colors either - the eye's response is also wavelength-dependent, varying across the visible spectrum and between observers,
changing with viewing conditions, age, and health.

Spectral colorimetry captures this by representing illuminants, surfaces, and filters as spectral
values, and using different observers' color-matching functions to predict how the eye perceives
color under any combination of light and material.

The result is a set of three, observer dependent, numbers — the **XYZ tristimulus values** — from which all other color
quantities are derived: chromaticity, CCT, CIELAB, color-rendering metrics, and more.

This approach is the foundation of all rigorous colorimetry work: lighting product development
and qualification, ICC color profile construction, paint-to-screen matching, gamut analysis,
and anything else where "what the camera sees" must match "what the eye sees."

# Usage

To use this library in a Rust application, run the command:

 ```bash
    cargo add colorimetry
```

or add this line to the dependencies in your Cargo.toml file:

```text
    colorimetry = "0.0.8"
```

## JavaScript and WebAssembly

The library also compiles to WebAssembly and ships a ready-to-use ES module, so the same
colorimetry calculations can run in a browser or in a JavaScript runtime such as Deno.

**Build the WASM package**

```bash
cargo xtask wasm          # requires wasm-pack and wasm-opt
```

This generates the `pkg/` directory containing `colorimetry.js` (ESM), `colorimetry_bg.wasm`,
and the TypeScript declarations `colorimetry.d.ts`.

**Browser**

Import from the [esm.sh](https://esm.sh) CDN inside a `<script type="module">` block.
The default export is an async `init()` function that compiles the `.wasm` binary; it
must be awaited before calling any library function.

> **Version note** — the published `colorimetry@0.0.8` WASM package predates enabling the
> `cri` export, so the `d65.cri().ra()` examples below require a newer package version
> (the first release after `0.0.8` that includes `cri`) or a locally built `pkg/` directory.

```html
<script type="module">
  import init, { Illuminant, CieIlluminant }
    from "https://esm.sh/colorimetry@0.0.8";

  await init();

  const d65 = Illuminant.illuminant(CieIlluminant.D65);
  console.log("D65 Ra:", d65.cri().ra());
</script>
```

**Deno**

Import from esm.sh using a URL import — no `npm:` prefix or local build required:

```typescript
import init, { Illuminant, CieIlluminant }
  from "https://esm.sh/colorimetry@0.0.8";

await init();

const d65 = Illuminant.illuminant(CieIlluminant.D65);
console.log("D65 Ra:", d65.cri().ra());
```

> **Note** — the WASM support is at an early stage.  Not all types and methods are
> exposed to JavaScript yet.  The TypeScript declarations in `colorimetry.d.ts`
> document the full current surface.

# Examples

<details>
<summary><strong>Calculate Tristimulus Values for Illuminants</strong></summary>

**XYZ tristimulus values** are the three numbers (X, Y, Z) that encode a color stimulus as seen
by a CIE standard observer. They are obtained by integrating the product of a spectrum with the
three color-matching functions x̄(λ), ȳ(λ), z̄(λ) of the observer. The Y value equals luminance
(or illuminance for an illuminant), and the x = X/(X+Y+Z), y = Y/(X+Y+Z) chromaticity
coordinates place the color in the familiar CIE 1931 diagram.

Different observers give slightly different numbers because their color-matching functions differ.
The **CIE 1931 2° observer** (based on a small 2° central foveal field) is the default used in
most product datasheets. The **CIE 2015 observers** are based on updated cone-fundamental
measurements and are more accurate, especially in the deep-blue region — important for LED
sources with strong short-wavelength components.

This example calculates the XYZ tristimulus values of the D65 illuminant for both the CIE 1931
2° standard observer and the CIE 2015 10° observer.

```
  use colorimetry::illuminant::D65;
# use approx::assert_abs_diff_eq as check;

  // D65 Tristimulus values, using the CIE 1931 standard observer by default.
  // set_illuminance(100.0) normalizes Y = 100, so X and Z scale accordingly.
  let xyz_d65 = D65.xyz(None).set_illuminance(100.0);

  let [x, y, z] = xyz_d65.to_array();
  // [95.04, 100.0, 108.86]  — the standard reference values for D65 / CIE 1931
# check!([x, y, z].as_ref(), [95.04, 100.0, 108.86].as_ref(),  epsilon = 5E-3);

  // D65 Tristimulus values using the CIE 2015 10° observer.
  // The 10° field is more appropriate when evaluating large colored surfaces
  // (walls, luminaires, architectural finishes viewed from a normal distance).
  use colorimetry::observer::Observer::Cie2015_10;
  let xyz_d65_10 = D65
    .xyz(Some(Cie2015_10)).set_illuminance(100.0);

  let [x_10, y_10, z_10] = xyz_d65_10.to_array();
  //[94.72, 100.0, 107.143]
# check!([x_10, y_10, z_10].as_ref(), [94.72, 100.0, 107.143].as_ref(), epsilon = 5E-3);
```

</details>

<details>
<summary><strong>Calculate Correlated Color Temperature and Tint</strong></summary>

**Correlated color temperature (CCT)** describes where a white light source sits on the
warm–cool axis, expressed in kelvin (K). It is defined as the temperature of the Planckian
(blackbody) radiator whose chromaticity is closest to that of the source, in the CIE 1960 UCS
u, v diagram (CIE 15:2018[^cie15] §9.4).

Typical CCT ranges used in practice:
- **2700–3000 K** — warm white (incandescent character), preferred for hospitality, residential
- **3500–4000 K** — neutral white, common in offices and retail
- **5000–6500 K** — cool white / daylight, used in task lighting, studios

The **Duv** value (Δuv) measures the perpendicular distance from the Planckian locus in the
u, v diagram (positive = above the locus, i.e. greenish tint; negative = below, i.e. pinkish/
magenta tint). ANSI C78.377 specifies a tolerance of ±0.006 for general-purpose lamps; tighter
specifications are used for professional and horticultural lighting.

This example calculates both CCT and Duv for the CIE standard illuminant A, which is the
reference for incandescent/halogen lamp color (a tungsten filament at ~2856 K):

```
# #[cfg(feature="cie-illuminants")]
  use colorimetry::illuminant::A;
# use approx::assert_abs_diff_eq as check;

  // Calculate CCT and Duv for the A illuminant.
  // Requires the `cct` and `cie-illuminants` features.
  // Result: ~2856 K on the Planckian locus (Duv = 0.0 — a blackbody by definition).
# #[cfg(all(feature="cct", feature="cie-illuminants"))]
  let [cct, duv] = A.cct().unwrap().to_array();
# #[cfg(all(feature="cct", feature="cie-illuminants"))]
# check!([cct, duv].as_ref(), [2855.4977, 0.0].as_ref(),  epsilon = 5E-4);
  // [2855.4977, 0.0]
```

</details>

<details>
<summary><strong>Calculate Color Fidelity Index for Illuminants</strong></summary>

**Color fidelity** answers the question: *"How faithfully does this light source reproduce the
appearance of object colors compared to an ideal reference?"* The reference is a blackbody
(Planckian) radiator for sources below 4000 K, or a CIE daylight (D-series) illuminant above
5000 K, with a smooth linear blend between 4000 K and 5000 K.

The **CIE 2017 Colour Fidelity Index** (**R<sub>f</sub>**, CIE 224:2017[^cie224]) is the modern scientific
measure. It evaluates 99 Color Evaluation Samples (CES) — a carefully chosen set of
surface reflectances spanning a wide range of hues, saturations, and object categories (skin
tones, foliage, food, textiles, and architectural materials). For each sample, the perceived color
under the test source is compared to the reference using the CIECAM02-UCS J′a′b′ perceptual
color space, which correctly accounts for chromatic adaptation (the eye's automatic white-balance).
The 99 color differences are averaged and converted to an **R<sub>f</sub>** score from 0 to 100:

| R<sub>f</sub> | Interpretation |
|---|---|
| 90–100 | Excellent fidelity — suitable for color-critical work (museums, print, surgery) |
| 80–89 | Good fidelity — high-end retail, hospitality, professional spaces |
| 70–79 | Acceptable — general commercial and office use |
| < 70 | Poor fidelity — industrial or utility use where color appearance is not critical |

The older **CRI** (R<sub>a</sub>, CIE 13.3-1995[^2]) uses only 8 pastel test colors and the simpler
CIE 1964 (U\*V\*W\*) space, which makes it less accurate for spectrally narrow sources like
LEDs and tri-phosphor fluorescent lamps — a source can score R<sub>a</sub> = 80 yet have strongly
distorted reds or greens that R<sub>f</sub> correctly penalizes.

This library implements the version harmonized with **ANSI/IES TM-30-20/24** (scaling constant CF = 6.73).

Below is an example calculation of the general Colour Fidelity Index for the CIE F2 illuminant
(a standard cool-white halophosphate fluorescent lamp, ≈ 4230 K):

```
# #[cfg(feature = "cie-illuminants")]
  use colorimetry::illuminant::F2;
# use approx::assert_abs_diff_eq as check;

# #[cfg(all(feature = "cfi", feature = "cie-illuminants"))]
# {
  // Calculate the Color Fidelity Index of the CIE F2 standard illuminant.
  // Requires the `cfi` and `cie-illuminants` features.
  let cf_f2 = F2.cfi().unwrap();
  let rf = cf_f2.color_fidelity_index();
  // ≈ 70 — F2's broad spectral gaps between mercury lines cause visible color distortion
  // for saturated reds and blues; acceptable for utility fluorescent lighting.
# check!(rf, 70.3,  epsilon = 1E-1);
# }
```

</details>

<details>

<summary><strong>Calculate the Spectral Locus for Chromaticity Diagrams</strong></summary>

The **spectral locus** is the horseshoe-shaped boundary of the CIE 1931 chromaticity diagram.
Each point on it corresponds to a pure monochromatic (single-wavelength) stimulus: 380 nm
(deep violet) at one end, 780 nm (deep red) at the other, with the straight-line closing segment
(the *line of purples*) connecting the two ends.

The interior of the horseshoe contains all colors that can be produced by mixtures of real light —
it is the boundary of the *object color solid* projected onto the chromaticity plane.
It appears on every LED and luminaire datasheet as the backdrop for color point tolerance
MacAdam ellipses and ANSI chromaticity bins.

Below, we compute the chromaticity coordinates that define the spectral locus, producing a
list of `[wavelength_nm, x, y]` triples suitable for plotting:

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
<summary><strong>Calculate XYZ/RGB Transformation Matrices for Color Profiles</strong></summary>

An **ICC color profile** encodes how to convert device-native RGB values (from a camera,
display, or printer) to and from the device-independent CIE XYZ space. The key ingredient is a
3×3 linear matrix whose columns are the XYZ coordinates of the device's red, green, and blue
primaries, scaled so the white point matches.

This library computes these matrices **spectrally** — it integrates the actual spectral power
distributions of the RGB primaries (defined by their chromaticity and the illuminant) with the
observer's color-matching functions, rather than using fixed tabulated XYZ values.
This matters when switching observers: the same display has slightly different matrix coefficients
for the CIE 1931 2° observer vs. the CIE 2015 10° observer, because their color-matching
functions differ.

Here, we compute forward and inverse matrices for the `DisplayP3` color space with both
`Cie1931` and `Cie2015` observers:

```
# use approx::assert_abs_diff_eq as check;
  use colorimetry::observer::Observer;
  use colorimetry::rgb::RgbSpace::DisplayP3;

  let xyz2rgb_31 = Observer::Cie1931.xyz2rgb_matrix(DisplayP3);
# let want31 = nalgebra::Matrix3::new(
#    2.4933, -0.9313, -0.4027,
#   -0.8298,  1.7629,  0.0236,
#    0.0355, -0.076,   0.9574
# );
#   check!(xyz2rgb_31, &want31, epsilon=5E-4);
  //  2.4933, -0.9313, -0.4027
  // -0.8298,  1.7629,  0.0236
  //  0.0355, -0.076,   0.9574

  let rgb2xyz_31 = Observer::Cie1931.rgb2xyz_matrix(DisplayP3);
#  let want31inv = nalgebra::Matrix3::new(
#     0.4866, 0.2656, 0.1981,
#     0.2291, 0.6917, 0.0792,
#     0.0001, 0.0451, 1.0433,
# );
# check!(rgb2xyz_31, &want31inv, epsilon=5E-4);
  // 0.4866, 0.2656, 0.1981
  // 0.2291, 0.6917, 0.0792
  // 0.0001, 0.0451, 1.0433

  use colorimetry::observer::Observer::Cie2015;

  let xyz2rgb_15 = Cie2015.xyz2rgb_matrix(DisplayP3).clone();
# let want15 = nalgebra::Matrix3::new(
#     2.5258,  -1.0009, -0.3649,
#    -0.9006,   1.8546, -0.0011,
#     0.0279,  -0.0574,  0.95874
# );
# check!(xyz2rgb_15, want15, epsilon=5E-4);
  //  2.5258,  -1.0009, -0.3649
  // -0.9006,   1.8546, -0.0011
  //  0.0279,  -0.0574,  0.95874
```

</details>

<details>
<summary><strong>Find the closest spectral match in the Munsell Color Book</strong></summary>

The **Munsell Color System** is a perceptually uniform atlas of surface colors organized along
three axes: *Hue* (color direction: R, YR, Y, GY, G, BG, B, PB, P, RP), *Value* (lightness, 0–10),
and *Chroma* (saturation distance from neutral grey, typically 0–20+).
With over 1,600 physical chips, it is widely used in architecture and interior design to specify
paint colors, and in colorimetry as a reference collection of real-world surface reflectances.

The notation **5R 4/14** means: hue 5 Red, Value 4 (mid-dark), Chroma 14 (vivid red).
A **ΔE** of 3 is at the threshold of visible difference for most observers under controlled
viewing — distances below ~1.5 are imperceptible, above ~5 are clearly different.

This example finds the best matching Munsell chip for the CRI R9 test color — a saturated
deep red — using the `CieCam16::de_ucs` color difference metric and the `Cie2015_10` observer,
which more accurately models perception of large colored surfaces:

```
# #[cfg(all(feature= "cri", feature = "munsell"))]
# {
# use approx::assert_abs_diff_eq as check;
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
# assert_eq!(key, "5R4/14");
# check!(delta_e, 2.85, epsilon = 5e-2);
  // ("5R4/14", 2.85)
# }
```

</details>

<details>
<summary><strong>Match a paint color to a <i>sRGB</i> display value under realistic viewing conditions</strong></summary>

Showing a realistic preview of a wall color on a screen is harder than it looks.
The challenge is that a display and a painted wall are fundamentally different stimuli — the
display emits light (additive RGB), while the paint absorbs and reflects it (subtractive).
Moreover, the observer's eye adapts differently to indoor tungsten and outdoor daylight, making
the same paint chip look warmer or cooler depending on the ambient light.

This example matches the Munsell paint chip *5 BG 5/8* — a teal/blue-green color —
to its nearest *sRGB* equivalent, mimicking real-world viewing conditions.

<span style="display: inline-block; width: 1em; height: 1em; background-color: rgb(0, 113, 138);
border-radius: 50%; vertical-align: middle; border: 1px solid #000;"></span>
<span>rgb(0, 113, 138)</span>

Instead of the traditional *CIE 1931 2°* observer, this match uses the *CIE 2015 10° observer*,
which more accurately reflects how large wall areas appear. The illumination uses the warm-white
*LED_B2* standard illuminant (≈ 3000 K), typical of residential and hospitality lighting.
The CIECAM16 color appearance model then accounts for chromatic adaptation, so the display
pixel is not a naïve XYZ conversion but a perceptually faithful rendering of what the eye
actually sees on a freshly painted surface.

```
# #[cfg(all(feature = "munsell", feature = "cie-illuminants"))]
# {
  // requires `cie-illuminants` and `munsell` features
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

  // Re-express the CIE 2015 pixel in CIE 1931 XYZ, then convert to sRGB.
  // sRGB is defined under the CIE 1931 observer, so this step is required for
  // correct encoding — a pixel rendered for CIE 2015 must be re-tagged for CIE 1931
  // before being sent to a display.
  let xyz_1931 = Cie1931.xyz(&rgb_2015, None);
  let rgb_1931 = xyz_1931.rgb(SRGB).compress();
  let [r, g, b]: [u8; 3] = rgb_1931.into();
  //  (0, 113, 138)
# assert!(r == 0 && g == 113 && b == 138);
# }
```

</details>

# Capabilities

## Spectral representation

All spectra in this library use a fixed grid from **380 to 780 nm** at **1 nm intervals** (401 samples),
matching the wavelength domain recommended by CIE 15:2018[^cie15] §7.2 for accurate numerical
integration. Spectral values are stored as [`f64`] components of a fixed-size [`nalgebra`] vector,
so integration against color-matching functions is a single dot product.

- [`Spectrum`] — the core spectral type (401 values, 380–780 nm)
  - [`Spectrum::linear_interpolate`] and [`Spectrum::sprague_interpolate`] — resample measured
    data from instrument-native intervals (5 nm, 10 nm, or irregular) onto the 1 nm grid.
    The Sprague method is recommended for data at 10 nm intervals (CIE 15:2018[^cie15] §7.2.3).
  - [`Spectrum::smooth`] — optional Gaussian smoothing to reduce measurement noise

## Analytical spectral models

Generate spectra from physical or empirical models, without measured data:

- [`Illuminant::planckian`] — Planck's law blackbody spectrum at any temperature. The spectrum is
  normalized to 1 W/m² total irradiance over the full range. Use this as the reference illuminant
  below 4000 K in color-rendering calculations, or to explore CCT-dependent color appearance.
- [`Illuminant::led`] — Ohno (2005) model for an LED chip's spectral peak, parameterized by
  center wavelength and full-width-at-half-maximum (FWHM). Useful for synthesizing LED
  spectra when measured SPD data are not available.
- [`Colorant::gaussian`] — Gaussian bandpass filter, useful as a synthetic test colorant or
  for modeling narrow-band interference filters.
- [`Stimulus::from_rgb`] — spectral reconstruction of an RGB pixel using Gaussian primaries,
  allowing a display color to be used anywhere a spectral stimulus is needed.

## CIE standard illuminants (optional, `cie-illuminants` feature)

The CIE defines standard illuminants (CIE 15:2018[^cie15] §4) to ensure reproducible,
observer-independent color specifications. The most commonly used are:

- **Daylight**: [`D65`] (6504 K, global average daylight — the default for sRGB, ICC profiles,
  and most colorimetric work), [`D50`] (5003 K, horizon light — used in print and graphic arts)
- **Incandescent**: [`A`] (2856 K — tungsten filament lamp, the CRI reference for warm-white sources)
- **Fluorescent** (F-series): [`F1`]–[`F12`] — standard fluorescent lamp types from the 1980s,
  used as test sources in lighting evaluation. F1–F6 are halophosphate types; F7–F9 are
  broad-band tri-phosphor; F10–F12 are narrow-band tri-phosphor (highest CRI among the set).
- **Extended fluorescent** (F3.x series): [`F3_1`]–[`F3_15`] — additional fluorescent SPDs used
  in color appearance research.
- **LED**: [`LED_B1`]–[`LED_B5`] (single blue-chip + phosphor), [`LED_BH1`] (blue-chip hybrid),
  [`LED_RGB1`] (red-green-blue multiband), [`LED_V1`] (violet-pumped) — the LED illuminants
  introduced in CIE 15:2018[^cie15] for testing color metrics on modern solid-state sources.

## Colorant collections (optional)

- [`MunsellCollection`] — the full Munsell renotation atlas: over 1,600 physical surface
  reflectances organized by hue, value, and chroma. Use for gamut visualization,
  perceptual uniformity tests, and paint-chip matching.
- [`CES`] — the 99 Color Evaluation Samples (CES) defined in CIE 224:2017[^cie224] Annex A.
  These are the test reflectances used in all R<sub>f</sub> and R<sub>g</sub> calculations.

## Light source quality metrics

These metrics quantify the photometric and colorimetric properties of a light source and are
the standard output required on luminaire datasheets, in specification documents, and in
regulatory compliance testing.

### Correlated color temperature (CCT and Duv)

[`CCT`] — requires `cct` feature. Computed under the **CIE 1931 observer** as specified by
CIE 15:2018[^cie15] §9.4, using a high-accuracy implementation based on Ohno (2014).

Returns both the **CCT** (in kelvin) and the **Duv** (signed distance from the Planckian locus
in the CIE 1960 UCS diagram):
- Positive Duv (> 0) → chromaticity is *above* the locus → slightly **greenish tint**
- Negative Duv (< 0) → chromaticity is *below* the locus → slightly **pinkish/magenta tint**
- |Duv| < 0.006 is the ANSI C78.377 tolerance for white-light sources[^1]

### Color Rendering Index (CRI)

[`CRI`] — requires `cri` feature. The legacy CIE 13.3-1995[^2] metric, still widely specified
in purchase documents and regulatory requirements. Provides the general index R<sub>a</sub>
(average over 8 pastel test colors) and the 14 special indices R<sub>1</sub>–R<sub>14</sub>
(including R<sub>9</sub>, the saturated red that legacy CRI is known to handle poorly).

Use CRI when backward compatibility with existing specifications is required. For new
work, prefer the CFI (below), which is more accurate for modern LED and multi-band sources.

### CIE 2017 Colour Fidelity Index (CFI / TM-30)

[`CFI`] — requires `cfi` feature. The current scientific standard for color rendering
evaluation, following CIE 224:2017[^cie224] and ANSI/IES TM-30-20/24.[^3]

- [`CFI::color_fidelity_index`] — overall **R<sub>f</sub>** (0–100): average fidelity across all
  99 CES under the test source vs. the reference illuminant. Higher is better.
- [`CFI::color_gamut_index`] — **R<sub>g</sub>** (ANSI/IES TM-30): area of the 16-bin averaged
  polygon in the CIECAM02-UCS a′b′ plane, relative to the reference polygon.
  R<sub>g</sub> = 100 → same gamut area as the reference;
  R<sub>g</sub> > 100 → colors appear more saturated on average (gamut expansion);
  R<sub>g</sub> < 100 → colors appear more muted (gamut compression).
- [`CFI::local_color_fidelity_indices`] — **R<sub>f,hj</sub>**: fidelity broken down by each of the 16 hue sectors (bins
  of 22.5° in the a′b′ plane, ordered from red through yellow, green, cyan, blue, purple, and
  back). Use these to diagnose *which* hue region is problematic.
- [`CFI::chroma_shift_indices`] — **R<sub>cs,hj</sub>**: fractional chroma shift per hue bin
  (positive = saturation boost, negative = desaturation).
- [`CFI::hue_shift_indices`] — **R<sub>hs,hj</sub>**: hue shift per bin in radians, wrapped to (−π, π].

## Color appearance models and color difference

Beyond tristimulus values, this library provides models that account for viewing context —
how bright the surround is, what the eye has adapted to, and the non-linear nature of human
contrast sensitivity. These are needed whenever color differences must be perceptually meaningful.

- [`CieLab`] — the CIE 1976 L\*a\*b\* space (CIE 15:2018[^cie15] §8.2.1). Approximately
  uniform for small color differences under a fixed illuminant and observer. Color differences are
  expressed as ΔE\*<sub>ab</sub> (via [`CieLab::ciede`]) or the more accurate
  ΔE\*<sub>00</sub> (via [`CieLab::ciede2000`], CIE 15:2018[^cie15] §8.3.1).
- [`CieCam02`] — CIECAM02 color appearance model (CIE 15:2018[^cie15] §10.1 / CIE 248).
  Models lightness (J), colorfulness (M), chroma (C), hue (H/h), and brightness (Q) as
  functions of the viewing conditions (adapting field luminance, background luminance,
  surround type). Used internally for CFI calculations. Color difference: [`CieCam02::de_ucs`].
- [`CieCam16`] — the successor to CIECAM02, with a simplified and more stable chromatic
  adaptation transform. Recommended for new work. Color difference: [`CieCam16::de_ucs`].

## RGB color spaces

Transforms between CIE XYZ and device-dependent RGB, computed spectrally for any observer.

- [`RgbSpace::SRGB`] — the standard for the web, photography, and most consumer displays
  (D65 white point, CIE 1931 observer, power-law gamma ≈ 2.2)
- [`RgbSpace::Adobe`] — Adobe RGB 1998, wider gamut than sRGB (same D65/CIE 1931 basis)
- [`RgbSpace::DisplayP3`] — Apple's Display P3, used on modern phones and monitors; wider
  red gamut than sRGB, D65 white point

All spaces support spectral-basis matrix computation so the matrices are consistent for any
observer, not just CIE 1931.

## CIE standard observers

Each [`Observer`] variant carries its full set of color-matching functions
(CIE 15:2018[^cie15] §3) and computes XYZ tristimulus values by numerical integration
against any spectrum.

- [`Observer::Cie1931`] — **CIE 1931 2° observer** (CIE 15:2018[^cie15] §3.1).
  The default for almost all colorimetric work: ICC profiles, sRGB, display characterization,
  CCT calculation. Based on a 2° central foveal field. Use for small samples viewed at
  normal reading distance (patches up to ~17 mm diameter at 0.5 m).
- [`Observer::Cie1964`] — **CIE 1964 10° observer** (CIE 15:2018[^cie15] §3.2).
  Based on a 10° field; better for large fields such as broad wall areas or immersive scenes.
  Required by CIE 224:2017[^cie224] for CFI calculations (§4.4).
- [`Observer::Cie2015`] — **CIE 2015 2° cone-fundamentals observer** (CIE 15:2018[^cie15] §3.3).
  More physiologically accurate, especially in the blue (short-wave) region. Use for high-
  accuracy spectral work or when the source has strong violet/deep-blue components.
- [`Observer::Cie2015_10`] — **CIE 2015 10° cone-fundamentals observer**.
  Best choice for predicting the appearance of large surfaces (walls, ceilings, artwork)
  under spectrally complex sources.

## Gamut analysis

Compute the limits of realizable color under a given observer:

- [`OptimalColors`] / [`RelXYZGamut`] — the *object color solid* (Rösch–MacAdam optimal
  colors): the maximum chroma achievable at each lightness level, computed by accumulating
  all possible binary reflectance functions. Used to calculate the gamut boundary for
  reflective surfaces and to normalize the R<sub>g</sub> computation.
- [`CieLChGamut`] — gamut boundary in the CIE L\*C\*h\* space; useful for hue-neutral
  gamut mapping and plotting the maximum chroma envelope on chromaticity diagrams.

[^cie15]: CIE 15:2018 *Colorimetry*, 4th edition. Commission Internationale de l'Éclairage. ISBN 978-3-902842-13-8. Available at <https://cie.co.at/publications/colorimetry-4th-edition>.
[^cie224]: CIE 224:2017 *Colour Fidelity Index for Accurate Scientific Use*. Commission Internationale de l'Éclairage. ISBN 978-3-902842-55-8. Available at <https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use>.
[^1]: McCamy, C. S. (1992). Correlated color temperature as an explicit function of chromaticity coordinates. *Color Research & Application*, 17(2), 142–144.
[^2]: CIE 13.3-1995 *Method of Measuring and Specifying Colour Rendering Properties of Light Sources*. R<sub>a</sub> is computed from 8 Munsell-based pastel test color samples in the CIE 1964 (U\*V\*W\*) uniform color space.
[^3]: CIE 224:2017[^cie224] defines R<sub>f</sub> from 99 CES using CIECAM02-UCS J′a′b′ color differences and a softplus formula (CF = 6.73). R<sub>g</sub> and the per-bin metrics R<sub>f,hj</sub>, R<sub>cs,hj</sub>, R<sub>hs,hj</sub> follow ANSI/IES TM-30-20/24.

# Standards

This library implements calculations defined in the following CIE and IES standards.
The documents are not bundled with the library — each developer must obtain their own copy from the issuing body.

| Standard | Topic | Purchase |
|---|---|---|
| **CIE 15:2018** | Colorimetry, 4th edition — tristimulus values, standard observers, chromaticity, CCT, CRI | <https://cie.co.at/publications/colorimetry-4th-edition> |
| **CIE 224:2017** | Colour Fidelity Index for accurate scientific use — R<sub>f</sub>, R<sub>g</sub>, CES, CIECAM02-UCS | <https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use> |
| **ANSI/IES TM-30-20/24** | IES method for evaluating light source color rendition — per-bin R<sub>f,hj</sub>, R<sub>cs,hj</sub>, R<sub>hs,hj</sub>, CVG | <https://www.ies.org/store/> |

# Features

- `cie-illuminants`
  Adds the full set of CIE standard illuminants: the A illuminant (incandescent), the F-series
  (fluorescent), the F3.x extended set, and the LED-series illuminants introduced in CIE 15:2018.
  D65 and D50 are always available without this feature.

- `munsell`
  Includes the reflectance spectra for the Munsell renotation color atlas (1,600+ chips).
  Enables [`MunsellCollection`] for nearest-neighbor color matching.

- `cct`
  Enables correlated color temperature (CCT) and Duv calculations for illuminants.
  Generates a 4 096-entry lookup table (three `f64` values per entry) at first use; the table
  is allocated at compile time but populated on demand.
  Automatically included when the `cri` feature is enabled.

- `cri`
  Enables Color Rendering Index (CRI / R<sub>a</sub>) calculations following CIE 13.3-1995.
  Provides R<sub>a</sub> and the 14 special indices R<sub>1</sub>–R<sub>14</sub>.
  Loads 14 Munsell-based test color sample spectra.

- `cfi`
  Enables Color Fidelity Index (CFI / R<sub>f</sub>) calculations following CIE 224:2017 / ANSI/IES TM-30.
  Provides the overall R<sub>f</sub>, per-sample R<sub>f,1</sub>–R<sub>f,99</sub>, gamut index R<sub>g</sub>,
  and the 16-bin metrics R<sub>f,hj</sub> / R<sub>cs,hj</sub> / R<sub>hs,hj</sub>.
  Loads the 99 CES reflectance spectra at compile time.

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

# Command Line Tool

This library has an associated command-line tool, named `color`, contained in the `colorimetry-cli` crate.
It provides a convenient way to perform colorimetric calculations, convert spectral data, and make color plots directly from the terminal.

To use it you have to install it using `cargo`, which in turn requires that you have Rust and Cargo installed.
Install them first by following the instructions at <https://www.rust-lang.org/tools/install>.
After having done this, run the following command:

```bash
cargo install colorimetry-cli
```

That should be it — you can now use the `color` command in your terminal.
Run the following command to see the available options:

```bash
color --help
```

# Color Plots

The `colorimetry` library includes a plotting module, in an associated `colorimetry-plot` crate
that can be used to generate chromaticity diagrams, spectral plots, and color rendering
visualizations. Output is in SVG format, viewable in any modern web browser or vector
graphics editor such as Inkscape.

```bash
cargo add colorimetry-plot
```

# Developer Tasks with `xtask`

This project uses a Rust-based `xtask` utility for common development tasks:

- `cargo xtask check` – run clippy, fmt, build check, and README sync verification
- `cargo xtask test` – run tests across all feature configurations
- `cargo xtask doc` – build rustdoc and fail on any warnings
- `cargo xtask wasm` – compile WebAssembly output in `pkg/` (requires `wasm-pack` and `wasm-opt`)

# License

All content &copy;2026 Harbers Bik LLC, and licensed under either of the

- Apache License, Version 2.0,
  ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>), or the
- MIT license
  [LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>,

at your option.

# Contribution

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
[`CFI::color_fidelity_index`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.color_fidelity_index
[`CFI::color_gamut_index`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.color_gamut_index
[`CFI::local_color_fidelity_indices`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.local_color_fidelity_indices
[`CFI::chroma_shift_indices`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.chroma_shift_indices
[`CFI::hue_shift_indices`]: https://docs.rs/colorimetry/latest/colorimetry/illuminant/struct.CFI.html#method.hue_shift_indices
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
pub mod rgb;
pub mod spectrum;
pub mod stimulus;
pub mod traits;
pub mod xyz;

pub use error::Error;
