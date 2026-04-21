/* tslint:disable */
/* eslint-disable */
/**
 * A **lightweight enum** representing the CIE standard illuminants from the CIE 15:2018 datasets
 * (downloaded August 2024). Each variant holds a zero-cost reference to its precompiled spectrum,
 * making it easy to include as a field in your own types without pulling in heavy data structures.
 *
 * This enum implements `IntoEnumIterator`, so you can **iterate through every standard illuminant**
 * (useful for testing, batch conversions, or validation).
 *
 * - Use `CieIlluminant::iter()` or `CieIlluminant::spectrum()` to list or retrieve any built-in illuminant.
 * - For a generic D-series illuminant at any correlated color temperature, use
 *   `Spectrum::cie_d_illuminant(cct: f64)`.
 *
 * By default, only **D65** and **D50** are included. To pull in the full set of fluorescent “F3_X”
 * series and other CIE illuminants, enable the `"cie-illuminants"` feature in `Cargo.toml`
 * (or build with `--features cie-illuminants`). Omit that feature (or use `--no-default-features`)
 * to keep your binary lean.
 *
 * In JavaScript/WebAssembly builds, the `colorimetry` package excludes these extra spectra by default
 * for faster load times. To include them, use the `colorimetry-all` bundle instead.
 *
 * For more background, see the Wikipedia article on
 * [Standard illuminant white points](https://en.wikipedia.org/wiki/Standard_illuminant#White_points_of_standard_illuminants).
 *
 * # Examples
 * ```rust
 * use colorimetry::illuminant::CieIlluminant;
 * use strum::IntoEnumIterator;
 *
 * // Iterate through and print all available CIE illuminants:
 * for illum in CieIlluminant::iter() {
 *     println!("{illum}");
 * }
 * ```
 */
export enum CieIlluminant {
  D65 = 0,
  D50 = 1,
  E = 2,
  A = 3,
  F1 = 4,
  F2 = 5,
  F3 = 6,
  F4 = 7,
  F5 = 8,
  F6 = 9,
  F7 = 10,
  F8 = 11,
  F9 = 12,
  F10 = 13,
  F11 = 14,
  F12 = 15,
  F3_1 = 16,
  F3_2 = 17,
  F3_3 = 18,
  F3_4 = 19,
  F3_5 = 20,
  F3_6 = 21,
  F3_7 = 22,
  F3_8 = 23,
  F3_9 = 24,
  F3_10 = 25,
  F3_11 = 26,
  F3_12 = 27,
  F3_13 = 28,
  F3_14 = 29,
  F3_15 = 30,
  LED_B1 = 31,
  LED_B2 = 32,
  LED_B3 = 33,
  LED_B4 = 34,
  LED_B5 = 35,
  LED_BH1 = 36,
  LED_RGB1 = 37,
  LED_V1 = 38,
  LED_V2 = 39,
}
/**
 * Selects a CIE standard colorimetric observer.
 *
 * The tag is embedded in every [`XYZ`] and [`Rgb`](crate::rgb::Rgb) value so that
 * operations across incompatible observers can be detected at runtime.  Each variant
 * is a lightweight index; the color-matching function tables are stored in
 * [`observer_data`].
 */
export enum Observer {
  /**
   * CIE 1931 2° standard observer — the default for most colorimetry.
   *
   * Used by sRGB, ICC profiles, and CIE CRI Ra.
   */
  Cie1931 = 0,
  /**
   * CIE 1964 10° supplementary standard observer.
   *
   * Preferred when the viewed area subtends more than ~4° at the eye.
   * Used by CIE 224:2017 / ANSI/IES TM-30 for colour fidelity calculations.
   */
  Cie1964 = 1,
  /**
   * CIE 2015 2° observer — CMFs constructed as linear transforms of the Stockman & Sharpe
   * (2000) cone fundamentals.
   *
   * More accurate than `Cie1931` in the short-wavelength (blue) region.
   */
  Cie2015 = 2,
  /**
   * CIE 2015 10° observer — CMFs constructed as linear transforms of the Stockman & Sharpe
   * (2000) cone fundamentals.
   *
   * Wide-field counterpart of [`Cie2015`](Observer::Cie2015).
   */
  Cie2015_10 = 3,
}
/**
 * Spectrally based color space, using spectral representations of the primaries and the
 * reference white.
 *
 * Using the CIE 1931 standard observer, using a wavelength domain from 380 top 780
 * nanometer with 1 nanometer steps, these result in their usual chromaticity
 * values.  The most common _sRGB_ color space is obtained using the
 * `RgbSpace::srgb()` constructor. For this instance, the blue and green primaries
 * are direct Gaussian-filtered D65 spectra. A mixture of the blue primary and a
 * G1aussian-filtered red component is used for the red primary. Similar
 * constructors are provided for the `Adobe` and `DisplayP3` color spaces.
 *
 * The benefit of spectral primaries is that color management and color profiles
 * can use updated Colorimetric Observers, such as the Cone-Fundamental based CIE
 * 2015 observers, which don't have the CIE 1931 deficiencies. For example, they
 * can also be optimized for special observers by considering an observer's age or
 * health conditions.
 */
export enum RgbSpace {
  SRGB = 0,
  Adobe = 1,
  DisplayP3 = 2,
  CieRGB = 3,
}
/**
 * Container for CIE 2017 Colour Fidelity Index (**R<sub>f</sub>**) calculations,
 * including both the general color fidelity **R<sub>f</sub>** score and the 99 special color fidelity indices (**R<sub>f,1</sub>** to **R<sub>f,99</sub>**)
 * as specified in [CIE 224:2017](https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use).
 *
 * # Requirements
 * - Requires the `cfi` feature to access color evaluation samples (CES) used for testing.
 *
 * # Overview
 * The CIE 2017 Colour Fidelity Index (CFI, or **R<sub>f</sub>**) is a modern metric for evaluating how accurately a light source renders colors.
 * It uses 99 Color Evaluation Samples (CES) that cover a broad range of real-world colors, providing a much more comprehensive assessment
 * than older metrics like the Color Rendering Index (CRI, or **R<sub>a</sub>**).
 * - The general index (**R<sub>f</sub>**) gives an overall measure of color fidelity.
 * - The special indices (**R<sub>f,1</sub>** to **R<sub>f,99</sub>**) show fidelity for each specific color sample.
 *
 * # Comparison with CRI
 * The traditional **CRI** metric (Ra) uses only 8 or 14 pastel color samples and is known to be limited, especially for modern light sources such as LEDs.
 * **CFI (Rf)** is a newer, more robust standard: it uses a much wider set of samples and is based on state-of-the-art color appearance models,
 * providing a more accurate and reliable prediction of real-world color rendering.
 * - **Use CRI** if you need compatibility with legacy systems or must comply with standards that specify CRI.
 * - **Use CFI (Rf)** for a more precise and scientifically up-to-date assessment of color fidelity, especially with modern or tunable light sources.
 *
 * # TM-30 version
 * This implementation follows **ANSI/IES TM-30-20** and **TM-30-24**, which are harmonised with
 * **CIE 224:2017**. It does **not** implement the earlier TM-30-15 or TM-30-18 editions.
 * The key difference from TM-30-15 is the scaling constant in the Rf formula: TM-30-15 used
 * `CF = 7.54`; all later editions (TM-30-18, TM-30-20, TM-30-24, CIE 224:2017) use `CF = 6.73`,
 * which this library implements.
 *
 * # Reference
 * [CIE 224:2017 – CIE 2017 Colour Fidelity Index for accurate scientific use](https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use)
 */
export class CFI {
  private constructor();
  free(): void;
  /**
   * General colour fidelity index Rf (0–100).
   *
   * A single overall score measuring how faithfully the test source renders the 99
   * Colour Evaluation Samples compared to the reference illuminant.
   */
  colorFidelityIndex(): number;
  /**
   * General colour gamut index Rg.
   *
   * Measures the area of the gamut polygon relative to the reference (100 = same area).
   * Values above 100 indicate a wider gamut than the reference; below 100 means narrower.
   */
  colorGamutIndex(): number;
  /**
   * Local colour fidelity index Rf,hj for each of the 16 hue bins (TM-30 / CIE 224:2017 §4.5).
   *
   * Returns a `Float64Array` of 16 values.  Bin 0 starts at 0° (red), progressing
   * counter-clockwise in 22.5° steps around the hue circle.
   */
  localColorFidelityIndices(): Float64Array;
  /**
   * Chroma shift index Rcs,hj for each of the 16 hue bins (TM-30 / CIE 224:2017 §4.6).
   *
   * Returns a `Float64Array` of 16 values.  Positive means the test source boosts
   * saturation in that hue direction; negative means desaturation.  Typical range ≈ −0.5…+0.5.
   */
  chromaShiftIndices(): Float64Array;
  /**
   * Hue shift index Rhs,hj for each of the 16 hue bins (TM-30 / CIE 224:2017 §4.7).
   *
   * Returns a `Float64Array` of 16 values in radians, wrapped to (−π, π].
   * Positive means a counter-clockwise hue shift; negative means clockwise.
   */
  hueShiftIndices(): Float64Array;
  /**
   * Special colour fidelity indices Rf,i for all 99 CES (CIE 224:2017 §7).
   *
   * Returns a `Float64Array` of 99 values, one per Colour Evaluation Sample.
   */
  specialColorFidelityIndices(): Float64Array;
}
/**
 * The **Color Rendering Index (CRI)** for a light source, computed according to CIE 13.3-1995.
 *
 * This struct holds the 14 individual rendering indices R₁…R₁₄ for the standard test color samples,
 * and provides the general CRI, Rₐ, which is the average of the first eight Rᵢ values.  
 *
 * # Calculation Method
 * 1. The test illuminant is scaled to 100 lx and converted to CIE XYZ under the CIE 1931 observer.  
 * 2. Each of the 14 standard Colorant test spectra (TCS) is measured under both the test and the  
 *    reference illuminant (black-body or D-series at the test’s correlated color temperature).  
 * 3. For each sample, the color difference ΔE in CIE UVW space is computed, and  
 *    Rᵢ = 100 − 4.6 · ΔE.  
 * 4. The general CRI Rₐ is then  
 *    ```text
 *    Rₐ = (R₁ + R₂ + … + R₈) / 8
 *    ```
 *
 * # Examples
 * ```rust
 * use colorimetry::illuminant::{Illuminant, CRI};
 *
 * // Compute CRI for the D65 illuminant:
 * let cri: CRI = (&Illuminant::d65()).try_into().unwrap();
 *
 * // General CRI:
 * let ra = cri.ra();
 * println!("General CRI Rₐ = {:.1}", ra);
 *
 * // All 14 individual Rᵢ values:
 * let ri_values = cri.to_array();
 * for (i, &ri) in ri_values.iter().enumerate() {
 *     println!("R{} = {:.1}", i + 1, ri);
 * }
 * ```
 *
 * # Notes
 * - This implementation uses the **CIE 1931** color space and requires the `"cri"` feature to be enabled in the crate.  
 * - The CRI-metric is now considered somewhat outdated; newer metrics (e.g., TM-30) are recommended for modern lighting applications.  
 *   However, CRI remains widely used and understood across the lighting industry.
 *
 * # Errors
 * Constructing a `CRI` can fail if the illuminant’s correlated color temperature is out of the
 * valid range (1000–25000 K) or its distance from the Planckian locus exceeds 0.05 Δuv.
 */
export class CRI {
  private constructor();
  free(): void;
  /**
   * Returns the general colour rendering index Rₐ (0–100),
   * the average of the first eight special rendering indices R₁…R₈.
   */
  ra(): number;
}
/**
 * A chromaticity coordinate with x and y values.
 */
export class Chromaticity {
  private constructor();
  free(): void;
}
export class CieLab {
  private constructor();
  free(): void;
}
/**
 * # Illuminant
 *
 * An illuminant is a spectral power distribution that represents the
 * spectral power density of a light source (sun, bulb, LED, etc.) in
 * W/m²/nm over 380–780 nm (401 samples).
 */
export class Illuminant {
  free(): void;
  /**
   * Create a new illuminant spectrum from the given data.
   *
   * The data must be the 401 values from 380 to 780 nm, with an interval size of 1 nanometer.
   */
  constructor(data: Float64Array);
  /**
   * Returns the spectral data values, as a Float64Array containing 401 data
   * points, over a wavelength domain from 380 t0 780 nanometer, with a
   * stepsize of 1 nanometer.
   */
  Values(): Float64Array;
  /**
   * Calculates the Color Rendering Index for this illuminant spectrum.
   *
   * Returns a `CRI` object with a `ra()` method that gives the general colour
   * rendering index Rₐ (average of R₁…R₈, scaled 0–100).
   */
  cri(): CRI;
  /**
   * Calculates the Colour Fidelity Index (CFI / Rf) for this illuminant spectrum.
   *
   * Returns a `CFI` object exposing `rf()`, `rg()`, `rfHj()`, `rcsHj()`, `rhsHj()`,
   * and `specialIndices()`.  Requires the `cfi` feature.
   */
  cfi(): CFI;
  /**
   * Get the CieIlluminant spectrum. Typically you don't need to use the Spectrum itself, as many
   * methods just accept the CieIlluminant directly.
   */
  static illuminant(stdill: CieIlluminant): Illuminant;
}
/**
 * # Related Tristimulus Values
 *
 * Tristimulus Values for a given sample and reference white,
 * used to represent related colors as used in various color
 * models. Typically the reference white is normalized to have
 * an Y-value of 100
 */
export class RelXYZ {
  private constructor();
  free(): void;
}
/**
 * Represents a color stimulus using Red, Green, and Blue (RGB) values constrained to the `[0.0, 1.0]` range.
 * Each component is a floating-point value representing the relative intensity of the respective primary color
 * within a defined RGB color space.
 *
 * Unlike the CIE XYZ tristimulus values, which use imaginary primaries, RGB values are defined using real primaries
 * based on a specific color space. These primaries typically form a triangular area within a CIE (x,y) chromaticity
 * diagram, representing the gamut of colors the device can reproduce.
 *
 * # Usage
 * The `Rgb` struct is used to encapsulate color information in a device-independent manner, allowing for accurate color
 * representation, conversion, and manipulation within defined RGB spaces. It is particularly useful for applications
 * involving color management, digital imaging, and rendering where strict adherence to gamut boundaries is required.
 *
 * # Example
 * ```rust
 * # use colorimetry::rgb::Rgb;
 * # use approx::assert_abs_diff_eq;
 *
 * // Create an sRGB color with normalized RGB values
 * let rgb = Rgb::new(0.5, 0.25, 0.75, None, None).unwrap();
 * assert_abs_diff_eq!(rgb.to_array().as_ref(), [0.5, 0.25, 0.75].as_ref(), epsilon = 1e-6);
 * ```
 *
 * # Notes
 * - The `Rgb` struct strictly enforces the `[0.0, 1.0]` range for each component. Any attempt to create values
 *   outside this range will result in an error.
 * - The `observer` field allows for color conversion accuracy under different lighting and viewing conditions,
 *   enhancing the reliability of transformations to other color spaces such as XYZ.
 */
export class Rgb {
  private constructor();
  free(): void;
}
/**
 *
 * This container holds spectral values within a wavelength domain ranging from 380
 * to 780 nanometers, with an interval size of 1 nanometer and a total of 401
 * values. It also includes a category tag and an optional 'total' value for the
 * aggregate value associated with the spectrum.
 *
 * A `Spectrum` can be constructed from data, but many other construction methods
 * are available in this library, such as standard illuminants A and D65, Planckian
 * (Black Body) illuminants, or a `Stimulus` spectrum for a pixel of an sRGB
 * display.
 * 
 */
export class Spectrum {
  free(): void;
  /**
   * Create a new spectrum from the given data.
   *
   * The data must be the 401 values from 380 to 780 nm, with an interval size of 1 nanometer.
   *
   * If the Spectral data you have uses another wavelength domain and/or a different
   * wavelength interval, use the linear interpolate constructor,
   * which takes a wavelength domain and spectral data as arguments.
   */
  constructor(data: Float64Array);
  /**
   * Returns the spectral data values, as a Float64Array containing 401 data
   * points, over a wavelength domain from 380 t0 780 nanometer, with a
   * stepsize of 1 nanometer.
   */
  Values(): Float64Array;
  /**
   * This function maps spectral data with irregular intervals or intervals different than 1
   * nanometer to the standard spectrum as used in this library.
   *
   * For domains with a regular interval, the wavelength slice should have a size of two, containing
   * the minimum and maximum wavelength values, both also in units of meters or nanometers.
   *
   * For irregular domains, this function requires a slice of wavelengths and a slice of spectral
   * data, both of the same size. The wavelengths can be specified in units of meters or nanometers.
   *
   * In case of duplicate wavelength values the last data values is used, so it is impossible to
   * define filters with vertical edges using this method.
   *
   * ```rust
   * // Creates a linear gradient filter, with a zero transmission at 380 nanometer, and full
   * // transmission at 780 nanometer. This is an example using a uniform wavelength domain as input.
   * use colorimetry::prelude::*;
   * use approx::assert_ulps_eq;
   * let data = [0.0, 1.0];
   * let wl = [380.0, 780.0];
   * let mut spd = Spectrum::linear_interpolate(&wl, &data).unwrap();
   * assert_ulps_eq!(spd[380], 0.);
   * assert_ulps_eq!(spd[380+100], 0.25);
   * assert_ulps_eq!(spd[380+200], 0.5);
   * assert_ulps_eq!(spd[380+300], 0.75);
   * assert_ulps_eq!(spd[380+400], 1.0);
   *
   * // Creates a top hat filter, with slanted angles, using an irregular
   * // wavelength domain.
   * let data = vec![0.0, 1.0, 1.0, 0.0];
   * let wl = vec![480.0, 490.0, 570.0, 580.0];
   * let spd = Spectrum::linear_interpolate(&wl, &data).unwrap();
   * assert_ulps_eq!(spd[380+0], 0.0);
   * assert_ulps_eq!(spd[380+100], 0.0);
   * assert_ulps_eq!(spd[380+110], 1.0);
   * assert_ulps_eq!(spd[380+190], 1.0);
   * assert_ulps_eq!(spd[380+200], 0.0);
   * assert_ulps_eq!(spd[380+300], 0.0);
   * assert_ulps_eq!(spd[380+400], 0.0);
   * ```
   */
  static linearInterpolate(wavelengths: Float64Array, data: Float64Array): Spectrum;
}
/**
 * CIECAM viewing conditions.
 *
 * The ViewConditions as recommended by CIE248:2022 are provided for various scenarios as constants, and are included as:
 * - [`CIE248_CABINET`] Viewing a surface in a cabinet
 * - [`CIE248_HOME_SCREEN`] Viewing a self-luminous display at home
 * - [`CIE248_PROJECTED_DARK`] Viewing projected images in a darkened room
 * - [`CIE248_OFFICE_SCREEN`] Viewing a self-luminous display under office illumination
 *  
 * The TM30 and Color Fidelity ViewConditions are provided as [`TM30VC`].
 */
export class ViewConditions {
  private constructor();
  free(): void;
}
/**
 * Represents a color stimulus using unconstrained Red, Green, and Blue (RGB) floating-point values
 * within a device's RGB color space. The values can extend beyond the typical 0.0 to 1.0 range,
 * allowing for out-of-gamut colors that cannot be accurately represented by the device.
 */
export class WideRgb {
  private constructor();
  free(): void;
}
/**
 * Represents a color by its tristimulus value XYZ color space.
 *
 * The `XYZ` struct represents the tristimulus values (X, Y, Z) and the associated observer.
 * The observer defines the color matching functions used for the conversion.
 */
export class XYZ {
  free(): void;
  /**
   * Create an XYZ Tristimuls Values object.
   *
   * Accepts as arguments
   *
   * - x and y chromaticity coordinates only , using the "Cie::Cie1931" observer as default
   * - x and y chromaticity coordinates, and standard observer ID as 3rd argument
   * - X, Y, and Z tristimulus values, using the "Cie::Cie1931" observer as default
   * - X, Y, and Z tristimulus values, and a standard Observer ID as 4th argument
   *
   * When only x and y chromaticity coordinates are specified, the luminous
   * value is set to 100.0 candela per square meter.
   *
   * ```javascript, ignore
   * // Create a new XYZ object using D65 CIE 1931 chromaticity coordinates
   * const xyz = new cmt.XYZ(0.31272, 0.32903);
   *
   * // Get and check the corresponding tristimulus values, with a luminous value
   * // of 100.0
   * const [x, y, z] = xyz.to_array();
   * assert.assertAlmostEquals(x, 95.047, 5E-3); // D65 wikipedia
   * assert.assertAlmostEquals(y, 100.0);
   * assert.assertAlmostEquals(z, 108.883, 5E-3);
   *
   * // and get back the orgiinal chromaticity coordinates:
   * const [xc, yc] = xyz.chromaticity();
   * assert.assertAlmostEquals(xc, 0.31272);
   * assert.assertAlmostEquals(yc, 0.32903);
   *
   * // to get the luminous value:
   * const l = xyz.luminousValue();
   * assert.assertAlmostEquals(l, 100.0);
   * // D65 CIE 1931 chromaticity coordinates
   * const xyz = new cmt.XYZ(0.31272, 0.32903);
   * ```
   */
  constructor(x: number, y: number, ...opt: Array<any>);
  /**
   * Get the XYZ tristimulus value as an array.
   */
  values(): Array<any>;
  /**
   * Get the chromaticity coordinates
   */
  chromaticity(): Array<any>;
  /**
   * Get the luminous value, Y.
   */
  y(): number;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
  readonly memory: WebAssembly.Memory;
  readonly spectrum_new_js: (a: number, b: number) => [number, number, number];
  readonly spectrum_Values: (a: number) => [number, number];
  readonly spectrum_linearInterpolate: (a: number, b: number, c: number, d: number) => [number, number, number];
  readonly __wbg_cfi_free: (a: number, b: number) => void;
  readonly __wbg_widergb_free: (a: number, b: number) => void;
  readonly xyz_new_js: (a: number, b: number, c: any) => [number, number, number];
  readonly xyz_values: (a: number) => any;
  readonly xyz_chromaticity: (a: number) => any;
  readonly xyz_y: (a: number) => number;
  readonly __wbg_xyz_free: (a: number, b: number) => void;
  readonly __wbg_cri_free: (a: number, b: number) => void;
  readonly illuminant_new_js: (a: number, b: number) => [number, number, number];
  readonly illuminant_Values: (a: number) => [number, number];
  readonly illuminant_cri: (a: number) => [number, number, number];
  readonly illuminant_cfi: (a: number) => [number, number, number];
  readonly illuminant_illuminant: (a: number) => number;
  readonly cri_ra: (a: number) => number;
  readonly cfi_colorFidelityIndex: (a: number) => number;
  readonly cfi_colorGamutIndex: (a: number) => number;
  readonly cfi_localColorFidelityIndices: (a: number) => [number, number];
  readonly cfi_chromaShiftIndices: (a: number) => [number, number];
  readonly cfi_hueShiftIndices: (a: number) => [number, number];
  readonly cfi_specialColorFidelityIndices: (a: number) => [number, number];
  readonly __wbg_spectrum_free: (a: number, b: number) => void;
  readonly __wbg_viewconditions_free: (a: number, b: number) => void;
  readonly __wbg_cielab_free: (a: number, b: number) => void;
  readonly __wbg_chromaticity_free: (a: number, b: number) => void;
  readonly __wbg_rgb_free: (a: number, b: number) => void;
  readonly __wbg_relxyz_free: (a: number, b: number) => void;
  readonly __wbg_illuminant_free: (a: number, b: number) => void;
  readonly __wbindgen_export_0: WebAssembly.Table;
  readonly __wbindgen_malloc: (a: number, b: number) => number;
  readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
  readonly __externref_table_dealloc: (a: number) => void;
  readonly __wbindgen_free: (a: number, b: number, c: number) => void;
  readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;
/**
* Instantiates the given `module`, which can either be bytes or
* a precompiled `WebAssembly.Module`.
*
* @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
*
* @returns {InitOutput}
*/
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
* If `module_or_path` is {RequestInfo} or {URL}, makes a request and
* for everything else, calls `WebAssembly.instantiate` directly.
*
* @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
*
* @returns {Promise<InitOutput>}
*/
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
