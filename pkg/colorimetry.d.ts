/* tslint:disable */
/* eslint-disable */
export enum CieIlluminant {
  D65 = 0,
  D50 = 1,
}
/**
 * Light-weight identifier added to the `XYZ` and `RGB` datasets,
 *    representing the colorimetric standard observer used.
 *
 *    No data included here, which would be the Rust way, but that does not work with wasm-bindgen.
 *    This can be directly used in JavaScript, and has the benefit to be just an index.
 */
export enum Observer {
  Cie1931 = 0,
}
/**
 * A Light Weight tag, representing an RGB color space.
 * Used for example in the RGB value set, to identify the color space being used.
 */
export enum RgbSpace {
  SRGB = 0,
  Adobe = 1,
  DisplayP3 = 2,
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
 * assert_abs_diff_eq!(rgb.values().as_ref(), [0.5, 0.25, 0.75].as_ref(), epsilon = 1e-6);
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
   * const [x, y, z] = xyz.values();
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
  readonly illuminant_new_js: (a: number, b: number) => [number, number, number];
  readonly illuminant_Values: (a: number) => [number, number];
  readonly illuminant_illuminant: (a: number) => number;
  readonly spectrum_new_js: (a: number, b: number) => [number, number, number];
  readonly spectrum_Values: (a: number) => [number, number];
  readonly spectrum_linearInterpolate: (a: number, b: number, c: number, d: number) => [number, number, number];
  readonly xyz_new_js: (a: number, b: number, c: any) => [number, number, number];
  readonly xyz_values: (a: number) => any;
  readonly xyz_chromaticity: (a: number) => any;
  readonly xyz_y: (a: number) => number;
  readonly __wbg_widergb_free: (a: number, b: number) => void;
  readonly __wbg_relxyz_free: (a: number, b: number) => void;
  readonly __wbg_chromaticity_free: (a: number, b: number) => void;
  readonly __wbg_spectrum_free: (a: number, b: number) => void;
  readonly __wbg_rgb_free: (a: number, b: number) => void;
  readonly __wbg_illuminant_free: (a: number, b: number) => void;
  readonly __wbg_xyz_free: (a: number, b: number) => void;
  readonly __wbg_cielab_free: (a: number, b: number) => void;
  readonly __wbg_viewconditions_free: (a: number, b: number) => void;
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
