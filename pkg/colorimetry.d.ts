/* tslint:disable */
/* eslint-disable */
/**
* Stefan Boltzmann law: Blackbody's radiant emittance (W m<sup>-2</sup>), as function of its absolute
* temperature (K).
* @param {number} temperature
* @returns {number}
*/
export function stefanBoltzmann(temperature: number): number;
/**
*/
export enum StdIlluminant {
  D65 = 0,
  D50 = 1,
  A = 2,
  F1 = 3,
  F2 = 4,
  F3 = 5,
  F4 = 6,
  F5 = 7,
  F6 = 8,
  F7 = 9,
  F8 = 10,
  F9 = 11,
  F10 = 12,
  F11 = 13,
  F12 = 14,
  F3_1 = 15,
  F3_2 = 16,
  F3_3 = 17,
  F3_4 = 18,
  F3_5 = 19,
  F3_6 = 20,
  F3_7 = 21,
  F3_8 = 22,
  F3_9 = 23,
  F3_10 = 24,
  F3_11 = 25,
  F3_12 = 26,
  F3_13 = 27,
  F3_14 = 28,
  F3_15 = 29,
  LED_B1 = 30,
  LED_B2 = 31,
  LED_B3 = 32,
  LED_B4 = 33,
  LED_B5 = 34,
  LED_BH1 = 35,
  LED_RGB1 = 36,
  LED_V1 = 37,
  LED_V2 = 38,
}
/**
*
*    Light-weight identifier added to the `XYZ` and `RGB` datasets,
*    representing the colorimetric standard observer used.
*
*    No data included here, which would be the Rust way, to maintain
*    compatibility with wasm-bindgen, and to allow this enum to be directly used
*    in JavaScript.
* 
*/
export enum Observer {
  Std1931 = 0,
  Std1976 = 1,
  Std2015 = 2,
  Std2015_10 = 3,
}
/**
*
*A Light Weight tag, representing an RGB color space.
*Used for example in the RGB value set, to identify the color space being used.  
* 
*/
export enum RgbSpace {
  SRGB = 0,
  ADOBE = 1,
  DisplayP3 = 2,
}
/**
* Encapcsulated Array of calculated Ri values, from a test light source.
*/
export class CRI {
  free(): void;
}
/**
*/
export class CieLab {
  free(): void;
}
/**
*/
export class Illuminant {
  free(): void;
/**
* Creates a new Spectrum object, using as input a `Category`, a
* Float64Array with exactly 401 datapoints, and an optional third
* parameter called total, representing the total irradiance, transmission,
* or reflectivity of the values, depending on the category of the
* spectrum. The spectral values should be associated with a wavelength
* domain from 380 to 480 nanometer, with an interval size of 1 nanometer.
*
* If the Spectral data you have uses another wavelength domain and/or a different
* wavelength interval, use the linear or sprague interpolate constructors,
* which takes a wavelength domain and spectral data as arguments.
* @param {Float64Array} data
*/
  constructor(data: Float64Array);
/**
* Returns the spectral data values, as a Float64Array containing 401 data
* points, over a wavelength domain from 380 t0 780 nanometer, with a
* stepsize of 1 nanometer.
* @returns {Float64Array}
*/
  Values(): Float64Array;
/**
* Calculates the Color Rendering Index values for illuminant spectrum.
* 
* To use this function, first use `await CRI.init()`, which downloads the
* Test Color Samples required for the calculation.  These are downloaded
* seperately to limit the size of the main web assembly library.
* @returns {CRI}
*/
  cri(): CRI;
/**
* Get the StdIlluminant spectrum. Typically you don't need to use the Spectrum itself, as many
* methods just accept the StdIlluminant directly.
* @param {StdIlluminant} stdill
* @returns {Illuminant}
*/
  static illuminant(stdill: StdIlluminant): Illuminant;
}
/**
*/
export class MunsellMatt {
  free(): void;
}
/**
*
*    A data structure to define Standard Observers, such as the CIE 1931 2º and
*    the CIE 2015 standard observers.
*    
*    These are defined in the form of the three color matching functions,
*    typically denoted by $\hat{x}(\lamda)$,$\hat{y}{\lambda}$, and $\hat{z}(\lambda)$.
*    Traditionally, the CIE1931 Colorimetric Standard Observer is used almost exclusively,
*    but is known to be not a good representation of human vision in the blue region of the
*    spectrum. We also know now that the way you see color varies with age, and your healty,
*    and that not everyone sees to same color.
*
*    In this library colors are represented by spectral distributions, to allow color modelling
*    with newer, and better standard observers, such as the CIE2015 Observer, derived from
*    the sensitivities of the cones in the retina of your eye, the biological color receptors
*    of light.
*
*    It's main purpose is to calculate `XYZ` tristimulus values for a general stimulus,
*    in from of a `Spectrum`.
*/
export class ObserverData {
  free(): void;
}
/**
* Representation of a color stimulus in a set of Red, Green, and Blue (RGB) values,
* representing its relative composition using standard primaries.
* 
* RGB values are commonly used in digital images, with the relative intensity
* of the primaries defined as three 8-bit values, with range from 0 to 255.
* As ooposed to CIE XYZ tristimulus values, which used imaginary primaries,
* displays use real primaries, typically defined in the CIE 1931 diagram.
* They cover a triangular area, referred to the _color gamut_ of a display.
*/
export class RGB {
  free(): void;
}
/**
*
*This container holds spectral values within a wavelength domain ranging from 380
*to 780 nanometers, with an interval size of 1 nanometer and a total of 401
*values. It also includes a category tag and an optional 'total' value for the
*aggregate value associated with the spectrum.
*
*A `Spectrum` can be constructed from data, but many other construction methods
*are available in this library, such as standard illuminants A and D65, Planckian
*(Black Body) illuminants, or a `Stimulus` spectrum for a pixel of an sRGB
*display.
* 
*/
export class Spectrum {
  free(): void;
/**
* Creates a new Spectrum object, using as input a `Category`, a
* Float64Array with exactly 401 datapoints, and an optional third
* parameter called total, representing the total irradiance, transmission,
* or reflectivity of the values, depending on the category of the
* spectrum. The spectral values should be associated with a wavelength
* domain from 380 to 480 nanometer, with an interval size of 1 nanometer.
*
* If the Spectral data you have uses another wavelength domain and/or a different
* wavelength interval, use the linear or sprague interpolate constructors,
* which takes a wavelength domain and spectral data as arguments.
* @param {Float64Array} data
*/
  constructor(data: Float64Array);
/**
* Returns the spectral data values, as a Float64Array containing 401 data
* points, over a wavelength domain from 380 t0 780 nanometer, with a
* stepsize of 1 nanometer.
* @returns {Float64Array}
*/
  Values(): Float64Array;
/**
*
*    This function maps spectral data with irregular intervals or intervals
*    different than 1 nanometer to the standard spectrum as used in this
*    library.
*
*    For domains with a regular interval, the wavelength slice should have a size
*    of two, containing the minimum and maximum wavelength values, both also in
*    units of meters or nanometers.
*
*    For irregular domains, this function requires a slice of wavelengths and
*    a slice of spectral data, both of the same size. The wavelengths can be
*    specified in units of meters or nanometers.
*
*    In case of duplicate wavelength values the last data values is used, so it
*    is impossible to define filters with vertical edges using this method.
*
*    ```ts, ignore
*    // Creates a linear gradient filter, with a zero transmission at 380
*    // nanometer, and full transmission at 780 nanometer. This is an example
*    // using a uniform wavelength domain as input.
*    use colorimetry as cmt;
*    # use approx::assert_ulps_eq;
*    let data = [0.0, 1.0];
*    let wl = [380.0, 780.0];
*    let mut spd = cmt::Spectrum::linear_interpolate(cmt::Category::Colorant, &wl, &data, None).unwrap().values();
*    assert_ulps_eq!(spd[0], 0.);
*    assert_ulps_eq!(spd[100], 0.25);
*    assert_ulps_eq!(spd[200], 0.5);
*    assert_ulps_eq!(spd[300], 0.75);
*    assert_ulps_eq!(spd[400], 1.0);
*
*    // Creates a top hat filter, with slanted angles, using an irregular
*    // wavelength domain.
*    let data = vec![0.0, 1.0, 1.0, 0.0];
*    let wl = vec![480.0, 490.0, 570.0, 580.0];
*    let spd = cmt::Spectrum::linear_interpolate(cmt::Category::Colorant, &wl, &data, None).unwrap().values();
*    assert_ulps_eq!(spd[0], 0.0);
*    assert_ulps_eq!(spd[100], 0.0);
*    assert_ulps_eq!(spd[110], 1.0);
*    assert_ulps_eq!(spd[190], 1.0);
*    assert_ulps_eq!(spd[200], 0.0);
*    assert_ulps_eq!(spd[300], 0.0);
*    assert_ulps_eq!(spd[400], 0.0);
*    ```
*    
* @param {Float64Array} wavelengths
* @param {Float64Array} data
* @param {any} total_js
* @returns {Spectrum}
*/
  static linearInterpolate(wavelengths: Float64Array, data: Float64Array, total_js: any): Spectrum;
/**
* Calculates the Color Rendering Index values for illuminant spectrum.
* 
* To use this function, first use `await CRI.init()`, which downloads the
* Test Color Samples required for the calculation.  These are downloaded
* seperately to limit the size of the main web assembly library.
* @returns {CRI}
*/
  cri(): CRI;
}
/**
* A set of two CIE XYZ Tristimulus values, for a Standard Observer.
* 
* One is associated with an illuminant or a reference white value, denoted by the fieldname `xyzn`, and
* a second, optional value, a set of tristimulus values for a stimulus, or 'a color' value.
* XYZ values are not often used directly, but form the basis for many colorimetric models, such as
* CIELAB and CIECAM.
*/
export class XYZ {
  free(): void;
/**
*
*    Create an XYZ Tristimuls Values object.
*
*    Accepts as arguments 
*
*    - x and y chromaticity coordinates only , using the "Cie::Std1931" observer as default
*    - x and y chromaticity coordinates, and standard observer ID as 3rd argument
*    - X, Y, and Z tristimulus values, using the "Cie::Std1931" observer as default
*    - X, Y, and Z tristimulus values, and a standard Observer ID as 4th argument
*
*    When only x and y chromaticity coordinates are specified, the luminous
*    value is set to 100.0 candela per square meter.
*
*    ```javascript, ignore
*    // Create a new XYZ object using D65 CIE 1931 chromaticity coordinates
*    const xyz = new cmt.XYZ(0.31272, 0.32903);
*    
*    // Get and check the corresponding tristimulus values, with a luminous value
*    // of 100.0
*    const [x, y, z] = xyz.values();
*    assert.assertAlmostEquals(x, 95.047, 5E-3); // D65 wikipedia
*    assert.assertAlmostEquals(y, 100.0);
*    assert.assertAlmostEquals(z, 108.883, 5E-3);
*
*    // and get back the orgiinal chromaticity coordinates:
*    const [xc, yc] = xyz.chromaticity();
*    assert.assertAlmostEquals(xc, 0.31272);
*    assert.assertAlmostEquals(yc, 0.32903);
*
*
*    // to get the luminous value:
*    const l = xyz.luminousValue();
*    assert.assertAlmostEquals(l, 100.0);
*    // D65 CIE 1931 chromaticity coordinates
*    const xyz = new cmt.XYZ(0.31272, 0.32903);
*    ```
*    
* @param {number} x
* @param {number} y
* @param {...Array<any>} opt
*/
  constructor(x: number, y: number, ...opt: Array<any>);
/**
* Get the XYZ tristimulus value as an array.
* Values of the stimulus, if present, else the illuminant.
* @returns {Array<any>}
*/
  values(): Array<any>;
/**
* Get the chromaticity coordinates
* @returns {Array<any>}
*/
  chromaticity(): Array<any>;
/**
* Get the luminous value
* @returns {number}
*/
  luminousValue(): number;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
  readonly memory: WebAssembly.Memory;
  readonly __wbg_cielab_free: (a: number, b: number) => void;
  readonly __wbg_rgb_free: (a: number, b: number) => void;
  readonly __wbg_xyz_free: (a: number, b: number) => void;
  readonly xyz_new_js: (a: number, b: number, c: number, d: number) => void;
  readonly xyz_values: (a: number) => number;
  readonly xyz_chromaticity: (a: number) => number;
  readonly xyz_luminousValue: (a: number) => number;
  readonly __wbg_observerdata_free: (a: number, b: number) => void;
  readonly stefanBoltzmann: (a: number) => number;
  readonly __wbg_spectrum_free: (a: number, b: number) => void;
  readonly spectrum_new_js: (a: number, b: number, c: number) => void;
  readonly spectrum_Values: (a: number, b: number) => void;
  readonly spectrum_linearInterpolate: (a: number, b: number, c: number, d: number, e: number, f: number) => void;
  readonly spectrum_cri: (a: number, b: number) => void;
  readonly __wbg_cri_free: (a: number, b: number) => void;
  readonly __wbg_munsellmatt_free: (a: number, b: number) => void;
  readonly __wbg_illuminant_free: (a: number, b: number) => void;
  readonly illuminant_new_js: (a: number, b: number, c: number) => void;
  readonly illuminant_Values: (a: number, b: number) => void;
  readonly illuminant_cri: (a: number, b: number) => void;
  readonly illuminant_illuminant: (a: number) => number;
  readonly __wbindgen_malloc: (a: number, b: number) => number;
  readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
  readonly __wbindgen_add_to_stack_pointer: (a: number) => number;
  readonly __wbindgen_free: (a: number, b: number, c: number) => void;
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
