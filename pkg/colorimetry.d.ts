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
export enum Category {
/**
* The spectral distribution of onne or more sources, illuminating a color sample
*/
  Illuminant = 0,
/**
* A Filter spectrum , such as a wratten or glass filter, which changes the properties of an illuminant.
*/
  Filter = 1,
/**
* The spectrum of a color patch, typically consisting of a paint or ink on a substrate, as measured with a spectrophotomteer.
*/
  ColorPatch = 2,
/**
* A ray of light from object we are looking at, typically an illuminated by an illuminant.
*/
  Stimulus = 3,
/**
* The type of spectrum is unknown.
*/
  Unknown = 4,
}
/**
*/
export enum RGBSpace {
  SRGB = 0,
  SRGBFloat = 1,
  SRGB16 = 2,
}
/**
*/
export enum ObsId {
  Std1931 = 0,
  Std1976 = 1,
  Std2015 = 2,
  Std2015_10 = 3,
}
/**
*/
export class Lab {
  free(): void;
}
/**
* A data structure to define Standard Observers, such as the CIE 1931 2º and the CIE 2015 standard observers.
* These are defined in the form of three discrete representations of color matching functions.
*/
export class Observer {
  free(): void;
}
/**
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
*The categories are:
*
*- `Illuminant`: a spectral irradiance distribution with values given in watts
*    per square meter per nanometer, and a `total` value given in watts per square
*    meter.
*- `Filter`: a spectral transmission function with unitless values ranging from
*    0.0 to 1.0, and the `total` value representing the total power transmission of
*    the filter.
*- `Substrate`: a spectral transmission function when combined with a `Filter`
*    and spectral reflectivity function combined with a `ColorPatch`.
*- `ColorPatch`: a spectral reflectivity function with unitless values ranging from
*    0.0 to 1.0.
*- `Stimulus`: a spectral radiance distribution of a beam of light entering
*    through the pupil of our eyes, on its way to be processed and triggering a
*    sensation of color in our mind. Spectral data of a stimulus have a unit of watt
*    per square meter per nanometer per steradian, and a total.
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
* @param {Category} cat
* @param {Float64Array} data
* @param {number | undefined} [total]
*/
  constructor(cat: Category, data: Float64Array, total?: number);
/**
* Returns the spectral data values, as a Float64Array containing 401 data
* points, over a wavelength domain from 380 t0 780 nanometer, with a
* stepsize of 1 nanometer.
* @returns {Float64Array}
*/
  Values(): Float64Array;
}
/**
* A set of CIE XYZ Tristimulus values, associated with a Standard Observer.
* They are calculated using the spectrum of a Stimulus, such as beam of light
* reflected from a color patch, or emitted from a pixel of a display.
* XYZ values are not often used directly, but form the basis for many colorimetric models,
* such as CIELAB and CIECAM.
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
  readonly __wbg_spectrum_free: (a: number) => void;
  readonly spectrum_new_js: (a: number, b: number, c: number, d: number, e: number, f: number) => void;
  readonly spectrum_Values: (a: number, b: number) => void;
  readonly __wbg_observer_free: (a: number) => void;
  readonly __wbg_rgb_free: (a: number) => void;
  readonly stefanBoltzmann: (a: number) => number;
  readonly __wbg_xyz_free: (a: number) => void;
  readonly xyz_new_js: (a: number, b: number, c: number, d: number) => void;
  readonly xyz_values: (a: number) => number;
  readonly xyz_chromaticity: (a: number) => number;
  readonly xyz_luminousValue: (a: number) => number;
  readonly __wbg_lab_free: (a: number) => void;
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
* @param {SyncInitInput} module
*
* @returns {InitOutput}
*/
export function initSync(module: SyncInitInput): InitOutput;

/**
* If `module_or_path` is {RequestInfo} or {URL}, makes a request and
* for everything else, calls `WebAssembly.instantiate` directly.
*
* @param {InitInput | Promise<InitInput>} module_or_path
*
* @returns {Promise<InitOutput>}
*/
export default function __wbg_init (module_or_path?: InitInput | Promise<InitInput>): Promise<InitOutput>;
