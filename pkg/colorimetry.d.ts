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
  Illuminant = 0,
  Filter = 1,
  Substrate = 2,
  Colorant = 3,
  Stimulus = 4,
  Unknown = 5,
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
*    and spectral reflectivity function combined with a `Colorant`.
*- `Colorant`: a spectral reflectivity function with unitless values ranging from
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
* XYZ Tristimuls Values JavaScript Constructor
* 
* Accepts as arguments 
* - x and y chromaticity coordinates only , using the "Cie::Std1931" observer as default
* - x and y chromaticity coordinates, and standard observer ID as 3rd argument
* - X, Y, and Z tristimulus values, using the "Cie::Std1931" observer as default
* - X, Y, and Z tristimulus values, and standard Observer ID as 4th argument
* @param {number} _x
* @param {number} _y
* @param {any} _z
* @param {any} _obs
*/
  constructor(_x: number, _y: number, _z: any, _obs: any);
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
  readonly xyz_new_js: (a: number, b: number, c: number, d: number) => number;
  readonly __wbg_lab_free: (a: number) => void;
  readonly __wbindgen_add_to_stack_pointer: (a: number) => number;
  readonly __wbindgen_malloc: (a: number, b: number) => number;
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
