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
*/
export class Spectrum {
  free(): void;
/**
* @returns {Spectrum}
*/
  static d65(): Spectrum;
/**
* @returns {Spectrum}
*/
  static d50(): Spectrum;
/**
* @returns {Spectrum}
*/
  static a(): Spectrum;
/**
* @param {number} gval
* @returns {Spectrum}
*/
  static grey(gval: number): Spectrum;
/**
* @param {number} gval
* @returns {Spectrum}
*/
  static grey_filter(gval: number): Spectrum;
/**
* @returns {Spectrum}
*/
  static white(): Spectrum;
/**
* @returns {Spectrum}
*/
  static black(): Spectrum;
/**
* @returns {Spectrum}
*/
  static equal_energy(): Spectrum;
/**
* @param {number} center
* @param {number} width
* @returns {Spectrum}
*/
  static led(center: number, width: number): Spectrum;
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
  readonly spectrum_d65: () => number;
  readonly spectrum_d50: () => number;
  readonly spectrum_a: () => number;
  readonly spectrum_grey: (a: number) => number;
  readonly spectrum_grey_filter: (a: number) => number;
  readonly spectrum_white: () => number;
  readonly spectrum_black: () => number;
  readonly spectrum_equal_energy: () => number;
  readonly spectrum_led: (a: number, b: number) => number;
  readonly __wbg_observer_free: (a: number) => void;
  readonly __wbg_xyz_free: (a: number) => void;
  readonly xyz_new_js: (a: number, b: number, c: number, d: number) => number;
  readonly __wbg_lab_free: (a: number) => void;
  readonly __wbg_rgb_free: (a: number) => void;
  readonly stefanBoltzmann: (a: number) => number;
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
