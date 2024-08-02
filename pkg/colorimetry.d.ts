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
export enum ObsId {
  Std1931 = 0,
  Std1976 = 1,
  Std2015 = 2,
  Std2015_10 = 3,
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
