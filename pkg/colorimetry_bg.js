let wasm;
export function __wbg_set_wasm(val) {
    wasm = val;
}


const lTextDecoder = typeof TextDecoder === 'undefined' ? (0, module.require)('util').TextDecoder : TextDecoder;

let cachedTextDecoder = new lTextDecoder('utf-8', { ignoreBOM: true, fatal: true });

cachedTextDecoder.decode();

let cachedUint8Memory0 = null;

function getUint8Memory0() {
    if (cachedUint8Memory0 === null || cachedUint8Memory0.byteLength === 0) {
        cachedUint8Memory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8Memory0;
}

function getStringFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return cachedTextDecoder.decode(getUint8Memory0().subarray(ptr, ptr + len));
}

const heap = new Array(128).fill(undefined);

heap.push(undefined, null, true, false);

let heap_next = heap.length;

function addHeapObject(obj) {
    if (heap_next === heap.length) heap.push(heap.length + 1);
    const idx = heap_next;
    heap_next = heap[idx];

    heap[idx] = obj;
    return idx;
}
/**
* Stefan Boltzmann law: Blackbody's radiant emittance (W m<sup>-2</sup>), as function of its absolute
* temperature (K).
* @param {number} temperature
* @returns {number}
*/
export function stefanBoltzmann(temperature) {
    const ret = wasm.stefanBoltzmann(temperature);
    return ret;
}

/**
*/
export const ObsId = Object.freeze({ Std1931:0,"0":"Std1931",Std1976:1,"1":"Std1976",Std2015:2,"2":"Std2015",Std2015_10:3,"3":"Std2015_10", });
/**
*/
export const RGBSpace = Object.freeze({ SRGB:0,"0":"SRGB",SRGBFloat:1,"1":"SRGBFloat",SRGB16:2,"2":"SRGB16", });

const LabFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_lab_free(ptr >>> 0));
/**
*/
export class Lab {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        LabFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_lab_free(ptr);
    }
}

const ObserverFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_observer_free(ptr >>> 0));
/**
* A data structure to define Standard Observers, such as the CIE 1931 2º and the CIE 2015 standard observers.
* These are defined in the form of three discrete representations of color matching functions.
*/
export class Observer {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ObserverFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_observer_free(ptr);
    }
}

const RGBFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_rgb_free(ptr >>> 0));
/**
*/
export class RGB {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        RGBFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_rgb_free(ptr);
    }
}

const SpectrumFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_spectrum_free(ptr >>> 0));
/**
*/
export class Spectrum {

    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(Spectrum.prototype);
        obj.__wbg_ptr = ptr;
        SpectrumFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        SpectrumFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_spectrum_free(ptr);
    }
    /**
    * @returns {Spectrum}
    */
    static d65() {
        const ret = wasm.spectrum_d65();
        return Spectrum.__wrap(ret);
    }
    /**
    * @returns {Spectrum}
    */
    static d50() {
        const ret = wasm.spectrum_d50();
        return Spectrum.__wrap(ret);
    }
    /**
    * @returns {Spectrum}
    */
    static a() {
        const ret = wasm.spectrum_a();
        return Spectrum.__wrap(ret);
    }
    /**
    * @param {number} gval
    * @returns {Spectrum}
    */
    static grey(gval) {
        const ret = wasm.spectrum_grey(gval);
        return Spectrum.__wrap(ret);
    }
    /**
    * @param {number} gval
    * @returns {Spectrum}
    */
    static grey_filter(gval) {
        const ret = wasm.spectrum_grey_filter(gval);
        return Spectrum.__wrap(ret);
    }
    /**
    * @returns {Spectrum}
    */
    static white() {
        const ret = wasm.spectrum_white();
        return Spectrum.__wrap(ret);
    }
    /**
    * @returns {Spectrum}
    */
    static black() {
        const ret = wasm.spectrum_black();
        return Spectrum.__wrap(ret);
    }
    /**
    * @returns {Spectrum}
    */
    static equal_energy() {
        const ret = wasm.spectrum_equal_energy();
        return Spectrum.__wrap(ret);
    }
    /**
    * @param {number} center
    * @param {number} width
    * @returns {Spectrum}
    */
    static led(center, width) {
        const ret = wasm.spectrum_led(center, width);
        return Spectrum.__wrap(ret);
    }
}

const XYZFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_xyz_free(ptr >>> 0));
/**
* A set of CIE XYZ Tristimulus values, associated with a Standard Observer.
* They are calculated using the spectrum of a Stimulus, such as beam of light
* reflected from a color patch, or emitted from a pixel of a display.
* XYZ values are not often used directly, but form the basis for many colorimetric models,
* such as CIELAB and CIECAM.
*/
export class XYZ {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        XYZFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_xyz_free(ptr);
    }
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
    constructor(_x, _y, _z, _obs) {
        const ret = wasm.xyz_new_js(_x, _y, addHeapObject(_z), addHeapObject(_obs));
        this.__wbg_ptr = ret >>> 0;
        return this;
    }
}

export function __wbindgen_throw(arg0, arg1) {
    throw new Error(getStringFromWasm0(arg0, arg1));
};

