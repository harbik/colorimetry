/// <reference types="./colorimetry.d.ts"/>
let wasm;

const heap = new Array(128).fill(undefined);

heap.push(undefined, null, true, false);

function getObject(idx) { return heap[idx]; }

let heap_next = heap.length;

function addHeapObject(obj) {
    if (heap_next === heap.length) heap.push(heap.length + 1);
    const idx = heap_next;
    heap_next = heap[idx];

    heap[idx] = obj;
    return idx;
}

function isLikeNone(x) {
    return x === undefined || x === null;
}

let cachedFloat64Memory0 = null;

function getFloat64Memory0() {
    if (cachedFloat64Memory0 === null || cachedFloat64Memory0.byteLength === 0) {
        cachedFloat64Memory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64Memory0;
}

let cachedInt32Memory0 = null;

function getInt32Memory0() {
    if (cachedInt32Memory0 === null || cachedInt32Memory0.byteLength === 0) {
        cachedInt32Memory0 = new Int32Array(wasm.memory.buffer);
    }
    return cachedInt32Memory0;
}

function dropObject(idx) {
    if (idx < 132) return;
    heap[idx] = heap_next;
    heap_next = idx;
}

function takeObject(idx) {
    const ret = getObject(idx);
    dropObject(idx);
    return ret;
}

let WASM_VECTOR_LEN = 0;

let cachedUint8Memory0 = null;

function getUint8Memory0() {
    if (cachedUint8Memory0 === null || cachedUint8Memory0.byteLength === 0) {
        cachedUint8Memory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8Memory0;
}

const cachedTextEncoder = (typeof TextEncoder !== 'undefined' ? new TextEncoder('utf-8') : { encode: () => { throw Error('TextEncoder not available') } } );

const encodeString = (typeof cachedTextEncoder.encodeInto === 'function'
    ? function (arg, view) {
    return cachedTextEncoder.encodeInto(arg, view);
}
    : function (arg, view) {
    const buf = cachedTextEncoder.encode(arg);
    view.set(buf);
    return {
        read: arg.length,
        written: buf.length
    };
});

function passStringToWasm0(arg, malloc, realloc) {

    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8Memory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8Memory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }

    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8Memory0().subarray(ptr + offset, ptr + len);
        const ret = encodeString(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

const cachedTextDecoder = (typeof TextDecoder !== 'undefined' ? new TextDecoder('utf-8', { ignoreBOM: true, fatal: true }) : { decode: () => { throw Error('TextDecoder not available') } } );

if (typeof TextDecoder !== 'undefined') { cachedTextDecoder.decode(); };

function getStringFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return cachedTextDecoder.decode(getUint8Memory0().subarray(ptr, ptr + len));
}

function passArrayF64ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 8, 8) >>> 0;
    getFloat64Memory0().set(arg, ptr / 8);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64Memory0().subarray(ptr / 8, ptr / 8 + len);
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

let stack_pointer = 128;

function addBorrowedObject(obj) {
    if (stack_pointer == 1) throw new Error('out of js stack');
    heap[--stack_pointer] = obj;
    return stack_pointer;
}
/**
*/
export const Category = Object.freeze({
/**
* The spectral distribution of onne or more sources, illuminating a color sample
*/
Illuminant:0,"0":"Illuminant",
/**
* A Filter spectrum , such as a wratten or glass filter, which changes the properties of an illuminant.
*/
Filter:1,"1":"Filter",
/**
* The spectrum of a color patch, typically consisting of a paint or ink on a substrate, as measured with a spectrophotomteer.
*/
ColorPatch:2,"2":"ColorPatch",
/**
* A ray of light from object we are looking at, typically an illuminated by an illuminant.
*/
Stimulus:3,"3":"Stimulus",
/**
* The type of spectrum is unknown.
*/
Unknown:4,"4":"Unknown", });
/**
*/
export const RGBSpace = Object.freeze({ SRGB:0,"0":"SRGB",SRGBFloat:1,"1":"SRGBFloat",SRGB16:2,"2":"SRGB16", });
/**
*/
export const ObsId = Object.freeze({ Std1931:0,"0":"Std1931",Std1976:1,"1":"Std1976",Std2015:2,"2":"Std2015",Std2015_10:3,"3":"Std2015_10", });

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
    constructor(cat, data, total) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArrayF64ToWasm0(data, wasm.__wbindgen_malloc);
            const len0 = WASM_VECTOR_LEN;
            wasm.spectrum_new_js(retptr, cat, ptr0, len0, !isLikeNone(total), isLikeNone(total) ? 0 : total);
            var r0 = getInt32Memory0()[retptr / 4 + 0];
            var r1 = getInt32Memory0()[retptr / 4 + 1];
            var r2 = getInt32Memory0()[retptr / 4 + 2];
            if (r2) {
                throw takeObject(r1);
            }
            this.__wbg_ptr = r0 >>> 0;
            return this;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
    * Returns the spectral data values, as a Float64Array containing 401 data
    * points, over a wavelength domain from 380 t0 780 nanometer, with a
    * stepsize of 1 nanometer.
    * @returns {Float64Array}
    */
    Values() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.spectrum_Values(retptr, this.__wbg_ptr);
            var r0 = getInt32Memory0()[retptr / 4 + 0];
            var r1 = getInt32Memory0()[retptr / 4 + 1];
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_free(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
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
    constructor(x, y, ...opt) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.xyz_new_js(retptr, x, y, addBorrowedObject(opt));
            var r0 = getInt32Memory0()[retptr / 4 + 0];
            var r1 = getInt32Memory0()[retptr / 4 + 1];
            var r2 = getInt32Memory0()[retptr / 4 + 2];
            if (r2) {
                throw takeObject(r1);
            }
            this.__wbg_ptr = r0 >>> 0;
            return this;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            heap[stack_pointer++] = undefined;
        }
    }
    /**
    * Get the XYZ tristimulus value as an array.
    * @returns {Array<any>}
    */
    values() {
        const ret = wasm.xyz_values(this.__wbg_ptr);
        return takeObject(ret);
    }
    /**
    * Get the chromaticity coordinates
    * @returns {Array<any>}
    */
    chromaticity() {
        const ret = wasm.xyz_chromaticity(this.__wbg_ptr);
        return takeObject(ret);
    }
    /**
    * Get the luminous value
    * @returns {number}
    */
    luminousValue() {
        const ret = wasm.xyz_luminousValue(this.__wbg_ptr);
        return ret;
    }
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);

            } catch (e) {
                if (module.headers.get('Content-Type') != 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else {
                    throw e;
                }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);

    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };

        } else {
            return instance;
        }
    }
}

function __wbg_get_imports() {
    const imports = {};
    imports.wbg = {};
    imports.wbg.__wbindgen_try_into_number = function(arg0) {
        let result;
    try { result = +getObject(arg0) } catch (e) { result = e }
    const ret = result;
    return addHeapObject(ret);
};
imports.wbg.__wbindgen_number_get = function(arg0, arg1) {
    const obj = getObject(arg1);
    const ret = typeof(obj) === 'number' ? obj : undefined;
    getFloat64Memory0()[arg0 / 8 + 1] = isLikeNone(ret) ? 0 : ret;
    getInt32Memory0()[arg0 / 4 + 0] = !isLikeNone(ret);
};
imports.wbg.__wbindgen_object_drop_ref = function(arg0) {
    takeObject(arg0);
};
imports.wbg.__wbindgen_number_new = function(arg0) {
    const ret = arg0;
    return addHeapObject(ret);
};
imports.wbg.__wbindgen_string_get = function(arg0, arg1) {
    const obj = getObject(arg1);
    const ret = typeof(obj) === 'string' ? obj : undefined;
    var ptr1 = isLikeNone(ret) ? 0 : passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    getInt32Memory0()[arg0 / 4 + 1] = len1;
    getInt32Memory0()[arg0 / 4 + 0] = ptr1;
};
imports.wbg.__wbindgen_error_new = function(arg0, arg1) {
    const ret = new Error(getStringFromWasm0(arg0, arg1));
    return addHeapObject(ret);
};
imports.wbg.__wbindgen_string_new = function(arg0, arg1) {
    const ret = getStringFromWasm0(arg0, arg1);
    return addHeapObject(ret);
};
imports.wbg.__wbg_get_bd8e338fbd5f5cc8 = function(arg0, arg1) {
    const ret = getObject(arg0)[arg1 >>> 0];
    return addHeapObject(ret);
};
imports.wbg.__wbg_length_cd7af8117672b8b8 = function(arg0) {
    const ret = getObject(arg0).length;
    return ret;
};
imports.wbg.__wbg_of_647f9238b4d5407a = function(arg0, arg1) {
    const ret = Array.of(getObject(arg0), getObject(arg1));
    return addHeapObject(ret);
};
imports.wbg.__wbg_of_6a70eed8d41f469c = function(arg0, arg1, arg2) {
    const ret = Array.of(getObject(arg0), getObject(arg1), getObject(arg2));
    return addHeapObject(ret);
};
imports.wbg.__wbindgen_throw = function(arg0, arg1) {
    throw new Error(getStringFromWasm0(arg0, arg1));
};

return imports;
}

function __wbg_init_memory(imports, maybe_memory) {

}

function __wbg_finalize_init(instance, module) {
    wasm = instance.exports;
    __wbg_init.__wbindgen_wasm_module = module;
    cachedFloat64Memory0 = null;
    cachedInt32Memory0 = null;
    cachedUint8Memory0 = null;


    return wasm;
}

function initSync(module) {
    if (wasm !== undefined) return wasm;

    const imports = __wbg_get_imports();

    __wbg_init_memory(imports);

    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }

    const instance = new WebAssembly.Instance(module, imports);

    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(input) {
    if (wasm !== undefined) return wasm;

    if (typeof input === 'undefined') {
        input = new URL('colorimetry_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof input === 'string' || (typeof Request === 'function' && input instanceof Request) || (typeof URL === 'function' && input instanceof URL)) {
        input = fetch(input);
    }

    __wbg_init_memory(imports);

    const { instance, module } = await __wbg_load(await input, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync }
export default __wbg_init;
