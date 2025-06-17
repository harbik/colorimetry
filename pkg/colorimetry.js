let wasm;

const cachedTextDecoder = (typeof TextDecoder !== 'undefined' ? new TextDecoder('utf-8', { ignoreBOM: true, fatal: true }) : { decode: () => { throw Error('TextDecoder not available') } } );

if (typeof TextDecoder !== 'undefined') { cachedTextDecoder.decode(); };

let cachedUint8ArrayMemory0 = null;

function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function getStringFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

function isLikeNone(x) {
    return x === undefined || x === null;
}

let cachedDataViewMemory0 = null;

function getDataViewMemory0() {
    if (cachedDataViewMemory0 === null || cachedDataViewMemory0.buffer.detached === true || (cachedDataViewMemory0.buffer.detached === undefined && cachedDataViewMemory0.buffer !== wasm.memory.buffer)) {
        cachedDataViewMemory0 = new DataView(wasm.memory.buffer);
    }
    return cachedDataViewMemory0;
}

let WASM_VECTOR_LEN = 0;

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
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8ArrayMemory0();

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
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = encodeString(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

let cachedFloat64ArrayMemory0 = null;

function getFloat64ArrayMemory0() {
    if (cachedFloat64ArrayMemory0 === null || cachedFloat64ArrayMemory0.byteLength === 0) {
        cachedFloat64ArrayMemory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64ArrayMemory0;
}

function passArrayF64ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 8, 8) >>> 0;
    getFloat64ArrayMemory0().set(arg, ptr / 8);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function takeFromExternrefTable0(idx) {
    const value = wasm.__wbindgen_export_0.get(idx);
    wasm.__externref_table_dealloc(idx);
    return value;
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
}
/**
 * @enum {0 | 1}
 */
export const CieIlluminant = Object.freeze({
    D65: 0, "0": "D65",
    D50: 1, "1": "D50",
});
/**
 * Light-weight identifier added to the `XYZ` and `RGB` datasets,
 *    representing the colorimetric standard observer used.
 *
 *    No data included here, which would be the Rust way, but that does not work with wasm-bindgen.
 *    This can be directly used in JavaScript, and has the benefit to be just an index.
 * @enum {0}
 */
export const Observer = Object.freeze({
    Cie1931: 0, "0": "Cie1931",
});
/**
 * A Light Weight tag, representing an RGB color space.
 * Used for example in the RGB value set, to identify the color space being used.
 * @enum {0 | 1 | 2}
 */
export const RgbSpace = Object.freeze({
    SRGB: 0, "0": "SRGB",
    Adobe: 1, "1": "Adobe",
    DisplayP3: 2, "2": "DisplayP3",
});

const ChromaticityFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_chromaticity_free(ptr >>> 0, 1));
/**
 * A chromaticity coordinate with x and y values.
 */
export class Chromaticity {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ChromaticityFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_chromaticity_free(ptr, 0);
    }
}

const CieLabFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_cielab_free(ptr >>> 0, 1));

export class CieLab {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        CieLabFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_cielab_free(ptr, 0);
    }
}

const IlluminantFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_illuminant_free(ptr >>> 0, 1));
/**
 * # Illuminant
 *
 * An illuminant is a spectral power distribution that represents the
 * spectral power density of a light source (sun, bulb, LED, etc.) in
 * W/m²/nm over 380–780 nm (401 samples).
 */
export class Illuminant {

    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(Illuminant.prototype);
        obj.__wbg_ptr = ptr;
        IlluminantFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        IlluminantFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_illuminant_free(ptr, 0);
    }
    /**
     * Create a new illuminant spectrum from the given data.
     *
     * The data must be the 401 values from 380 to 780 nm, with an interval size of 1 nanometer.
     * @param {Float64Array} data
     */
    constructor(data) {
        const ptr0 = passArrayF64ToWasm0(data, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.illuminant_new_js(ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0] >>> 0;
        IlluminantFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Returns the spectral data values, as a Float64Array containing 401 data
     * points, over a wavelength domain from 380 t0 780 nanometer, with a
     * stepsize of 1 nanometer.
     * @returns {Float64Array}
     */
    Values() {
        const ret = wasm.illuminant_Values(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Get the CieIlluminant spectrum. Typically you don't need to use the Spectrum itself, as many
     * methods just accept the CieIlluminant directly.
     * @param {CieIlluminant} stdill
     * @returns {Illuminant}
     */
    static illuminant(stdill) {
        const ret = wasm.illuminant_illuminant(stdill);
        return Illuminant.__wrap(ret);
    }
}

const RelXYZFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_relxyz_free(ptr >>> 0, 1));
/**
 * # Related Tristimulus Values
 *
 * Tristimulus Values for a given sample and reference white,
 * used to represent related colors as used in various color
 * models. Typically the reference white is normalized to have
 * an Y-value of 100
 */
export class RelXYZ {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        RelXYZFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_relxyz_free(ptr, 0);
    }
}

const RgbFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_rgb_free(ptr >>> 0, 1));
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

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        RgbFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_rgb_free(ptr, 0);
    }
}

const SpectrumFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_spectrum_free(ptr >>> 0, 1));
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
        wasm.__wbg_spectrum_free(ptr, 0);
    }
    /**
     * Create a new spectrum from the given data.
     *
     * The data must be the 401 values from 380 to 780 nm, with an interval size of 1 nanometer.
     *
     * If the Spectral data you have uses another wavelength domain and/or a different
     * wavelength interval, use the linear interpolate constructor,
     * which takes a wavelength domain and spectral data as arguments.
     * @param {Float64Array} data
     */
    constructor(data) {
        const ptr0 = passArrayF64ToWasm0(data, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.spectrum_new_js(ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0] >>> 0;
        SpectrumFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Returns the spectral data values, as a Float64Array containing 401 data
     * points, over a wavelength domain from 380 t0 780 nanometer, with a
     * stepsize of 1 nanometer.
     * @returns {Float64Array}
     */
    Values() {
        const ret = wasm.spectrum_Values(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
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
     * @param {Float64Array} wavelengths
     * @param {Float64Array} data
     * @returns {Spectrum}
     */
    static linearInterpolate(wavelengths, data) {
        const ptr0 = passArrayF64ToWasm0(wavelengths, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF64ToWasm0(data, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.spectrum_linearInterpolate(ptr0, len0, ptr1, len1);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return Spectrum.__wrap(ret[0]);
    }
}

const ViewConditionsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_viewconditions_free(ptr >>> 0, 1));
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

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ViewConditionsFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_viewconditions_free(ptr, 0);
    }
}

const WideRgbFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_widergb_free(ptr >>> 0, 1));
/**
 * Represents a color stimulus using unconstrained Red, Green, and Blue (RGB) floating-point values
 * within a device's RGB color space. The values can extend beyond the typical 0.0 to 1.0 range,
 * allowing for out-of-gamut colors that cannot be accurately represented by the device.
 */
export class WideRgb {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        WideRgbFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_widergb_free(ptr, 0);
    }
}

const XYZFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_xyz_free(ptr >>> 0, 1));
/**
 * Represents a color by its tristimulus value XYZ color space.
 *
 * The `XYZ` struct represents the tristimulus values (X, Y, Z) and the associated observer.
 * The observer defines the color matching functions used for the conversion.
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
        wasm.__wbg_xyz_free(ptr, 0);
    }
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
     * @param {number} x
     * @param {number} y
     * @param {...Array<any>} opt
     */
    constructor(x, y, ...opt) {
        const ret = wasm.xyz_new_js(x, y, opt);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0] >>> 0;
        XYZFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Get the XYZ tristimulus value as an array.
     * @returns {Array<any>}
     */
    values() {
        const ret = wasm.xyz_values(this.__wbg_ptr);
        return ret;
    }
    /**
     * Get the chromaticity coordinates
     * @returns {Array<any>}
     */
    chromaticity() {
        const ret = wasm.xyz_chromaticity(this.__wbg_ptr);
        return ret;
    }
    /**
     * Get the luminous value, Y.
     * @returns {number}
     */
    y() {
        const ret = wasm.xyz_y(this.__wbg_ptr);
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
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

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
    imports.wbg.__wbg_get_0c3cc364764a0b98 = function(arg0, arg1) {
        const ret = arg0[arg1 >>> 0];
        return ret;
    };
    imports.wbg.__wbg_length_12246a78d2f65d3a = function(arg0) {
        const ret = arg0.length;
        return ret;
    };
    imports.wbg.__wbg_of_894f51209b51cdf4 = function(arg0, arg1) {
        const ret = Array.of(arg0, arg1);
        return ret;
    };
    imports.wbg.__wbg_of_968f342775805f3a = function(arg0, arg1, arg2) {
        const ret = Array.of(arg0, arg1, arg2);
        return ret;
    };
    imports.wbg.__wbindgen_error_new = function(arg0, arg1) {
        const ret = new Error(getStringFromWasm0(arg0, arg1));
        return ret;
    };
    imports.wbg.__wbindgen_init_externref_table = function() {
        const table = wasm.__wbindgen_export_0;
        const offset = table.grow(4);
        table.set(0, undefined);
        table.set(offset + 0, undefined);
        table.set(offset + 1, null);
        table.set(offset + 2, true);
        table.set(offset + 3, false);
        ;
    };
    imports.wbg.__wbindgen_number_get = function(arg0, arg1) {
        const obj = arg1;
        const ret = typeof(obj) === 'number' ? obj : undefined;
        getDataViewMemory0().setFloat64(arg0 + 8 * 1, isLikeNone(ret) ? 0 : ret, true);
        getDataViewMemory0().setInt32(arg0 + 4 * 0, !isLikeNone(ret), true);
    };
    imports.wbg.__wbindgen_number_new = function(arg0) {
        const ret = arg0;
        return ret;
    };
    imports.wbg.__wbindgen_string_get = function(arg0, arg1) {
        const obj = arg1;
        const ret = typeof(obj) === 'string' ? obj : undefined;
        var ptr1 = isLikeNone(ret) ? 0 : passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len1 = WASM_VECTOR_LEN;
        getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
        getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
    };
    imports.wbg.__wbindgen_string_new = function(arg0, arg1) {
        const ret = getStringFromWasm0(arg0, arg1);
        return ret;
    };
    imports.wbg.__wbindgen_throw = function(arg0, arg1) {
        throw new Error(getStringFromWasm0(arg0, arg1));
    };
    imports.wbg.__wbindgen_try_into_number = function(arg0) {
        let result;
        try { result = +arg0 } catch (e) { result = e }
        const ret = result;
        return ret;
    };

    return imports;
}

function __wbg_init_memory(imports, memory) {

}

function __wbg_finalize_init(instance, module) {
    wasm = instance.exports;
    __wbg_init.__wbindgen_wasm_module = module;
    cachedDataViewMemory0 = null;
    cachedFloat64ArrayMemory0 = null;
    cachedUint8ArrayMemory0 = null;


    wasm.__wbindgen_start();
    return wasm;
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (typeof module !== 'undefined') {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();

    __wbg_init_memory(imports);

    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }

    const instance = new WebAssembly.Instance(module, imports);

    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (typeof module_or_path !== 'undefined') {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (typeof module_or_path === 'undefined') {
        module_or_path = new URL('colorimetry_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    __wbg_init_memory(imports);

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync };
export default __wbg_init;
