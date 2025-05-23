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
 * @enum {0 | 1 | 2 | 3}
 */
export const Observer = Object.freeze({
    Std1931: 0, "0": "Std1931",
    Std1964: 1, "1": "Std1964",
    Std2015: 2, "2": "Std2015",
    Std2015_10: 3, "3": "Std2015_10",
});
/**
 *
 * A Light Weight tag, representing an RGB color space.
 * Used for example in the RGB value set, to identify the color space being used.
 *
 * @enum {0 | 1 | 2}
 */
export const RgbSpace = Object.freeze({
    SRGB: 0, "0": "SRGB",
    ADOBE: 1, "1": "ADOBE",
    DisplayP3: 2, "2": "DisplayP3",
});
/**
 * @enum {0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 | 33 | 34 | 35 | 36 | 37 | 38}
 */
export const StdIlluminant = Object.freeze({
    D65: 0, "0": "D65",
    D50: 1, "1": "D50",
    A: 2, "2": "A",
    F1: 3, "3": "F1",
    F2: 4, "4": "F2",
    F3: 5, "5": "F3",
    F4: 6, "6": "F4",
    F5: 7, "7": "F5",
    F6: 8, "8": "F6",
    F7: 9, "9": "F7",
    F8: 10, "10": "F8",
    F9: 11, "11": "F9",
    F10: 12, "12": "F10",
    F11: 13, "13": "F11",
    F12: 14, "14": "F12",
    F3_1: 15, "15": "F3_1",
    F3_2: 16, "16": "F3_2",
    F3_3: 17, "17": "F3_3",
    F3_4: 18, "18": "F3_4",
    F3_5: 19, "19": "F3_5",
    F3_6: 20, "20": "F3_6",
    F3_7: 21, "21": "F3_7",
    F3_8: 22, "22": "F3_8",
    F3_9: 23, "23": "F3_9",
    F3_10: 24, "24": "F3_10",
    F3_11: 25, "25": "F3_11",
    F3_12: 26, "26": "F3_12",
    F3_13: 27, "27": "F3_13",
    F3_14: 28, "28": "F3_14",
    F3_15: 29, "29": "F3_15",
    LED_B1: 30, "30": "LED_B1",
    LED_B2: 31, "31": "LED_B2",
    LED_B3: 32, "32": "LED_B3",
    LED_B4: 33, "33": "LED_B4",
    LED_B5: 34, "34": "LED_B5",
    LED_BH1: 35, "35": "LED_BH1",
    LED_RGB1: 36, "36": "LED_RGB1",
    LED_V1: 37, "37": "LED_V1",
    LED_V2: 38, "38": "LED_V2",
});

const CRIFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_cri_free(ptr >>> 0, 1));
/**
 * Encapcsulated Array of calculated Ri values, from a test light source.
 */
export class CRI {

    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(CRI.prototype);
        obj.__wbg_ptr = ptr;
        CRIFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        CRIFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_cri_free(ptr, 0);
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
     * Calculates the Color Rendering Index values for illuminant spectrum.
     *
     * To use this function, first use `await CRI.init()`, which downloads the
     * Test Color Samples required for the calculation.  These are downloaded
     * seperately to limit the size of the main web assembly library.
     * @returns {CRI}
     */
    cri() {
        const ret = wasm.illuminant_cri(this.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return CRI.__wrap(ret[0]);
    }
    /**
     * Get the StdIlluminant spectrum. Typically you don't need to use the Spectrum itself, as many
     * methods just accept the StdIlluminant directly.
     * @param {StdIlluminant} stdill
     * @returns {Illuminant}
     */
    static illuminant(stdill) {
        const ret = wasm.illuminant_illuminant(stdill);
        return Illuminant.__wrap(ret);
    }
}

const MunsellMattFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_munsellmatt_free(ptr >>> 0, 1));

export class MunsellMatt {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MunsellMattFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_munsellmatt_free(ptr, 0);
    }
}

const ObserverDataFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_observerdata_free(ptr >>> 0, 1));
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

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ObserverDataFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_observerdata_free(ptr, 0);
    }
}

const RGBFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_rgb_free(ptr >>> 0, 1));
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

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        RGBFinalization.unregister(this);
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
    static linearInterpolate(wavelengths, data, total_js) {
        const ptr0 = passArrayF64ToWasm0(wavelengths, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF64ToWasm0(data, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.spectrum_linearInterpolate(ptr0, len0, ptr1, len1, total_js);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return Spectrum.__wrap(ret[0]);
    }
    /**
     * Calculates the Color Rendering Index values for illuminant spectrum.
     *
     * To use this function, first use `await CRI.init()`, which downloads the
     * Test Color Samples required for the calculation.  These are downloaded
     * seperately to limit the size of the main web assembly library.
     * @returns {CRI}
     */
    cri() {
        const ret = wasm.spectrum_cri(this.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return CRI.__wrap(ret[0]);
    }
}

const ViewConditionsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_viewconditions_free(ptr >>> 0, 1));

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
    /**
     * Degree of Adaptation, if omitted, formula 4.3 of CIE248:2022 is used.``
     * @returns {number | undefined}
     */
    get dopt() {
        const ret = wasm.__wbg_get_viewconditions_dopt(this.__wbg_ptr);
        return ret[0] === 0 ? undefined : ret[1];
    }
    /**
     * Degree of Adaptation, if omitted, formula 4.3 of CIE248:2022 is used.``
     * @param {number | null} [arg0]
     */
    set dopt(arg0) {
        wasm.__wbg_set_viewconditions_dopt(this.__wbg_ptr, !isLikeNone(arg0), isLikeNone(arg0) ? 0 : arg0);
    }
    /**
     * @returns {number}
     */
    get f() {
        const ret = wasm.__wbg_get_viewconditions_f(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set f(arg0) {
        wasm.__wbg_set_viewconditions_f(this.__wbg_ptr, arg0);
    }
    /**
     * Adaptation Luminance, in cd/m2.
     * La = Lw/5, with Lw: luminance of a perfect white object
     * @returns {number}
     */
    get la() {
        const ret = wasm.__wbg_get_viewconditions_la(this.__wbg_ptr);
        return ret;
    }
    /**
     * Adaptation Luminance, in cd/m2.
     * La = Lw/5, with Lw: luminance of a perfect white object
     * @param {number} arg0
     */
    set la(arg0) {
        wasm.__wbg_set_viewconditions_la(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get nc() {
        const ret = wasm.__wbg_get_viewconditions_nc(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set nc(arg0) {
        wasm.__wbg_set_viewconditions_nc(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get yb() {
        const ret = wasm.__wbg_get_viewconditions_yb(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set yb(arg0) {
        wasm.__wbg_set_viewconditions_yb(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get c() {
        const ret = wasm.__wbg_get_viewconditions_c(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set c(arg0) {
        wasm.__wbg_set_viewconditions_c(this.__wbg_ptr, arg0);
    }
}

const XYZFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_xyz_free(ptr >>> 0, 1));
/**
 * A set of two CIE XYZ Tristimulus values, for a Standard Observer.
 *
 * One is associated with an illuminant or a reference white value, denoted by the fieldname `xyzn`, and
 * a second, optional value, a set of tristimulus values for a stimulus, or 'a color' value.
 * XYZ values are not often used directly, but form the basis for many colorimetric models, such as
 * CIELAB and CIECAM.
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
     * Values of the stimulus, if present, else the illuminant.
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
