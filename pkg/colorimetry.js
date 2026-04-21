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
 * A **lightweight enum** representing the CIE standard illuminants from the CIE 15:2018 datasets
 * (downloaded August 2024). Each variant holds a zero-cost reference to its precompiled spectrum,
 * making it easy to include as a field in your own types without pulling in heavy data structures.
 *
 * This enum implements `IntoEnumIterator`, so you can **iterate through every standard illuminant**
 * (useful for testing, batch conversions, or validation).
 *
 * - Use `CieIlluminant::iter()` or `CieIlluminant::spectrum()` to list or retrieve any built-in illuminant.
 * - For a generic D-series illuminant at any correlated color temperature, use
 *   `Spectrum::cie_d_illuminant(cct: f64)`.
 *
 * By default, only **D65** and **D50** are included. To pull in the full set of fluorescent “F3_X”
 * series and other CIE illuminants, enable the `"cie-illuminants"` feature in `Cargo.toml`
 * (or build with `--features cie-illuminants`). Omit that feature (or use `--no-default-features`)
 * to keep your binary lean.
 *
 * In JavaScript/WebAssembly builds, the `colorimetry` package excludes these extra spectra by default
 * for faster load times. To include them, use the `colorimetry-all` bundle instead.
 *
 * For more background, see the Wikipedia article on
 * [Standard illuminant white points](https://en.wikipedia.org/wiki/Standard_illuminant#White_points_of_standard_illuminants).
 *
 * # Examples
 * ```rust
 * use colorimetry::illuminant::CieIlluminant;
 * use strum::IntoEnumIterator;
 *
 * // Iterate through and print all available CIE illuminants:
 * for illum in CieIlluminant::iter() {
 *     println!("{illum}");
 * }
 * ```
 * @enum {0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 | 33 | 34 | 35 | 36 | 37 | 38 | 39}
 */
export const CieIlluminant = Object.freeze({
    D65: 0, "0": "D65",
    D50: 1, "1": "D50",
    E: 2, "2": "E",
    A: 3, "3": "A",
    F1: 4, "4": "F1",
    F2: 5, "5": "F2",
    F3: 6, "6": "F3",
    F4: 7, "7": "F4",
    F5: 8, "8": "F5",
    F6: 9, "9": "F6",
    F7: 10, "10": "F7",
    F8: 11, "11": "F8",
    F9: 12, "12": "F9",
    F10: 13, "13": "F10",
    F11: 14, "14": "F11",
    F12: 15, "15": "F12",
    F3_1: 16, "16": "F3_1",
    F3_2: 17, "17": "F3_2",
    F3_3: 18, "18": "F3_3",
    F3_4: 19, "19": "F3_4",
    F3_5: 20, "20": "F3_5",
    F3_6: 21, "21": "F3_6",
    F3_7: 22, "22": "F3_7",
    F3_8: 23, "23": "F3_8",
    F3_9: 24, "24": "F3_9",
    F3_10: 25, "25": "F3_10",
    F3_11: 26, "26": "F3_11",
    F3_12: 27, "27": "F3_12",
    F3_13: 28, "28": "F3_13",
    F3_14: 29, "29": "F3_14",
    F3_15: 30, "30": "F3_15",
    LED_B1: 31, "31": "LED_B1",
    LED_B2: 32, "32": "LED_B2",
    LED_B3: 33, "33": "LED_B3",
    LED_B4: 34, "34": "LED_B4",
    LED_B5: 35, "35": "LED_B5",
    LED_BH1: 36, "36": "LED_BH1",
    LED_RGB1: 37, "37": "LED_RGB1",
    LED_V1: 38, "38": "LED_V1",
    LED_V2: 39, "39": "LED_V2",
});
/**
 * Selects a CIE standard colorimetric observer.
 *
 * The tag is embedded in every [`XYZ`] and [`Rgb`](crate::rgb::Rgb) value so that
 * operations across incompatible observers can be detected at runtime.  Each variant
 * is a lightweight index; the color-matching function tables are stored in
 * [`observer_data`].
 * @enum {0 | 1 | 2 | 3}
 */
export const Observer = Object.freeze({
    /**
     * CIE 1931 2° standard observer — the default for most colorimetry.
     *
     * Used by sRGB, ICC profiles, and CIE CRI Ra.
     */
    Cie1931: 0, "0": "Cie1931",
    /**
     * CIE 1964 10° supplementary standard observer.
     *
     * Preferred when the viewed area subtends more than ~4° at the eye.
     * Used by CIE 224:2017 / ANSI/IES TM-30 for colour fidelity calculations.
     */
    Cie1964: 1, "1": "Cie1964",
    /**
     * CIE 2015 2° observer — CMFs constructed as linear transforms of the Stockman & Sharpe
     * (2000) cone fundamentals.
     *
     * More accurate than `Cie1931` in the short-wavelength (blue) region.
     */
    Cie2015: 2, "2": "Cie2015",
    /**
     * CIE 2015 10° observer — CMFs constructed as linear transforms of the Stockman & Sharpe
     * (2000) cone fundamentals.
     *
     * Wide-field counterpart of [`Cie2015`](Observer::Cie2015).
     */
    Cie2015_10: 3, "3": "Cie2015_10",
});
/**
 * Spectrally based color space, using spectral representations of the primaries and the
 * reference white.
 *
 * Using the CIE 1931 standard observer, using a wavelength domain from 380 top 780
 * nanometer with 1 nanometer steps, these result in their usual chromaticity
 * values.  The most common _sRGB_ color space is obtained using the
 * `RgbSpace::srgb()` constructor. For this instance, the blue and green primaries
 * are direct Gaussian-filtered D65 spectra. A mixture of the blue primary and a
 * G1aussian-filtered red component is used for the red primary. Similar
 * constructors are provided for the `Adobe` and `DisplayP3` color spaces.
 *
 * The benefit of spectral primaries is that color management and color profiles
 * can use updated Colorimetric Observers, such as the Cone-Fundamental based CIE
 * 2015 observers, which don't have the CIE 1931 deficiencies. For example, they
 * can also be optimized for special observers by considering an observer's age or
 * health conditions.
 * @enum {0 | 1 | 2 | 3}
 */
export const RgbSpace = Object.freeze({
    SRGB: 0, "0": "SRGB",
    Adobe: 1, "1": "Adobe",
    DisplayP3: 2, "2": "DisplayP3",
    CieRGB: 3, "3": "CieRGB",
});

const CFIFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_cfi_free(ptr >>> 0, 1));
/**
 * Container for CIE 2017 Colour Fidelity Index (**R<sub>f</sub>**) calculations,
 * including both the general color fidelity **R<sub>f</sub>** score and the 99 special color fidelity indices (**R<sub>f,1</sub>** to **R<sub>f,99</sub>**)
 * as specified in [CIE 224:2017](https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use).
 *
 * # Requirements
 * - Requires the `cfi` feature to access color evaluation samples (CES) used for testing.
 *
 * # Overview
 * The CIE 2017 Colour Fidelity Index (CFI, or **R<sub>f</sub>**) is a modern metric for evaluating how accurately a light source renders colors.
 * It uses 99 Color Evaluation Samples (CES) that cover a broad range of real-world colors, providing a much more comprehensive assessment
 * than older metrics like the Color Rendering Index (CRI, or **R<sub>a</sub>**).
 * - The general index (**R<sub>f</sub>**) gives an overall measure of color fidelity.
 * - The special indices (**R<sub>f,1</sub>** to **R<sub>f,99</sub>**) show fidelity for each specific color sample.
 *
 * # Comparison with CRI
 * The traditional **CRI** metric (Ra) uses only 8 or 14 pastel color samples and is known to be limited, especially for modern light sources such as LEDs.
 * **CFI (Rf)** is a newer, more robust standard: it uses a much wider set of samples and is based on state-of-the-art color appearance models,
 * providing a more accurate and reliable prediction of real-world color rendering.
 * - **Use CRI** if you need compatibility with legacy systems or must comply with standards that specify CRI.
 * - **Use CFI (Rf)** for a more precise and scientifically up-to-date assessment of color fidelity, especially with modern or tunable light sources.
 *
 * # TM-30 version
 * This implementation follows **ANSI/IES TM-30-20** and **TM-30-24**, which are harmonised with
 * **CIE 224:2017**. It does **not** implement the earlier TM-30-15 or TM-30-18 editions.
 * The key difference from TM-30-15 is the scaling constant in the Rf formula: TM-30-15 used
 * `CF = 7.54`; all later editions (TM-30-18, TM-30-20, TM-30-24, CIE 224:2017) use `CF = 6.73`,
 * which this library implements.
 *
 * # Reference
 * [CIE 224:2017 – CIE 2017 Colour Fidelity Index for accurate scientific use](https://cie.co.at/publications/colour-fidelity-index-accurate-scientific-use)
 */
export class CFI {

    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(CFI.prototype);
        obj.__wbg_ptr = ptr;
        CFIFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        CFIFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_cfi_free(ptr, 0);
    }
    /**
     * General colour fidelity index Rf (0–100).
     *
     * A single overall score measuring how faithfully the test source renders the 99
     * Colour Evaluation Samples compared to the reference illuminant.
     * @returns {number}
     */
    colorFidelityIndex() {
        const ret = wasm.cfi_colorFidelityIndex(this.__wbg_ptr);
        return ret;
    }
    /**
     * General colour gamut index Rg.
     *
     * Measures the area of the gamut polygon relative to the reference (100 = same area).
     * Values above 100 indicate a wider gamut than the reference; below 100 means narrower.
     * @returns {number}
     */
    colorGamutIndex() {
        const ret = wasm.cfi_colorGamutIndex(this.__wbg_ptr);
        return ret;
    }
    /**
     * Local colour fidelity index Rf,hj for each of the 16 hue bins (TM-30 / CIE 224:2017 §4.5).
     *
     * Returns a `Float64Array` of 16 values.  Bin 0 starts at 0° (red), progressing
     * counter-clockwise in 22.5° steps around the hue circle.
     * @returns {Float64Array}
     */
    localColorFidelityIndices() {
        const ret = wasm.cfi_localColorFidelityIndices(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Chroma shift index Rcs,hj for each of the 16 hue bins (TM-30 / CIE 224:2017 §4.6).
     *
     * Returns a `Float64Array` of 16 values.  Positive means the test source boosts
     * saturation in that hue direction; negative means desaturation.  Typical range ≈ −0.5…+0.5.
     * @returns {Float64Array}
     */
    chromaShiftIndices() {
        const ret = wasm.cfi_chromaShiftIndices(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Hue shift index Rhs,hj for each of the 16 hue bins (TM-30 / CIE 224:2017 §4.7).
     *
     * Returns a `Float64Array` of 16 values in radians, wrapped to (−π, π].
     * Positive means a counter-clockwise hue shift; negative means clockwise.
     * @returns {Float64Array}
     */
    hueShiftIndices() {
        const ret = wasm.cfi_hueShiftIndices(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Special colour fidelity indices Rf,i for all 99 CES (CIE 224:2017 §7).
     *
     * Returns a `Float64Array` of 99 values, one per Colour Evaluation Sample.
     * @returns {Float64Array}
     */
    specialColorFidelityIndices() {
        const ret = wasm.cfi_specialColorFidelityIndices(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
}

const CRIFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_cri_free(ptr >>> 0, 1));
/**
 * The **Color Rendering Index (CRI)** for a light source, computed according to CIE 13.3-1995.
 *
 * This struct holds the 14 individual rendering indices R₁…R₁₄ for the standard test color samples,
 * and provides the general CRI, Rₐ, which is the average of the first eight Rᵢ values.
 *
 * # Calculation Method
 * 1. The test illuminant is scaled to 100 lx and converted to CIE XYZ under the CIE 1931 observer.
 * 2. Each of the 14 standard Colorant test spectra (TCS) is measured under both the test and the
 *    reference illuminant (black-body or D-series at the test’s correlated color temperature).
 * 3. For each sample, the color difference ΔE in CIE UVW space is computed, and
 *    Rᵢ = 100 − 4.6 · ΔE.
 * 4. The general CRI Rₐ is then
 *    ```text
 *    Rₐ = (R₁ + R₂ + … + R₈) / 8
 *    ```
 *
 * # Examples
 * ```rust
 * use colorimetry::illuminant::{Illuminant, CRI};
 *
 * // Compute CRI for the D65 illuminant:
 * let cri: CRI = (&Illuminant::d65()).try_into().unwrap();
 *
 * // General CRI:
 * let ra = cri.ra();
 * println!("General CRI Rₐ = {:.1}", ra);
 *
 * // All 14 individual Rᵢ values:
 * let ri_values = cri.to_array();
 * for (i, &ri) in ri_values.iter().enumerate() {
 *     println!("R{} = {:.1}", i + 1, ri);
 * }
 * ```
 *
 * # Notes
 * - This implementation uses the **CIE 1931** color space and requires the `"cri"` feature to be enabled in the crate.
 * - The CRI-metric is now considered somewhat outdated; newer metrics (e.g., TM-30) are recommended for modern lighting applications.
 *   However, CRI remains widely used and understood across the lighting industry.
 *
 * # Errors
 * Constructing a `CRI` can fail if the illuminant’s correlated color temperature is out of the
 * valid range (1000–25000 K) or its distance from the Planckian locus exceeds 0.05 Δuv.
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
    /**
     * Returns the general colour rendering index Rₐ (0–100),
     * the average of the first eight special rendering indices R₁…R₈.
     * @returns {number}
     */
    ra() {
        const ret = wasm.cri_ra(this.__wbg_ptr);
        return ret;
    }
}

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
     * Calculates the Color Rendering Index for this illuminant spectrum.
     *
     * Returns a `CRI` object with a `ra()` method that gives the general colour
     * rendering index Rₐ (average of R₁…R₈, scaled 0–100).
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
     * Calculates the Colour Fidelity Index (CFI / Rf) for this illuminant spectrum.
     *
     * Returns a `CFI` object exposing `rf()`, `rg()`, `rfHj()`, `rcsHj()`, `rhsHj()`,
     * and `specialIndices()`.  Requires the `cfi` feature.
     * @returns {CFI}
     */
    cfi() {
        const ret = wasm.illuminant_cfi(this.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return CFI.__wrap(ret[0]);
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
 * assert_abs_diff_eq!(rgb.to_array().as_ref(), [0.5, 0.25, 0.75].as_ref(), epsilon = 1e-6);
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
     * const [x, y, z] = xyz.to_array();
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
