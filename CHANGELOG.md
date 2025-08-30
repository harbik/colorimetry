# Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## Categories each change fall into

* **Added**: for new features.
* **Changed**: for changes in existing functionality.
* **Deprecated**: for soon-to-be removed features.
* **Removed**: for now removed features.
* **Fixed**: for any bug fixes.
* **Security**: in case of vulnerabilities.

## Unreleased

### Added

* `colorimetry-cli` as a member package, in the `cli` folder (see `Colorimetry-cli` crate), implemening the `color` binary tool
* `colorimetry-plot` as a member package, in the `plot` folder (see `Colorimetry-plot` crate)
* `CieLChGamut`, to calculate the Chroma limit for a particular Lightness and hue value.
* `GammaCurve` to API, and add zero argument vector to indicate color spaces without a gamma correction (such as CieRGB), which only use floating point representations
* `CIE RGB` observer, using the CIE E Illuminant
* `CIE E` Standard Illuminant
* `Spectrum::from_wavelength_map` to create a spectrum from a set of wavelength and value pairs.
* `cielab-hues` and `cielab-areas` examples for plotting CieLab iso-hue lines and lightness lines in xy-chromaticity diagrams
* `WideRgb::is_within_gamut` with tolerance parameter
* `CieLab::is_black` to test if a CieLab value is black
* `OptimalColors`, `RelXYZGamut to build gamut hashmaps for the envelopes of the XYZ and CieLab color spaces.
* `xtask` package as a sub-package, as a build and check tool, to use:

  * `cargo xtask gen` to generate source files for `rgb-spaces` and `rgb-transforms`
  * `cargo xtask check` to check the files using `clippy`, `fmt`, and `rdme`,
  * `cargo xtask test` to run test with various feature settings,
  * `cargo xtask doc` to generate the documentation, failing on warnings, and
  * `cargo xtask wasm` to generate the web-assembly files in `pkg` (requires `wasm-pack` and `wasm-opt`)

### Changed

* Rename `Colorant::from_xyz` to `from_rxyz`
* Rename `CieLab::xyz` to `rxyz`
* Rename `CieLab::xyzn` to `white_point`
* Rename `Light::xyz100` to `white_point`
* Rename `XYZ::from_vecs` to `XYZ::from_vec`
* RGB <-> XYZ conversion matrices for the various color spaces hard coded through `xtask gen` `gen_rgbxyz.rs` function, using a handlebars template

### Fixed

* Fix `Observer::calc_xyz2rgb` to use `calc_rgb2xyz_matrix` to bootstrap transformation matrices, instead of using the xtask generated version

### Removed

* `supplemental-observers` feature, all are now always included.

## [0.0.7] - 2025-6-17

### Added [0.0.7]

* `Gaussian` struct, to represent normal distributions as used in this library as spectral
* distributions and used for filtering.
* `Planck` struct, with methods for calculating blackbody emission spectra.
* Add `CieCam02` implementation, with same methods as CieCam16.
* Add `CFI` Color Fidelity Index calculations with `general_color_fidelity_index` and `special_color_fidelity_indices`.
* Add standard CIECAM viewing conditions as recommended by CIE248 as constants.
* A Spectral RGB matching use case to convert spectral reflectance into sRGB values using the CIE 2015 10° observer.
* New helper functions in the `cam` module (and submodules `cam02`, `cam16`, `viewconditions`) for spectral color appearance calculations and surround presets.
* A spectral color management use case, demonstrating matching the Munsell paint chip _5 BG 5/8_ to its sRGB equivalent (LED_B2 illuminant, CIE 2015 10° observer).
* Add `WideRgb::to_rgb` that returns an `Rgb` instance with the same channel values, if the wide
  RGB was not out-of-gamut.
* Add `WideRgb::is_in_gamut` to check if the RGB values are within the RGB gamut.
* Add `RelXYZ` to represent related colors, sharing a single white reference

### Changed [0.0.7]

* Refactored `physics.rs` with the planck functions and LED functions now moved to the `Illuminant` module,
  and the remaining functions combined with teh `geometry.rs` functions into `math.rs`.
* Changed `stefan_boltzmann` method name to `total_radiance`.
* Renamed `ObserverData::xyz_fn_illuminant` to `xyz_from_fn`.
* Renamed `ObserverData::xyz_from_std_illuminant_x_fn` to `xyz_from_colorant_fn`.
* Replaced all instances of `&ObserverData` in the public API with `Observer`, which is an `enum`, to avoid confusion between the use of `ObserverData` and `Observer`.
* Renamed `CieCam16::ciede2016` to `de_ucs`.
* Replaced all instances of `&RgbSpace` in the API with `RgbSpace`, as both were used inconsistently.
* Updated inline code comments for clarity in color conversion routines and view‐condition constructors.
* Refined README section layout and wording to better guide users through the new spectral examples.
* Mark all enums that add variants when cargo features are enabled as `#[non_exhaustive]`.
  This includes `Observer` and `CieIlluminant`.
* Mark all enums that might get new variants in future non-API breaking releases as
  `#[non_exhaustive]`. This includes `RgbSpace`, `Cam` and `Error`.
* Renamed `CES_DATA` to `CES`
* `CieCam02`, `CieCam03` and `CieLab` constructors now taking a single `RelXYZ` argument, instead of two `XYZ` values.

### Removed [0.0.7]

* Various normal distribution (gaussian) helper functions, now all collected as methods of the `Gaussian` struct.
* Various Planck's law helper functions, now all collected as methods of the `Planck` struct.

## [0.0.6] - 2025-05-27

### Added [0.0.6]

* Add `r()`, `g()` and `b()` methods to `WideRgb` for easy access to each channel value.
* Add `Chromaticity` struct and use it instead of `[f64; 2]` to represent chromaticity coordinates.
* Add `jch` method to `CieCam16` to get JCh values as a `[f64;3]`-array.
* Add `ciede2016` method to `CieCam16` to get the CIECAM16-UCS color difference between two `CieCam16` values.
* Add `cam` and `CieCam16` documentation.
* Add `ciede2000` method to `CieLab`.
* Add `AsRef<Spectrum> for Illuminant`, to replace the implicit `Deref`.

### Removed [0.0.6]

* Remove undocumented `XYZ::srgb` method that both clamped out-of-gamut values and converted
  directly to a gamma encoded `[u8; 3]`. Obtain the same result with the more explicit
  `xyz.rgb(Some(RgbSpace::SRGB)).clamp().values()`.
* Remove conversion directly from `WideRgb` to clamped and gamma encoded `[u8; 3]`.
  Prefer being more explicit by converting to the `Rgb` type in between with one of the provided conversion methods.
* Remove `DeRef<Spectrum> for Illuminant` to be replaced with the more explcicit `AsRef<Spectrum>`.
* Remove `XYZ_D65`, using `xyz_d65` instead.

### Changed [0.0.6]

* `XYZ::new` now only uses a single array of tristimulus values; the white reference value was dropped in the XYZ representation.
* Make indexing into a `Spectrum` with out of bounds wavelengths cause a panic, instead of returning `NaN`s or modifying the first or last values.
* Exposed `CieCam16::jab` methods, which was private before.
* Renamed private `CieCam16::jabp` to `CieCam16::jab_prime` method, and made it public.
* Renamed private `CieCam16::jchp` to `CieCam16::jch_prime` method, and made it public.
* Made `CieCam16` conversion matrices private, as being implementation details.
* Added expclicit reference white tristimulus values to the `CieLab` constructor, which replaces the reference white tristimulus values which were previously included in `XYZ`.
* Renamed `delta_e` to `ciede` for CieLab, to align with the common name for this in Colorimetry.
* Renamed the `CieLab::new` method `CieLab::from_xyz`, as it takes `XYZ` values as arguments.
* `CieLab::new` takes now a CIE L*a*b* [f64;3]-array and a reference white `XYZ` value.
* Renamed the `CieCam16::new` method `CieCam16::from_xyz`, as it takes `XYZ` values as arguments.
* `CieCam16` new takes now an _JCh_ `[f64;3]`-array, a reference white `XYZ` value, and a `ViewConditions` instance.
* Changed library module structure combining related modules in sub-modules, make the structure less flat, and re-export of most sub-module entities.
* Renamed `StdIlluminant`-`enum` to `CieIlluminant`
* Rename `CmtError` to `Error`.

### Fixed [0.0.6]

* Fix bug `XYZ::set_illuminance`. Avoid divide be zero when luminous values is zero, or negative.
* Fix bug in `WideRgb::compress`. Previously `WideRgb` instances with only positive channel values would have its lowest channel value invalidly scaled down to 0.0. And instances with only values below 1.0 would have its highest channel value invalidly scaled up to 1.0.

## [0.0.5] - 2025-05-14

### Added [0.0.5]

* `WideRgb::clamp` and `WideRgb::compress` methods to create `Rgb` from `WideRgb` values, to pull out-of-gamut values within the colorspace,
* `WideRgb` type allowing unconstrained, out-of-gamut RGB values
* Implement strums `EnumIter` on `Observer`. Allows easy iteration over all available observers.
* Add `r()`, `g()` and `b()` methods to `RGB` for easy access to each channel value.
* Add `x()`, `y()` and `z()` methods to `XYZ` for easy access to each channel value.

### Changed [0.0.5]

* Constrain `Rgb` type to in-gamut values only, i.e. all R, G, and B values are required to be the
  range of [0..=1.0].
* Renamed `RGB` type to `Rgb`
* Stop normalizing XYZ values to illuminance = 100 in `XYZ::values()`.
* Replace `spectral_locus_nm_min` and `spectral_locus_nm_max` with a single
  `spectral_locus_wavelength_range` method that return both values as a typed range.
* Rename `Observer::spectral_locus_by_nm` to `xyz_at_wavelength` and relax the constraints on
  the allowed wavelength range. This method can now sample the color matching functions in the
  full range from 380 - 780 nm.

### Removed [0.0.5]

* Make `spectral_locus_index_min`, `spectral_locus_index_max` and `spectral_locus_by_index`
  private.
* Remove `RGB::from_xyz` method, which requires XYZ values to be in the range from 0.0 to 1.0;
    use `XYZ::rgb` instead, as that uses the reference illuminance for scaling.

### Fixed [0.0.5]

* Fix `CIE2015_10` data error (X and Y CMF's were identical)
* Fix caching bug in `Observer::spectral_locus_index_min` and
  `Observer::spectral_locus_index_max`. Previously only the values for the first
  observer it was called on was returned for all subsequent calls.
* Fix caching bug in `RgbSpaceData::primaries_as_colorants`. Previously only the values
  for the first colorspace it was called on was returned for all subsequent calls.
* Fix `Observer::spectral_locus_index_min` and `Observer::spectral_locus_index_max` to
  not panic for the `Std2015` observer.
* Slightly relax constraints on x-y values in `XYZ::from_chromaticity`. This allows converting
  the spectral locus positions of all included observers to chromaticity coordinates and back
  to XYZ values again.
* Fix caching bug in `Observer::rgb2xyz` and `Observer::xyz2rgb`. If multiple observers are used,
  only the computed matrixes for the first one to call into these methods would be returned in
  subsequent invocations.

## [0.0.4] - 2025-05-06

### Added [0.0.4]

* Add `Stimulus::new(Spectrum)` to create a stimulus from spectral values
* Add `Illuminant::cct()` method to calculate correlated color temperature and Planckian distance.
* Add `CieLab::values()` method for easy access to the CIELAB L*, a*, and b* values.
* Add `RGB::values()` method for easy access to the raw red, green and blue values.
* Enable all features when building documentation for docs.rs.
* Add `Spectrum::values() -> &[f64; 401]` for direct access to the underlying spectral data.
* Add `Illuminant::new(Spectrum)` and `From<Spectrum> for Illuminant` for easy creation of
  custom illuminants from spectral data.
* Add `Colorant::new(Spectrum)` and `TryFrom<Spectrum> for Colorant` for easy creation of
  custom colorants from spectral data.
* Add `Illuminant::xyz` to compute XYZ tristimulus values directly from an illuminant.
  convenience shorthand for calling `xyz_from_spectrum` on the observer.
* Add `Colorant::cielab` method for calculating the CIELAB values of a colorant, given an
  illuminant and an observer.

### Changed [0.0.4]

* Changed `Stimulus::srgb` to `Stimulus::from_srgb`
* Changed `Stimulus::rgb` to `Stimulus::from_rgb`
* Changed `XYZ::try_from_luv60` to `XYZ::from_luv60`
* Changed `XYZ::try_from_chromaticity` to `XYZ::from_chromaticity`
* Changed `Triangle::try_new` to `Triangle::new`
* Changed `LineAB::try_new` to `LineAB::new`
* Changed `CCT::try_from_xyz` to `CCT::from_xyz`
* Changed `CCT::try_new_with_tint` to `CCT::new_with_tint`
* Changed `CCT::try_new` to `CCT::new`
* Take fixed size arrays instead of slices in `XYZ::new`. Makes the function signature clearer,
  and removes a possible panic case.
* Replace the public global static `XY_PRIMARIES` with the methods `name`,
  `primaries_chromaticity` and `white` on the `RgbSpace` type.

### Removed [0.0.4]

* Remove the `paste` and `url` dependencies, since they were unused. And move the
  `colored` dependency to `dev-dependencies` since it was only used in examples.
* Remove `mul` and `mul_f64` methods from `Spectrum`. Multiplication can still be done via the
  `Mul<Spectrum>` and `Mul<f64>` trait implementations.
* Remove `Deref` and `DerefMut` implementations on `Colorant` that allowed dereference to
  `Spectrum`s. The underlying spectrum can still be obtained via the `Filter::spectrum` method.
* Remove `TryFrom<&[f64]> for Illuminant`. Instead use `Illuminant::new(Spectrum::try_from(slice))`
  or the infallible `From<Spectrum> for Illuminant` if you already have a spectrum.
* Remove `TryFrom<&[f64]> for Colorant`
* Remove all error variants in `CmtError` that were unused.

### Fixed [0.0.4]

* Ensure the spectral values of a `Colorant` stays within the range 0.0 - 1.0. It was previously
  possible to move outside the allowed range due to the `DerefMut<Target=Spectrum>` implementation,
  and some missing clamping in other arithmetic operations on the colorant.
* Relax bounds checking in `XYZ::try_from_chromaticity` to allow the sum of `x` and `y` to be
  exactly `1.0` instead of just strictly smaller.

## [0.0.3] - 2025-04-16

## [0.0.2] - 2024-08-09

## [0.0.1] - 2024-08-09
