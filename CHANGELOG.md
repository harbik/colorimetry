# Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

### Categories each change fall into

* **Added**: for new features.
* **Changed**: for changes in existing functionality.
* **Deprecated**: for soon-to-be removed features.
* **Removed**: for now removed features.
* **Fixed**: for any bug fixes.
* **Security**: in case of vulnerabilities.

## Unreleased

## [0.0.4] - 2025-05-01

### Added
- Add `Stimulus::new(Spectrum)` to create a stimulus from spectral values
- Add `Illuminant::cct()` method to calculate correlated color temperature and Planckian distance.
- Add `CieLab::values()` method for easy access to the CIELAB L*, a*, and b* values.
- Add `RGB::values()` method for easy access to the raw red, green and blue values.
- Enable all features when building documentation for docs.rs.
- Add `Spectrum::values() -> &[f64; 401]` for direct access to the underlying spectral data.
- Add `Illuminant::new(Spectrum)` and `From<Spectrum> for Illuminant` for easy creation of
  custom illuminants from spectral data.
- Add `Colorant::new(Spectrum)` and `TryFrom<Spectrum> for Colorant` for easy creation of
  custom colorants from spectral data.
- Add `Illuminant::xyz` to compute XYZ tristimulus values directly from an illuminant.
  convenience shorthand for calling `xyz_from_spectrum` on the observer.
- Add `Colorant::cielab` method for calculating the CIELAB values of a colorant, given an
  illuminant and an observer.

### Changed
- Changed `Stimulus::srgb` to `Stimulus::from_srgb`
- Changed `Stimulus::rgb` to `Stimulus::from_rgb`
- Changed `XYZ::try_from_luv60` to `XYZ::from_luv60`
- Changed `XYZ::try_from_chromaticity` to `XYZ::from_chromaticity`
- Changed `Triangle::try_new` to `Triangle::new`
- Changed `LineAB::try_new` to `LineAB::new`
- Changed `CCT::try_from_xyz` to `CCT::from_xyz`
- Changed `CCT::try_new_with_tint` to `CCT::new_with_tint`
- Changed `CCT::try_new` to `CCT::new`
- Take fixed size arrays instead of slices in `XYZ::new`. Makes the function signature clearer,
  and removes a possible panic case.
- Replace the public global static `XY_PRIMARIES` with the methods `name`,
  `primaries_chromaticity` and `white` on the `RgbSpace` type.

### Removed
- Remove the `paste` and `url` dependencies, since they were unused. And move the
  `colored` dependency to `dev-dependencies` since it was only used in examples.
- Remove `mul` and `mul_f64` methods from `Spectrum`. Multiplication can still be done via the
  `Mul<Spectrum>` and `Mul<f64>` trait implementations.
- Remove `Deref` and `DerefMut` implementations on `Colorant` that allowed dereference to
  `Spectrum`s. The underlying spectrum can still be obtained via the `Filter::spectrum` method.
- Remove `TryFrom<&[f64]> for Illuminant`. Instead use `Illuminant::new(Spectrum::try_from(slice))`
  or the infallible `From<Spectrum> for Illuminant` if you already have a spectrum.
- Remove `TryFrom<&[f64]> for Colorant`
- Remove all error variants in `CmtError` that were unused.

### Fixed
- Ensure the spectral values of a `Colorant` stays within the range 0.0 - 1.0. It was previously
  possible to move outside the allowed range due to the `DerefMut<Target=Spectrum>` implementation,
  and some missing clamping in other arithmetic operations on the colorant.
- Relax bounds checking in `XYZ::try_from_chromaticity` to allow the sum of `x` and `y` to be
  exactly `1.0` instead of just strictly smaller.


## [0.0.3] - 2025-04-16


## [0.0.2] - 2024-08-09


## [0.0.1] - 2024-08-09



