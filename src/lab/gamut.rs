// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

use crate::lab::CieLCh;
use crate::observer::Observer;
use crate::rgb::RgbSpace;
use crate::traits::Light;
use crate::xyz::XYZ;

const CONVERGENCE_THRESHOLD: f64 = 1e-5;
const CHROMA_HIGH_MAX: f64 = 500.0;

/// Number of bisection steps needed to narrow `[0, CHROMA_HIGH_MAX]` down to
/// `CONVERGENCE_THRESHOLD`: 500 / 2^26 is about 7.5e-6.
const MAX_ITERATIONS: usize = 26;

pub struct CieLChGamut {
    rgb_space: RgbSpace,
    white_point: XYZ,
}

impl CieLChGamut {
    pub fn new(observer: Observer, rgb_space: RgbSpace) -> Self {
        let white_point = rgb_space.white().white_point(observer);
        CieLChGamut {
            white_point,
            rgb_space,
        }
    }

    pub fn oberver(&self) -> Observer {
        self.white_point.observer()
    }

    /// Determines the maximum chroma for a given lightness (`l`) and hue (`h`).
    ///
    /// This method uses a binary search to find the maximum valid chroma value (`c`).
    /// It starts with an initial guess for chroma and iteratively narrows down the range
    /// until it finds the maximum chroma that valid color.
    /// Validity is determined by checking if the all RGB values of the resulting color
    /// are within the RGB gamut defined by the `rgb_space`, which means all three
    /// channels are in the range [0.0, 1.0].
    ///
    /// The chroma (`c`) parameter is expected to be in the range [0.0, 500.0].
    ///
    /// # Parameters
    /// - `l`: The lightness value (0.0 to 100.0).
    /// - `h`: The hue angle in degrees (0.0 to 360.0).
    ///
    /// # Returns
    /// A `CieLCh` color with the specified lightness and hue, and the maximum chroma within the
    /// gamut, or `None` if no color at this lightness and hue is realizable at all — either because
    /// the lightness itself lies outside the RGB gamut, or because the result falls outside the
    /// spectral locus.
    ///
    /// # Notes
    ///
    /// - This method checks the validity of a CieLCh color by converting it to RGB and ensuring all
    ///   RGB values are in the range [0.0, 1.0], and that it is located within the spectral locus area.
    ///   A channel value above 1.0 or below 0.0 both mean the color lies outside the RGB gamut of the
    ///   color space.
    /// - The method performs a binary search for chroma, starting from 0.0 to `CHROMA_HIGH_MAX`.
    ///   The value of `CHROMA_HIGH_MAX` (500.0) is chosen as a practical upper limit based on
    ///   empirical observations and theoretical considerations of typical chroma ranges in color spaces.
    /// - It iterates up to `MAX_ITERATIONS` times, which is enough for the bisection to reach
    ///   `CONVERGENCE_THRESHOLD` over the full `CHROMA_HIGH_MAX` range.
    /// - The bisection assumes that the in-gamut chroma values at a given lightness and hue form a
    ///   single interval starting at `c = 0`. That holds almost everywhere, but not quite
    ///   everywhere: because a line of constant lightness and hue in CIELAB is a curve in XYZ, it
    ///   can leave the (convex) RGB gamut and re-enter it near a corner of the RGB cube. Where that
    ///   happens the first crossing is returned rather than the outermost one, and the chroma is
    ///   under-reported. For sRGB this affects roughly 0.01% of the lightness/hue plane, confined
    ///   to L\* above about 96 near the yellow corner.
    pub fn max_chroma(&self, l: f64, h: f64) -> Option<CieLCh> {
        let mut c_low = 0.0;
        let mut c_high = CHROMA_HIGH_MAX;

        // The bisection below only ever raises `c_low` to a chroma it has found to be in gamut, so
        // it can only return an in-gamut color if its starting point is one. That is not a given:
        // for an RGB space whose primaries do not reproduce its own white point, the neutral axis
        // leaves the gamut above some lightness, and no chroma at that lightness is realizable.
        if !CieLCh::new([l, c_low, h], self.white_point)
            .rgb(self.rgb_space)
            .is_in_gamut()
        {
            return None;
        }
        let mut c = CHROMA_HIGH_MAX / 2.0; // Initial guess for chroma
        for _ in 0..MAX_ITERATIONS {
            let cielch = CieLCh::new([l, c, h], self.white_point);
            let rgb = cielch.rgb(self.rgb_space);
            if rgb.is_in_gamut() {
                c_low = c; // Found a valid chroma, increase lower bound
            } else {
                c_high = c; // Not in gamut, decrease upper bound
            }
            if (c_high - c_low).abs() < CONVERGENCE_THRESHOLD {
                // Convergence threshold
                break;
            }
            c = (c_low + c_high) / 2.0; // Update guess
        }

        // Use c_low as the chroma value because it represents the largest chroma
        // that is still within the RGB gamut after the binary search.
        let cielch = CieLCh::new([l, c_low, h], self.white_point);

        // Ensure the resulting CieLCh color is within the spectral locus area
        let xy = cielch.rxyz().xyz().chromaticity();
        if cielch.observer().spectral_locus().contains(xy.to_array()) {
            Some(cielch)
        } else {
            None
        }
    }

    /// Determines the maximum chroma for a given lightness (`l`) and hue (`h`) that is within the RGB gamut.
    /// This method is similar to `max_chroma`, but it checks if the resulting RGB values are within the gamut explicitly.
    /// # Parameters
    /// - `l`: The lightness value (0.0 to 100.0).
    /// - `h`: The hue angle in degrees (0.0 to 360.0).
    /// # Returns
    /// An `Option<CieLCh>` color with the specified lightness and hue, and the maximum chroma that is within the RGB gamut.
    /// # Notes
    /// - This method checks the validity of a CieLCh color by converting it to RGB and ensuring all RGB values are in the range [0.0, 1.0].
    pub fn max_chroma_in_gamut(&self, l: f64, h: f64) -> Option<CieLCh> {
        let cielch = self.max_chroma(l, h)?;
        let rgb = cielch.rgb(self.rgb_space);
        // Check if RGB values are within the gamut explicitly
        if rgb.to_array().iter().all(|&v| (0.0..=1.0).contains(&v)) {
            Some(cielch)
        } else {
            None
        }
    }

    pub fn rgb_space(&self) -> RgbSpace {
        self.rgb_space
    }

    pub fn white_point(&self) -> XYZ {
        self.white_point
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lab::CieLCh;
    use crate::rgb::WideRgb;
    use crate::xyz::RelXYZ;

    /// Every lightness and hue inside the RGB gamut has a realizable maximum chroma, and the
    /// returned color must itself be in gamut.
    ///
    /// Regression test: the in-gamut check used to reject only channel values above 1.0, never
    /// negative ones. The search then ran away to `CHROMA_HIGH_MAX` for any hue where no channel
    /// exceeds 1.0 before the spectral locus is left, and the locus check turned that into `None`
    /// — for every hue below L\* 30, among others.
    #[test]
    fn max_chroma_is_some_and_in_gamut_over_the_whole_lightness_range() {
        for space in [
            RgbSpace::SRGB,
            RgbSpace::Adobe,
            RgbSpace::DisplayP3,
            RgbSpace::CieRGB,
        ] {
            let gamut = CieLChGamut::new(Observer::Cie1931, space);
            for li in 1..100 {
                let l = li as f64;
                for hi in 0..72 {
                    let h = hi as f64 * 5.0;
                    // `CieRGB`'s primaries do not reproduce its own white point, so its neutral
                    // axis leaves the gamut above L* 90 and nothing at that lightness is
                    // realizable. Everywhere the neutral axis is in gamut, a maximum chroma exists.
                    let neutral_in_gamut = CieLCh::new([l, 0.0, h], gamut.white_point())
                        .rgb(space)
                        .is_in_gamut();
                    let Some(lch) = gamut.max_chroma(l, h) else {
                        assert!(
                            !neutral_in_gamut,
                            "{space:?}: no max chroma at L*={l}, h={h}"
                        );
                        continue;
                    };
                    assert!(neutral_in_gamut);
                    assert!(lch.c() > 0.0, "{space:?}: zero chroma at L*={l}, h={h}");
                    assert!(
                        lch.rgb(space).is_in_gamut(),
                        "{space:?}: out of gamut at L*={l}, h={h}: {:?}",
                        lch.rgb(space).to_array()
                    );
                }
            }
        }
    }

    /// The chroma just outside the reported maximum must be out of gamut, so the value returned is
    /// the boundary and not merely some in-gamut chroma.
    #[test]
    fn max_chroma_is_on_the_gamut_boundary() {
        let gamut = CieLChGamut::new(Observer::Cie1931, RgbSpace::SRGB);
        for li in 1..100 {
            let l = li as f64;
            for hi in 0..72 {
                let h = hi as f64 * 5.0;
                let c = gamut.max_chroma(l, h).unwrap().c();
                let outside = CieLCh::new([l, c + 0.01, h], gamut.white_point());
                assert!(
                    !outside.rgb(RgbSpace::SRGB).is_in_gamut(),
                    "L*={l}, h={h}: c={c} is not the boundary"
                );
            }
        }
    }

    /// The saturated sRGB primaries and secondaries lie exactly on the gamut surface, so the
    /// maximum chroma at their own lightness and hue must be their own chroma.
    ///
    /// Yellow is excluded: it sits in the re-entrant region documented on `max_chroma`, where the
    /// bisection returns the first gamut crossing instead of the outermost one.
    #[test]
    fn max_chroma_reproduces_the_srgb_primaries() {
        let gamut = CieLChGamut::new(Observer::Cie1931, RgbSpace::SRGB);
        for rgb in [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
        ] {
            let xyz = WideRgb::new(
                rgb[0],
                rgb[1],
                rgb[2],
                Some(Observer::Cie1931),
                Some(RgbSpace::SRGB),
            )
            .xyz();
            let want = CieLCh::from_xyz(RelXYZ::new(xyz.into(), gamut.white_point()));
            let got = gamut.max_chroma(want.l(), want.h()).unwrap();
            assert!(
                approx::abs_diff_eq!(got.c(), want.c(), epsilon = 1e-3),
                "rgb {rgb:?}: got {}, want {}",
                got.c(),
                want.c()
            );
        }
    }

    /// `max_chroma` only ever returns colors that are already in gamut, so `max_chroma_in_gamut`
    /// agrees with it everywhere.
    #[test]
    fn max_chroma_in_gamut_agrees_with_max_chroma() {
        let gamut = CieLChGamut::new(Observer::Cie1931, RgbSpace::SRGB);
        for li in 1..100 {
            let l = li as f64;
            for hi in 0..24 {
                let h = hi as f64 * 15.0;
                assert_eq!(gamut.max_chroma(l, h), gamut.max_chroma_in_gamut(l, h));
            }
        }
    }
}
