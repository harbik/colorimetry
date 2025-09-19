// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

use crate::lab::CieLCh;
use crate::observer::Observer;
use crate::rgb::RgbSpace;
use crate::traits::Light;
use crate::xyz::XYZ;

const CONVERGENCE_THRESHOLD: f64 = 1e-5;
const CHROMA_HIGH_MAX: f64 = 500.0;

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
    /// are within the RGB gamut defined by the `rgb_space`.
    ///
    /// The chroma (`c`) parameter is expected to be in the range [0.0, 500.0].
    ///
    /// # Parameters
    /// - `l`: The lightness value (0.0 to 100.0).
    /// - `h`: The hue angle in degrees (0.0 to 360.0).
    ///
    /// # Returns
    /// A `CieLCh` color with the specified lightness and hue, and the maximum chroma within the gamut.
    ///
    /// # Notes
    ///
    /// - This method checks the validity of a CieLCh color by converting it to RGB and ensuring all RGB values are less than 1.0 and is located within the spectral locus area.
    /// - Negative RGB values indicate that the CieLCh color lies outside the RGB gamut of the color space, but are further constrained by a spectral boundary check.
    /// - The method performs a binary search for chroma, starting from 0.0 to `CHROMA_HIGH_MAX`.
    ///   The value of `CHROMA_HIGH_MAX` (500.0) is chosen as a practical upper limit based on
    ///   empirical observations and theoretical considerations of typical chroma ranges in color spaces.
    /// - It iterates up to 20 times to ensure convergence. This limit was chosen as a balance between computational efficiency and convergence accuracy. It can be adjusted if higher precision or faster computation is required.
    /// - The method performs a binary search for chroma, starting from 0.0 to `CHROMA_HIGH_MAX`.
    pub fn max_chroma(&self, l: f64, h: f64) -> Option<CieLCh> {
        let mut c_low = 0.0;
        let mut c_high = CHROMA_HIGH_MAX;
        let mut c = CHROMA_HIGH_MAX / 2.0; // Initial guess for chroma
        for _ in 0..20 {
            // 20 iterations for convergence
            let cielch = CieLCh::new([l, c, h], self.white_point);
            let rgb = cielch.rgb(self.rgb_space);
            if rgb.to_array().iter().all(|&v| v < 1.0) {
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
mod tests {}
