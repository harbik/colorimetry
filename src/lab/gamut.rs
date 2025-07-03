const CONVERGENCE_THRESHOLD: f64 = 1e-5;
use crate::lab::CieLCh;
use crate::observer::Observer;
use crate::rgb::RgbSpace;
use crate::traits::Light;
use crate::xyz::XYZ;

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

    /// Determines the maximum chroma for a given lightness (`l`) and hue (`h`).
    ///
    /// This method uses a binary search to find the maximum valid chroma value (`c`).
    /// It starts with an initial guess for chroma and iteratively narrows down the range
    /// until it finds the maximum chroma that valid color.
    /// Validity is determined by checking if the all RGB values of the resulting color
    /// are within the RGB gamut defined by the `rgb_space`.
    ///
    /// The chroma (`c`) parameter is expected to be in the range [0.0, 200.0].
    ///
    /// # Parameters
    /// - `l`: The lightness value (0.0 to 100.0).
    /// - `h`: The hue angle in degrees (0.0 to 360.0).
    ///
    /// # Returns
    /// A `CieLCh` color with the specified lightness and hue, and the maximum chroma within the gamut.
    pub fn full(&self, l: f64, h: f64) -> CieLCh {
        let mut c_low = 0.0;
        let mut c_high = 200.0;
        let mut c = 50.0; // Initial guess, can be adjusted
        for _ in 0..20 {
            // 20 iterations for convergence
            let cielch = CieLCh::new([l, c, h], self.white_point);
            let rgb = cielch.rgb(self.rgb_space);
            if rgb.values().iter().all(|&v| v < 1.0) {
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
        CieLCh::new([l, c_low, h], self.white_point)
    }

    pub fn rgb_gamut(&self, l: f64, h: f64) -> Option<CieLCh> {
        let cielch = self.full(l, h);
        let rgb = cielch.rgb(self.rgb_space);
        // Check if RGB values are within the gamut explicitly
        if rgb.values().iter().all(|&v| (0.0..=1.0).contains(&v)) {
            // If all RGB values are within the gamut, return the CieLCh color
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
