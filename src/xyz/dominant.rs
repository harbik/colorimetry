use crate::{
    error::Error,
    math::{LineAB, Orientation},
};

use super::XYZ;

impl XYZ {
    /// The Dominant Wavelength of a color point is the wavelength of spectral
    /// color, obtained from the intersection of a line through a white point
    /// and itself, with the spectral locus.  Points on this line were
    /// historically thought off as having the same hue, but that has been
    /// proven wrong.  This value has limited practical use, but is sometimes
    /// used for color definition of LEDs.
    ///
    /// The spectral locus, being the boundary of all possible colors in the CIE
    /// 1931 diagram, collapses to one point beyond a wavelength of 699nm. As a
    /// result, the maxium range of dominant wavelengths which can be obtained
    /// is from 380 to 699 nanometer;
    ///
    pub fn dominant_wavelength(&self, white: XYZ) -> Result<f64, Error> {
        let mut sign = 1.0;
        let wavelength_range = self.observer.data().spectral_locus_wavelength_range();
        let mut low = *wavelength_range.start();
        let mut high = *wavelength_range.end();
        let mut mid = 540usize; // 200 fails, as its tail overlaps into the blue region
        if white.observer != self.observer {
            Err(Error::RequireSameObserver)
        } else {
            let chromaticity = self.chromaticity();
            let [mut x, mut y] = [chromaticity.x(), chromaticity.y()];
            let white_chromaticity = white.chromaticity();
            // if color point is in the purple rotate it around the white point by 180º, and give wavelength a negative value
            let blue_edge = LineAB::new(
                white_chromaticity.to_array(),
                self.observer
                    .data()
                    .xyz_at_wavelength(low)
                    .unwrap()
                    .chromaticity()
                    .to_array(),
            )
            .unwrap();
            let red_edge = LineAB::new(
                white_chromaticity.to_array(),
                self.observer
                    .data()
                    .xyz_at_wavelength(high)
                    .unwrap()
                    .chromaticity()
                    .to_array(),
            )
            .unwrap();
            match (blue_edge.orientation(x, y), red_edge.orientation(x, y)) {
                (Orientation::Colinear, _) => return Ok(380.0),
                (_, Orientation::Colinear) => return Ok(699.0),
                (Orientation::Left, Orientation::Right) => {
                    // mirror point into non-purple region
                    sign = -1.0;
                    x = 2.0 * white_chromaticity.x() - x;
                    y = 2.0 * white_chromaticity.y() - y;
                }
                _ => {} // do nothing
            }
            // start bisectional search
            while high - low > 1 {
                let bisect = LineAB::new(
                    white_chromaticity.to_array(),
                    self.observer
                        .data()
                        .xyz_at_wavelength(mid)
                        .unwrap()
                        .chromaticity()
                        .to_array(),
                )
                .unwrap();
                //   let a = bisect.angle_deg();
                match bisect.orientation(x, y) {
                    Orientation::Left => high = mid,
                    Orientation::Right => low = mid,
                    Orientation::Colinear => {
                        low = mid;
                        high = mid;
                    }
                }
                mid = (low + high) / 2;
            }
            if low == high {
                Ok(sign * low as f64)
            } else {
                let low_ab = LineAB::new(
                    white.chromaticity().to_array(),
                    self.observer
                        .data()
                        .xyz_at_wavelength(low)
                        .unwrap()
                        .chromaticity()
                        .to_array(),
                )
                .unwrap();
                let dlow = low_ab.distance_with_sign(x, y);
                let high_ab = LineAB::new(
                    white.chromaticity().to_array(),
                    self.observer
                        .data()
                        .xyz_at_wavelength(high)
                        .unwrap()
                        .chromaticity()
                        .to_array(),
                )
                .unwrap();
                let dhigh = high_ab.distance_with_sign(x, y);
                if dlow < 0.0 || dhigh > 0.0 {
                    // not ended up between two lines
                    let s = format!("bisection error in dominant wavelength search:  {dlow} {low} {dhigh} {high}");
                    return Err(Error::ErrorString(s));
                }
                let dl = (dlow.abs() * high as f64 + dhigh.abs() * low as f64)
                    / (dlow.abs() + dhigh.abs());
                Ok(sign * dl)
            }
        }
    }
}

#[cfg(test)]
mod xyz_test {
    use crate::{math::LineAB, prelude::*};
    use approx::assert_ulps_eq;

    #[test]
    fn dominant_wavelength_test() {
        let d65 = CIE1931.xyz_d65().set_illuminance(50.0);

        // 550 nm
        let sl = CIE1931
            .xyz_at_wavelength(550)
            .unwrap()
            .set_illuminance(50.0);
        let t = d65.try_add(sl).unwrap();
        let dl = t.dominant_wavelength(d65).unwrap();
        assert_ulps_eq!(dl, 550.0);

        for wl in 380..=699usize {
            let sl2 = CIE1931.xyz_at_wavelength(wl).unwrap();
            //let [slx, sly] = sl2.chromaticity();
            //println!("sl xy: {slx} {sly}");
            let dl = sl2.dominant_wavelength(d65).unwrap();
            assert_ulps_eq!(dl, wl as f64, epsilon = 1E-10);
        }
    }

    #[test]
    fn dominant_wavelength_purple_test() {
        let d65 = CIE1931.xyz_d65();
        let white_chromaticity = d65.chromaticity();

        // get purple line
        let xyzb = CIE1931.xyz_at_wavelength(380).unwrap();
        let [xb, yb] = xyzb.chromaticity().to_array();
        let xyzr = CIE1931.xyz_at_wavelength(699).unwrap();
        let [xr, yr] = xyzr.chromaticity().to_array();
        let line_t = LineAB::new([xb, yb], [xr, yr]).unwrap();
        for wl in 380..=699usize {
            let sl = CIE1931.xyz_at_wavelength(wl).unwrap();
            let chromaticity = sl.chromaticity();
            let line_u =
                LineAB::new(chromaticity.to_array(), white_chromaticity.to_array()).unwrap();
            let ([_xi, yi], t, _) = line_t.intersect(&line_u).unwrap();
            if t > 0.0 && t < 1.0 {
                // see https://en.wikipedia.org/wiki/CIE_1931_color_space#Mixing_colors_specified_with_the_CIE_xy_chromaticity_diagram
                let b = xyzb.set_illuminance(100.0 * (yb * (yr - yi)));
                let r = xyzr.set_illuminance(100.0 * (yr * (yi - yb)));
                let s = b.try_add(r).unwrap();
                let dl = s.dominant_wavelength(d65).unwrap();
                assert_ulps_eq!(dl, -(wl as f64));
            }
        }
    }
}
