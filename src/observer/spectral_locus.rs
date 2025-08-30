// SPDX-License-Identifier: Apache-2.0 OR MIT
// Copyright (c) 2024-2025, Harbers Bik LLC

//! Spectral Locus module provides functionality to work with the spectral locus in colorimetry.
//!
//! The spectral locus represents the boundary of visible colors in a chromaticity diagram,
//! formed by plotting monochromatic light sources of different wavelengths.
//!
//! The module provides:
//! - `SpectralLocus`: Main struct representing the spectral locus
//! - `SpectralLocusIterator`: Iterator for traversing spectral locus points
//!
//! The spectral locus is represented as a closed polygon, where each point corresponds
//! to the chromaticity coordinates of monochromatic light at different wavelengths.

use std::ops::RangeBounds;

use geo::{Contains, Polygon};

use crate::{
    math::QuadraticCurve,
    observer::Observer,
    spectrum::{NS, SPECTRUM_WAVELENGTH_RANGE},
    xyz::XYZ,
};

/// Represents the spectral locus, which includes:
/// - An observer for spectral data.
/// - A vector of tuples where each tuple contains:
///   - A wavelength (usize).
///   - A chromaticity coordinate ([f64; 2]).
pub struct SpectralLocus {
    observer: Observer,
    polygon: Polygon,
}

impl SpectralLocus {
    pub fn new(observer: Observer) -> Self {
        let obs_data = observer.data().data;
        let mut v = Vec::with_capacity(NS + 1);
        for i in 0..NS {
            let xyz = obs_data.column(i).into();
            let chromaticity = XYZ::from_vec(xyz, observer).chromaticity();
            v.push(chromaticity.to_array());
        }
        v.push(v[0]); // close the polygon
        let polygon = vec_to_polygon(v);
        SpectralLocus { observer, polygon }
    }

    pub fn contains(&self, point: [f64; 2]) -> bool {
        self.polygon.contains(&geo::Point::new(point[0], point[1]))
    }

    pub fn observer(&self) -> Observer {
        self.observer
    }

    /// Iterate over the spectral locus points over a given wavelength range and step size producing
    /// the spectral locus coordinates (f64, f64) and slope angles f64 in radians.
    /// The slope is calculated using a quadratic curve fit for interior points, and a simple
    /// difference for the endpoints.
    ///
    /// This iterator is useful for applications that require not only the locus points,
    /// but also the direction of the locus at each point, for example when plotting tick marks along the spectral locus
    /// in plots, to find its boundaries.
    ///
    /// # Item
    /// Each iteration yields a tuple:
    /// - `((f64, f64), f64)`
    ///   - The chromaticity coordinate (x, y) of the locus point.
    ///   - The slope angle (in radians) at that point.
    pub fn iter_range_with_slope<'a>(
        &'a self,
        range: impl RangeBounds<usize>,
        step: usize,
    ) -> SpectralLocusSlopeIterator<'a> {
        let min_wavelength = *SPECTRUM_WAVELENGTH_RANGE.start();
        let clip = |s, lim| {
            if s > min_wavelength {
                s - min_wavelength
            } else {
                lim
            }
        };
        let start = match range.start_bound() {
            std::ops::Bound::Included(&s) => clip(s, 0),
            std::ops::Bound::Excluded(&s) => clip(s + 1, 0),
            std::ops::Bound::Unbounded => 0,
        };
        let end = match range.end_bound() {
            std::ops::Bound::Included(&e) => clip(e + 1, NS - 1),
            std::ops::Bound::Excluded(&e) => clip(e, NS - 1),
            std::ops::Bound::Unbounded => NS - 1,
        };

        if end <= start {
            panic!("End of range must be greater than start");
        }

        SpectralLocusSlopeIterator {
            locus: self,
            index: start,
            step,
            end,
        }
    }

    pub fn polygon(&self) -> &Polygon {
        &self.polygon
    }
}

pub struct SpectralLocusIterator<'a> {
    locus: &'a SpectralLocus,
    index: usize,
}

impl<'a> IntoIterator for &'a SpectralLocus {
    type Item = (f64, f64);
    type IntoIter = SpectralLocusIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        SpectralLocusIterator {
            locus: self,
            index: 0,
        }
    }
}

impl<'a> Iterator for SpectralLocusIterator<'a> {
    type Item = (f64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.locus.polygon.exterior().points().count() {
            let coord = self.locus.polygon.exterior()[self.index];
            self.index += 1;
            Some((coord.x, coord.y))
        } else {
            None
        }
    }
}

pub struct SpectralLocusSlopeIterator<'a> {
    locus: &'a SpectralLocus,
    index: usize,
    step: usize,
    end: usize,
}

impl Iterator for SpectralLocusSlopeIterator<'_> {
    type Item = ((f64, f64), f64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.end {
            let coord = self.locus.polygon.exterior()[self.index];
            let angle = match self.index {
                0 => {
                    let coord_plus = self.locus.polygon.exterior()[1];
                    (coord_plus.y - coord.y).atan2(coord_plus.x - coord.x)
                }
                i if i < NS - 1 => {
                    let coord_min = self.locus.polygon.exterior()[i - 1];
                    let coord_plus = self.locus.polygon.exterior()[i + 1];
                    let qurve = QuadraticCurve::new(
                        (coord_min.x, coord_min.y),
                        (coord.x, coord.y),
                        (coord_plus.x, coord_plus.y),
                    )
                    .unwrap();
                    qurve.slope_angle(0.5)
                }
                n if n == NS - 1 => {
                    let coord_min = self.locus.polygon.exterior()[n - 1];
                    (coord.y - coord_min.y).atan2(coord.x - coord_min.x)
                }
                _ => 0.0, // Should not be reached given the if self.index < self.end condition
            };
            self.index += self.step;
            Some(((coord.x, coord.y), angle))
        } else {
            None
        }
    }
}

impl<'a> DoubleEndedIterator for SpectralLocusSlopeIterator<'a> {
    fn next_back(&mut self) -> Option<<Self as Iterator>::Item> {
        if self.index < self.end && self.end > 0 {
            self.end = self.end.saturating_sub(self.step);
            let coord = self.locus.polygon.exterior()[self.end];
            let angle = match self.end {
                0 => {
                    let coord_plus = self.locus.polygon.exterior()[1];
                    (coord_plus.y - coord.y).atan2(coord_plus.x - coord.x)
                }
                i if i < NS - 1 => {
                    let coord_min = self.locus.polygon.exterior()[i - 1];
                    let coord_plus = self.locus.polygon.exterior()[i + 1];
                    let qurve = QuadraticCurve::new(
                        (coord_min.x, coord_min.y),
                        (coord.x, coord.y),
                        (coord_plus.x, coord_plus.y),
                    )
                    .unwrap();
                    qurve.slope_angle(0.5)
                }
                n if n == NS - 1 => {
                    let coord_min = self.locus.polygon.exterior()[n - 1];
                    (coord.y - coord_min.y).atan2(coord.x - coord_min.x)
                }
                _ => 0.0,
            };
            Some(((coord.x, coord.y), angle))
        } else {
            None
        }
    }
}

// make sure Vec is closed
use geo::{Coord, LineString};
fn vec_to_polygon(coords: Vec<[f64; 2]>) -> Polygon {
    let linestring = LineString::from(
        coords
            .into_iter()
            .map(|[x, y]| Coord { x, y })
            .collect::<Vec<_>>(),
    );

    Polygon::new(linestring, vec![]) // no interior rings
}

#[cfg(test)]
mod test {
    use approx::assert_abs_diff_eq;

    use super::*;
    use crate::observer::Observer::Cie1931;

    #[test]
    fn test_spectral_locus() {
        let locus = SpectralLocus::new(Cie1931);
        assert_eq!(locus.observer(), Cie1931);
        assert!(locus.contains([0.3, 0.3]));
        assert!(!locus.contains([0.05, 0.05]));
        assert!(!locus.contains([0.7, 0.7]));

        let points: Vec<(f64, f64)> = locus.into_iter().collect();
        assert_eq!(points.len(), NS + 1); // closed polygon
    }

    #[test]
    fn test_spectral_locus_iterator() {
        let locus = SpectralLocus::new(Cie1931);
        let mut iter = locus.into_iter();
        let point = iter.next().unwrap();
        assert_abs_diff_eq!(point.0, 0.17411, epsilon = 0.00005);
        assert_abs_diff_eq!(point.1, 0.00496, epsilon = 0.00005);
        let last = iter.last().unwrap();
        assert_abs_diff_eq!(last.0, 0.17411, epsilon = 0.00005);
        assert_abs_diff_eq!(last.1, 0.00496, epsilon = 0.00005);
    }

    #[test]
    fn test_spectral_locus_slope_iterator() {
        let locus = SpectralLocus::new(Cie1931);
        let mut iter = locus.iter_range_with_slope(400..700, 10);

        // Check the first point
        let (point, angle) = iter.next().unwrap();
        assert_abs_diff_eq!(point.0, 0.173336, epsilon = 0.00005);
        assert_abs_diff_eq!(point.1, 0.0047967, epsilon = 0.00005);
        assert_abs_diff_eq!(angle, -2.84106, epsilon = 0.00005); // slope at first point

        // Check the last point
        let (point, angle) = iter.next_back().unwrap();
        assert_abs_diff_eq!(point.0, 0.734390, epsilon = 0.00005);
        assert_abs_diff_eq!(point.1, 0.265609, epsilon = 0.00005);
        assert_abs_diff_eq!(angle, -0.785398, epsilon = 0.00005); // slope at last point
    }
}
