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
use geo::{Contains, Polygon};

use crate::{observer::Observer, spectrum::NS, xyz::XYZ};

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
}
