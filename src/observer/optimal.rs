use std::collections::HashMap;
use nalgebra::{DMatrix, Vector3};

use super::Observer;
use crate::{lab::CieLab, prelude::CieIlluminant, spectrum::NS, xyz::{RelXYZ, XYZ}};

pub struct OptimalColors(XYZ, DMatrix<Vector3<f64>>);

impl Observer {
    pub fn optimal_colors(&self, ref_white: CieIlluminant) -> OptimalColors {
        let mut optcol: DMatrix<Vector3<f64>> = DMatrix::from_fn(NS, NS, |_, _| Vector3::zeros());
        let spectral_locus = self.spectral_locus(ref_white);
        let white_point = spectral_locus[0].1.white_point();
        let white_point_vec = white_point.xyz;

        // fill top left, including diagonal, with direct colors
        for r in 0..NS {
            for c in 0..NS - r {
                for w in 0..=r {
                    optcol[(r, c)] += spectral_locus[c + w].1.xyz().xyz;
                }
            }
        }

        // fill bottom right triangle (inverse colors)
        for r in 0..NS {
            for c in 0..NS {
                if r + c <= NS - 1 {
                    let rr = NS - 1 - r;
                    let cc = NS - 1 - c;
                    optcol[(rr, cc)] = white_point_vec - optcol[(r, c)];
                }
            }
        }
        OptimalColors(white_point, optcol)
    }
}

impl OptimalColors {
    pub fn white_point(&self) -> XYZ {
        self.0
    }

    pub fn colors(&self) -> &DMatrix<Vector3<f64>> {
        &self.1
    }

    pub fn observer(&self) -> Observer {
        self.0.observer
    }

    /// Returns a map of maximum chroma values for each (L, H) pair in the CIE LCh color space.
    /// The keys are tuples of (L, H) where L is the lightness (0-100) and H is the hue angle divided by 2 (0-180).
    /// The values are the maximum chroma values (C) for those (L, H) pairs.
    /// This is useful for determining the gamut of colors that can be represented in the CIE LCh color space.
    pub fn max_chroma(&self) -> HashMap<(u8, u8), f32 > {
        let mut map = HashMap::new();
        let white_point = self.0;
        for r in 0..NS {
            for c in 0..NS {
                let xyz = self.1[(r, c)];
                let rel_xyz= RelXYZ::new(xyz.into(), white_point);
                let lab = CieLab::from_xyz(rel_xyz);
                let [l, chroma, h] = lab.lch();
                let l = l.round() as u8;
                let h = (h/2.0).round() as u8;
                let c_f32 = chroma as f32;
                map.entry((l,h)).and_modify(|existing| {
                    if *existing < c_f32 {
                        *existing = c_f32;
                    }
                })
                .or_insert(c_f32);
            }
        }
        map
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observer::Observer::Cie1931;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_optimum_colors() {
        let observer = Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = observer.optimal_colors(ref_white);
        let xyz_ref_white = observer.xyz_d65();

        let first: Vec<XYZ> = observer
            .spectral_locus(ref_white)
            .iter()
            .map(|(_w, xyz)| xyz.xyz())
            .collect();

        first
            .iter()
            .map(|xyz| xyz.xyz)
            .zip(opt_colors.colors().row(0).iter())
            .for_each(|(x, y)| {
                assert_abs_diff_eq!(x, y, epsilon = 1e-10);
            });
        
        let last: Vec<XYZ> = observer
            .spectral_locus(ref_white)
            .iter()
            .rev()
            .map(|(_w, xyz)| xyz_ref_white - xyz.xyz())
            .collect();

        last 
            .iter()
            .map(|xyz| xyz.xyz)
            .zip(opt_colors.colors().row(NS-1).iter())
            .for_each(|(x, y)| {
                println!("{:.5} {:.5} {:.5} {:.5} {:.5} {:.5}", x.x, x.y, x.z, y.x, y.y, y.z);
                assert_abs_diff_eq!(x, y, epsilon = 1e-10);
            });

    }

    #[test]
    fn test_cielch() {
        let observer = Cie1931;
        let ref_white = CieIlluminant::D65;
        let opt_colors = observer.optimal_colors(ref_white);
        let cielch = opt_colors.max_chroma();
        for l in 0..=100 {
            for h in 0..=180 {
                if let Some(c) = cielch.get(&(l, h)) {
                    println!("{}, {}, {}", l, h, c);
                }
            }
        }   
    }
}



